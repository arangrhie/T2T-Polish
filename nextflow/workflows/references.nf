/*
 * workflows/references.nf
 *
 * BUILD_REFS sub-workflow:
 *   - Concatenates input pieces into hap1, hap2, and dip FA.GZ references.
 *   - Runs MERYL_REPETITIVE and BWA_INDEX in parallel on all three refs,
 *     so that both long-read and short-read mapping can start as soon as
 *     each reference is ready.
 *
 * BUILD_REFS_FROM_FILES sub-workflow (polishing rounds 2+):
 *   - Takes already-polished plain FASTA files ([hap, fa, gzi, fai]) from
 *     PREPARE_NEXT_ROUND — no BUILD_HAP_REFERENCES step needed.
 *   - Runs MERYL_REPETITIVE and BWA_INDEX in parallel, identical to BUILD_REFS.
 *
 * Emits:
 *   wm_refs  — [ hap, ref_fa_gz, ref_fai, repetitive_k15.txt ]
 *   bwa_refs — [ hap, ref_fa_gz, ref_fai, amb, ann, bwt, pac, sa ]
 */

include { BUILD_HAP_REFERENCES; BUILD_DIP_REFERENCE } from '../modules/references'
include { MERYL_REPETITIVE                } from '../modules/winnowmap'
include { BWA_INDEX                       } from '../modules/bwa'

workflow BUILD_REFS {

    main:

    def no_file = file('NO_FILE')
    def ebv     = params.ebv_fasta_gz ? file(params.ebv_fasta_gz) : no_file

    // Re-entry: if assemblies_dir is set, infer the pre-built FA.GZ paths from
    // asm_name + asm_ver and load them directly, skipping BUILD_HAP_REFERENCES
    // and BUILD_DIP_REFERENCE.  The .gzi and .fai files must exist alongside.
    def allRefs   // [ hap, ref_fa_gz, ref_gzi, ref_fai ]
    if ( params.assemblies_dir ) {
        def d    = params.assemblies_dir.replaceAll('/$', '')
        def pfx  = "${d}/${params.asm_name}_${params.asm_ver}"
        allRefs = Channel.of(
            tuple('hap1', file("${pfx}.hap1.fa.gz"), file("${pfx}.hap1.fa.gz.gzi"), file("${pfx}.hap1.fa.gz.fai")),
            tuple('hap2', file("${pfx}.hap2.fa.gz"), file("${pfx}.hap2.fa.gz.gzi"), file("${pfx}.hap2.fa.gz.fai")),
            tuple('dip',  file("${pfx}.dip.fa.gz"),  file("${pfx}.dip.fa.gz.gzi"),  file("${pfx}.dip.fa.gz.fai"))
        )
    } else {
        // ---- hap1 + hap2 ---------------------------------------------------
        def refInputs = Channel.of(
            tuple('hap1', file(params.hap1_fasta_gz),
                  file(params.mito_exemplar_fasta_gz), ebv,
                  file(params.rdna_exemplar_fasta_gz)),
            tuple('hap2', file(params.hap2_fasta_gz),
                  file(params.mito_exemplar_fasta_gz), ebv,
                  file(params.rdna_exemplar_fasta_gz))
        )
        def builtHaps = BUILD_HAP_REFERENCES(refInputs)

        // ---- dip (hap1 + hap2 + mito + [ebv] + rdna in one step) ----------
        def dipRef = BUILD_DIP_REFERENCE(
            file(params.hap1_fasta_gz),
            file(params.hap2_fasta_gz),
            file(params.mito_exemplar_fasta_gz),
            ebv,
            file(params.rdna_exemplar_fasta_gz)
        )

        allRefs = builtHaps.mix(dipRef)
    }

    // allRefs: [ hap, ref_fa_gz, ref_gzi, ref_fai ]  (hap1, hap2, dip)

    // Drop the .gzi — MERYL_REPETITIVE and BWA_INDEX only need [ hap, fa, fai ]
    // Inject params.asm_ver as ver_from so downstream processes can name outputs correctly.
    def allRefs_4 = allRefs.map { hap, fa, gzi, fai -> tuple(hap, params.asm_ver, fa, fai) }

    // ---- prepare for long-read mapping: build meryl repetitive kmers -------
    // ---- prepare for short-read mapping: build BWA index -------------------
    // If assemblies_dir is set, the index and repetitive files already exist
    // alongside the FAs — construct the output channels statically, no process needed.
    def wm_refs_ch
    def bwa_refs_ch
    if ( params.assemblies_dir ) {
        def d   = params.assemblies_dir.replaceAll('/$', '')
        def pfx = "${d}/${params.asm_name}_${params.asm_ver}"
        wm_refs_ch = allRefs_4.map { hap, ver, fa, fai ->
            def rep = file("${fa.parent}/${params.asm_name}_${ver}.${hap}.repetitive_k15.txt")
            if ( !rep.exists() ) error "assemblies_dir: missing repetitive file: ${rep}"
            tuple(hap, ver, fa, fai, rep)
        }
        bwa_refs_ch = allRefs_4.map { hap, ver, fa, fai ->
            def idxFiles = ['amb','ann','bwt','pac','sa'].collect { ext -> file("${fa}.${ext}") }
            idxFiles.each { f -> if ( !f.exists() ) error "assemblies_dir: missing BWA index file: ${f}" }
            tuple(hap, ver, fa, fai, idxFiles[0], idxFiles[1], idxFiles[2], idxFiles[3], idxFiles[4])
        }
    } else {
        // allRefs_4 is consumed by two downstream helpers — fork it first.
        def allRefs_split = allRefs_4.multiMap { it -> meryl: it; bwa: it }
        wm_refs_ch  = _meryl_repetitive_or_reuse(allRefs_split.meryl).refs
        bwa_refs_ch = _bwa_index_or_reuse(allRefs_split.bwa).refs
    }

    emit:
    wm_refs  = wm_refs_ch   // [ hap, ver_from, ref_fa_gz, ref_fai, repetitive_k15.txt ]
    bwa_refs = bwa_refs_ch  // [ hap, ver_from, ref_fa_gz, ref_fai, amb, ann, bwt, pac, sa ]
    dv_refs  = allRefs.map { hap, fa, gzi, fai -> tuple(hap, params.asm_ver, fa, gzi, fai) }
              // [ hap, ver_from, ref_fa_gz, ref_fa_gzi, ref_fai ] — .gzi required by make_examples
}

/*
 * BUILD_REFS_FROM_FILES — used for polishing rounds 2, 3, … N.
 *
 * Takes a channel of already-polished [ hap, fa, gzi, fai ] tuples produced by
 * PREPARE_NEXT_ROUND and builds the Winnowmap (meryl) and BWA indexes in
 * parallel, emitting the same shapes as BUILD_REFS.
 *
 * No BUILD_HAP_REFERENCES step is needed: the polished FASTAs are already
 * complete assemblies (mito/EBV/rDNA were embedded in round 1).
 */
workflow BUILD_REFS_FROM_FILES {

    take:
    hap_refs  // [ hap, fa, gzi, fai ] — polished hap1, hap2, dip from PREPARE_NEXT_ROUND
    ver_from  // String — version label for this round (e.g. 'v0.2')

    main:

    // Inject ver_from; drop the .gzi — MERYL_REPETITIVE only needs [ hap, ver, fa, fai ]
    def hap_refs_4 = hap_refs.map { hap, fa, gzi, fai -> tuple(hap, ver_from, fa, fai) }

    // hap_refs_4 is consumed by two downstream helpers — fork it first.
    def hap_refs_split = hap_refs_4.multiMap { it -> meryl: it; bwa: it }

    def wm_refs_ch  = _meryl_repetitive_or_reuse(hap_refs_split.meryl).refs

    // Pass ver_from through _bwa_index_or_reuse — no re-injection combine needed.
    def bwa_refs_ch = _bwa_index_or_reuse(hap_refs_split.bwa).refs

    emit:
    wm_refs  = wm_refs_ch
    bwa_refs = bwa_refs_ch
    dv_refs  = hap_refs.map { hap, fa, gzi, fai -> tuple(hap, ver_from, fa, gzi, fai) }
              // [ hap, ver_from, ref_fa_gz, ref_fa_gzi, ref_fai ]
}

/*
 * Helper: skip BWA_INDEX if all five index files already exist alongside
 * the reference FA.GZ (whether in assemblies_dir or outdir/assemblies/).
 *
 * Takes:  [ hap, ver_from, ref_fa_gz, ref_fai ]
 * Emits:  [ hap, ver_from, ref_fa_gz, ref_fai, amb, ann, bwt, pac, sa ]
 */
workflow _bwa_index_or_reuse {

    take:
    refs_4  // [ hap, ver_from, ref_fa_gz, ref_fai ]

    main:

    // Split into two branches based on whether the index already exists.
    def branched = refs_4.branch { hap, ver, fa, fai ->
        // Check next to the FA itself (works for both assemblies_dir and outdir/assemblies)
        def idxFiles = ['amb','ann','bwt','pac','sa'].collect { ext -> file("${fa}.${ext}") }
        ready: fa.exists() && idxFiles.every { it.exists() }
            return tuple(hap, ver, fa, fai, idxFiles[0], idxFiles[1], idxFiles[2], idxFiles[3], idxFiles[4])
        index: true
    }

    def bwa_out = BWA_INDEX( branched.index ).mix( branched.ready )

    emit:
    refs = bwa_out
}

/*
 * Helper: skip MERYL_REPETITIVE if the repetitive_k15 file already exists
 * alongside the reference FA.GZ (whether in assemblies_dir or outdir/assemblies/).
 * Expected filename: {asm_name}_{ver_from}.{hap}.repetitive_k15.txt
 *
 * Takes:  [ hap, ver_from, ref_fa_gz, ref_fai ]
 * Emits:  [ hap, ver_from, ref_fa_gz, ref_fai, repetitive_k15.txt ]
 */
workflow _meryl_repetitive_or_reuse {

    take:
    refs_4  // [ hap, ver_from, ref_fa_gz, ref_fai ]

    main:

    def branched = refs_4.branch { hap, ver_from, fa, fai ->
        // Check next to the FA itself (works for both assemblies_dir and outdir/assemblies)
        def rep = file("${fa.parent}/${params.asm_name}_${ver_from}.${hap}.repetitive_k15.txt")
        ready: rep.exists()
            return tuple(hap, ver_from, fa, fai, rep)
        build: true
    }

    def meryl_out = MERYL_REPETITIVE( branched.build ).mix( branched.ready )

    emit:
    refs = meryl_out
}
