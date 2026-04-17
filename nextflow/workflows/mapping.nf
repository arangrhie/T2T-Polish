/*
 * workflows/mapping.nf
 *
 * MAPPING sub-workflow.
 *
 * Each process (WINNOWMAP_MAP, WINNOWMAP_MERGE, etc.) is called exactly once.
 * All platforms' reads are combined into a single channel before each process
 * call so NF 25.x sees only one invocation per process.
 *
 * Inputs (from BUILD_REFS):
 *   wm_refs  — [ hap, ver_from, ref_fa_gz, ref_fai, repetitive_k15.txt ]
 *   bwa_refs — [ hap, ver_from, ref_fa_gz, ref_fai, amb, ann, bwt, pac, sa ]
 */

include { WINNOWMAP_MAP; WINNOWMAP_MERGE;
          WINNOWMAP_FILTER; SAM2PAF     } from '../modules/winnowmap'
include { BWA_MAP; BWA_MERGE            } from '../modules/bwa'

// ---------------------------------------------------------------------------
def requireParam(String name, Object value) {
    if ( value == null || value.toString().trim().isEmpty() )
        error "Missing required parameter: params.${name}"
}

// ---------------------------------------------------------------------------
workflow MAPPING {
// ---------------------------------------------------------------------------
    take:
    wm_refs   // [ hap, ver_from, ref_fa_gz, ref_fai, repetitive_k15.txt ]
    bwa_refs  // [ hap, ver_from, ref_fa_gz, ref_fai, amb, ann, bwt, pac, sa ]

    main:

    // ---- Winnowmap (HiFi + ONT) --------------------------------------------
    // Build one unified reads channel across all long-read platforms, then
    // call each process exactly once.  The platform val in the tuple routes
    // outputs to the right publishDir/filename at runtime.

    def wmPlatforms = ['hifi', 'ont'].findAll { params.platforms.contains(it) }

    // Collect all [plat, reads_file] pairs into one channel
    def allWmReads = Channel.empty()
    wmPlatforms.each { plat ->
        def glob = (plat == 'hifi') ? params.read_glob_hifi : params.read_glob_ont
        requireParam("read_glob_${plat}", glob)
        allWmReads = allWmReads.mix(
            Channel.fromPath(glob).map { f -> tuple(plat, f) }
        )
    }

    def wm_pri_bams = Channel.empty()

    if ( wmPlatforms ) {
        // For ONT (when ont_map_haps=false), only use dip refs; for HiFi use all.
        // We split by platform inside the combine to apply the filter per-plat.
        def hifiReads = allWmReads.filter { plat, f -> plat == 'hifi' }
        def ontReads  = allWmReads.filter { plat, f -> plat == 'ont'  }

        def hifiInput = wm_refs.combine(hifiReads)
        def ontRefs   = params.ont_map_haps ? wm_refs
                                            : wm_refs.filter { hap, ver, ref, fai, rep -> hap == 'dip' }
        def ontInput  = ontRefs.combine(ontReads)

        // Single WINNOWMAP_MAP call — all platforms combined
        def mapped   = WINNOWMAP_MAP( hifiInput.mix(ontInput) )
        def merged   = WINNOWMAP_MERGE( mapped.groupTuple(by: [0, 1, 2]) )
        def filtered = WINNOWMAP_FILTER(merged)
        SAM2PAF(filtered)
        wm_pri_bams = filtered
    }

    // ---- BWA (Illumina + Element) ------------------------------------------
    def bwaPlatforms = ['illumina', 'element'].findAll { params.platforms.contains(it) }

    def allBwaReads = Channel.empty()
    bwaPlatforms.each { plat ->
        def glob = (plat == 'illumina') ? params.read_glob_illumina : params.read_glob_element
        requireParam("read_glob_${plat}", glob)
        allBwaReads = allBwaReads.mix(
            Channel.fromFilePairs(glob).map { sid, files -> tuple(plat, files[0], files[1]) }
        )
    }

    def bwa_mrg_bams = Channel.empty()

    if ( bwaPlatforms ) {
        // Single BWA_MAP call — all short-read platforms combined
        def mapped = BWA_MAP( bwa_refs.combine(allBwaReads) )
        def merged = BWA_MERGE( mapped.groupTuple(by: [0, 1, 2]) )
        bwa_mrg_bams = merged
    }

    emit:
    wm_pri_bams  = wm_pri_bams   // [ hap, ver_from, platform, bam, bai ]
    bwa_mrg_bams = bwa_mrg_bams  // [ hap, ver_from, platform, bam, csi ]
}
