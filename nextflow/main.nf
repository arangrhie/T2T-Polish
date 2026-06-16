#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * T2T-Polish mapping pipeline — entry point
 *
 * Usage:
 *   nextflow run main.nf -c user.config [-profile slurm] [--platforms hifi,ont]
 *
 * All process definitions live under modules/.
 * The mapping orchestration logic lives in workflows/mapping.nf.
 * Resource labels and Slurm settings are in resources.config (loaded by nextflow.config).
 */

// ------------ parameters (defaults; override in user.config) ---------------
params.outdir          = params.outdir    ?: 'results'
params.asm_name        = params.asm_name  ?: 'assembly'
params.asm_ver         = params.asm_ver   ?: 'v0.1'
if ( !params.asm_name ) error "params.asm_name must be set (e.g. 'bTaeGut7')"
if ( !params.asm_ver  ) error "params.asm_ver must be set (e.g. 'v0.1')"
// Full filename prefix used throughout: {asm_name}_{asm_ver}  (e.g. 'bTaeGut7_v0.1')
params.platforms       = (params.platforms ?: 'hifi,ont,illumina,element')
                           .tokenize(',').collect { it.trim() }.findAll { it }

params.hap1_fasta_gz          = params.hap1_fasta_gz          ?: null
params.hap2_fasta_gz          = params.hap2_fasta_gz          ?: null
params.ebv_fasta_gz           = params.ebv_fasta_gz           ?: null
params.mito_exemplar_fasta_gz = params.mito_exemplar_fasta_gz ?: null
params.rdna_exemplar_fasta_gz = params.rdna_exemplar_fasta_gz ?: null

// Normalise globs: treat the string 'null' (common config mistake) as real null
params.read_glob_hifi     = (params.read_glob_hifi     in [null,'null','']) ? null : params.read_glob_hifi
params.read_glob_ont      = (params.read_glob_ont      in [null,'null','']) ? null : params.read_glob_ont
params.read_glob_illumina = (params.read_glob_illumina in [null,'null','']) ? null : params.read_glob_illumina
params.read_glob_element  = (params.read_glob_element  in [null,'null','']) ? null : params.read_glob_element

params.samtools          = params.samtools          ?: 'samtools'
params.k8                = params.k8                ?: 'k8'
params.mapping_outdir    = params.mapping_outdir    ?: "${params.outdir}/mapping"
params.keep_intermediates = params.keep_intermediates ?: false
params.ont_map_haps       = params.ont_map_haps       ?: false  // map ONT to hap1/hap2 in addition to dip

// Re-entry directories
params.assemblies_dir   = params.assemblies_dir   ?: null
params.mapping_dir      = params.mapping_dir      ?: null
params.deepvariant_dir  = params.deepvariant_dir  ?: null

// DeepVariant parameters
params.dv_outdir          = params.dv_outdir          ?: "${params.outdir}/deepvariant"
params.dv_sample          = params.dv_sample          ?: null   // sample name in VCF; defaults to hap tag
params.dv_mq_hap          = params.dv_mq_hap          ?: 0      // MQ filter for hap1/hap2 merged BAMs
params.dv_mq_dip          = params.dv_mq_dip          ?: 5      // MQ filter for dip   merged BAMs
params.ont_chemistry      = params.ont_chemistry      ?: 'r10'  // 'r10' or 'r9'
params.dv_n_shard         = params.dv_n_shard         ?: 12
params.dv_long_platforms  = params.dv_long_platforms  ?: 'hifi' // long-read platforms to use in hybrid merge
params.dv_short_platforms = params.dv_short_platforms ?: 'illumina,element' // short-read platforms

// SNV candidates / polishing parameters
// Disable with --run_snv_candidates false (and run_dv_ont = true for the ONT track).
params.run_snv_candidates    = params.containsKey('run_snv_candidates') ? params.run_snv_candidates : true
params.snv_outdir            = params.snv_outdir            ?: "${params.outdir}/snv_candidates"
params.hybrid_meryl          = params.hybrid_meryl          ?: null   // path to hybrid read k-mer meryl DB
params.merfin_peak           = params.merfin_peak           ?: null   // integer peak coverage for Merfin
params.merfin                = params.merfin                ?: 'merfin'
params.asm_ver_next          = params.asm_ver_next          ?: null   // first polished version; auto-bumped if null
params.keep_snv_intermediates = params.keep_snv_intermediates ?: false
// Number of SNV polishing rounds to run (default 2).
// Each round maps reads to the polished assembly from the previous round,
// calls variants with DeepVariant, and applies Merfin-validated SNVs.
params.polish_rounds         = params.polish_rounds         ?: 2

// ------------ startup validation --------------------------------------------
// Resolve effective short-read platform set (union of --platforms and --dv_short_platforms).
def _plat_list   = params.platforms instanceof List ? params.platforms
                       : params.platforms.tokenize(',').collect { it.trim() }.findAll { it }
def _short_plat  = params.dv_short_platforms.tokenize(',').collect { it.trim() }.findAll { it }
def _active_short = (_plat_list.intersect(['illumina', 'element']) + _short_plat.findAll { it in ['illumina','element'] }).unique()

if ( !params.platforms.contains('hifi') || !params.read_glob_hifi )
    error "HiFi reads are required.\n" +
          "Add 'hifi' to --platforms and set params.read_glob_hifi in user.config."

if ( !params.platforms.contains('ont') || !params.read_glob_ont )
    error "ONT reads are required.\n" +
          "Add 'ont' to --platforms and set params.read_glob_ont in user.config."

if ( _active_short && !params.read_glob_illumina && !params.read_glob_element )
    error "At least one short-read glob must be provided when 'illumina' or 'element' is active.\n" +
          "Set params.read_glob_illumina and/or params.read_glob_element in user.config,\n" +
          "or exclude both with --platforms hifi,ont --dv_short_platforms ''."

// ------------ imports -------------------------------------------------------
include { BUILD_REFS                                  } from './workflows/references'
include { BUILD_REFS_FROM_FILES as BUILD_REFS_R2      } from './workflows/references'
include { BUILD_REFS_FROM_FILES as BUILD_REFS_R3      } from './workflows/references'
include { BUILD_REFS_FROM_FILES as BUILD_REFS_R4      } from './workflows/references'
include { BUILD_REFS_FROM_FILES as BUILD_REFS_R5      } from './workflows/references'

include { MAPPING_R1  } from './workflows/mapping_r1'
include { MAPPING_R2  } from './workflows/mapping_r2'
include { MAPPING_R3  } from './workflows/mapping_r3'
include { MAPPING_R4  } from './workflows/mapping_r4'
include { MAPPING_R5  } from './workflows/mapping_r5'

include { DEEPVARIANT_R1 } from './workflows/deepvariant_r1'
include { DEEPVARIANT_R2 } from './workflows/deepvariant_r2'
include { DEEPVARIANT_R3 } from './workflows/deepvariant_r3'
include { DEEPVARIANT_R4 } from './workflows/deepvariant_r4'
include { DEEPVARIANT_R5 } from './workflows/deepvariant_r5'

include { SNV_CANDIDATES as SNV_CANDIDATES_R1 } from './workflows/snv_candidates'
include { SNV_CANDIDATES as SNV_CANDIDATES_R2 } from './workflows/snv_candidates'
include { SNV_CANDIDATES as SNV_CANDIDATES_R3 } from './workflows/snv_candidates'
include { SNV_CANDIDATES as SNV_CANDIDATES_R4 } from './workflows/snv_candidates'
include { SNV_CANDIDATES as SNV_CANDIDATES_R5 } from './workflows/snv_candidates'

// ------------ helpers -------------------------------------------------------
def requireParam(String name, Object value) {
    if ( value == null || value.toString().trim().isEmpty() )
        error "Missing required parameter: params.${name}"
}

/**
 * Increment the trailing integer in a version string.
 *   'v0.1'           → 'v0.2'
 *   'assembly_v1.3'  → 'assembly_v1.4'
 *   'myasm'          → 'myasm_polished'  (no trailing digit — fallback)
 */
def bumpVersion(String ver) {
    def m = (ver =~ /^(.*\D)(\d+)$/)
    m ? "${m[0][1]}${(m[0][2] as int) + 1}" : "${ver}_polished"
}

// ------------ entry workflow ------------------------------------------------
workflow {

    // Verkko FASTA inputs are only required when pre-built refs are not provided
    def _refs_provided = params.assemblies_dir as boolean
    if ( !_refs_provided ) {
        requireParam('hap1_fasta_gz',          params.hap1_fasta_gz)
        requireParam('hap2_fasta_gz',          params.hap2_fasta_gz)
    }

    def rounds = (params.polish_rounds as int)
    if ( rounds < 1 ) error "params.polish_rounds must be ≥ 1 (got ${params.polish_rounds})"
    if ( rounds > 5 ) error "params.polish_rounds > 5 is not supported (got ${params.polish_rounds})"
    if ( ! params.run_dv ) {
        log.warn "DeepVariant stages are disabled (--run_dv false); overriding params.polish_rounds to 1."
        rounds = 1
    }

    def ver = []
    ver[0] = params.asm_ver ?: 'v0.1'
    ver[1] = params.asm_ver_next ?: bumpVersion(ver[0])
    (2..5).each { i -> ver[i] = bumpVersion(ver[i - 1]) }

    // =========================================================================
    // Each round is fully sequential: BUILD_REFS → MAPPING → DEEPVARIANT →
    // SNV_CANDIDATES.  Per-round workflow files carry aliased process names
    // (_R1 … _R5) so NF 25.x sees each process bound exactly once.
    // =========================================================================

    // ---- Round 1 ------------------------------------------------------------
    BUILD_REFS()
    def r1_refs    = BUILD_REFS.out.wm_refs.map { hap, ver_from, ref, fai, _rep -> tuple(hap, ver_from, ref, fai) }
    def r1_dv_refs = BUILD_REFS.out.dv_refs  // [ hap, ver_from, ref_fa_gz, ref_fa_gzi, ref_fai ]

    // Re-entry: supply mapping_dir to skip mapping stages, or deepvariant_dir
    // to skip make_examples.  Both use the {asm_name}_{asm_ver} prefix to
    // discover files automatically.
    // When only deepvariant_dir is set, MAPPING will be skipped.
    def _pfx = "${params.asm_name}_${params.asm_ver}"

    def r1_wm_bams
    def r1_bwa_bams
    if ( !params.deepvariant_dir && params.mapping_dir ) {
        def mdir = params.mapping_dir.replaceAll('/$', '')

        // ---- Detect which (hap, platform) merged BAMs already exist on disk ----
        // Done at launch time (plain Groovy) so we can gate the ref channels
        // without needing any extra params.
        // NOTE: use new File() (java.io.File), NOT file() — the Nextflow file()
        // helper returns a java.nio.file.Path which lacks listFiles().
        // Scan only subdirs matching this version's prefix to avoid picking up
        // BAMs from other versions that happen to share the same hap/platform.
        def _mdir_file = new File(mdir)
        def _pfx_subdirs = (_mdir_file.listFiles() ?: [])
            .findAll { it.isDirectory() && it.name.startsWith(_pfx) }

        def wm_on_disk = _pfx_subdirs
            .collectMany { dir -> dir.listFiles()?.toList() ?: [] }
            .findAll { it.name ==~ /.*\.([^.]+)\.([^.]+)\.pri\.bam$/ &&
                      !(it.name ==~ /.*\.dedup\.pri\.bam$/) }
            .collect { f ->
                def m = (f.name =~ /\.([^.]+)\.([^.]+)\.pri\.bam$/)[0]
                [m[1], m[2]]
            } as Set

        def bwa_on_disk = _pfx_subdirs
            .collectMany { dir -> dir.listFiles()?.toList() ?: [] }
            .findAll { it.name ==~ /.*\.([^.]+)\.([^.]+)\.dedup\.pri\.bam$/ }
            .collect { f ->
                def m = (f.name =~ /\.([^.]+)\.([^.]+)\.dedup\.pri\.bam$/)[0]
                [m[1], m[2]]
            } as Set

        log.info "mapping_dir scan: found ${_pfx_subdirs.size()} subdirs for prefix '${_pfx}'"
        log.info "mapping_dir scan: wm_on_disk  = ${wm_on_disk}"
        log.info "mapping_dir scan: bwa_on_disk = ${bwa_on_disk}"

        // A hap is "done" only when ALL expected platforms for that hap are present.
        // Derive the expected platform sets from what's actually on disk (union across
        // all haps), rather than from params.platforms — this way the check is correct
        // even when the current invocation doesn't list all previously-run platforms.
        // ONT exception: ONT is only expected for hap1/hap2 if params.ont_map_haps=true.
        def _found_long_plats  = wm_on_disk.collect { it[1] }.unique()
        def _found_short_plats = bwa_on_disk.collect { it[1] }.unique()

        log.info "mapping_dir scan: found long-read platforms on disk:  ${_found_long_plats}"
        log.info "mapping_dir scan: found short-read platforms on disk: ${_found_short_plats}"

        def wm_haps_done = ['hap1', 'hap2', 'dip'].findAll { hap ->
            def expected = _found_long_plats.findAll { plat ->
                plat == 'hifi' || hap == 'dip' || params.ont_map_haps
            }
            expected && expected.every { plat -> wm_on_disk.contains([hap, plat]) }
        }
        def bwa_haps_done = ['hap1', 'hap2', 'dip'].findAll { hap ->
            _found_short_plats && _found_short_plats.every { plat -> bwa_on_disk.contains([hap, plat]) }
        }

        log.info "mapping_dir re-entry: wm complete for haps ${wm_haps_done}, bwa complete for haps ${bwa_haps_done}"

        // Pass only haps that need (re-)mapping to MAPPING_R1.
        def r1_wm_refs_needed  = BUILD_REFS.out.wm_refs
            .filter { hap, v, ref, fai, rep -> !wm_haps_done.contains(hap) }
        def r1_bwa_refs_needed = BUILD_REFS.out.bwa_refs
            .filter { hap, v, ref, fai, amb, ann, bwt, pac, sa -> !bwa_haps_done.contains(hap) }

        MAPPING_R1( r1_wm_refs_needed, r1_bwa_refs_needed )

        // Load completed BAMs from disk (only haps that are fully done).
        // BAMs live in per-hap-platform subdirs: mdir/{_pfx}.{hap}.{plat}/{_pfx}.{hap}.{plat}.pri.bam
        // Pick whichever sidecar index exists (.bai or .csi); prefer .bai.
        def pickIdx = { bam ->
            def b = file("${bam}.bai")
            b.exists() ? b : file("${bam}.csi")
        }

        def r1_wm_from_disk = Channel.fromPath("${mdir}/${_pfx}.*.*/${_pfx}.*.*.pri.bam")
            .filter { !(it.name ==~ /.*\.dedup\.pri\.bam$/) }
            .map { f ->
                def m = (f.name =~ /\.([^.]+)\.([^.]+)\.pri\.bam$/)[0]
                tuple(m[1], params.asm_ver, m[2], f, pickIdx(f))
            }
            .filter { hap, v, plat, bam, idx -> wm_haps_done.contains(hap) }

        def r1_bwa_from_disk = Channel.fromPath("${mdir}/${_pfx}.*.*/${_pfx}.*.*.dedup.pri.bam")
            .map { f ->
                def m = (f.name =~ /\.([^.]+)\.([^.]+)\.dedup\.pri\.bam$/)[0]
                tuple(m[1], params.asm_ver, m[2], f, pickIdx(f))
            }
            .filter { hap, v, plat, bam, idx -> bwa_haps_done.contains(hap) }

        r1_wm_bams  = r1_wm_from_disk.mix( MAPPING_R1.out.wm_pri_bams )
        r1_bwa_bams = r1_bwa_from_disk.mix( MAPPING_R1.out.bwa_mrg_bams )
    } else if ( !params.deepvariant_dir ) {
        MAPPING_R1( BUILD_REFS.out.wm_refs, BUILD_REFS.out.bwa_refs )
        r1_wm_bams  = MAPPING_R1.out.wm_pri_bams
        r1_bwa_bams = MAPPING_R1.out.bwa_mrg_bams
    } else {
        log.info "deepvariant_dir is set; skipping MAPPING_R1 and mapping_dir re-entry for Round-1."
        r1_wm_bams  = Channel.empty()
        r1_bwa_bams = Channel.empty()
    }

    // Optional re-entry inventory for Round-1 DV/SNV outputs.
    def r1_dv_records      = []
    def r1_dv_missing_keys = []
    def r1_all_dv_present  = false
    def r1_snv_present     = false
    def r1_snv_next_refs_present = false

    if ( params.deepvariant_dir ) {
        def dvdir = params.deepvariant_dir.replaceAll('/$', '')
        def snvdir = params.snv_outdir.replaceAll('/$', '')
        def asmout = "${params.outdir}/assemblies".replaceAll('/$', '')

        def short_for_hybrid = _plat_list.contains('illumina') ? 'illumina' : (_plat_list.contains('element') ? 'element' : (_active_short ? _active_short[0] : 'illumina'))
        def hybrid_combo = "hifi_${short_for_hybrid}"
        def ont_mode = (params.ont_chemistry == 'r9') ? 'ONT_R9' : 'ONT_R104'

        def expectedDv = [
            [hap: 'hap1', combo: hybrid_combo, mode: 'HYBRID_PACBIO_ILLUMINA', mq: params.dv_mq_hap.toInteger()],
            [hap: 'hap2', combo: hybrid_combo, mode: 'HYBRID_PACBIO_ILLUMINA', mq: params.dv_mq_hap.toInteger()],
            [hap: 'dip',  combo: hybrid_combo, mode: 'HYBRID_PACBIO_ILLUMINA', mq: params.dv_mq_dip.toInteger()],
            [hap: 'dip',  combo: 'ont',        mode: ont_mode,                  mq: params.dv_mq_dip.toInteger()]
        ]

        log.info "deepvariant_dir inventory: scanning ${dvdir} for Round-1 expected DV outputs"
        expectedDv.each { e ->
            def runDir = "${dvdir}/${_pfx}.${e.hap}.${e.combo}.MQ${e.mq}"
            def vcfPath = "${runDir}/dv_${e.mode}_MQ${e.mq}.${e.hap}.vcf.gz"
            def tbiPath = "${vcfPath}.tbi"
            def gvcfPath = "${runDir}/dv_${e.mode}_MQ${e.mq}.${e.hap}.gvcf.gz"
            def gtbiPath = "${gvcfPath}.tbi"

            def vcfOk  = new File(vcfPath).exists()
            def tbiOk  = new File(tbiPath).exists()
            def gvcfOk = new File(gvcfPath).exists()
            def gtbiOk = new File(gtbiPath).exists()
            def key = "${e.hap}|${e.combo}|${e.mode}|MQ${e.mq}"

            if (vcfOk && tbiOk && gvcfOk && gtbiOk) {
                log.info "deepvariant_dir inventory: FOUND    ${key}"
                r1_dv_records << tuple(e.hap, params.asm_ver, e.combo, e.mode, e.mq,
                                       file(vcfPath), file(tbiPath), file(gvcfPath), file(gtbiPath))
            } else {
                def missingParts = []
                if (!vcfOk)  missingParts << 'vcf'
                if (!tbiOk)  missingParts << 'vcf.tbi'
                if (!gvcfOk) missingParts << 'gvcf'
                if (!gtbiOk) missingParts << 'gvcf.tbi'
                log.info "deepvariant_dir inventory: MISSING  ${key} -> ${missingParts.join(',')}"
                r1_dv_missing_keys << key
            }
        }

        r1_all_dv_present = (r1_dv_missing_keys.isEmpty())
        log.info "deepvariant_dir inventory: Round-1 DV complete=${r1_all_dv_present} (${r1_dv_records.size()}/${expectedDv.size()} combos found)"

        def snvBase = "${snvdir}/${ver[0]}_to_${ver[1]}"
        def snvLooseVcf = "${snvBase}/snv_candidates.merfin-loose.vcf.gz"
        def snvLooseTbi = "${snvLooseVcf}.tbi"
        def snvCandidatesVcf = "${snvBase}/intermediates/snv_candidates.vcf.gz"

        def hap1Fa = "${asmout}/${params.asm_name}_${ver[1]}.hap1.fa.gz"
        def hap2Fa = "${asmout}/${params.asm_name}_${ver[1]}.hap2.fa.gz"
        def dipFa  = "${asmout}/${params.asm_name}_${ver[1]}.dip.fa.gz"

        def hap1Ok = new File(hap1Fa).exists() && new File("${hap1Fa}.gzi").exists() && new File("${hap1Fa}.fai").exists()
        def hap2Ok = new File(hap2Fa).exists() && new File("${hap2Fa}.gzi").exists() && new File("${hap2Fa}.fai").exists()
        def dipOk  = new File(dipFa).exists()  && new File("${dipFa}.gzi").exists()  && new File("${dipFa}.fai").exists()

        r1_snv_next_refs_present = hap1Ok && hap2Ok && dipOk
        r1_snv_present = new File(snvLooseVcf).exists() && new File(snvLooseTbi).exists() && r1_snv_next_refs_present

        log.info "SNV_CANDIDATES inventory: loose_vcf=${new File(snvLooseVcf).exists()} tbi=${new File(snvLooseTbi).exists()} candidates_vcf=${new File(snvCandidatesVcf).exists()}"
        log.info "SNV_CANDIDATES inventory: next_refs hap1=${hap1Ok} hap2=${hap2Ok} dip=${dipOk}"
        log.info "SNV_CANDIDATES inventory: Round-1 complete=${r1_snv_present}"
    }

    if ( params.run_dv ) {
        def r1_dv_vcfs
        if ( params.deepvariant_dir && r1_all_dv_present ) {
            log.info "Round-1 DeepVariant: all expected outputs found in deepvariant_dir; skipping DEEPVARIANT_R1"
            r1_dv_vcfs = Channel.fromList(r1_dv_records)
        } else {
            if (params.deepvariant_dir) {
                log.info "Round-1 DeepVariant: missing outputs detected (${r1_dv_missing_keys.size()} combos); running DEEPVARIANT_R1 for missing items only"
            }
            DEEPVARIANT_R1( r1_dv_refs, r1_wm_bams, r1_bwa_bams )
            r1_dv_vcfs = DEEPVARIANT_R1.out.dv_vcfs
        }

        def r1_next_refs

        if ( params.run_snv_candidates ) {
            if ( params.deepvariant_dir && r1_snv_present ) {
                log.info "Round-1 SNV_CANDIDATES: all expected outputs found; skipping SNV_CANDIDATES_R1"
                r1_next_refs = Channel.of(
                    tuple('hap1', file("${params.outdir}/assemblies/${params.asm_name}_${ver[1]}.hap1.fa.gz"), file("${params.outdir}/assemblies/${params.asm_name}_${ver[1]}.hap1.fa.gz.gzi"), file("${params.outdir}/assemblies/${params.asm_name}_${ver[1]}.hap1.fa.gz.fai")),
                    tuple('hap2', file("${params.outdir}/assemblies/${params.asm_name}_${ver[1]}.hap2.fa.gz"), file("${params.outdir}/assemblies/${params.asm_name}_${ver[1]}.hap2.fa.gz.gzi"), file("${params.outdir}/assemblies/${params.asm_name}_${ver[1]}.hap2.fa.gz.fai")),
                    tuple('dip',  file("${params.outdir}/assemblies/${params.asm_name}_${ver[1]}.dip.fa.gz"),  file("${params.outdir}/assemblies/${params.asm_name}_${ver[1]}.dip.fa.gz.gzi"),  file("${params.outdir}/assemblies/${params.asm_name}_${ver[1]}.dip.fa.gz.fai"))
                )
            } else {
                requireParam('hybrid_meryl', params.hybrid_meryl)
                requireParam('merfin_peak',  params.merfin_peak)
                SNV_CANDIDATES_R1( r1_dv_vcfs, r1_refs, ver[0], ver[1] )
                r1_next_refs = SNV_CANDIDATES_R1.out.next_refs
            }
        }

        // ---- Round 2 ------------------------------------------------------------
        if ( rounds >= 2 && params.run_snv_candidates ) {
            BUILD_REFS_R2( r1_next_refs, ver[1] )
            def r2_refs = BUILD_REFS_R2.out.wm_refs.map { hap, ver_from, ref, fai, _rep -> tuple(hap, ver_from, ref, fai) }

            MAPPING_R2( BUILD_REFS_R2.out.wm_refs, BUILD_REFS_R2.out.bwa_refs )
            DEEPVARIANT_R2( BUILD_REFS_R2.out.dv_refs, MAPPING_R2.out.wm_pri_bams, MAPPING_R2.out.bwa_mrg_bams )
            SNV_CANDIDATES_R2( DEEPVARIANT_R2.out.dv_vcfs, r2_refs, ver[1], ver[2] )
        }

        // ---- Round 3 ------------------------------------------------------------
        if ( rounds >= 3 && params.run_snv_candidates ) {
            BUILD_REFS_R3( SNV_CANDIDATES_R2.out.next_refs, ver[2] )
            def r3_refs = BUILD_REFS_R3.out.wm_refs.map { hap, ver_from, ref, fai, _rep -> tuple(hap, ver_from, ref, fai) }

            MAPPING_R3( BUILD_REFS_R3.out.wm_refs, BUILD_REFS_R3.out.bwa_refs )
            DEEPVARIANT_R3( BUILD_REFS_R3.out.dv_refs, MAPPING_R3.out.wm_pri_bams, MAPPING_R3.out.bwa_mrg_bams )
            SNV_CANDIDATES_R3( DEEPVARIANT_R3.out.dv_vcfs, r3_refs, ver[2], ver[3] )
        }

        // ---- Round 4 ------------------------------------------------------------
        if ( rounds >= 4 && params.run_snv_candidates ) {
            BUILD_REFS_R4( SNV_CANDIDATES_R3.out.next_refs, ver[3] )
            def r4_refs = BUILD_REFS_R4.out.wm_refs.map { hap, ver_from, ref, fai, _rep -> tuple(hap, ver_from, ref, fai) }

            MAPPING_R4( BUILD_REFS_R4.out.wm_refs, BUILD_REFS_R4.out.bwa_refs )
            DEEPVARIANT_R4( BUILD_REFS_R4.out.dv_refs, MAPPING_R4.out.wm_pri_bams, MAPPING_R4.out.bwa_mrg_bams )
            SNV_CANDIDATES_R4( DEEPVARIANT_R4.out.dv_vcfs, r4_refs, ver[3], ver[4] )
        }

        // ---- Round 5 ------------------------------------------------------------
        if ( rounds >= 5 && params.run_snv_candidates ) {
            BUILD_REFS_R5( SNV_CANDIDATES_R4.out.next_refs, ver[4] )
            def r5_refs = BUILD_REFS_R5.out.wm_refs.map { hap, ver_from, ref, fai, _rep -> tuple(hap, ver_from, ref, fai) }

            MAPPING_R5( BUILD_REFS_R5.out.wm_refs, BUILD_REFS_R5.out.bwa_refs )
            DEEPVARIANT_R5( BUILD_REFS_R5.out.dv_refs, MAPPING_R5.out.wm_pri_bams, MAPPING_R5.out.bwa_mrg_bams )
            SNV_CANDIDATES_R5( DEEPVARIANT_R5.out.dv_vcfs, r5_refs, ver[4], ver[5] )
        }
    } else {
        log.info "Skipping DeepVariant and SNV candidate stages (--run_dv false)."
    }
}
