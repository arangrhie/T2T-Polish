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
params.evaluate           = params.evaluate           ?: false  // evaluation-only mode (skip DeepVariant/SNV)

// Backward compatibility for legacy run_dv flag:
// - --evaluate true always enables evaluation-only mode.
// - Otherwise, honor explicit run_dv if provided.
// - Default behavior remains polishing mode (DeepVariant enabled).
def run_dv_enabled = params.evaluate ? false
                    : (params.containsKey('run_dv') ? (params.run_dv as boolean) : true)

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

include { MAPPING_R1  } from './workflows/mapping_r1'
include { MAPPING_R2  } from './workflows/mapping_r2'

include { DEEPVARIANT_R1 } from './workflows/deepvariant_r1'
include { DEEPVARIANT_R2 } from './workflows/deepvariant_r2'

include { SNV_CANDIDATES as SNV_CANDIDATES_R1 } from './workflows/snv_candidates'
include { SNV_CANDIDATES as SNV_CANDIDATES_R2 } from './workflows/snv_candidates'

include { EVALUATION } from './workflows/evaluation'

// Standalone SAM2PAF alias used to backfill PAFs in mapping_dir re-entry when
// a hap's BAM is on disk but its sibling .pri.paf is missing. Aliased so it
// does not collide with the SAM2PAF_R1 used inside MAPPING_R1 (NF 25.x
// requires unique process names per pipeline run).
include { SAM2PAF as SAM2PAF_BACKFILL } from './modules/winnowmap'

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

// ------------ disk-scan re-entry helpers -----------------------------------
// Inspect the published output tree for a given ${asm_name}_${ver} prefix and
// report which BAMs / VCFs / polished FASTAs already exist.  Used to gate
// MAPPING_R{n}, DEEPVARIANT_R{n} and SNV_CANDIDATES_R{n} so a partial run is
// picked up automatically even when Nextflow's -resume cache is unavailable
// (e.g. the work/ directory was cleared or the run was launched fresh).
//
// These helpers are pure Groovy and return plain Maps; the workflow body is
// responsible for building Channels from the reported files.

def _pickBamIdx(bam) {
    def b = file("${bam}.bai")
    b.exists() ? b : file("${bam}.csi")
}

def scanMappingDir(String mdir, String pfx, platforms, boolean ontMapHaps) {
    // Compute the expected (hap, platform) BAM combos from params.platforms:
    //   hifi         → hap1, hap2, dip
    //   ont          → dip     (also hap1, hap2 iff ont_map_haps)
    //   illumina|element (short) → hap1, hap2, dip
    //   hybrid (synthetic; needed by DEEPVARIANT_R{n} hybrid track)
    //                → present iff both hifi bam AND short-read bam exist for that hap
    // Only combos matching platforms in `platforms` are ever "expected".
    // Accept either a List<String> or a comma-separated String — Nextflow may
    // hand us the raw config value, and `String as Set` in Groovy silently
    // char-splits, so normalise here.
    def platList = (platforms instanceof List)
        ? platforms
        : (platforms == null ? [] : platforms.toString().tokenize(',').collect { it.trim() }.findAll { it })
    def platSet    = platList as Set
    def haps       = ['hap1', 'hap2', 'dip']
    def shortPlats = ['illumina', 'element'].findAll { platSet.contains(it) }
    // The DV hybrid track combines HiFi with a single short-read platform;
    // in practice illumina wins if both are provided (matches _short_for_hybrid).
    def usedShort  = shortPlats ? (shortPlats.contains('illumina') ? 'illumina' : shortPlats[0]) : null

    def expectedHifi  = platSet.contains('hifi')
        ? haps.collect { [it, 'hifi'] } : []
    def expectedOnt   = platSet.contains('ont')
        ? (ontMapHaps ? haps.collect { [it, 'ont'] } : [['dip', 'ont']]) : []
    def expectedLong  = (expectedHifi + expectedOnt)
    def expectedShort = shortPlats.collectMany { plat -> haps.collect { [it, plat] } }
    def expectedHybridHaps = (platSet.contains('hifi') && usedShort) ? haps : []

    def bamPath   = { hap, plat, suffix ->
        "${mdir}/${pfx}.${hap}.${plat}/${pfx}.${hap}.${plat}${suffix}"
    }
    def existsBam = { hap, plat, suffix -> new File(bamPath(hap, plat, suffix)).exists() }

    def emptyResult = [
        expected_long   : expectedLong  as Set,
        expected_short  : expectedShort as Set,
        expected_hybrid_haps : expectedHybridHaps as Set,
        present_long    : [] as Set,
        present_short   : [] as Set,
        missing_long    : expectedLong  as Set,
        missing_short   : expectedShort as Set,
        hybrid_ready_haps: [] as Set,
        missing_hybrid_haps: expectedHybridHaps as Set,
        wm_haps_done    : [],
        bwa_haps_done   : [],
        used_short      : usedShort,
        paf_on_disk     : [] as Set
    ]
    if ( !mdir || !(new File(mdir).isDirectory()) ) return emptyResult

    // Per-combo presence check (constructed paths — no globs, no regex).
    def presentLong  = (expectedLong .findAll { hap, plat -> existsBam(hap, plat, '.pri.bam')       }) as Set
    def presentShort = (expectedShort.findAll { hap, plat -> existsBam(hap, plat, '.dedup.pri.bam') }) as Set
    def missingLong  = ((expectedLong  as Set) - presentLong)  as Set
    def missingShort = ((expectedShort as Set) - presentShort) as Set

    // Hybrid readiness (per hap): hifi bam + short-read bam both present.
    def hybridReady = expectedHybridHaps.findAll { hap ->
        presentLong .contains([hap, 'hifi']) &&
        presentShort.contains([hap, usedShort])
    } as Set
    def missingHybrid = ((expectedHybridHaps as Set) - hybridReady) as Set

    // Per-hap mapping-stage completeness (all its expected combos present).
    def wmDone = haps.findAll { hap ->
        def exp = expectedLong.findAll { it[0] == hap }
        exp && exp.every { presentLong.contains(it) }
    }
    def bwaDone = haps.findAll { hap ->
        def exp = expectedShort.findAll { it[0] == hap }
        exp && exp.every { presentShort.contains(it) }
    }

    // Sibling PAFs (used by R1 evaluation; only for winnowmap combos).
    def pafOn = (expectedLong.findAll { hap, plat ->
        new File("${mdir}/${pfx}.${hap}.${plat}/${pfx}.${hap}.${plat}.pri.paf").exists()
    }) as Set

    return [
        expected_long   : expectedLong  as Set,
        expected_short  : expectedShort as Set,
        expected_hybrid_haps : expectedHybridHaps as Set,
        present_long    : presentLong,
        present_short   : presentShort,
        missing_long    : missingLong,
        missing_short   : missingShort,
        hybrid_ready_haps: hybridReady,
        missing_hybrid_haps: missingHybrid,
        wm_haps_done    : wmDone,
        bwa_haps_done   : bwaDone,
        used_short      : usedShort,
        paf_on_disk     : pafOn
    ]
}

def scanDeepVariantDir(String dvdir, String pfx, String ver,
                        int mqHap, int mqDip, String hybridCombo, String ontMode) {
    def records = []
    def missing = []
    def hybridNeeded = [] as Set   // haps with a missing hybrid combo
    def ontDipNeeded = false        // true iff dip ONT combo is missing
    if ( !dvdir || !(new File(dvdir).isDirectory()) ) {
        return [ records: records, missing_keys: ['<no dvdir>'], all_present: false,
                 hybrid_needed_haps: ['hap1','hap2','dip'] as Set, ont_dip_needed: true ]
    }
    def expected = [
        [hap: 'hap1', combo: hybridCombo, mode: 'HYBRID_PACBIO_ILLUMINA', mq: mqHap],
        [hap: 'hap2', combo: hybridCombo, mode: 'HYBRID_PACBIO_ILLUMINA', mq: mqHap],
        [hap: 'dip',  combo: hybridCombo, mode: 'HYBRID_PACBIO_ILLUMINA', mq: mqDip],
        [hap: 'dip',  combo: 'ont',       mode: ontMode,                  mq: mqDip]
    ]
    expected.each { e ->
        def runDir = "${dvdir}/${pfx}.${e.hap}.${e.combo}.MQ${e.mq}"
        def vcf   = "${runDir}/dv_${e.mode}_MQ${e.mq}.${e.hap}.vcf.gz"
        def tbi   = "${vcf}.tbi"
        def gvcf  = "${runDir}/dv_${e.mode}_MQ${e.mq}.${e.hap}.gvcf.gz"
        def gtbi  = "${gvcf}.tbi"
        def ok = [vcf,tbi,gvcf,gtbi].every { new File(it).exists() }
        def key = "${e.hap}|${e.combo}|${e.mode}|MQ${e.mq}"
        if (ok) {
            records << tuple(e.hap, ver, e.combo, e.mode, e.mq,
                             file(vcf), file(tbi), file(gvcf), file(gtbi))
        } else {
            missing << key
            if (e.combo == 'ont') ontDipNeeded = true
            else                  hybridNeeded << e.hap
        }
    }
    return [ records: records, missing_keys: missing, all_present: missing.isEmpty(),
             hybrid_needed_haps: hybridNeeded, ont_dip_needed: ontDipNeeded ]
}

/*
 * For each DV combo whose final VCF is missing, check whether upstream
 * intermediate outputs are already on disk so DEEPVARIANT_R{n} can start
 * from the latest available stage:
 *   - MERGE_HYBRID:      ${mdir}/${pfx}.{hap}.${hybridCombo}/${pfx}.{hap}.${hybridCombo}.bam(+.bai)
 *   - DV_MAKE_EXAMPLES:  ${dvdir}/${pfx}.{hap}.{combo}.MQ{mq}/examples/  (directory)
 * Only combos in neededHybridHaps / ontDipNeeded are scanned — combos whose
 * VCF is already present are handled by scanDeepVariantDir and never reach
 * DEEPVARIANT_R{n}.
 */
def scanDvIntermediates(String dvdir, String mdir, String pfx, String ver,
                        int mqHap, int mqDip, String hybridCombo, String ontMode,
                        neededHybridHaps, boolean ontDipNeeded) {
    def hybHaps = neededHybridHaps ? (neededHybridHaps as List) : []
    def mergedRecords       = []
    def hybExamplesRecords  = []
    def ontExamplesRecords  = []

    hybHaps.each { hap ->
        def mq = (hap == 'dip') ? mqDip : mqHap
        def mergedBam = "${mdir}/${pfx}.${hap}.${hybridCombo}/${pfx}.${hap}.${hybridCombo}.bam"
        def mergedBai = "${mergedBam}.bai"
        if (new File(mergedBam).exists() && new File(mergedBai).exists()) {
            mergedRecords << tuple(hap, ver, hybridCombo, file(mergedBam), file(mergedBai))
        }
        def hybEx = "${dvdir}/${pfx}.${hap}.${hybridCombo}.MQ${mq}/examples"
        if (new File(hybEx).isDirectory()) {
            hybExamplesRecords << tuple(hap, ver, hybridCombo, 'HYBRID_PACBIO_ILLUMINA', mq, file(hybEx))
        }
    }

    if (ontDipNeeded) {
        def hap = 'dip'
        def mq  = mqDip
        def ontEx = "${dvdir}/${pfx}.${hap}.ont.MQ${mq}/examples"
        if (new File(ontEx).isDirectory()) {
            ontExamplesRecords << tuple(hap, ver, 'ont', ontMode, mq, file(ontEx))
        }
    }

    return [
        merged_records       : mergedRecords,
        hyb_examples_records : hybExamplesRecords,
        ont_examples_records : ontExamplesRecords
    ]
}

def scanSnvOutputs(String snvdir, String asmout, String asmName,
                    String verFrom, String verTo) {
    def snvBase  = "${snvdir}/${verFrom}_to_${verTo}"
    def loose    = "${snvBase}/snv_candidates.merfin-loose.vcf.gz"
    def looseTbi = "${loose}.tbi"

    def hap1Fa = "${asmout}/${asmName}_${verTo}.hap1.fa.gz"
    def hap2Fa = "${asmout}/${asmName}_${verTo}.hap2.fa.gz"
    def dipFa  = "${asmout}/${asmName}_${verTo}.dip.fa.gz"

    def hap1Ok = [hap1Fa, "${hap1Fa}.gzi", "${hap1Fa}.fai"].every { new File(it).exists() }
    def hap2Ok = [hap2Fa, "${hap2Fa}.gzi", "${hap2Fa}.fai"].every { new File(it).exists() }
    def dipOk  = [dipFa,  "${dipFa}.gzi",  "${dipFa}.fai" ].every { new File(it).exists() }

    def nextRefsPresent = hap1Ok && hap2Ok && dipOk
    def snvPresent = new File(loose).exists() && new File(looseTbi).exists() && nextRefsPresent
    def tuples = nextRefsPresent ? [
        tuple('hap1', file(hap1Fa), file("${hap1Fa}.gzi"), file("${hap1Fa}.fai")),
        tuple('hap2', file(hap2Fa), file("${hap2Fa}.gzi"), file("${hap2Fa}.fai")),
        tuple('dip',  file(dipFa),  file("${dipFa}.gzi"),  file("${dipFa}.fai"))
    ] : []

    return [ snv_present: snvPresent, next_refs_present: nextRefsPresent,
             next_refs_tuples: tuples ]
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
    if ( rounds > 2 ) error "params.polish_rounds > 2 is not supported (got ${params.polish_rounds})"
    if ( params.containsKey('run_dv') ) {
        log.warn "params.run_dv is deprecated; use --evaluate true for evaluation-only mode."
    }
    if ( !run_dv_enabled ) {
        log.warn "Evaluation-only mode is enabled (--evaluate true); overriding params.polish_rounds to 1."
        rounds = 1
    }

    def ver = []
    ver[0] = params.asm_ver ?: 'v0.1'
    ver[1] = params.asm_ver_next ?: bumpVersion(ver[0])
    ver[2] = bumpVersion(ver[1])

    // =========================================================================
    // Each round is fully sequential: BUILD_REFS → MAPPING → DEEPVARIANT →
    // SNV_CANDIDATES.  Per-round workflow files carry aliased process names
    // (_R1 / _R2) so NF 25.x sees each process bound exactly once.
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
    def r1_wm_pafs
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

        // PAFs are published alongside BAMs by SAM2PAF; only mix in PAFs for haps
        // whose mapping is fully complete on disk so EVALUATION doesn't fire on
        // missing files.
        def r1_wm_pafs_from_disk = Channel.fromPath("${mdir}/${_pfx}.*.*/${_pfx}.*.*.pri.paf")
            .map { f ->
                def m = (f.name =~ /\.([^.]+)\.([^.]+)\.pri\.paf$/)[0]
                tuple(m[1], params.asm_ver, m[2], f)
            }
            .filter { hap, v, plat, paf -> wm_haps_done.contains(hap) }

        // ---- Backfill: BAMs on disk whose sibling .pri.paf is missing ----
        // Without this, EVALUATION silently no-ops (the channel filter for the
        // missing dip hifi/ont PAF emits nothing, so combine() produces zero
        // items and FULL_EVALUATION is never instantiated).
        def paf_on_disk = _pfx_subdirs
            .collectMany { dir -> dir.listFiles()?.toList() ?: [] }
            .findAll { it.name ==~ /.*\.([^.]+)\.([^.]+)\.pri\.paf$/ }
            .collect { f ->
                def m = (f.name =~ /\.([^.]+)\.([^.]+)\.pri\.paf$/)[0]
                [m[1], m[2]]
            } as Set

        def wm_needs_paf = wm_on_disk.findAll { !paf_on_disk.contains(it) && wm_haps_done.contains(it[0]) }
        log.info "mapping_dir scan: paf_on_disk = ${paf_on_disk}"
        if (wm_needs_paf) {
            log.info "mapping_dir backfill: regenerating PAF via SAM2PAF for BAMs with no sibling .pri.paf: ${wm_needs_paf}"
        }

        def r1_bams_need_paf = r1_wm_from_disk
            .filter { hap, v, plat, bam, idx -> wm_needs_paf.contains([hap, plat]) }
        def r1_wm_pafs_backfilled = SAM2PAF_BACKFILL(r1_bams_need_paf)

        r1_wm_bams  = r1_wm_from_disk.mix( MAPPING_R1.out.wm_pri_bams )
        r1_wm_pafs  = r1_wm_pafs_from_disk
            .mix( MAPPING_R1.out.wm_pri_pafs )
            .mix( r1_wm_pafs_backfilled )
        r1_bwa_bams = r1_bwa_from_disk.mix( MAPPING_R1.out.bwa_mrg_bams )
    } else if ( !params.deepvariant_dir ) {
        MAPPING_R1( BUILD_REFS.out.wm_refs, BUILD_REFS.out.bwa_refs )
        r1_wm_bams  = MAPPING_R1.out.wm_pri_bams
        r1_wm_pafs  = MAPPING_R1.out.wm_pri_pafs
        r1_bwa_bams = MAPPING_R1.out.bwa_mrg_bams
    } else {
        log.info "deepvariant_dir is set; skipping MAPPING_R1 and mapping_dir re-entry for Round-1."
        r1_wm_bams  = Channel.empty()
        r1_wm_pafs  = Channel.empty()
        r1_bwa_bams = Channel.empty()
    }

    if ( run_dv_enabled ) {
        // ---- shared context for disk-scan re-entry (used by all rounds) ----
        def _short_for_hybrid = _plat_list.contains('illumina') ? 'illumina' :
                                (_plat_list.contains('element') ? 'element' :
                                 (_active_short ? _active_short[0] : 'illumina'))
        def hybrid_combo = "hifi_${_short_for_hybrid}"
        def ont_mode     = (params.ont_chemistry == 'r9') ? 'ONT_R9' : 'ONT_R104'

        def _mdir_all     = params.mapping_outdir.replaceAll('/$', '')
        def _dvdir_all    = params.dv_outdir.replaceAll('/$', '')
        def _snvdir_all   = params.snv_outdir.replaceAll('/$', '')
        def _asmout_all   = "${params.outdir}/assemblies".replaceAll('/$', '')
        def _ont_map_haps = params.ont_map_haps as boolean
        def _mq_hap       = params.dv_mq_hap.toInteger()
        def _mq_dip       = params.dv_mq_dip.toInteger()

        // ---- per-round context builder --------------------------------------
        // Runs all three disk scans for a given round and returns a Map of
        // channels + predicates the round body uses to gate MAPPING_R{n} /
        // DEEPVARIANT_R{n} / SNV_CANDIDATES_R{n}.  Process invocations stay
        // inline in the per-round if-blocks because Nextflow's AST transformer
        // only rewrites literal identifiers as process/workflow calls.
        //
        // dvdirOverride / mdirOverride — used by Round 1 to point the DV and
        //                  mapping scans at user-supplied re-entry directories
        //                  (params.deepvariant_dir / params.mapping_dir).
        //                  Round 2 leaves them null and scans the default
        //                  publish dirs (params.dv_outdir / params.mapping_outdir).
        def roundContext = { int n, String ver_from, String ver_to,
                              String dvdirOverride = null, String mdirOverride = null ->
            def _mdir  = mdirOverride  ?: _mdir_all
            def _dvdir = dvdirOverride ?: _dvdir_all
            def pfx    = "${params.asm_name}_${ver_from}"
            // Use _plat_list (guaranteed List<String>) rather than
            // params.platforms directly — the top-of-script reassignment
            // of params.platforms may not persist into this closure scope,
            // and `String as Set` silently char-splits.
            def mscn = scanMappingDir(_mdir, pfx, _plat_list, _ont_map_haps)
            log.info "Round-${n} mapping scan (${_mdir}, prefix '${pfx}'):"
            log.info "  expected long  = ${mscn.expected_long  as List}"
            log.info "  expected short = ${mscn.expected_short as List}"
            log.info "  present  long  = ${mscn.present_long   as List}"
            log.info "  present  short = ${mscn.present_short  as List}"
            log.info "  missing  long  = ${mscn.missing_long   as List}"
            log.info "  missing  short = ${mscn.missing_short  as List}"
            log.info "  hybrid ready (hifi+${mscn.used_short}) haps = ${mscn.hybrid_ready_haps as List}"
            log.info "  wm_done=${mscn.wm_haps_done} bwa_done=${mscn.bwa_haps_done}"

            // Build disk BAM channels directly from the presence set — no globs.
            // Only include BAMs for haps that are FULLY done for their track;
            // partially-done haps get fully re-mapped by MAPPING_R{n} (mapping
            // does not support skipping individual platforms per hap).
            def wm_disk_records = mscn.present_long
                .findAll { hap, plat -> mscn.wm_haps_done.contains(hap) }
                .collect { hap, plat ->
                    def bam = file("${_mdir}/${pfx}.${hap}.${plat}/${pfx}.${hap}.${plat}.pri.bam")
                    tuple(hap, ver_from, plat, bam, _pickBamIdx(bam))
                }
            def bwa_disk_records = mscn.present_short
                .findAll { hap, plat -> mscn.bwa_haps_done.contains(hap) }
                .collect { hap, plat ->
                    def bam = file("${_mdir}/${pfx}.${hap}.${plat}/${pfx}.${hap}.${plat}.dedup.pri.bam")
                    tuple(hap, ver_from, plat, bam, _pickBamIdx(bam))
                }
            def wm_disk  = wm_disk_records  ? Channel.fromList(wm_disk_records)  : Channel.empty()
            def bwa_disk = bwa_disk_records ? Channel.fromList(bwa_disk_records) : Channel.empty()

            def dvscn = scanDeepVariantDir(_dvdir, pfx, ver_from,
                _mq_hap, _mq_dip, hybrid_combo, ont_mode)
            log.info "Round-${n} deepvariant scan (${_dvdir}): all_present=${dvscn.all_present} missing=${dvscn.missing_keys}"
            log.info "Round-${n} deepvariant: hybrid_needed_haps=${dvscn.hybrid_needed_haps as List} ont_dip_needed=${dvscn.ont_dip_needed}"
            if (dvscn.all_present)
                log.info "Round-${n} DeepVariant: all expected outputs found; skipping DEEPVARIANT_R${n}"

            // Intermediate DV output scan — for combos still needing work,
            // check for merged hybrid BAMs (MERGE_HYBRID output) and DV
            // examples/ dirs (DV_MAKE_EXAMPLES output) so those stages can be
            // skipped when their outputs are already on disk.
            def dv_int = scanDvIntermediates(_dvdir, _mdir, pfx, ver_from,
                _mq_hap, _mq_dip, hybrid_combo, ont_mode,
                dvscn.hybrid_needed_haps, dvscn.ont_dip_needed)
            log.info "Round-${n} DV intermediates on disk:"
            log.info "  merged hybrid BAMs = ${dv_int.merged_records.collect { it[0] }}"
            log.info "  hyb examples dirs  = ${dv_int.hyb_examples_records.collect { it[0] }}"
            log.info "  ont examples dirs  = ${dv_int.ont_examples_records.collect { it[0] }}"

            def pre_merged_hyb  = dv_int.merged_records       ? Channel.fromList(dv_int.merged_records)       : Channel.empty()
            def pre_examples_hyb = dv_int.hyb_examples_records ? Channel.fromList(dv_int.hyb_examples_records) : Channel.empty()
            def pre_examples_ont = dv_int.ont_examples_records ? Channel.fromList(dv_int.ont_examples_records) : Channel.empty()

            def snvscn = scanSnvOutputs(_snvdir_all, _asmout_all, params.asm_name, ver_from, ver_to)
            log.info "Round-${n} snv scan: snv_present=${snvscn.snv_present} next_refs_present=${snvscn.next_refs_present}"
            if (snvscn.snv_present)
                log.info "Round-${n} SNV_CANDIDATES: all expected outputs found; skipping SNV_CANDIDATES_R${n}"

            // Haps whose ref/BAMs are still needed by DEEPVARIANT_R{n}:
            //   - any hap with a missing hybrid combo (needs hifi + short-read BAMs)
            //   - dip if ONT combo is missing (needs ont BAM)
            def dv_needed_haps = (dvscn.hybrid_needed_haps as List) +
                                 (dvscn.ont_dip_needed ? ['dip'] : [])
            dv_needed_haps = dv_needed_haps.unique()

            return [
                wm_haps_done         : mscn.wm_haps_done,
                bwa_haps_done        : mscn.bwa_haps_done,
                wm_disk              : wm_disk,
                bwa_disk             : bwa_disk,
                dv_all_present       : dvscn.all_present,
                dv_records           : dvscn.records,
                dv_hybrid_needed_haps: dvscn.hybrid_needed_haps,
                dv_ont_dip_needed    : dvscn.ont_dip_needed,
                dv_needed_haps       : dv_needed_haps,
                pre_merged_hyb       : pre_merged_hyb,
                pre_examples_hyb     : pre_examples_hyb,
                pre_examples_ont     : pre_examples_ont,
                snv_present          : snvscn.snv_present,
                snv_next_refs_tuples : snvscn.next_refs_tuples
            ]
        }

        // ---- Round 1 --------------------------------------------------------
        // R1 uses params.deepvariant_dir (when set) as the DV scan location;
        // otherwise falls back to the default dv_outdir.  R1's mapping BAM
        // channels (r1_wm_bams / r1_bwa_bams) come from the mapping_dir
        // re-entry block above — roundContext's wm_disk / bwa_disk / wm_haps_done
        // are not consumed here.
        def r1_ctx = roundContext(1, ver[0], ver[1], params.deepvariant_dir, params.mapping_dir)

        def r1_dv_vcfs
        if ( r1_ctx.dv_all_present ) {
            r1_dv_vcfs = Channel.fromList(r1_ctx.dv_records)
        } else {
            def r1_dv_refs_needed = r1_dv_refs
                .filter { hap, v, ref, gzi, fai -> r1_ctx.dv_needed_haps.contains(hap) }
            def r1_wm_needed = r1_wm_bams.filter { hap, v, plat, bam, idx ->
                (plat != 'ont' && r1_ctx.dv_hybrid_needed_haps.contains(hap)) ||
                (plat == 'ont' && hap == 'dip' && r1_ctx.dv_ont_dip_needed)
            }
            def r1_bwa_needed = r1_bwa_bams.filter { hap, v, plat, bam, idx ->
                r1_ctx.dv_hybrid_needed_haps.contains(hap)
            }
            DEEPVARIANT_R1( r1_dv_refs_needed, r1_wm_needed, r1_bwa_needed,
                            r1_ctx.pre_merged_hyb, r1_ctx.pre_examples_hyb, r1_ctx.pre_examples_ont )
            r1_dv_vcfs = DEEPVARIANT_R1.out.dv_vcfs.mix( Channel.fromList(r1_ctx.dv_records) )
        }

        def r1_next_refs
        if ( params.run_snv_candidates ) {
            if ( r1_ctx.snv_present ) {
                r1_next_refs = Channel.fromList(r1_ctx.snv_next_refs_tuples)
            } else {
                requireParam('hybrid_meryl', params.hybrid_meryl)
                requireParam('merfin_peak',  params.merfin_peak)
                SNV_CANDIDATES_R1( r1_dv_vcfs, r1_refs, ver[0], ver[1] )
                r1_next_refs = SNV_CANDIDATES_R1.out.next_refs
            }
        }

        // ---- Round 2 --------------------------------------------------------
        if ( rounds >= 2 && params.run_snv_candidates ) {
            BUILD_REFS_R2( r1_next_refs, ver[1] )
            def r2_refs = BUILD_REFS_R2.out.wm_refs.map { hap, v, ref, fai, _rep -> tuple(hap, v, ref, fai) }
            def ctx = roundContext(2, ver[1], ver[2])

            MAPPING_R2(
                BUILD_REFS_R2.out.wm_refs .filter { hap, v, ref, fai, rep                -> !ctx.wm_haps_done .contains(hap) },
                BUILD_REFS_R2.out.bwa_refs.filter { hap, v, ref, fai, amb, ann, bwt, pac, sa -> !ctx.bwa_haps_done.contains(hap) }
            )
            def r2_wm_bams  = ctx.wm_disk .mix( MAPPING_R2.out.wm_pri_bams  )
            def r2_bwa_bams = ctx.bwa_disk.mix( MAPPING_R2.out.bwa_mrg_bams )

            def r2_dv_vcfs
            if ( ctx.dv_all_present ) {
                r2_dv_vcfs = Channel.fromList(ctx.dv_records)
            } else {
                // Feed only refs / BAMs needed by the missing DV combos:
                //   hybrid combo (hap): hifi BAM + short-read BAM for that hap
                //   ont    combo (dip): ont BAM for dip
                // Already-done DV records are mixed back in below so
                // SNV_CANDIDATES still receives all four combos.
                def r2_dv_refs_needed = BUILD_REFS_R2.out.dv_refs
                    .filter { hap, v, ref, gzi, fai -> ctx.dv_needed_haps.contains(hap) }
                def r2_wm_needed = r2_wm_bams.filter { hap, v, plat, bam, idx ->
                    (plat != 'ont' && ctx.dv_hybrid_needed_haps.contains(hap)) ||
                    (plat == 'ont' && hap == 'dip' && ctx.dv_ont_dip_needed)
                }
                def r2_bwa_needed = r2_bwa_bams.filter { hap, v, plat, bam, idx ->
                    ctx.dv_hybrid_needed_haps.contains(hap)
                }
                DEEPVARIANT_R2( r2_dv_refs_needed, r2_wm_needed, r2_bwa_needed,
                                ctx.pre_merged_hyb, ctx.pre_examples_hyb, ctx.pre_examples_ont )
                r2_dv_vcfs = DEEPVARIANT_R2.out.dv_vcfs.mix( Channel.fromList(ctx.dv_records) )
            }

            if ( !ctx.snv_present ) {
                SNV_CANDIDATES_R2( r2_dv_vcfs, r2_refs, ver[1], ver[2] )
            }
        }
    } else {
        log.info "Skipping DeepVariant and SNV candidate stages (--evaluate true)."
        
        // ---- Evaluation Stage (when evaluate = true) ----
        // Run post-assembly evaluation: Merqury QV, pattern analysis, issue detection
        
        if ( !params.hybrid_meryl ) {
            log.warn "params.hybrid_meryl not set; skipping EVALUATION workflow."
        } else {
            def asm_dip_input = Channel.of(
                tuple("${params.asm_name}_${params.asm_ver}", 
                      file("${params.outdir}/assemblies/${params.asm_name}_${params.asm_ver}.dip.fa.gz"),
                      file("${params.outdir}/assemblies/${params.asm_name}_${params.asm_ver}.dip.fa.gz.fai"),
                      params.asm_ver)
            )
            
            EVALUATION( asm_dip_input, r1_wm_pafs )
        }
    }
}
