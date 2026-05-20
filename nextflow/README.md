# T2T-Polish

Nextflow DSL2 pipeline for polishing a T2T diploid assembly.


## Overview

The pipeline takes a [Verkko](https://github.com/marbl/verkko) assembly and maps HiFi, ONT, Illumina or Element reads against three reference builds — `hap1`, `hap2`, and `dip` — all in parallel. Slurm submission is handled natively by Nextflow; no external submit scripts are needed.


```
nextflow/                    # Entry point: params + include + workflow {}
├── nextflow.config          # Global defaults; loads resources.config
├── resources.config         # Executor + per-label CPU/mem/time (Biowulf Slurm defaults)
├── user.config.example      # Copy → user.config and fill in paths/globs
├── modules/
│   ├── references.nf        # BUILD_HAP_REFERENCES, BUILD_DIP_REFERENCE
│   ├── winnowmap.nf         # MERYL_REPETITIVE, WINNOWMAP_MAP, WINNOWMAP_MERGE,
│   │                        #   WINNOWMAP_FILTER, SAM2PAF
│   ├── bwa.nf               # BWA_INDEX, BWA_MAP, BWA_MERGE
│   ├── deepvariant.nf       # MERGE_HYBRID, DV_MAKE_EXAMPLES,
│   │                        #   DV_CALL_VARIANTS, DV_POSTPROCESS,
│   │                        #   PEPPER_MARGIN_DV, DV_MERGE_CHR_VCFS
│   └── snv_candidates.nf    # SNV_FILTER_INTERSECT, SNV_MERFIN,
│                            #   SNV_APPLY_CONSENSUS, PREPARE_NEXT_ROUND
└── workflows/
    ├── references.nf        # BUILD_REFS (round 1 — builds from raw FASTAs)
    │                        # BUILD_REFS_FROM_FILES (rounds 2+ — takes polished FASTAs)
    ├── mapping_r{1..5}.nf   # MAPPING_R{N} sub-workflows (per-round aliases)
    ├── deepvariant_r{1..5}.nf  # DEEPVARIANT_R{N} sub-workflows (per-round aliases)
    └── snv_candidates.nf    # SNV_CANDIDATES sub-workflow (bcftools + Merfin + consensus)
```

## Setup

### 1. Copy and edit the user config

```sh
ln -s /path/to/T2T-Polish/nextflow
cp nextflow/user.config.example user.config
# Edit user.config with your paths/globs
```

Key parameters in `user.config`:

```groovy
params {
  // Point verkko_asm at the Verkko output directory.
  // The assembly.*.fasta.gz files are expected inside it.
  verkko_asm = '/path/to/verkko-output'
  hap1_fasta_gz          = "${params.verkko_asm}/assembly.haplotype1.fasta.gz" // or pri for evaluation (ex. primary + sex chromosomes + MT + EBV)
  hap2_fasta_gz          = "${params.verkko_asm}/assembly.haplotype2.fasta.gz" // or alt for evaluation (ex. alternate + unassigned + unloc + rDNA_morphs)
  mito_exemplar_fasta_gz = "${params.verkko_asm}/assembly.mito.exemplar.fasta.gz"  // optional; omit if absent
  rdna_exemplar_fasta_gz = "${params.verkko_asm}/assembly.rdna.exemplar.fasta.gz"  // optional; omit if absent
  ebv_fasta_gz           = "${params.verkko_asm}/assembly.ebv.fasta.gz"  // optional; omit if absent

  // Glob patterns for reads. HiFi and ONT are required. Paired-end globs
  // (illumina/element) must contain a {1,2} or R{1,2} wildcard so Nextflow
  // can form R1/R2 pairs.
  // At least one of read_glob_illumina or read_glob_element MUST be set
  // when either platform is active (the pipeline will error otherwise).
  read_glob_hifi     = '/path/to/hifi/*.fastq.gz'              // required
  read_glob_ont      = '/path/to/ont/*.fastq.gz'               // required
  read_glob_illumina = '/path/to/illumina/*.R{1,2}.fastq.gz'   // required if 'illumina' in platforms
  read_glob_element  = '/path/to/element/*_R{1,2}.fastq.gz'    // required if 'element'  in platforms
                                                               // (element preferred over illumina)
  // Assembly name and version — combined as {asm_name}_{asm_ver} for all output
  // filenames (BAMs, VCFs, FASTAs).  asm_ver is auto-incremented each polishing
  // round (v0.1 → v0.2 → v0.3, …).  Both values MUST be quoted strings.
  asm_name = 'assembly'   // e.g. 'bTaeGut7'; make sure to quote
  asm_ver  = 'v0.1'       // e.g. 'v0.6'  — just the version tag, no name prefix

  // Root of the shared tools installation.  Used to locate binaries such as
  // merfin when they are not on PATH.  Set to the directory that contains
  // subdirectories like merfin/1.1/bin/, k8/, etc.
  tools = '/path/to/tools'

  // k8 binary for SAM→PAF conversion (paftools.js sam2paf).
  // https://github.com/attractivechaos/k8/releases
  k8 = '/path/to/k8/k8'
}
```
#### Tips

* Mixed R9/R10 data: R10 data is preferred over R9.
* ONT mapping: ONT is only mapping to `dip` by default whilest HiFi is always mapped to all 3 haps (`dip`, `hap1` and `hap2`). Set `--ont_map_haps true` in case ONT mappings to each haplotype is needed.
* Short-read mapping: `element` is preferred over `illumina`.

### 2. Run

* Submit as a job to slurm (default: 2 rounds, output in `results` and intermediates in `work`)

```sh
sbatch nextflow/run.sh user.config
```

* Resume
If it crashed for some reason, edit user.config to point to the intermediate paths of the results folder
```sh
sbatch nextflow/run.sh user.config -resume
```

* Run only 1 round of polishing
```sh
sbatch nextflow/run.sh user.config --polish_rounds 1
```

* Post-polish mode
Stop after mapping for evaluation purposes. Skip DeepVariant and SNV_Candidates
```sh
sbatch nextflow/run.sh user.config --run_dv false
```

* Keep per-read / per-pair intermediate BAMs in the results directory
```sh
sbatch nextflow/run.sh -c user.config --keep_intermediates true
```

* Use R9 ONT chemistry
```sh
sbatch nextflow/run.sh -c user.config --ont_chemistry r9
```

* For debugging
On an interactive node, it is possible to monitor the status of jobs in terminal. Not recommended for routine polishing as some jobes take more than a day.

```sh
module load nextflow

nextflow run nextflow/main.nf -c user.config
# or
nextflow run nextflow/main.nf -c user.config -resume
# Dry-run — validate params and print all registered process slots (no jobs submitted)
nextflow run nextflow/main.nf -c user.config -preview
# Stub run — executes the full dataflow graph; stub: blocks replace real commands
# with touch/mkdir stubs so no real work is done.  NOTE: publishDir still fires,
# so stub outputs ARE written to your results directory.  Use a throw-away outdir:
nextflow run nextflow/main.nf -c user.config -stub --outdir stub_test
```
**`-preview` shows registered process slots, not jobs that will run.**
A process appearing in `-preview` means it *could* run; it does
not mean it *will* run.  To confirm what actually executed, check
`nextflow log last -f name,status` after a real run.

### 3. Monitoring and logs

Nextflow writes logs in several places relative to the directory where you run `nextflow run`:

| Path | Contents |
|------|----------|
| `.nextflow.log` | Mains engine log — config/startup errors, process submissions, retries. Rotated to `.nextflow.log.1`, `.nextflow.log.2`, … on subsequent runs. |
| `work/<xx>/<hash>/` | One subdirectory per process execution. Contains `.command.sh` (exact script), `.command.log` (stdout + stderr), `.command.out`, `.command.err`, `.exitcode`. |

**Graceful shutdown on Slurm:**
When the pipeline head job is submitted via `run.sh`, SIGINT will be passed on to Nextflow for a graceful shutdown. This lets currently running jobs finish before exiting (equivalent to pressing Ctrl+C once interactively). Cancel child jobs or interactive head jobs with plain `scancel` to send SIGTERM/SIGKILL, which will immediately kill and stop the running step.

| Command |	Signal sent |	Effect on Nextflow |
| ------- | ----------- | ------------------ |
| scancel --signal=INT <jobid> | SIGINT (Ctrl+C) | Waits for running child jobs to finish, then exits cleanly |
| scancel <jobid> (default) | SIGTERM | Nextflow may not handle it gracefully — running child jobs can be orphaned |
| scancel --signal=KILL <jobid> | SIGKILL | Immediate kill, uncatchable — definitely orphans child jobs |

**Inspecting a failed run:**

```sh
# List all past runs:
nextflow log

# Show work directory, exit code, and status for the most recent run:
nextflow log last -f name,status,exit,work_dir

# Jump straight to the error output of any failed task:
nextflow log last -f name,status,exit,work_dir | grep FAILED
# then:
cat work/xx/hash.../.command.log
```


### 4. Outputs

```
results/
├── assemblies/
│   ├── bTaeGut7_v0.1.hap1.fa.gz  (+.fai, .gzi, .repetitive_k15.txt, BWA index files)   ← initial (from Verkko)
│   ├── bTaeGut7_v0.1.hap2.fa.gz
│   ├── bTaeGut7_v0.1.dip.fa.gz
│   ├── bTaeGut7_v0.2.hap1.fa.gz  (+.fai, .gzi, .repetitive_k15.txt)              ← Round 1 polished
│   ├── bTaeGut7_v0.2.hap2.fa.gz
│   ├── bTaeGut7_v0.2.dip.fa.gz
│   ├── bTaeGut7_v0.3.hap1.fa.gz  (+.fai, .gzi, .repetitive_k15.txt)              ← Round 2 polished
│   ├── bTaeGut7_v0.3.hap2.fa.gz
│   └── bTaeGut7_v0.3.dip.fa.gz
├── mapping/
│   ├── bTaeGut7_v0.1.hap1.hifi/
│   │   ├── *.sort.bam                          (only with --keep_intermediates true)
│   │   ├── bTaeGut7_v0.1.hap1.hifi.bam         (merged)
│   │   ├── bTaeGut7_v0.1.hap1.hifi.pri.bam     (filtered: -F0x104)
│   │   └── bTaeGut7_v0.1.hap1.hifi.pri.paf
│   ├── bTaeGut7_v0.1.hap1.illumina/
│   │   ├── *.dedup.pri.bam                     (only with --keep_intermediates true)
│   │   └── bTaeGut7_v0.1.hap1.illumina.dedup.pri.bam
│   ├── bTaeGut7_v0.1.hap1.hifi_illumina/
│   │   └── bTaeGut7_v0.1.hap1.hifi_illumina.bam  (hybrid merged BAM for DV Track A)
│   ├── bTaeGut7_v0.1.hap2.{hifi,ont,illumina,element,hifi_illumina}/
│   │   └── (same structure)
│   └── bTaeGut7_v0.1.dip.{hifi,ont,illumina,element,hifi_illumina}/
│       └── (same structure)
├── deepvariant/                  (disable with --run_dv false)
│   │
│   │  Track A — Hybrid (all three haps)
│   ├── bTaeGut7_v0.1.hap1.hifi_illumina.MQ5/   (hap1/hap2 → MQ 5)
│   │   ├── examples/                            (make_examples output; always published)
│   │   ├── dv_HYBRID_PACBIO_ILLUMINA_MQ5.hap1.vcf.gz   (+.tbi)
│   │   └── dv_HYBRID_PACBIO_ILLUMINA_MQ5.hap1.gvcf.gz  (+.tbi)
│   ├── bTaeGut7_v0.1.hap2.hifi_illumina.MQ5/
│   │   └── (same structure)
│   ├── bTaeGut7_v0.1.dip.hifi_illumina.MQ0/    (dip → MQ 0)
│   │   ├── examples/
│   │   ├── dv_HYBRID_PACBIO_ILLUMINA_MQ0.dip.vcf.gz    (+.tbi)
│   │   └── dv_HYBRID_PACBIO_ILLUMINA_MQ0.dip.gvcf.gz   (+.tbi)
│   │
│   │  Track B — ONT, dip only
│   │  R10 (default):
│   ├── bTaeGut7_v0.1.dip.ont.MQ0/
│   │   ├── examples/
│   │   ├── dv_ONT_R104_MQ0.dip.vcf.gz   (+.tbi)
│   │   └── dv_ONT_R104_MQ0.dip.gvcf.gz  (+.tbi)
│   │  R9 (--ont_chemistry r9):
│   └── bTaeGut7_v0.1.dip.ont.MQ0/
│       └── dv_ONT_R9_MQ0.dip.vcf.gz     (+.tbi)   (no gVCF)
└── snv_candidates/                  (disable with --run_snv_candidates false)
    │
    │  Round 1 outputs (ver_from=v0.1, ver_to=v0.2):
    ├── v0.1_to_v0.2/
    │   ├── snv_candidates.merfin-loose.vcf.gz (+.tbi)  — Merfin-validated SNV candidates
    │   ├── v0.2.dip.fa.gz   (+.fai, .gzi)  — polished dip assembly
    │   ├── v0.2.hap1.fa.gz  (+.fai, .gzi)  — hap1 extracted from polished dip
    │   ├── v0.2.hap2.fa.gz  (+.fai, .gzi)  — hap2 extracted from polished dip
    │   ├── v0.1_to_v0.2.dip.chain          — liftover chain (dip)
    │   └── intermediates/                  — only with --keep_snv_intermediates true
    │       ├── hybrid_to_dip.PASS.merfin-strict.vcf.gz (+.tbi)
    │       ├── hybrid_to_dip.PASS.vcf.gz (+.tbi)
    │       ├── hybrid_to_hap.PASS.vcf.gz (+.tbi)
    │       ├── ont_to_dip.PASS.vcf.gz (+.tbi)
    │       └── snv_pre_merfin.vcf.gz (+.tbi)
    │
    │  Round 2 outputs (ver_from=v0.2, ver_to=v0.3):
    └── v0.2_to_v0.3/
        ├── snv_candidates.merfin-loose.vcf.gz (+.tbi)
        ├── v0.3.dip.fa.gz   (+.fai, .gzi)
        ├── v0.3.hap1.fa.gz  (+.fai, .gzi)
        ├── v0.3.hap2.fa.gz  (+.fai, .gzi)
        ├── v0.2_to_v0.3.dip.chain
        └── intermediates/                  — only with --keep_snv_intermediates true
    (additional rounds follow the same pattern: v0.3_to_v0.4/, etc.)
```


### 5. Re-enter from existing results

By default, the pipeline retries a failed job one more time to avoid stoppings from server issues.

In anycases, if a run was corrupted or `work` doesn't exists, you can skip finished stages
by pointing the pipeline at the published `results`.
This avoids re-running from scratch. All re-entry params default to `null` (disabled).

| Param | Skips | Example value |
|---|---|---|
| `assemblies_dir` | `BUILD_HAP_REFERENCES` + `BUILD_DIP_REFERENCE` | `'results/assemblies'` |
| `mapping_dir` | WM mapping, BWA mapping, and MERGE_HYBRID (per hap and per platform combination that is already present) | `'results/mapping'` |
| `deepvariant_dir` | DV `make_examples` step 1 (per hap — only missing haps are re-run) | `'results/deepvariant'` |

Set `mapping_dir` along with `deepvariant_dir`.

```groovy
// In user.config — example: re-enter at call_variants (step 2), all examples done
params {
  // Pre-built refs — skips BUILD_HAP_REFERENCES + BUILD_DIP_REFERENCE
  // Filenames inferred: bTaeGut7_v0.6.{hap1,hap2,dip}.fa.gz (+.gzi, .fai)
  assemblies_dir  = 'results/assemblies'
  // Pre-built BAMs — skips WM + BWA mapping and MERGE_HYBRID (per hap)
  mapping_dir     = 'results/mapping'
  // Pre-built examples — skips make_examples (step 1)
  deepvariant_dir = 'results/deepvariant'
}
```

### 6. Adjust resources for your cluster (optional)

Resource allocations (CPUs, memory, wall-time, queue names, `--gres`) are in
`resources.config`.  The defaults target the **NIH Biowulf** Slurm cluster.
If you are on a different system:

```sh
cp resources.config my_resources.config
# Edit queue names, memory limits, clusterOptions, executor, etc.
```

Then tell Nextflow to use your copy instead of the default one.
Either add this line at the **top** of your `user.config`:

```groovy
includeConfig '/absolute/path/to/my_resources.config'
```

Or pass it on the command line (it overrides `resources.config`):

```sh
nextflow run main.nf -c user.config -c my_resources.config
```

For non-Slurm schedulers change `executor = 'slurm'` to the appropriate
Nextflow executor name (`'sge'`, `'lsf'`, `'pbspro'`, `'local'`, …).
See the [Nextflow executor docs](https://www.nextflow.io/docs/latest/executor.html).

### 7. Full parameters

| Parameter               | Default                       | Description                         |
|-------------------------|-------------------------------|-------------------------------------|
| `params.outdir`         | `results`                     | Output directory               |
| `params.assemblies_dir` | `null`                        | Path to a directory containing pre-built `{asm_name}_{asm_ver}.{hap}.fa.gz` files (+ `.gzi`, `.fai`). When set, `BUILD_HAP_REFERENCES` and `BUILD_DIP_REFERENCE` are skipped and the Verkko FASTA params are not required. |
| `params.mapping_dir`    | `null`                        | Path to the `mapping/` output directory. BAMs are in per-hap subdirectories (`{pfx}.{hap}.{plat}/` for WM/BWA, `{pfx}.{hap}.{combo}/` for hybrid). All are discovered automatically from the `{asm_name}_{asm_ver}` prefix. Only missing combinations are re-submitted. |
| `params.deepvariant_dir`| `null`                        | Path to the `deepvariant/` output directory. DV `examples/` dirs are discovered automatically. Only missing haps are re-submitted. |
| `params.asm_name`       | `assembly`                    | Assembly name prefix (e.g. `'bTaeGut7'`). Combined with `asm_ver` to form filenames: `{asm_name}_{asm_ver}.{hap}.…` Must be a quoted string e.g. `''` in the config. |
| `params.asm_ver`        | `v0.1`                        | Version tag for the initial assembly (e.g. `'v0.6'`). Auto-bumped each polishing round. Must be a quoted string e.g. `''` — without quotes Groovy parses `v0.6` as a method call and throws an error. |
| `params.platforms`      | `hifi,ont,illumina`   | Comma-separated list to run         |
| `params.samtools`       | `samtools`                    | Path to samtools executable         |
| `params.mapping_outdir` | `${outdir}/mapping`           | Output directory for BAMs/PAFs and hybrid merged BAMs |
| `params.keep_intermediates` | `false`                   | Keep per-read/per-pair intermediate BAMs in results (see below) |
| `params.ont_map_haps`       | `false`                   | Map ONT reads to hap1 and hap2 in addition to dip. ONT→dip is sufficient for polishing. |
| `params.dv_outdir`          | `${outdir}/deepvariant`   | Output directory for DV VCFs        |
| `params.dv_sample`          | _(hap tag)_               | Sample name written into VCF header |
| `params.dv_mq_hap`          | `5`                       | MQ filter applied to hap1/hap2 merged BAMs |
| `params.dv_mq_dip`          | `0`                       | MQ filter applied to dip merged BAMs |
| `params.dv_n_shard`         | `12`                      | Number of `make_examples` shards (= step-1 CPUs) |
| `params.dv_long_platforms`  | `hifi`                    | Long-read platform(s) used in hybrid merge |
| `params.dv_short_platforms` | `illumina,element`        | Short-read platform(s) used in hybrid merge |
| `params.ont_chemistry`      | `r10`                     | ONT chemistry for Track B dip calling: `r10` (DeepVariant `ONT_R104`) or `r9` (PEPPER-Margin-DV). **If you have mixed R9/R10 data, use r10 only.** |
| `params.run_snv_candidates`  | `true`                    | Run SNV candidate collection + Merfin after DeepVariant. Disable with `--run_snv_candidates false`. |
| `params.snv_outdir`          | `${outdir}/snv_candidates`| Output directory for SNV candidate VCFs        |
| `params.hybrid_meryl`        | _(required)_              | Path to hybrid (HiFi + Illumina/Element) read k-mer meryl database (`*.k31.meryl` dir) |
| `params.merfin_peak`         | _(required)_              | Integer peak coverage value for Merfin (from GenomeScope / meryl histogram) |
| `params.merfin`              | `merfin`                  | Path or command name for the `merfin` binary    |
| `params.keep_snv_intermediates` | `false`               | Publish intermediate VCFs from `SNV_FILTER_INTERSECT` |
| `params.asm_ver_next`        | _(auto-bumped)_           | Version tag for the first polished assembly (e.g. `v0.2`); if unset, the trailing integer of `asm_ver` is incremented |
| `params.polish_rounds`       | `2`                       | Number of SNV polishing rounds to run (1–5). |

## Execution graph

```
main.nf
│
│  ╔══════════════════════════════════════════════════════════════════════════════════╗
├──╢  Round 1                                                                         ║
│  ╚══════════════════════════════════════════════════════════════════════════════════╝
│
│  BUILD_REFS()                                            (workflows/references.nf)
│  │
│  ├── BUILD_HAP_REFERENCES(hap1+mito+ebv+rdna) ────────────────┐
│  ├── BUILD_HAP_REFERENCES(hap2+mito+ebv+rdna) ────────────────┤
│  └── BUILD_DIP_REFERENCE(hap1+hap2+mito+ebv+rdna) ────────────┤
│                                                               │ allRefs (hap1, hap2, dip)
│                                           ┌───────────────────┤
│                                           │                   │
│                              ┌────────────▼────────────┐  ┌───▼───────────────┐
│                              │  MERYL_REPETITIVE  (×3) │  │  BWA_INDEX  (×3)  │
│                              │  [quick, 12c, 24g, 30m] │  │  [quick, 4c, 10g] │
│                              └────────────┬────────────┘  └───┬───────────────┘
│                                           │ wm_refs           │ bwa_refs
│                                           │                   │
│   ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ - ┼ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┤
│  MAPPING(wm_refs, bwa_refs)               │                   │  (workflows/mapping_r1.nf)
│  hifi: ×3 haps   ont: ×1 dip (×3 with --ont_map_haps true)    │
│  re-entry: mapping_dir skips WM + BWA mapping entirely        │
│                                           │                   │
│  wm_refs.combine(reads)                   │                   │
│     WINNOWMAP_MAP    (×haps×N)  [norm, 24c, 120g, 2d]         │
│     WINNOWMAP_MERGE  (×4-6)     [norm, 48c,  60g, 1d]         │
│     WINNOWMAP_FILTER (×4-6)     [norm, 12c,   8g, 1d] ──── wm_pri_bams
│     SAM2PAF          (×4-6)     [norm, 12c,   8g, 1d]         │
│     (×4: hifi→all3 + ont→dip; ×6: + ont_map_haps)             |
│                                                               │
│  bwa_refs.combine(read pairs)                                 │
│     BWA_MAP   (×3×M)  [norm, 24c, 120g, 2d] ──────────────────┤
│     BWA_MERGE (×3–6)  [norm, 24c,  48g, 1d] ───────── bwa_mrg_bams
│     (×3 per short-read platform × 3 haps; ×6 with illumina+element)
│                                                               │
│   ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ┤
│  DEEPVARIANT(refs, wm_pri_bams, bwa_mrg_bams)                 │  (workflows/deepvariant_r1.nf)
│                                                               │
│  Track A — Hybrid (HiFi + short-read), all haps               │
│     MERGE_HYBRID (×3)           [norm, 12c,  24g,  4h] ←──────┘  refs.combine(by: hap)
│       output → mapping/  (re-entry: mapping_dir skips this step per hap)
│     MQ: hap1/hap2 → 5,  dip → 0
│     DV_MAKE_EXAMPLES (×3)       [norm, 12c,  36g,  3d, 1000g scratch]
│       output → deepvariant/*/examples/  (always published)
│       re-entry: deepvariant_dir skips per-hap items that already exist
│     DV_CALL_VARIANTS (×3)       [gpu,   4c,  48g, 12h]
│     DV_POSTPROCESS   (×3)       [norm, 12c, 120g, 12h]
│        └─→ hybrid_{hap1,hap2,dip}.vcf.gz ────────────────────────┐
│                                                                  │
│  Track B — ONT, dip only  (runs in parallel with Track A)        │
│     R10:  DV_MAKE_EXAMPLES → DV_CALL_VARIANTS → DV_POSTPROCESS   │
│     R9:   PEPPER_MARGIN_DV (×N_chr) → DV_MERGE_CHR_VCFS          │
│        └─→ ont_dip.vcf.gz ───────────────────────────────────────┤
│                                                                  │
│   ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ── ┤
│  SNV_CANDIDATES(dv_vcfs, refs, ver[0]→ver[1])                    │  (workflows/snv_candidates.nf)
│  [disable with --run_snv_candidates false]                       │
│                                                                  │
│     hybrid_hap1 + hybrid_hap2 + hybrid_dip + ont_dip VCFs ←──────┘
│        │
│     SNV_FILTER_INTERSECT  [quick, 12c, 12g, 1h]
│        reheader → PASS filter → hap concat
│        → isec (hybrid_hap ∩ hybrid_dip, hom only)
│        → GT/GQ/AF filter → isec (hap_het ∩ ont_hom ∩ dip_alt)
│        │
│     SNV_MERFIN  [norm, 12c, 120g, 12h]
│        meryl count k=31 (dip → seqmer)
│        → merfin -strict → concat → GT→1/1 → merfin -loose
│        │
│     SNV_APPLY_CONSENSUS  [quick, 12c, 12g, 1h]
│        bcftools consensus -H1 --chain  (dip only)
│        | bgzip --index  →  ver[1].dip.fa.gz (+.gzi)  +  ver[0]_to_ver[1].dip.chain
│        │
│     PREPARE_NEXT_ROUND  [quick, 12c, 12g, 1h]
│        samtools faidx -r <hap1_names> | bgzip --index  →  ver[1].hap1.fa.gz (+.gzi)
│        samtools faidx -r <hap2_names> | bgzip --index  →  ver[1].hap2.fa.gz (+.gzi)
│        cp dip.fa.gz                            →  ver[1].dip.fa.gz (+.gzi)
│        └─→ next_refs: [hap, fa.gz, gzi, fai] ×3
│
│  ╔══════════════════════════════════════════════════════════════════════════════════╗
├──╢  Round 2  (if params.polish_rounds ≥ 2; same structure as Round 1)               ║
│  ╚══════════════════════════════════════════════════════════════════════════════════╝
│
│  BUILD_REFS_FROM_FILES(next_refs)                        (workflows/references.nf)
│     MERYL_REPETITIVE + BWA_INDEX on polished FASTAs  (skips BUILD_HAP_REFERENCES)
│     └─→ wm_refs, bwa_refs
│  MAPPING_R2      → wm_pri_bams, bwa_mrg_bams
│  DEEPVARIANT_R2  → hybrid_{hap1,hap2,dip}.vcf.gz + ont_dip.vcf.gz
│  SNV_CANDIDATES_R2(dv_vcfs, refs, ver[1]→ver[2])
│     └─→ ver[2].{hap1,hap2,dip}.fa.gz  +  ver[1]_to_ver[2].dip.chain
│         next_refs  →  Round 3
│
│  ╔══════════════════════════════════════════════════════════════════════════════════╗
└──╢  Rounds 3–5 follow the same pattern  (params.polish_rounds, default: 2)          ║
   ╚══════════════════════════════════════════════════════════════════════════════════╝

Final outputs per round  (snv_candidates/):
  ver[0]_to_ver[1]/  ver[1].{hap1,hap2,dip}.fa.gz  +  ver[0]_to_ver[1].dip.chain   ← Round 1
  ver[1]_to_ver[2]/  ver[2].{hap1,hap2,dip}.fa.gz  +  ver[1]_to_ver[2].dip.chain   ← Round 2
  …
```

**Concurrency notes:**
- `BUILD_HAP_REFERENCES` builds references for hap1 and hap2 alongside `BUILD_DIP_REFERENCE` which builds the dip reference. When `assemblies_dir` is set, all three refs are loaded from disk and these builders are skipped.
- `MERYL_REPETITIVE` and `BWA_INDEX` are submitted independently as soon as each ref is available.
- `mapping_dir` re-entry is only implemented for round 1. The workflow scans the published `results/mapping/` tree, remaps only the hap/platform combinations that are missing, and mixes those new outputs with the files already on disk.
- `deepvariant_dir` re-entry is also round-1 only. Existing `examples/` dirs and `call_variants_output` shards are reused per hap/combo/MQ item; anything missing is re-run, and `DV_POSTPROCESS` still needs the examples dir staged for each item.
- `WINNOWMAP_MAP` (HiFi) and `BWA_MAP` submits across the active haplotypes and read groups. `WINNOWMAP_MAP` for ONT maps only to dip by default. Within each haplotype, every read file / read pair gets its own Slurm job.
- Track A (hybrid) and Track B (ONT) run in parallel. The hybrid merge/call/postprocess chain runs independently per hap, while the ONT branch runs only on dip.
- `PEPPER_MARGIN_DV` and `DV_MERGE_CHR_VCFS` are always registered, but the r9 path receives input only when `--ont_chemistry r9`; the r10 path is the active one otherwise. R9 scatters one GPU job per chromosome, then merges the per-chromosome VCFs.
- `SNV_FILTER_INTERSECT` waits for all four DV VCFs (hybrid hap1, hybrid hap2, hybrid dip, ONT dip) to be ready before starting. `SNV_MERFIN`, `SNV_APPLY_CONSENSUS`, and `PREPARE_NEXT_ROUND` run sequentially within each round.
- `PREPARE_NEXT_ROUND` feeds the next round's `BUILD_REFS_FROM_FILES`, using the polished FASTAs.


## Process → original script mapping

| Original Slurm script      | Nextflow process     | Label / resources (as in resources.config)      |
|----------------------------|----------------------|-------------------------------------------------|
| _(ref build — hap)_        | `BUILD_HAP_REFERENCES` | `norm_build_ref` (4 CPU / 16g / 2h)           |
| _(ref build — dip)_        | `BUILD_DIP_REFERENCE`  | `norm_build_ref` (4 CPU / 16g / 2h)           |
| `winnowmap/init.sh`        | `MERYL_REPETITIVE`   | `quick_meryl` (12 CPU / 24g / 30min)            |
| `winnowmap/map.sh`         | `WINNOWMAP_MAP`      | `norm_map` (24 CPU / 120g / 2d / 900g scratch)  |
| `winnowmap/merge.sh`       | `WINNOWMAP_MERGE`    | `norm_merge_wm` (48 CPU / 60g / 1d)             |
| `winnowmap/filt.sh`        | `WINNOWMAP_FILTER`   | `norm_filter` (12 CPU / 8g / 1d)                |
| `coverage/sam2paf.sh`      | `SAM2PAF`            | `norm_filter` (12 CPU / 8g / 1d)                |
| `bwa/bwa_index.sh`         | `BWA_INDEX`          | `quick_small` (4 CPU / 10g / 4h)                |
| `bwa/bwa.sh`               | `BWA_MAP`            | `norm_bwa_map` (24 CPU / 120g / 5d / 2000g scratch) |
| `bwa/merge.sh`             | `BWA_MERGE`          | `norm_merge_bwa` (24 CPU / 48g / 1d)            |
| `deepvariant/merge_hybrid.sh` | `MERGE_HYBRID`    | `norm_merge_hybrid` (12 CPU / 24g / 4h)         |
| `deepvariant/step1_with_minqual.sh` | `DV_MAKE_EXAMPLES` | `norm_dv_make_examples` (12 CPU / 48g / 10d / 1000g scratch) |
| `deepvariant/step2_with_minqual.sh` | `DV_CALL_VARIANTS` | `norm_dv_call_variants` (12 CPU / 48g / 12h / 1 GPU) |
| `deepvariant/step3_with_minqual.sh` | `DV_POSTPROCESS`   | `norm_dv_postprocess` (12 CPU / 160g / 12h)    |
| `deepvariant/ont_r9_pepper_margin_dv.sh` | `PEPPER_MARGIN_DV` | `norm_pepper_margin_dv` (24 CPU / 48g / 2h / k80×4 GPU, per-chr) |
| `deepvariant/merge_per_chr_vcfs.sh` | `DV_MERGE_CHR_VCFS` | `quick_small` (4 CPU / 10g / 4h)             |
| `variant_call/snv_candidates.sh` (bcftools steps) | `SNV_FILTER_INTERSECT` | `quick_snv_filter` (12 CPU / 12g / 1h) |
| `variant_call/snv_candidates.sh` (merfin steps)   | `SNV_MERFIN`           | `norm_snv_merfin` (12 CPU / 120g / 12h) |
| _(apply candidates to each hap assembly)_         | `SNV_APPLY_CONSENSUS`  | `quick_snv_filter` (12 CPU / 12g / 1h) |
| _(rebuild hap1/hap2/dip for the next round)_      | `PREPARE_NEXT_ROUND`   | `quick_snv_filter` (12 CPU / 12g / 1h) |

## Disk space

### Where files are written

Nextflow writes **two separate copies** of every output file:

| Location | Purpose | Read by `-resume`? |
|---|---|---|
| `work/<hash>/` | Canonical task outputs; source of truth for the cache | ✅ yes |
| `results/` (publishDir) | Human-readable final outputs; hard-links or copies from `work/` | ❌ no |

`results/` is a **one-way publication sink**.  Nextflow never reads it back by default.
On `-resume`, the cache is checked against `work/` only unless [re-entry is attempted](### 5. Re-enter from existing results).
If `work/` is deleted, every task re-runs from scratch even if `results/` is fully intact.

### Intermediates and `work/`

Per-read `*.sort.bam` (Winnowmap) and per-pair `*.dedup.pri.bam` (BWA) files
are always written to `work/` — the `--keep_intermediates` flag only controls
whether they are also **copied into `results/`**.  With many read files these
intermediates can be the largest consumers of `work/` space.

### Cleaning up

**After the pipeline completes successfully**, `work/` can be removed to
reclaim disk space.  This permanently breaks `-resume` for that run.

```sh
# List all recorded runs (shows run names and IDs)
nextflow log

# Remove work/ directories for a single completed run (by run name or ID)
nextflow clean -f <run-name>

# Remove work/ directories for all runs recorded in .nextflow/
nextflow clean -f

# Dry-run: show what would be deleted without actually deleting
nextflow clean -n
nextflow clean -n <run-name>
```

**Tip — selective cleaning with `-but-less-recent`:**
If you have a long run that you want to keep resumable but want to free
space from earlier failed attempts, use:
```sh
nextflow clean -f -but-less-recent <run-name>
```
This keeps the most recent run's `work/` intact and deletes all older ones.

### Recommended workflow for large assemblies

1. Run the pipeline to completion.
2. Verify `results/` looks correct.
3. Archive `results/` to long-term storage (e.g. `rsync` to a data store).
4. Run `nextflow clean -f <run-name>` to delete `work/`.
5. Delete intermediate BAMs from `results/` if they were published
