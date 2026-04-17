# T2T-Polish

Nextflow DSL2 pipeline for polishing a T2T diploid assembly.

---

## Overview

The pipeline takes a [Verkko](https://github.com/marbl/verkko) assembly and maps HiFi, ONT, Illumina, and/or Element reads against three reference builds — `hap1`, `hap2`, and `dip` — all in parallel. Slurm submission is handled natively by Nextflow; no external submit scripts are needed.


```
nextflow/
├── ma> - `PEPPER_MARGIN_DV` and `DV_MERGE_CHR_VCFS` are always listed even when
>   `--ont_chemistry r10` is set.  The R9/R10 branch is evaluated at runtime
>   from channel data; `-preview` sees both branches as registered slots.                  # Entry point: params + include + workflow {}
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

---

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
│  re-entry: mapping_dir skips WM + BWA mapping entirely          │
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
- `BUILD_HAP_REFERENCES` for hap1 and hap2 runs simultaneously with `BUILD_DIP_REFERENCE` for dip.
- `MERYL_REPETITIVE` and `BWA_INDEX` fire in parallel the moment each ref is ready — they do not wait for each other.
- `WINNOWMAP_MAP` (HiFi) and `BWA_MAP` fan out across all three haplotypes simultaneously. `WINNOWMAP_MAP` for ONT maps only to dip by default; set `--ont_map_haps true` to also map to hap1 and hap2. Within each haplotype, every read file / read pair gets its own Slurm job.
- Track A (hybrid) and Track B (ONT) run fully in parallel — the hybrid jobs and the ONT-only dip jobs are independent.
- `MERGE_HYBRID` → DV three-step runs independently for each of hap1, hap2, dip in parallel, each with its own reference.
- R9 Track B scatters one `PEPPER_MARGIN_DV` GPU job per chromosome; all chromosomes run simultaneously then merge.
- `SNV_FILTER_INTERSECT` waits for all four DV VCFs (hybrid hap1, hybrid hap2, hybrid dip, ONT dip) to be ready before starting. `SNV_MERFIN`, `SNV_APPLY_CONSENSUS`, and `PREPARE_NEXT_ROUND` run sequentially within each round.
- Each polishing round runs fully sequentially after the previous one — Round N's `PREPARE_NEXT_ROUND` output feeds Round N+1's `BUILD_REFS_FROM_FILES`. Rounds 2+ skip `BUILD_HAP_REFERENCES`/`BUILD_DIP_REFERENCE` and run only `MERYL_REPETITIVE` + `BWA_INDEX` on the polished FASTAs.
- **Mixed R9/R10 data:** set `--ont_chemistry r10`. R10 is preferred — it uses DeepVariant's native `ONT_R104` model and produces a gVCF alongside the VCF. R9 requires the separate PEPPER-Margin-DeepVariant toolchain and does not emit a gVCF.
- `×3` = one job per hap (hap1, hap2, dip); `×N` / `×M` = one job per read file / pair; `×N_chr` = one job per chromosome/contig.
- **Winnowmap merge/filter/PAF**: `×4` by default with HiFi mapped to all three haps + ONT→dip (HiFi always maps all haps); `×6` with `--ont_map_haps true` (ONT also maps to hap1 and hap2).
- **BWA merge**: `×3` per active short-read platform (always all three haps); `×6` when both `illumina` and `element` are in `--platforms`.

---

## Process → original script mapping

| Original Slurm script      | Nextflow process     | Label / resources                               |
|----------------------------|----------------------|-------------------------------------------------|
| _(ref build — hap)_        | `BUILD_HAP_REFERENCES` | `norm_build_ref` (4 CPU / 16g / 2h)            |
| _(ref build — dip)_        | `BUILD_DIP_REFERENCE`  | `norm_build_ref` (4 CPU / 16g / 2h)            |
| `winnowmap/init.sh`        | `MERYL_REPETITIVE`   | `quick_meryl` (12 CPU / 24g / 30min)            |
| `winnowmap/map.sh`         | `WINNOWMAP_MAP`      | `norm_map` (24 CPU / 120g / 2d / 900g scratch)  |
| `winnowmap/merge.sh`       | `WINNOWMAP_MERGE`    | `norm_merge_wm` (48 CPU / 60g / 1d)             |
| `winnowmap/filt.sh`        | `WINNOWMAP_FILTER`   | `norm_filter` (12 CPU / 8g / 1d)                |
| `coverage/sam2paf.sh`      | `SAM2PAF`            | `norm_filter` (12 CPU / 8g / 1d)                |
| `bwa/bwa_index.sh`         | `BWA_INDEX`          | `quick_small` (4 CPU / 10g / 4h)                |
| `bwa/bwa.sh`               | `BWA_MAP`            | `norm_bwa_map` (24 CPU / 120g / 2d / 2000g scratch) |
| `bwa/merge.sh`             | `BWA_MERGE`          | `norm_merge_bwa` (24 CPU / 48g / 1d)            |
| `deepvariant/merge_hybrid.sh` | `MERGE_HYBRID`    | `norm_merge_hybrid` (12 CPU / 24g / 4h)         |
| `deepvariant/step1_with_minqual.sh` | `DV_MAKE_EXAMPLES` | `norm_dv_make_examples` (12 CPU / 36g / 3d / 1000g scratch) |
| `deepvariant/step2_with_minqual.sh` | `DV_CALL_VARIANTS` | `norm_dv_call_variants` (12 CPU / 48g / 12h / 1 GPU) |
| `deepvariant/step3_with_minqual.sh` | `DV_POSTPROCESS`   | `norm_dv_postprocess` (12 CPU / 120g / 12h)    |
| `deepvariant/ont_r9_pepper_margin_dv.sh` | `PEPPER_MARGIN_DV` | `norm_pepper_margin_dv` (24 CPU / 48g / 2h / k80×4 GPU, per-chr) |
| `deepvariant/merge_per_chr_vcfs.sh` | `DV_MERGE_CHR_VCFS` | `quick_small` (4 CPU / 10g / 4h)             |
| `variant_call/snv_candidates.sh` (bcftools steps) | `SNV_FILTER_INTERSECT` | `quick_snv_filter` (12 CPU / 12g / 1h) |
| `variant_call/snv_candidates.sh` (merfin steps)   | `SNV_MERFIN`           | `norm_snv_merfin` (12 CPU / 120g / 12h) |
| _(apply candidates to each hap assembly)_         | `SNV_APPLY_CONSENSUS`  | `quick_snv_filter` (12 CPU / 12g / 1h) |
| _(rebuild hap1/hap2/dip for the next round)_      | `PREPARE_NEXT_ROUND`   | `quick_snv_filter` (12 CPU / 12g / 1h) |

---

## Setup

### 1. Copy and edit the user config

```sh
cp user.config.example user.config
# Edit user.config with your paths/globs
```

Key parameters in `user.config`:

```groovy
params {
  // Point verkko_asm at the Verkko output directory.
  // The assembly.*.fasta.gz files are expected inside it.
  verkko_asm = '/path/to/verkko-output'

  hap1_fasta_gz          = "${params.verkko_asm}/assembly.haplotype1.fasta.gz"
  hap2_fasta_gz          = "${params.verkko_asm}/assembly.haplotype2.fasta.gz"
  ebv_fasta_gz           = "${params.verkko_asm}/assembly.ebv.fasta.gz"  // optional; omit if absent
  mito_exemplar_fasta_gz = "${params.verkko_asm}/assembly.mito.exemplar.fasta.gz"
  rdna_exemplar_fasta_gz = "${params.verkko_asm}/assembly.rdna.exemplar.fasta.gz"

  // Glob patterns for reads. HiFi and ONT are required. Paired-end globs
  // (illumina/element) must contain a {1,2} or R{1,2} wildcard so Nextflow
  // can form R1/R2 pairs.
  // At least one of read_glob_illumina or read_glob_element MUST be set
  // when either platform is active (the pipeline will error otherwise).
  read_glob_hifi     = '/path/to/hifi/*.fastq.gz'              // required
  read_glob_ont      = '/path/to/ont/*.fastq.gz'               // required
  read_glob_illumina = '/path/to/illumina/*.R{1,2}.fastq.gz'   // required if 'illumina' in platforms
  read_glob_element  = '/path/to/element/*_R{1,2}.fastq.gz'    // required if 'element'  in platforms
                                                               // (at least one of the two must be set)

  // k8 binary for SAM→PAF conversion (paftools.js sam2paf).
  // https://github.com/attractivechaos/k8/releases
  k8 = '/path/to/k8/k8'
}
```

Optional parameters (with defaults):

| Parameter               | Default                       | Description                         |
|-------------------------|-------------------------------|-------------------------------------|
| `params.outdir`         | `results`                     | Root output directory               |
| `params.assemblies_dir` | `null`                        | Path to a directory containing pre-built `{asm_name}_{asm_ver}.{hap}.fa.gz` files (+ `.gzi`, `.fai`). When set, `BUILD_HAP_REFERENCES` and `BUILD_DIP_REFERENCE` are skipped and the Verkko FASTA params are not required. |
| `params.mapping_dir`    | `null`                        | Path to the `mapping/` output directory. BAMs are in per-hap subdirectories (`{pfx}.{hap}.{plat}/` for WM/BWA, `{pfx}.{hap}.{combo}/` for hybrid); all are discovered automatically from the `{asm_name}_{asm_ver}` prefix. WM and BWA mapping are skipped; hybrid merge uses a per-item join (only missing haps are re-merged). |
| `params.deepvariant_dir`| `null`                        | Path to the `deepvariant/` output directory. DV `examples/` dirs are discovered automatically. `make_examples` uses a per-item join (only missing haps are re-run). |
| `params.asm_name`       | `assembly`                    | Assembly name prefix (e.g. `bTaeGut7`). Combined with `asm_ver` to form filenames: `{asm_name}_{asm_ver}.{hap}.…` Must be a **quoted string** in the config. |
| `params.asm_ver`        | `v0.1`                        | Version tag for the initial assembly (e.g. `v0.6`). Auto-bumped each polishing round. Must be a **quoted string** — without quotes Groovy parses `v0.6` as a method call and throws an error. |
| `params.platforms`      | `hifi,ont,illumina,element`   | Comma-separated list to run         |
| `params.samtools`       | `samtools`                    | Path to samtools executable         |
| `params.mapping_outdir` | `${outdir}/mapping`           | Output directory for BAMs/PAFs and hybrid merged BAMs |
| `params.keep_intermediates` | `false`                   | Keep per-read/per-pair intermediate BAMs in results (see below) |
| `params.ont_map_haps`       | `false`                   | Also map ONT reads to hap1 and hap2 (in addition to dip). Off by default — ONT→dip is sufficient for DeepVariant. |
| `params.dv_outdir`          | `${outdir}/deepvariant`   | Output directory for DV VCFs        |
| `params.dv_sample`          | _(hap tag)_               | Sample name written into VCF header |
| `params.dv_mq_hap`          | `5`                       | MQ filter applied to hap1/hap2 merged BAMs |
| `params.dv_mq_dip`          | `0`                       | MQ filter applied to dip merged BAMs |
| `params.dv_n_shard`         | `12`                      | Number of `make_examples` shards (= step-1 CPUs) |
| `params.dv_long_platforms`  | `hifi`                    | Long-read platform(s) used in hybrid merge |
| `params.dv_short_platforms` | `illumina,element`        | Short-read platform(s) used in hybrid merge |
| `params.ont_chemistry`      | `r10`                     | ONT chemistry for Track B dip calling: `r10` (DeepVariant `ONT_R104`) or `r9` (PEPPER-Margin-DV). **If you have mixed R9/R10 data, use `r10`.** |
| `params.keep_dv_intermediates` | `false`                | Publish `make_examples` tfrecords and `call_variants` output. The `examples/` directory is always published regardless. |
| `params.run_snv_candidates`  | `true`                    | Run SNV candidate collection + Merfin after DeepVariant. Disable with `--run_snv_candidates false`. |
| `params.snv_outdir`          | `${outdir}/snv_candidates`| Output directory for SNV candidate VCFs        |
| `params.hybrid_meryl`        | _(required)_              | Path to hybrid (HiFi + Illumina/Element) read k-mer meryl database (`*.k31.meryl` dir) |
| `params.merfin_peak`         | _(required)_              | Integer peak coverage value for Merfin (from GenomeScope / meryl histogram) |
| `params.merfin`              | `merfin`                  | Path or command name for the `merfin` binary    |
| `params.keep_snv_intermediates` | `false`               | Publish intermediate VCFs from `SNV_FILTER_INTERSECT` |
| `params.asm_ver_next`        | _(auto-bumped)_           | Version tag for the first polished assembly (e.g. `v0.2`); if unset, the trailing integer of `asm_ver` is incremented |
| `params.polish_rounds`       | `2`                       | Number of SNV polishing rounds to run (1–5). |

> **Intermediate BAMs** — by default the per-read `*.sort.bam` files produced
> by `WINNOWMAP_MAP` and the per-pair `*.dedup.pri.bam` files produced by
> `BWA_MAP` are *not* copied to the results directory; they remain available
> in Nextflow's `work/` cache for re-use with `-resume`.  Set
> `--keep_intermediates true` on the command line (or
> `params.keep_intermediates = true` in `user.config`) to also publish them
> alongside the merged BAMs.

### 3. Re-enter from existing results (optional)

If a run was interrupted after some stages completed, you can skip those stages
by pointing the pipeline at the already-published files instead of re-running
from scratch. All re-entry params default to `null` (disabled).

| Param | Skips | Example value |
|---|---|---|
| `assemblies_dir` | `BUILD_HAP_REFERENCES` + `BUILD_DIP_REFERENCE` | `'/path/to/results/assemblies'` |
| `mapping_dir` | WM mapping, BWA mapping, and MERGE_HYBRID (per hap — only missing hybrid BAMs are re-merged) | `'/path/to/results/mapping'` |
| `deepvariant_dir` | DV `make_examples` step 1 (per hap — only missing haps are re-run) | `'/path/to/results/deepvariant'` |

> **Use literal paths.** Cross-file `params.*` references in GStrings are not
> resolved at config parse time and cause `Unknown config attribute` errors.
>
> `assemblies_dir` must point to a directory containing
> `{asm_name}_{asm_ver}.{hap1,hap2,dip}.fa.gz` with matching `.gzi` and `.fai`
> files alongside — exactly as published by `BUILD_HAP_REFERENCES` /
> `BUILD_DIP_REFERENCE`.  The filenames are inferred from `asm_name` and
> `asm_ver`; you do not need to specify them individually.  When set, the
> Verkko FASTA params (`hap1_fasta_gz`, etc.) are not required.
> `MERYL_REPETITIVE` and `BWA_INDEX` still run unless their outputs
> (`{asm_name}_{asm_ver}.{hap}.repetitive_k15.txt` and the five BWA index
> files) are also present in that directory.
>
> `mapping_dir` must point to the `mapping/` output directory.  Winnowmap pri
> BAMs (`{pfx}.{hap}.{plat}/{pfx}.{hap}.{plat}.pri.bam`), BWA BAMs
> (`{pfx}.{hap}.{plat}/{pfx}.{hap}.{plat}.dedup.pri.bam`), and hybrid merged
> BAMs (`{pfx}.{hap}.{combo}/{pfx}.{hap}.{combo}.bam`) are all discovered
> automatically from the `{asm_name}_{asm_ver}` prefix.  WM mapping, BWA
> mapping, and `MERGE_HYBRID` are all skipped entirely at graph-build time —
> the process is never wired into the execution graph.
>
> `deepvariant_dir` must point to the `deepvariant/` output directory.
> Examples dirs (`{pfx}.{hap}.{combo}.MQ{mq}/examples`) are discovered
> automatically.  `DV_MAKE_EXAMPLES` is skipped entirely at graph-build time —
> the process is never wired into the execution graph.  If a hap is missing
> from `deepvariant_dir`, the pipeline will silently emit no items for that
> hap rather than re-running `make_examples`; if you need partial re-entry,
> unset `deepvariant_dir` and provide `mapping_dir` instead so BAMs are
> available for the missing hap.  MAPPING still runs when only
> `deepvariant_dir` is set (no `mapping_dir` provided).

```groovy
// In user.config — example: re-enter at call_variants (step 2), all examples done
params {
  asm_name        = 'bTaeGut7'              // must be a quoted string
  asm_ver         = 'v0.6'                  // must be a quoted string
  // Pre-built refs — skips BUILD_HAP_REFERENCES + BUILD_DIP_REFERENCE
  // Filenames inferred: bTaeGut7_v0.6.{hap1,hap2,dip}.fa.gz (+.gzi, .fai)
  assemblies_dir  = '/path/to/results/assemblies'
  // Pre-built BAMs — skips WM + BWA mapping and MERGE_HYBRID (per hap)
  mapping_dir     = '/path/to/results/mapping'
  // Pre-built examples — skips make_examples (step 1)
  deepvariant_dir = '/path/to/results/deepvariant'
}
```

### 4. Adjust resources for your cluster (optional)

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

### 5. Run

```sh
# Submit as a job to slurm

# load nextflow and link the nextflow path
module load nextflow
ln -s /data/Phillippy/tools/T2T-Polish/nextflow

cp nextflow/user.config.example user.config ## edit user.config

# submit (default: 2 rounds, output in results/)
sbatch nextflow/run.sh user.config

# If it crashed for some reason, edit user.config to point to the intermediate paths of the results folder
sbatch nextflow/run.sh user.config -resume

# Standard run
nextflow run main.nf -c user.config -ansi-log false \
  -with-trace logs/trace.txt -with-report logs/report.html -overwrite

# Resume after a partial run (requires intact work/ directory)
nextflow run main.nf -c user.config -resume

# Map ONT reads to hap1 and hap2 as well (default: dip only)
nextflow run main.nf -c user.config --ont_map_haps true

# Keep per-read / per-pair intermediate BAMs in the results directory
nextflow run main.nf -c user.config --keep_intermediates true

# Disable SNV candidate collection (mapping + DeepVariant only)
nextflow run main.nf -c user.config --run_snv_candidates false

# Use R9 ONT chemistry for Track B (PEPPER-Margin-DV on dip, per-chromosome scatter)
nextflow run main.nf -c user.config --ont_chemistry r9

# Run only 1 polishing round
nextflow run main.nf -c user.config --polish_rounds 1

# Dry-run — validate params and print all registered process slots (no jobs submitted)
nextflow run main.nf -c user.config -preview

# Stub run — executes the full dataflow graph; stub: blocks replace real commands
# with touch/mkdir stubs so no real work is done.  NOTE: publishDir still fires,
# so stub outputs ARE written to your results directory.  Use a throw-away outdir:
nextflow run main.nf -c user.config -stub --outdir stub_test
```

> **`-preview` shows registered process slots, not jobs that will run.**
> Every process that is wired into the workflow graph appears in the list with
> status `[-  ]` ("no tasks yet"), regardless of whether it will actually
> execute.  In particular:
>
> - `BUILD_HAP_REFERENCES` and `BUILD_DIP_REFERENCE` do **not** appear when
>   `assemblies_dir` is set, because they are guarded by a plain Groovy `if`
>   block that is evaluated at graph-build time — `-preview` can resolve it.
> - `MERYL_REPETITIVE` and `BWA_INDEX` **always** appear even when
>   `assemblies_dir` is set.  Their skip logic is a runtime `channel.branch`
>   that checks whether the output file already exists on disk — `-preview`
>   cannot evaluate file-system state, so both slots are always registered.
> - `MERGE_HYBRID_R1` is always listed because it is registered in
>   `DEEPVARIANT_R1`.  When `mapping_dir` is set it only fires for hap/combo
>   items whose hybrid BAM is absent from that directory; present items bypass
>   it via the per-item join.  `-preview` cannot evaluate this — it has no data.
> - `MAPPING_R1:WINNOWMAP_*` / `BWA_*` do **not** appear in `-preview` when
>   `mapping_dir` is set, because in that case `MAPPING_R1` is never called.
>   `MAPPING_R2:*` (and later rounds) **always** appear because per-round
>   re-entry via `mapping_dir` is only implemented for Round 1.
> - `PEPPER_MARGIN_DV` and `DV_MERGE_CHR_VCFS` are always listed even when
>   `--ont_chemistry r10` is set.  The R9/R10 branch is evaluated at runtime
>   from channel data; `-preview` sees both branches as registered slots.
>
> In short: a process appearing in `-preview` means it *could* run; it does
> not mean it *will* run.  To confirm what actually executed, check
> `nextflow log last -f name,status` after a real run.

**Graceful shutdown on Slurm:**  
When the pipeline head job is submitted via `sbatch`, use `scancel --signal=INT <jobid>` rather than plain `scancel`. This sends SIGINT to Nextflow, which lets currently running jobs finish before exiting (equivalent to pressing Ctrl+C once interactively). Plain `scancel` sends SIGTERM/SIGKILL and may leave orphaned Slurm jobs.

| Command |	Signal sent |	Effect on Nextflow |
| ------- | ----------- | ------------------ |
| scancel --signal=INT <jobid> | SIGINT (Ctrl+C) | Nextflow catches it, waits for running child jobs to finish, then exits cleanly |
| scancel <jobid> (default) | SIGTERM | Nextflow may not handle it gracefully — running child jobs can be orphaned |
| scancel --signal=KILL <jobid> | SIGKILL | Immediate kill, uncatchable — definitely orphans child jobs |


---

## Monitoring and logs

Nextflow writes logs in several places relative to the directory where you run `nextflow run`:

| Path | Contents |
|------|----------|
| `.nextflow.log` | Main engine log — config/startup errors, process submissions, retries. Rotated to `.nextflow.log.1`, `.nextflow.log.2`, … on subsequent runs. |
| `work/<xx>/<hash>/` | One subdirectory per process execution. Contains `.command.sh` (exact script), `.command.log` (stdout + stderr), `.command.out`, `.command.err`, `.exitcode`. |

**Useful run flags:**

```sh
# Plain line-by-line output (better for log files / non-interactive sessions):
nextflow run main.nf -c user.config -ansi-log false

# Write a TSV of task timings, CPU, memory, and I/O:
nextflow run main.nf -c user.config -with-trace trace.txt

# Write an HTML summary report:
nextflow run main.nf -c user.config -with-report report.html

# Draw the workflow DAG:
nextflow run main.nf -c user.config -with-dag dag.svg
```

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

---

## Outputs

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
├── deepvariant/
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

---

## Notes

- **No external submit scripts required.** The original `winnowmap/_submit.sh` and `bwa/_submit_bwa.sh` are fully replaced by native Nextflow processes. Slurm job dependencies, array logic, and `--gres=lscratch` are all handled via `resources.config` resource labels.
- **`-resume` is safe.** Nextflow caches each process by its inputs. Re-running with `-resume` skips any already-completed step, equivalent to the `if [[ -s $out.bam ]]` guards in the original scripts.
- **`module load` calls** inside process scripts rely on the Biowulf (NIH HPC) module system. Adjust them if running on a different cluster.
- **Cluster portability.** Copy `resources.config`, change `executor`, queue names, and `clusterOptions` to suit your scheduler, then load your copy via `includeConfig` in `user.config` or with `-c my_resources.config` on the command line.
- **DeepVariant always runs.** The hybrid-merge + variant-calling sub-workflow fires automatically once the relevant BAMs for each haplotype are ready. No flag is needed to enable it.
- **Two parallel DV tracks.** The pipeline runs two independent tracks simultaneously: (A) hybrid HiFi + short-read calling on all three haps, and (B) ONT-only calling on the dip reference.
- **Multiple short-read platforms.** When both `illumina` and `element` are in `--platforms`, `MERGE_HYBRID` merges all active short-read BAMs together with the HiFi BAM into a single hybrid BAM before calling variants — so DeepVariant still produces exactly one VCF per haplotype regardless of how many short-read platforms were used. The merged BAM combo tag (e.g. `HYBRID_PACBIO_ILLUMINA_ELEMENT`) is reflected in the output directory name but does not affect downstream processing.
- **ONT chemistry selection.** Track B uses DeepVariant's `ONT_R104` model by default (`--ont_chemistry r10`). For R9 data set `--ont_chemistry r9`, which switches to PEPPER-Margin-DeepVariant (module `pepper_deepvariant/0.8`) with a per-chromosome GPU scatter. **If you have mixed R9 and R10 data, use `--ont_chemistry r10`** — the R10 model is preferred, it runs as a single pipeline step, and it produces a gVCF. R9 mode does not emit a gVCF.
- **Per-hap MQ thresholds.** DeepVariant applies MQ 5 to hap1/hap2 and MQ 0 to dip — stricter for haplotypes, permissive for the diploid reference where multi-mapper signal is informative. Override with `--dv_mq_hap` and `--dv_mq_dip` respectively.
- **Assembly naming.** Set `params.asm_name` to your assembly prefix (e.g. `bTaeGut7`) and `params.asm_ver` to the starting version tag (e.g. `v0.1`). All output paths use `{asm_name}_{asm_ver}` as a prefix. Polished assemblies are written with auto-incremented version tags (`v0.2`, `v0.3`, …).
- **SNV candidates always run.** Pass `--run_snv_candidates false` to skip the SNV candidate collection sub-workflow. When enabled (the default), requires `params.hybrid_meryl` (path to a pre-built hybrid read k-mer meryl DB, e.g. `hybrid.k31.meryl`) and `params.merfin_peak` (integer peak coverage from `meryl histogram` / GenomeScope).
- **Two Merfin resource tiers.** The bcftools reheader / filter / isec steps (`SNV_FILTER_INTERSECT`) run on a 12 CPU / 12 GB node. Both Merfin runs (`-strict` on the hybrid→dip VCF and `-loose` on the final candidate set) run on a 12 CPU / 120 GB node (`norm_snv_merfin`) — the meryl k-mer DB lookup is the memory bottleneck. Adjust these labels in `resources.config` for your cluster.
- **SNV candidate pipeline mirrors `variant_call/snv_candidates.sh`.** The bcftools logic is identical to the original script: reheader → PASS filter → hap concat → isec consensus errors → GT/GQ/AF filters → merfin-strict → final concat → GT→1/1 → merfin-loose. The only difference is that the seqmer `meryl count k=31` step is folded into the `SNV_MERFIN` process so it shares the same high-memory allocation as the k-mer lookups.
- **SNV consensus is applied to dip only.** `SNV_APPLY_CONSENSUS` applies the merfin-loose VCF to the dip reference with `bcftools consensus -H 1`, producing a single polished `<ver_to>.dip.fa`. hap1 and hap2 are then extracted from the polished dip in `PREPARE_NEXT_ROUND` using `samtools faidx -r <names>` — where the sequence names come from the current-round hap1/hap2 FAI files. This means MT, EBV, rDNA, and any other shared sequences appear exactly once in the dip and are distributed correctly to each haplotype without duplication.
- **Multi-round polishing.** After each round `PREPARE_NEXT_ROUND` produces `[hap, fa, fai]` for hap1, hap2, and dip, which feed directly into `BUILD_REFS_FROM_FILES` → `MAPPING` → `DEEPVARIANT` → `SNV_CANDIDATES` for round N+1. Version tags are auto-bumped each round (`v0.1 → v0.2 → v0.3`, etc.). Set `params.polish_rounds` (default `2`) to control how many rounds run; the valid range is 1–5.
- **ONT DV VCF feeds SNV candidates.** The `ont_to_dip` VCF (from Track B) is used by `SNV_CANDIDATES` to corroborate hap-het sites with ONT homozygous support. If Track B produced no VCF, the `ont_dip_vcf` sub-channel will be empty and SNV candidates will silently skip. `PEPPER_MARGIN_DV` (R9 mode) requests four K80 GPUs per chromosome (`--gres=gpu:k80:4`). Adjust the `norm_dv_call_variants` and `norm_pepper_margin_dv` labels in `resources.config` for your cluster's GPU partition.
- **call_variants retry.** The original `step2_with_minqual.sh` had an explicit retry loop for zero-size GPU output files. This is replaced by `maxRetries = 1` on the `norm_dv_call_variants` label in `resources.config` — Nextflow will automatically re-submit the job once if the task exits with a non-zero code.
- **Paired-end globs** for `read_glob_illumina` and `read_glob_element` must be compatible with Nextflow's `Channel.fromFilePairs()`, i.e. contain a `{1,2}` or `R{1,2}` wildcard so R1 and R2 are grouped together. At least one of the two must be set when `illumina` or `element` appears in `--platforms` or `--dv_short_platforms` — the pipeline will error at startup if neither is provided. To run without any short-read data, exclude both with `--platforms hifi,ont --dv_short_platforms ''`.
- **Startup validation.** The pipeline checks three conditions before submitting any jobs and errors immediately if any fail: (1) `hifi` must be in `--platforms` and `params.read_glob_hifi` must be set; (2) `ont` must be in `--platforms` and `params.read_glob_ont` must be set; (3) when `illumina` or `element` is active (via `--platforms` or `--dv_short_platforms`), at least one of `read_glob_illumina` / `read_glob_element` must be non-null.

---

## Disk space

### Where files are written

Nextflow writes **two separate copies** of every output file:

| Location | Purpose | Read by `-resume`? |
|---|---|---|
| `work/<hash>/` | Canonical task outputs; source of truth for the cache | ✅ yes |
| `results/` (publishDir) | Human-readable final outputs; hard-links or copies from `work/` | ❌ no |

`results/` is a **one-way publication sink**.  Nextflow never reads it back.
On `-resume`, the cache is checked against `work/` only — so if `work/` is
deleted, every task re-runs from scratch even if `results/` is fully intact.

### Intermediate BAMs and `work/`

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

> **Tip — selective cleaning with `-but-less-recent`:**
> If you have a long run that you want to keep resumable but want to free
> space from earlier failed attempts, use:
> ```sh
> nextflow clean -f -but-less-recent <run-name>
> ```
> This keeps the most recent run's `work/` intact and deletes all older ones.

### Recommended workflow for large assemblies

1. Run the pipeline to completion.
2. Verify `results/` looks correct.
3. Archive `results/` to long-term storage (e.g. `rsync` to a data store).
4. Run `nextflow clean -f <run-name>` to delete `work/`.
5. Delete intermediate BAMs from `results/` if they were published
   (`--keep_intermediates true`) and are no longer needed.


