# T2T-Polish Nextflow — Dev Session Notes
_Last updated: 2026-03-26_

This file records bugs discovered and fixed during active development/debugging
sessions so that context is not lost between conversations.

---

## Session: SNV candidate pipeline bring-up (Mar 2026)

### Pipeline stage reached
The pipeline ran successfully through DeepVariant (all rounds), then entered
the SNV candidate / polishing stage (`SNV_CANDIDATES_R1`).  The following bugs
were found and fixed sequentially.

---

### Bug 1 — `bcftools index` produces `.csi`, not `.tbi`
**File:** `modules/snv_candidates.nf` — `SNV_FILTER_INTERSECT` and `SNV_MERFIN`

**Symptom:**
```
cp: cannot stat 'hybrid_to_dip.PASS.vcf.gz.tbi': No such file or directory
```

**Root cause:** `bcftools index` (without flags) writes a `.csi` index.
All `output:` declarations and the final `cp` commands expected `.tbi` (tabix format).

**Fix:** Added `-t` to every `bcftools index` call in both processes:
```bash
bcftools index -t <file>.vcf.gz
# and inside filter_PASS():
bcftools index -t -f $OUT
```

---

### Bug 2 — `merfin: command not found`
**File:** `modules/snv_candidates.nf` — `SNV_MERFIN`; `nextflow.config`; `user.config.example`

**Symptom:**
```
.command.sh: line 15: merfin: command not found
```

**Root cause:** `params.merfin` defaulted to `'merfin'` (bare name), but merfin
is not on `$PATH` on Biowulf — it lives at `$tools/merfin/1.1/bin/merfin`.

**Fix:**
- Added `params.tools = null` to `nextflow.config` (root of the shared tools tree).
- Changed `params.merfin` default to `null` in `nextflow.config`.
- Resolution chain in `SNV_MERFIN` script block:
  ```groovy
  def merfin = params.merfin ?: (params.tools ? "${params.tools}/merfin/1.1/bin/merfin" : 'merfin')
  ```
- Added `tools = '/data/Phillippy/tools'` to `user.config.example`.
- `k8` path in `user.config.example` updated to use `"${params.tools}/k8/k8"`.

**User action required:** Set `tools = '/path/to/tools'` in `user.config`.

---

### Bug 3 — `bgzip --index` fails when writing to stdout (SNV_APPLY_CONSENSUS)
**File:** `modules/snv_candidates.nf` — `SNV_APPLY_CONSENSUS`

**Symptom:**
```
[bgzip] Index file name expected when writing to stdout
```

**Root cause:** `bgzip --index` requires a named output file; it cannot write
the index when the compressed data is piped to stdout via `>`.

**Fix:** Added `-I <filename>.gzi` to name the index file explicitly:
```bash
| bgzip --index -I ${ver_to}.dip.fa.gz.gzi -@ ${task.cpus} > ${ver_to}.dip.fa.gz
```

---

### Bug 4 — Same `bgzip --index` stdout error in `PREPARE_NEXT_ROUND`
**File:** `modules/snv_candidates.nf` — `PREPARE_NEXT_ROUND`

**Symptom:** Same as Bug 3, for hap1 and hap2 extractions.

**Fix:** Added `-I` for both:
```bash
| bgzip --index -I ${ver_to}.hap1.fa.gz.gzi -@ ${task.cpus} > ${ver_to}.hap1.fa.gz
| bgzip --index -I ${ver_to}.hap2.fa.gz.gzi -@ ${task.cpus} > ${ver_to}.hap2.fa.gz
```

---

### Bug 5 — `cp: 'v0.7.dip.fa.gz' and 'v0.7.dip.fa.gz' are the same file`
**File:** `modules/snv_candidates.nf` — `PREPARE_NEXT_ROUND`
**File:** `workflows/snv_candidates.nf`

**Root cause:** `PREPARE_NEXT_ROUND` received the dip FA from
`SNV_APPLY_CONSENSUS` and tried to `cp` it to the same name in the same work
directory. The dip is already published by `SNV_APPLY_CONSENSUS`; there is no
need to re-copy it.

**Fix (module):**
- Removed `dip_ref` output tuple from `PREPARE_NEXT_ROUND`.
- Removed the `cp ${dip_fa_gz} …` block from the script.
- Removed `dip` entry from the stub.
- Kept `path(dip_fai)` in `input:` — it must still be staged so that
  `samtools faidx` can find the index alongside `dip_fa_gz`.

**Fix (workflow):**
- `dip` ref in `next_refs_ch` now comes directly from
  `SNV_APPLY_CONSENSUS.out.polished_dip_fa` mapped to `tuple('dip', fa_gz, gzi, fai)`.
- `PREPARE_NEXT_ROUND(...)` call still passes `dip_fai` as the fourth argument
  (required for `samtools faidx` index staging).

```groovy
// workflows/snv_candidates.nf
next_refs_ch = PREPARE_NEXT_ROUND.out.hap1_ref
    .mix(PREPARE_NEXT_ROUND.out.hap2_ref)
    .mix(SNV_APPLY_CONSENSUS.out.polished_dip_fa
        .map { ver, fa_gz, gzi, fai -> tuple('dip', fa_gz, gzi, fai) })
```

---

## Current state (as of 2026-03-26)

| Process | Status |
|---|---|
| `SNV_FILTER_INTERSECT` | ✅ Fixed (`.tbi` indexing) |
| `SNV_MERFIN` | ✅ Fixed (merfin path via `params.tools`) |
| `SNV_APPLY_CONSENSUS` | ✅ Fixed (`bgzip -I`) |
| `PREPARE_NEXT_ROUND` | ✅ Fixed (`bgzip -I`; removed redundant dip copy) |

Pipeline has not yet been confirmed to complete end-to-end.
Next expected stage after `PREPARE_NEXT_ROUND`: Round 2 mapping → DeepVariant → SNV candidates.

---

## Key configuration notes

- Set `tools = '/data/Phillippy/tools'` in `user.config` to resolve `merfin` and `k8`.
- `merfin` can be overridden explicitly: `merfin = '/path/to/merfin'`.
- `bgzip --index` always needs `-I <file>.gzi` when output goes to stdout.
- `bcftools index` always needs `-t` to produce `.tbi` (tabix) instead of `.csi`.
- `PREPARE_NEXT_ROUND` requires `dip_fai` staged as input so `samtools faidx`
  can find the index — even though `dip_fai` is not referenced by name in the script.
