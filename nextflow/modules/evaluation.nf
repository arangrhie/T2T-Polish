/*
 * modules/evaluation.nf
 *
 * Consolidated process for post-assembly evaluation: Merqury QV, pattern analysis, and issue detection.
 * Replaces: coverage/submit_evaluation.sh workflow
 * 
 * All steps run in a single process to minimize scheduling overhead (steps take ~seconds each).
 */

/*
 * Unified evaluation process: Merqury → Pattern → Issues → Intersection
 * Runs all quality metrics and issue detection in a single job.
 */
process FULL_EVALUATION {
    label 'evaluation'
    tag "${name}"
    publishDir "${params.outdir}", mode: 'link', overwrite: false

    input:
    tuple val(name), val(ver), 
          path(asm_fa_gz), path(asm_fai), 
          path(meryl_db),
          path(hifi_paf), path(ont_paf)

    output:
    tuple val(name),
          path("merqury/${name}_only.bed"),
          path("merqury/${name}_only.wig"),
          path("merqury/${name}.qv"),
          path("merqury/${name}.${name}.dip.qv"),
          path("pattern/${name}.bed"),
          path("pattern/${name}.telo.bed"),
          path("pattern/${name}.exclude.bed"),
          path("pattern/${name}.error.bed"),
          path("pattern/${name}/${name}.microsatellite.{GC,AT,GA,TC}.128.bw"),
          path("issues/hifi/${name}.issues.bed"),
          path("issues/hifi/${name}.issues.fm.bed"),
          path("issues/hifi/${name}.cov.wig"),
          path("issues/hifi/${name}.cov.bw"),
          path("issues/hifi/${name}.issues.bb"),
          path("issues/hifi/cov.mean.txt"),
          path("issues/hifi/cov.med.txt"),
          path("issues/hifi/cov.sd.txt"),
          path("issues/hifi/cov.theta.txt"),
          path("issues/hifi/high_cutoff.txt"),
          path("issues/hifi/low_cutoff.txt"),
          path("issues/ont/${name}.issues.bed"),
          path("issues/ont/${name}.issues.fm.bed"),
          path("issues/ont/${name}.cov.wig"),
          path("issues/ont/${name}.cov.bw"),
          path("issues/ont/${name}.issues.bb"),
          path("issues/ont/cov.mean.txt"),
          path("issues/ont/cov.med.txt"),
          path("issues/ont/cov.sd.txt"),
          path("issues/ont/cov.theta.txt"),
          path("issues/ont/high_cutoff.txt"),
          path("issues/ont/low_cutoff.txt"),
          path("issues/${name}.issues.bed"),
          path("issues/${name}.issues.bb")

    script:
    def tools_path = params.tools ?: '/unknown/tools'
    """
    set -euo pipefail
    module load meryl samtools seqtk bedtools ucsc || true
    
    export tools=${tools_path}
    
    # ---- Step 1: Merqury QV Evaluation ----
    echo "Step 1: Running Merqury QV evaluation..."
    mkdir -p merqury
    cd merqury
    if [[ ! -s ${name}_only.bed ]]; then
        \$MERQURY/eval/qv.sh ../${meryl_db} ../${asm_fa_gz} ${name}
    fi
    cd ../

    # ---- Step 2: Initialize Pattern Files ----
    echo "Step 2: Initializing pattern files..."
    mkdir -p pattern
    cd pattern
    
    # Decompress if needed
    if [[ "${asm_fa_gz}" == *.gz ]]; then
        pigz -dc ../${asm_fa_gz} > ${name}.fa
        samtools faidx -@12 ${name}.fa
    else
        ln -sf ../${asm_fa_gz} ${name}.fa
        ln -sf ../${asm_fai} ${name}.fa.fai
    fi

    # Generate bed files
    awk '{print \$1"\\t0\\t"\$2}' ${name}.fa.fai > ${name}.bed
    seqtk telo -d 50000 ${name}.fa > ${name}.telo.bed
    seqtk gap -l0 ${name}.fa > ${name}.exclude.bed

    # Generate pattern/microsatellite files
    \${tools}/T2T-Polish/pattern/microsatellites.sh ${name}.fa
    cd ../

    # ---- Step 3: Merge Merqury Error Regions ----
    echo "Step 3: Merging Merqury error regions with pattern..."
    cd pattern
    bedtools merge -i ../merqury/${name}_only.bed > ${name}.error.bed
    cd ../

    # ---- Step 4: Collect Issues from PAF Files ----
    echo "Step 4: Collecting issues from HiFi and ONT mappings..."
    mkdir -p issues/hifi issues/ont
    
    for platform in hifi ont; do
        paf_file=\${platform}_paf
        if [[ \$platform == "hifi" ]]; then
            paf_file=../../${hifi_paf}
        else
            paf_file=../../${ont_paf}
        fi
        
        echo "  Processing \$platform alignment..."
        cd issues/\$platform
        
        # Link pattern directory
        ln -sf ../../pattern pattern
        
        # Verify pattern files exist
        if [[ ! -d pattern/${name} || ! -s pattern/${name}.exclude.bed ]]; then
            echo "Missing pattern files for \$platform. Exiting."
            exit 1
        fi

        # Link PAF
        paf_base=${name}
        ln -sf \$paf_file ${name}.paf || true

        # Collect coverage statistics
        \${tools}/T2T-Polish/coverage/collect_summary.sh ${name}.paf \${platform}
        \${tools}/T2T-Polish/coverage/coverage_stat.sh ${name}.cov.wig pattern/${name}.exclude.bed

        # Collect low support regions
        \${tools}/T2T-Polish/coverage/low_support.sh ${name}.paf ${name} \${platform} pattern

        # Fix mapping limits and convert to bigBed
        awk -v OFS='\\t' -v FS='\\t' 'NR==FNR { mapping[\$1] = \$3; next } \\
          {if (\$3 > mapping[\$1]) { \$3=mapping[\$1]; \$8=mapping[\$1]; } print}' \\
          pattern/${name}.bed ${name}.issues.bed > ${name}.issues.fm.bed

        wigToBigWig -clip ${name}.cov.wig pattern/${name}.fa.fai ${name}.cov.bw
        bedToBigBed -type=bed9 ${name}.issues.fm.bed pattern/${name}.fa.fai ${name}.issues.bb

        cd ../../
    done

    # ---- Step 5: Intersect HiFi and ONT Issues ----
    echo "Step 5: Intersecting HiFi and ONT issues..."
    cd issues
    
    # Find regions supported by both HiFi and ONT
    bedtools intersect -u -a hifi/${name}.issues.fm.bed -b ont/${name}.issues.fm.bed > ${name}.issues.bed
    bedToBigBed -type=bed9 ${name}.issues.bed ../pattern/${name}.fa.fai ${name}.issues.bb
    
    cd ../
    echo "Evaluation complete!"
    """

    stub:
    """
    mkdir -p merqury pattern issues/hifi issues/ont
    mkdir -p pattern/${name}
    touch merqury/${name}_only.bed merqury/${name}_only.wig merqury/${name}.qv merqury/${name}.${name}.dip.qv
    touch pattern/${name}.{bed,telo.bed,exclude.bed,error.bed}
    touch pattern/${name}/${name}.microsatellite.{GC,AT,GA,TC}.128.bw
    touch issues/hifi/${name}.{issues.bed,issues.fm.bed,cov.wig,cov.bw,issues.bb}
    touch issues/hifi/cov.{mean,med,sd,theta}.txt issues/hifi/{high,low}_cutoff.txt
    touch issues/ont/${name}.{issues.bed,issues.fm.bed,cov.wig,cov.bw,issues.bb}
    touch issues/ont/cov.{mean,med,sd,theta}.txt issues/ont/{high,low}_cutoff.txt
    touch issues/${name}.{issues.bed,issues.bb}
    """
}
