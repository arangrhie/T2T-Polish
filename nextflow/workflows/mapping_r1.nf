/* Round-1 — per-round process aliases for NF 25.x uniqueness constraint */
include { WINNOWMAP_MAP    as WINNOWMAP_MAP_R1         } from '../modules/winnowmap'
include { WINNOWMAP_MERGE  as WINNOWMAP_MERGE_HIFI_R1  } from '../modules/winnowmap'
include { WINNOWMAP_MERGE  as WINNOWMAP_MERGE_ONT_R1   } from '../modules/winnowmap'
include { WINNOWMAP_FILTER as WINNOWMAP_FILTER_R1      } from '../modules/winnowmap'
include { SAM2PAF          as SAM2PAF_R1               } from '../modules/winnowmap'
include { BWA_MAP          as BWA_MAP_R1               } from '../modules/bwa'
include { BWA_MERGE        as BWA_MERGE_R1             } from '../modules/bwa'

workflow MAPPING_R1 {
    take:
    wm_refs   // [ hap, ver_from, ref_fa_gz, ref_fai, repetitive_k15.txt ]
    bwa_refs  // [ hap, ver_from, ref_fa_gz, ref_fai, amb, ann, bwt, pac, sa ]

    main:
    // ---- Winnowmap (HiFi + ONT) ----
    def hifiGlob = params.platforms.contains('hifi') ? params.read_glob_hifi : null
    def ontGlob  = params.platforms.contains('ont')  ? params.read_glob_ont  : null

    def hifiReads = hifiGlob ? Channel.fromPath(hifiGlob).map { f -> tuple('hifi', f) }
                             : Channel.empty()
    def ontReads  = ontGlob  ? Channel.fromPath(ontGlob ).map { f -> tuple('ont',  f) }
                             : Channel.empty()

    // wm_refs is a hot channel — multiMap to avoid consuming it twice
    def wm = wm_refs.multiMap { it -> hifi: it; ont: it }
    def hifiInput = wm.hifi.combine(hifiReads)
    def ontBase   = params.ont_map_haps ? wm.ont
                                        : wm.ont.filter { hap, ver, ref, fai, rep -> hap == 'dip' }
    def ontInput  = ontBase.combine(ontReads)

    def mapped   = WINNOWMAP_MAP_R1( hifiInput.mix(ontInput) )
    // Branch by platform before groupTuple so each platform's sub-channel
    // closes independently — FILTER/SAM2PAF start per (hap,platform) as
    // soon as that platform's merges are done, without waiting for the other.
    def mappedByPlat = mapped.branch { hap, ver, plat, bam, bai ->
        hifi: plat == 'hifi'
        ont:  true
    }
    def mergedHifi = WINNOWMAP_MERGE_HIFI_R1( mappedByPlat.hifi.groupTuple(by: [0, 1, 2]) )
    def mergedOnt  = WINNOWMAP_MERGE_ONT_R1(  mappedByPlat.ont .groupTuple(by: [0, 1, 2]) )
    def filtered = WINNOWMAP_FILTER_R1( mergedHifi.mix(mergedOnt) )
    def filteredPafs = SAM2PAF_R1(filtered)

    // ---- BWA (Illumina + Element) ----
    def illGlob  = params.platforms.contains('illumina') ? params.read_glob_illumina : null
    def elemGlob = params.platforms.contains('element')  ? params.read_glob_element  : null

    def illReads  = illGlob  ? Channel.fromFilePairs(illGlob ).map { sid, f -> tuple('illumina', f[0], f[1]) }
                             : Channel.empty()
    def elemReads = elemGlob ? Channel.fromFilePairs(elemGlob).map { sid, f -> tuple('element',  f[0], f[1]) }
                             : Channel.empty()

    def allBwaReads = illReads.mix(elemReads)
    def bwaMapped = BWA_MAP_R1( bwa_refs.combine(allBwaReads) )
    def bwaMerged = BWA_MERGE_R1( bwaMapped.groupTuple(by: [0, 1, 2]) )

    emit:
    wm_pri_bams  = filtered         // [ hap, ver_from, platform, bam, bai ]
    wm_pri_pafs  = filteredPafs     // [ hap, ver_from, platform, paf ]
    bwa_mrg_bams = bwaMerged        // [ hap, ver_from, platform, bam, csi ]
}
