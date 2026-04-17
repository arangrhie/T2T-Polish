/* Round-4 — per-round process aliases for NF 25.x uniqueness constraint */
include { WINNOWMAP_MAP    as WINNOWMAP_MAP_R4         } from '../modules/winnowmap'
include { WINNOWMAP_MERGE  as WINNOWMAP_MERGE_HIFI_R4  } from '../modules/winnowmap'
include { WINNOWMAP_MERGE  as WINNOWMAP_MERGE_ONT_R4   } from '../modules/winnowmap'
include { WINNOWMAP_FILTER as WINNOWMAP_FILTER_R4      } from '../modules/winnowmap'
include { SAM2PAF          as SAM2PAF_R4               } from '../modules/winnowmap'
include { BWA_MAP          as BWA_MAP_R4               } from '../modules/bwa'
include { BWA_MERGE        as BWA_MERGE_R4             } from '../modules/bwa'

workflow MAPPING_R4 {
    take:
    wm_refs   // [ hap, ver_from, ref_fa_gz, ref_fai, repetitive_k15.txt ]
    bwa_refs  // [ hap, ver_from, ref_fa_gz, ref_fai, amb, ann, bwt, pac, sa ]

    main:
    // ---- Winnowmap (HiFi + ONT) ----
    // Build unified reads channel; skip any platform whose glob is null.
    def hifiGlob = params.platforms.contains('hifi') ? params.read_glob_hifi : null
    def ontGlob  = params.platforms.contains('ont')  ? params.read_glob_ont  : null

    def hifiReads = hifiGlob ? Channel.fromPath(hifiGlob).map { f -> tuple('hifi', f) }
                             : Channel.empty()
    def ontReads  = ontGlob  ? Channel.fromPath(ontGlob ).map { f -> tuple('ont',  f) }
                             : Channel.empty()

    def nHifi = hifiGlob ? file(hifiGlob).size() : 0
    def nOnt  = ontGlob  ? file(ontGlob ).size() : 0

    def wm = wm_refs.multiMap { it -> hifi: it; ont: it }
    def hifiInput = wm.hifi.combine(hifiReads)
    def ontBase   = params.ont_map_haps ? wm.ont
                                        : wm.ont.filter { hap, ver, ref, fai, rep -> hap == 'dip' }
    def ontInput  = ontBase.combine(ontReads)

    def mapped   = WINNOWMAP_MAP_R4( hifiInput.mix(ontInput) )
    def mappedByPlat = mapped.branch { hap, ver, plat, bam, bai ->
        hifi: plat == 'hifi'
        ont:  true
    }
    def mergedHifi = WINNOWMAP_MERGE_HIFI_R4( mappedByPlat.hifi.groupTuple(by: [0, 1, 2], size: nHifi) )
    def mergedOnt  = WINNOWMAP_MERGE_ONT_R4(  mappedByPlat.ont .groupTuple(by: [0, 1, 2], size: nOnt ) )
    def filtered = WINNOWMAP_FILTER_R4( mergedHifi.mix(mergedOnt) )
    SAM2PAF_R4(filtered)

    // ---- BWA (Illumina + Element) ----
    def illGlob  = params.platforms.contains('illumina') ? params.read_glob_illumina : null
    def elemGlob = params.platforms.contains('element')  ? params.read_glob_element  : null

    def illReads  = illGlob  ? Channel.fromFilePairs(illGlob ).map { sid, f -> tuple('illumina', f[0], f[1]) }
                             : Channel.empty()
    def elemReads = elemGlob ? Channel.fromFilePairs(elemGlob).map { sid, f -> tuple('element',  f[0], f[1]) }
                             : Channel.empty()

    def allBwaReads = illReads.mix(elemReads)
    def bwaMapped = BWA_MAP_R4( bwa_refs.combine(allBwaReads) )
    def bwaMerged = BWA_MERGE_R4( bwaMapped.groupTuple(by: [0, 1, 2]) )

    emit:
    wm_pri_bams  = filtered         // [ hap, ver_from, platform, bam, bai ]
    bwa_mrg_bams = bwaMerged        // [ hap, ver_from, platform, bam, csi ]
}
