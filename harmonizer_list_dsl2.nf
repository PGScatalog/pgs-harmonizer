#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { HmPOS as HM_POS_37 } from './nf_modules/hmpos'
include { HmPOS as HM_POS_38 } from './nf_modules/hmpos'
include { HmVCF as HM_VCF_37 } from './nf_modules/hmvcf'
include { HmVCF as HM_VCF_38 } from './nf_modules/hmvcf'
include { Finalise as HM_FINALISE_37 } from './nf_modules/finalise'
include { Finalise as HM_FINALISE_38 } from './nf_modules/finalise'

workflow hm_37 {
    take:
        pgs_ids
    main:
        genebuild='37'
        // Data analysis
        HM_POS_37(pgs_ids, genebuild)
        HM_VCF_37(HM_POS_37.out, genebuild)

        // Data post-processing
        HM_FINALISE_37(HM_VCF_37.out, genebuild)
}

workflow hm_38 {
    take:
        pgs_ids
    main:
        genebuild='38'
        // Data analysis
        HM_POS_38(pgs_ids, genebuild)
        HM_VCF_38(HM_POS_38.out, genebuild)

        // Data post-processing
        HM_FINALISE_38(HM_VCF_38.out, genebuild)
}

workflow {
    // Get list of PGS IDs
    pgs_ids_list = channel.from(params.pgs)

    // Run sub-workflows
    hm_37(pgs_ids_list)
    hm_38(pgs_ids_list)
}
