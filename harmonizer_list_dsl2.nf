#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { HmPOS as HM_POS } from './nf_modules/hmpos'
include { HmVCF as HM_VCF } from './nf_modules/hmvcf'
include { Finalise as HM_FINALISE } from './nf_modules/finalise'


workflow {
    // Get list of PGS IDs
    pgs_ids_list = channel.from(params.pgs)

    // Data analysis
    HM_POS(pgs_ids_list)
    HM_VCF(HM_POS.out)

    // Data post-processing
    HM_FINALISE(HM_VCF.out)
}
