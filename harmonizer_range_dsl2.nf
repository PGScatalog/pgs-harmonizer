#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { HmPOS as HM_POS_37 } from './nf_modules/hmpos'
include { HmPOS as HM_POS_38 } from './nf_modules/hmpos'
include { HmVCF as HM_VCF_37 } from './nf_modules/hmvcf'
include { HmVCF as HM_VCF_38 } from './nf_modules/hmvcf'
include { Finalise as HM_FINALISE_37 } from './nf_modules/finalise'
include { Finalise as HM_FINALISE_38 } from './nf_modules/finalise'


process get_pgs_ids_list {
  input:
    val from_id
    val to_id

  output:
    path params.pgs_ids_file, emit: pgs_ids_list_file

  script:
  """
  python $params.loc_scripts/pgs_harmonizer/get_pgs_ids_list.py \
    --num_from $from_id \
    --num_to $to_id \
    --output $params.pgs_ids_file \
    --rest_server $params.rest_api_url \
    --loc_files $params.loc_files
  """
}


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
    // Channels
    pgs_from = Channel.from(params.pgs_num_from)
    pgs_to = Channel.from(params.pgs_num_to)

    // Get list of PGS IDs
    get_pgs_ids_list(pgs_from,pgs_to)
    pgs_ids_list = get_pgs_ids_list.out.pgs_ids_list_file.splitText{it.strip()}

    // Run sub-workflows
    hm_37(pgs_ids_list)
    hm_38(pgs_ids_list)
}
