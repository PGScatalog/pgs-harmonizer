#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { HmPOS as HM_POS } from './nf_modules/hmpos'
include { HmVCF as HM_VCF } from './nf_modules/hmvcf'
include { Finalise as HM_FINALISE } from './nf_modules/finalise'


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


workflow {
    // Channels
    pgs_from = Channel.from(params.pgs_num_from)
    pgs_to = Channel.from(params.pgs_num_to)

    // Get list of PGS IDs
    get_pgs_ids_list(pgs_from,pgs_to)
    pgs_ids_list = get_pgs_ids_list.out.pgs_ids_list_file.splitText{it.strip()}

    // Data analysis
    HM_POS(pgs_ids_list)
    HM_VCF(HM_POS.out)

    // Data post-processing
    HM_FINALISE(HM_VCF.out)
}
