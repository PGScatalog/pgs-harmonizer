#!/usr/bin/env nextflow

nextflow.enable.dsl=2


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
    --rest_server $params.rest_api_url
  """
}


process HmPOS {
    label 'retry_increasing_mem'
    input:
      val pgs_id

    output:
      val pgs_id

    script:
    """
    python $params.loc_scripts/Harmonize.py HmPOS -loc_files $params.loc_files -loc_hmoutput $params.loc_hmoutput -var2location $params.loc_var2location --gzip $pgs_id $params.genebuild_grc
    """
}

process HmVCF {
    label 'retry_increasing_mem'
    input:
      val pgs_id

    output:
      val pgs_id

    script:
    """
    python $params.loc_scripts/Harmonize.py HmVCF -loc_files $params.loc_hmoutput -loc_hmoutput $params.loc_hmoutput_vcf -loc_vcfs $params.loc_vcfs --gzip $pgs_id $params.genebuild_grc
    """
}


process Finalise {
    input:
      val pgs_id

    script:
    """
    python $params.loc_scripts/pgs_harmonizer/finalise_harmonized_file.py --score_id $pgs_id --input_dir $params.loc_hmoutput --staged_dir $params.loc_staged --genebuild $params.genebuild --sqlite_file $params.hm_version_sqlite_file_path
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
    HmPOS(pgs_ids_list)
    HmVCF(HmPOS.out)

    // Data post-processing
    Finalise(HmVCF.out)
}
