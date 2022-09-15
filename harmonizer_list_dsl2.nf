#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process HmPOS {
    label 'retry_increasing_mem'
    input:
      val pgs_id

    output:
      val pgs_id

    script:
    """
    python $params.loc_scripts/Harmonize.py HmPOS -loc_files $params.loc_files -loc_hmoutput $params.loc_hmoutput -var2location $params.loc_var2location --gzip $pgs_id $params.genomebuild_grch
    """
}

process HmVCF {
    label 'retry_increasing_big_mem'
    input:
      val pgs_id

    output:
      val pgs_id

    script:
    """
    python $params.loc_scripts/Harmonize.py HmVCF -loc_files $params.loc_hmoutput -loc_hmoutput $params.loc_hmoutput_vcf -loc_vcfs $params.loc_vcfs --gzip $pgs_id $params.genomebuild_grch
    """
}


process Finalise {
    label 'retry_increasing_mem'
    input:
      val pgs_id

    script:
    """
    python $params.loc_scripts/pgs_harmonizer/finalise_harmonized_file.py --score_id $pgs_id --input_dir $params.loc_hmoutput_vcf --staged_dir $params.loc_staged --genomebuild $params.genomebuild --sqlite_file $params.hm_version_sqlite_file_path
    """
}

workflow {
    // Get list of PGS IDs
    pgs_ids_list = channel.from(params.pgs)

    // Data analysis
    HmPOS(pgs_ids_list)
    HmVCF(HmPOS.out)

    // Data post-processing
    Finalise(HmVCF.out)
}
