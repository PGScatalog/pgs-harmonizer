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