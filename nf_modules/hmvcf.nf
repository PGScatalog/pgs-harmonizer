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