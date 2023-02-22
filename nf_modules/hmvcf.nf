process HmVCF {
    label 'retry_increasing_big_mem'
    input:
      val pgs_id
      val genomebuild
    output:
      val pgs_id

    script:
    """
    python $params.loc_scripts/Harmonize.py HmVCF -loc_files $params.loc_hmoutput/$genomebuild -loc_hmoutput $params.loc_hmoutput_vcf/$genomebuild -loc_vcfs $params.loc_vcfs/$genomebuild --gzip $pgs_id GRCh$genomebuild
    """
}