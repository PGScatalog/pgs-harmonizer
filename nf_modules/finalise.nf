process Finalise {
    label 'retry_increasing_mem'
    input:
      val pgs_id
      val genomebuild
    script:
    """
    python $params.loc_scripts/pgs_harmonizer/finalise_harmonized_file.py --score_id $pgs_id --input_dir $params.loc_hmoutput_vcf/$genomebuild --staged_dir $params.loc_staged --genomebuild $genomebuild --sqlite_file $params.hm_version_sqlite_file_path
    """
}