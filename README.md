# pgs-harmonizer
A pipeline to format and harmonize [Polygenic Score (PGS) Catalog Scoring Files](http://www.pgscatalog.org/downloads/#dl_ftp) 
within and between different genome builds. 

_External users NB: the pipeline will run very slowly if only rsIDs are given, because it looks up each rsID in 200 
variant batches using the Ensembl API. There are ways to speed this up by looking up the total set of rsIDs once using the [pgs_variants_coords](https://github.com/PGScatalog/pgs_variants_coords) tool which uses Ensembl VCF files._

---


## Requirements
* Python packages: requests, pandas, pyliftover, tqdm, cyvcf2
  ```pip install requirements.txt```
* PGS Catalog rsID -> position mappings (see [pgs_variants_coords](https://github.com/PGScatalog/pgs_variants_coords) to generate the variants position knowledge base).
* VCFs from Ensembl (in `/map/vcf_ref/`) or a specific cohort (in `/map/cohort_ref/`). Ensembl VCFs can be downloaded 
by running the `DownloadMappings.py` script in the project directory - **only needed for the  Harmonize HmVCF script and the pipelines**.


## Run single script
### Current options
* Harmonize script
```
$ python ./Harmonize.py -h
usage: Harmonize.py [-h] {HmPOS,HmVCF} ...

Harmonize a PGS Catalog Scoring file (PGS######.txt.gz) to a specific genome
build.

positional arguments:
  {HmPOS,HmVCF}  Harmonization Commands help
    HmPOS        HmPOS - Harmonizing position information (adding/updating
                 chr/pos information)
    HmVCF        HmVCF - Checking positional information and/or adding
                 other_alleles

optional arguments:
  -h, --help     show this help message and exit
```

* HmPOS option
```
python ./Harmonize.py HmPOS -h
usage: Harmonize.py HmPOS [-h] [-loc_files DIR] [-source_build GENOMEBUILD] 
                          [-loc_hmoutput DIR] [-var2location DIR] [--useAPI] 
                          [--catchmissingAPI] [--silent_tqdm] [--ignore_rsid] [--gzip]
                          PGS###### GRCh3#

positional arguments:
  PGS######             PGS Catalog Score ID
  GRCh3#                Target genome build choices: 'GRCh37'or GRCh38'

options:
  -h, --help            show this help message and exit
  -loc_files DIR        Root directory where the PGS files are located, 
                        otherwise assumed to be in: ../pgs_ScoringFiles/
  -source_build GENOMEBUILD
                        Source genome build [overwrites information in the 
                        scoring file header]
  -loc_hmoutput DIR     Directory where the harmonization output will be saved 
                        (default: PGS_HmPOS/)
  -var2location DIR     Root directory where DB of PGS Catalog rsID to chr/pos 
                        mappings is stored (default: ./map/ENSEMBL/)
  --useAPI              Uses the ENSEMBL API (not tractable for scores >1000 variants)
  --catchmissingAPI     Query the ENSEMBL API for variants missing from the 
                        PGS Catalog var2location DB
  --silent_tqdm         Disables tqdm progress bar
  --ignore_rsid         Ignores rsID mappings and harmonizes variants using 
                        only liftover
  --gzip                Writes gzipped harmonized output
```


* HmVCF option
```
python ./Harmonize.py HmVCF -h
usage: Harmonize.py HmVCF [-h] [-loc_files DIR] [-loc_hmoutput DIR] 
                          [-loc_vcfs DIR] [-cohort_vcf COHORT] 
                          [--addOtherAllele] [--addVariantID] 
                          [--skip_strandflips] [--split_unmappable] 
                          [--keep_duplicates] [--gzip] [--silent_tqdm]
                          PGS###### GRCh3#

positional arguments:
  PGS######           PGS Catalog Score ID
  GRCh3#              Target genome build choices: 'GRCh37'or GRCh38'

options:
  -h, --help          show this help message and exit
  -loc_files DIR      Root directory where the PGS files are located, 
                      otherwise assumed to be in: PGS_HmPOS/
  -loc_hmoutput DIR   Directory where the harmonization output will be saved 
                      (default: PGS_HmVCF/)
  -loc_vcfs DIR       Directory where the VCF files are located, otherwise 
                      assumed to be in: map/vcf_ref/
  -cohort_vcf COHORT  Cohort VCF: Used to check if a variant is present in the 
                      genotyped/imputed variants for a cohort and add other 
                      allele when the information from ENSEMBL is ambiguous 
                      (multiple potential alleles)
  --addOtherAllele    Adds a other_allele(s) column for PGS that only have 
                      a recorded effect_allele
  --addVariantID      Returns a column with the ID from the VCF corresponding 
                      to the matched variant/allele(s)
  --skip_strandflips  This flag will stop the harmonizer from correcting 
                      strand flips
  --split_unmappable  This flag will write unmapped & uncorrected variants 
                      (hm_code < 0) to separate files (suffixes: [.mapped, 
                      .unmatched])
  --keep_duplicates   This flag will allows duplicate variants to be present 
                      in the mapped variant file. The default behaviour is to 
                      drop them.
  --gzip              Writes gzipped harmonized output
  --silent_tqdm       Disables tqdm progress bar
```

### Examples
To test that the pipeline can run try these commands on the `tests/data` in the project directory:
    
    # PGS000015
    ## GRCh37
    python Harmonize.py HmPOS PGS000015 GRCh37 -loc_files ./tests/data/ --gzip
    python Harmonize.py HmVCF PGS000015 GRCh37 --gzip
    ## GRCh38
    python Harmonize.py HmPOS PGS000015 GRCh38 -loc_files ./tests/data/ --gzip
    python Harmonize.py HmVCF PGS000015 GRCh38 --gzip
    
    # PGS000065
    ## GRCh37
    python Harmonize.py HmPOS PGS000065 GRCh37 -loc_files ./tests/data/ --gzip
    python Harmonize.py HmVCF PGS000065 GRCh37 --gzip
    ## GRCh38
    python Harmonize.py HmPOS PGS000065 GRCh38 -loc_files ./tests/data/ --gzip
    python Harmonize.py HmVCF PGS000065 GRCh38 --gzip
   

## Run via Nextflow pipeline

The pipeline generates the HmPOS, HmVCF and the finalised version of each Scoring file.

### Additional Requirements

There are additional requirements to run **pgs-harmonizer** as a pipeline:
* [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation)
* [SQLite](https://www.sqlite.org/index.html): The pipeline generates a SQLite PGS Harmonized Knowledge Base (or populate an existing one).

Example of command to create a new Knowledge Base:
```
sqlite3 pgs_harmonizated_kb.db < harmonize_KB_schema.sql
```


### Run via a list of PGS IDs

This pipeline take a list of PGS IDs and run them in parallel.
It generates, for the given genome build:
 * 1 Harmonized Position file
 * 1 Harmonized VCF file
 * 1 Harmonized "Final" file: starting from the HmVCF file it updates the metadata, change the file name and store information in the SQLite Knowledge Base.

#### Nextflow pipeline configuration (nextflow_list.config)
```groovy
gb_version = <target genome build - numeric version, e.g. '38'>
root_dir = <path to the working directory (other than Nextflow's)>

params {
    pgs = <list of PGS IDs to harmonize, e.g. ['PGS000020','PGS000140']>
    genomebuild_grch = "GRCh$gb_version" <target genome build - string version, e.g. 'GRCh38'>
    genomebuild = "$gb_version" <target genome build - numeric version, e.g. '38'>
    pgs_ids_file = <name of the locally generated PGS IDs list file, e.g. "pgs_ids.txt">
    rest_api_url = <URL to the REST API server>
    loc_scripts = <path to the pgs-harmonizer directory, e.g. /Users/my_account/pgs-harmonizer/ (absolute) or ../../../ (relative)>
    loc_files = <path to the Scoring files directory>
    loc_hmoutput = <path to the HmPOS output directory, e.g. "$root_dir/HmPOS/$gb_version/">
    loc_hmoutput_vcf = <path to the HmVCF output directory, e.g. "$root_dir/HmVCF/$gb_version/">
    loc_var2location = <path to the directory containing the variant coordinates Knowledge Bases>
    loc_vcfs = <path to the Ensembl VCFs - optional if they are located in 'map/vcf_ref/'>
    loc_staged = <path to the directory containing the final version of the harmonized files, e.g. "$root_dir/HmFinal/$gb_version/">
    hm_version_sqlite_file_path = <path to the Harmonization Knowledge Base, e.g. 'pgs_harmonizated_kb.db'>
}
```

If it runs on LSF, here is an example about what to add to nextflow.config:
```groovy
// Only to be used if the script is running on LSF, on the cloud or in a container.
// More information is available here: https://www.nextflow.io/docs/latest/executor.html
executor {
  name = 'lsf'
  queueSize = 25
}

process {
  queue = 'short'
  withLabel: medium_mem {
    memory = '1 GB'
  }
  withLabel: retry_increasing_mem {
    errorStrategy = 'retry'
    memory = {2.GB * task.attempt}
    maxRetries = { task.exitStatus in [130,140] ? 10 : 2 }
  }
  withLabel: retry_increasing_big_mem {
    errorStrategy = 'retry'
    memory = {3.GB * task.attempt}
    maxRetries = { task.exitStatus in [130,140] ? 10 : 2 }
  }
}
```

#### Run Nextflow pipeline
Example of command to run the pipeline:
```
nextflow harmonizer_list_dsl2.nf -c nextflow_list.config
```

### Run via a range of PGS IDs

This pipeline take a range of PGS numeric IDs and run the rang of PGS IDs in parallel.
For instance `pgs_num_from=1` and `pgs_num_to=5`: it will generates harmonized files for PGS000001,PGS000002,PGS000003,PGS000004,PGS000005, for the given genome build:
 * 1 Harmonized Position file
 * 1 Harmonized VCF file
 * 1 Finalized file: it updates the metadata, change the file name and store information in the SQLite Knowledge Base 

#### Nextflow pipeline configuration (nextflow_range.config)

The only change from `nextflow_list.config` is that it replace **pgs** by **pgs_num_from** and **pgs_num_to**.

```groovy
gb_version = <target genome build - numeric version, e.g. '38'>
root_dir = <path to the working directory (other than Nextflow's)>

params {
    pgs_num_from = <PGS ID start number, e.g. '20' for PGS000020>
    pgs_num_to = <PGS ID end number, e.g. '140' for PGS000140>
    genomebuild_grch = "GRCh$gb_version" <target genome build - string version, e.g. 'GRCh38'>
    genomebuild = "$gb_version" <target genome build - numeric version, e.g. '38'>
    pgs_ids_file = <name of the locally generated PGS IDs list file, e.g. "pgs_ids.txt">
    rest_api_url = <URL to the REST API server>
    loc_scripts = <path to the pgs-harmonizer directory, e.g. /Users/my_account/pgs-harmonizer/ (absolute) or ../../../ (relative)>
    loc_files = <path to the Scoring files directory>
    loc_hmoutput = <path to the HmPOS output directory, e.g. "$root_dir/HmPOS/$gb_version/">
    loc_hmoutput_vcf = <path to the HmVCF output directory, e.g. "$root_dir/HmVCF/$gb_version/">
    loc_var2location = <path to the directory containing the variant coordinates Knowledge Bases>
    loc_vcfs = <path to the Ensembl VCFs - optional if they are located in 'map/vcf_ref/'>
    loc_staged = <path to the directory containing the final version of the harmonized files, e.g. "$root_dir/HmFinal/$gb_version/">
    hm_version_sqlite_file_path = <path to the Harmonization Knowledge Base, e.g. 'pgs_harmonizated_kb.db'>
}
```

If it runs on LSF, here is an example about what to add to nextflow.config:
```groovy
// Only to be used if the script is running on LSF, on the cloud or in a container.
// More information is available here: https://www.nextflow.io/docs/latest/executor.html
executor {
  name = 'lsf'
  queueSize = 25
}

process {
  queue = 'short'
  withLabel: medium_mem {
    memory = '1 GB'
  }
  withLabel: retry_increasing_mem {
    errorStrategy = 'retry'
    memory = {2.GB * task.attempt}
    maxRetries = { task.exitStatus in [130,140] ? 10 : 2 }
  }
  withLabel: retry_increasing_big_mem {
    errorStrategy = 'retry'
    memory = {3.GB * task.attempt}
    maxRetries = { task.exitStatus in [130,140] ? 10 : 2 }
  }
}
```

#### Run Nextflow pipeline
Example of command to run the pipeline:
```
nextflow harmonizer_range_dsl2.nf -c nextflow_range.config
```


## Description of Harmonization Codes (`hm_code`)
Once the HmVCF function is run each variant is assigned a value in the harmonization code (`hm_code`) column that 
reflects how the variant appears in the target variant data.  

     +----+--------------------------------------------------------------+
     |Code|Description of harmonisation process                          |
     +----+--------------------------------------------------------------+
     | 5  | Mapped [variant exists in REFERENCE VCF]                     |
     | 4  | Mapped [variant exists in REFERENCE VCF, other allele(s) with| 
     |    |         ambiguous orientations exist at the locus]           |
     | 3  | Mapped [variant in REFERENCE VCF with ambiguous orientation  |
     |    |         (e.g. A/T, C/G variants) ]                           |
     | 1  | Duplicated harmonized variant (multiple lines map to 1 ID    |
     | 0  | Author-reported variant information (not found in VCF)       |
     |-1  | Unable to map the variant                                    |
     |-4  | Strands flipped? [reverse complement alleles exist in VCF,   |
     |    |                  default behaviour is to correct by flipping]|
     |-5  | Variant doesn't exist in the REFERENCE VCF                   |
     +----+--------------------------------------------------------------+

## Related resources
### OpenTargets [GWAS summary statistics harmoniser](https://github.com/opentargets/genetics-sumstat-harmoniser)
- Harmonization flowchart : https://github.com/opentargets/genetics-sumstat-harmoniser/blob/master/flowchart_v3.svg
### GWAS Catalog [Summary Statistics harmonisation](https://github.com/EBISPOT/gwas-sumstats-harmoniser)
- Slightly different but uses the OpenTargets code/pipeline
### Liftover tools
- UCSC liftover tools 
- [pyliftover](https://pypi.org/project/pyliftover/) - code: https://github.com/konstantint/pyliftover

### ENSEMBL APIs
#### Variation POST (used to find location by rsID):
- GRCh38 (Current release): https://rest.ensembl.org/documentation/info/variation_post
- GRCh37 (Past release): http://grch37.rest.ensembl.org/documentation/info/variation_post

