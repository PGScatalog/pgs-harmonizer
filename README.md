# pgs-harmonizer
A pipeline to format and harmonize [Polygenic Score (PGS) Catalog Scoring Files](http://www.pgscatalog.org/downloads/#dl_ftp) 
within and between different genome builds. 

_External users NB: the pipeline will run very slowly if only rsIDs are given, because it looks up each rsID in 200 
variant batches using the Ensembl API. There are ways to speed this up by looking up the total set of rsIDs once using the 
`EnsemblMappings/var2location_3738.pl` script and the public MySQL connection._


## Requirements
- Python packages: requests, pandas, pyliftover, tqdm, cyvcf2
- PGS Catalog rsID -> position mappings (downloadable from ...)
- VCFs from Ensembl (in `/map/vcf_ref/`) or a specific cohort (in `/map/cohort_ref/`). Ensembl VCFs can be downloaded 
by running the `DownloadMappings.py` script in the project directory.

## Current options
<pre>$ python ./Harmonize.py -h
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
  -h, --help     show this help message and exit</pre>
  
<pre>$ python ./Harmonize.py HmPOS -h
usage: Harmonize.py HmPOS [-h] [-loc_files DIR] [-source_build GENOMEBUILD]
                          [-loc_hmoutput DIR] [--var2location] [--silent_tqdm]
                          [--ignore_rsid] [--gzip]
                          PGS###### GRCh3#

positional arguments:
  PGS######             PGS Catalog Score ID
  GRCh3#                Target genome build choices: 'GRCh37'or GRCh38'

optional arguments:
  -h, --help            show this help message and exit
  -loc_files DIR        Root directory where the PGS files are located,
                        otherwise assumed to be in: ../pgs_ScoringFiles/
  -source_build GENOMEBUILD
                        Source genome build [overwrites information in the
                        scoring file header]
  -loc_hmoutput DIR     Directory where the harmonization output will be saved
                        (default: PGS_HmPOS/)
  -var2location DIR     Root directory where DB of PGS Catalog rsID to chr/pos
                        mappings is stored (default: ./EnsemblMappings)
  --useAPI              Uses the ENSEMBL API (not tractable for scores >1000
                        variants)
  --catchmissingAPI     Query the ENSEMBL API for variants missing from the
                        PGS Catalog var2location DB
  --silent_tqdm         Disables tqdm progress bar
  --ignore_rsid         Ignores rsID mappings and harmonizes variants using
                        only liftover
  --gzip                Writes gzipped harmonized output</pre>

<pre>$ python ./Harmonize.py HmVCF -h
usage: Harmonize.py HmVCF [-h] [-loc_files DIR] [-loc_hmoutput DIR]
                          [-cohort_vcf COHORT] [--addOtherAllele]
                          [--addVariantID] [--author_reported]
                          [--skip_strandflips] [--silent_tqdm] [--gzip]
                          PGS###### GRCh3#

positional arguments:
  PGS######           PGS Catalog Score ID
  GRCh3#              Target genome build choices: 'GRCh37'or GRCh38'

optional arguments:
  -h, --help          show this help message and exit
  -loc_files DIR      Root directory where the PGS files are located,
                      otherwise assumed to be in: PGS_HmPOS/
  -loc_hmoutput DIR   Directory where the harmonization output will be saved
                      (default: PGS_HmPOS/)
  -loc_vcfs DIR       Directory where the VCF files are located, otherwise
                      assumed to be in: map/vcf_ref/
  -cohort_vcf COHORT  Cohort VCF: Used to check if a variant is present in the
                      genotyped/imputed variants for a cohort and add other
                      allele when the information from ENSEMBL is ambiguous
                      (multiple potential alleles)
  --addOtherAllele    Adds a other_allele(s) column for PGS that only have a
                      recorded effect_allele
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
</pre>

## Examples
To test that the pipeline can run try these commands on the test data in the project directory:
    
    # PGS000015
    ## GRCh37
    python Harmonize.py HmPOS PGS000015 GRCh37 -loc_files ./test_data/ --gzip
    python Harmonize.py HmVCF PGS000015 GRCh37 --gzip
    ## GRCh38
    python Harmonize.py HmPOS PGS000015 GRCh38 -loc_files ./test_data/ --gzip
    python Harmonize.py HmVCF PGS000015 GRCh38 --gzip
    
    # PGS000065
    ## GRCh37
    python Harmonize.py HmPOS PGS000065 GRCh37 -loc_files ./test_data/ --gzip
    python Harmonize.py HmVCF PGS000065 GRCh37 --gzip
    ## GRCh38
    python Harmonize.py HmPOS PGS000065 GRCh38 -loc_files ./test_data/ --gzip
    python Harmonize.py HmVCF PGS000065 GRCh38 --gzip
   
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
### Variation POST (used to find location by rsID):
- GRCh38 (Current release): https://rest.ensembl.org/documentation/info/variation_post
- GRCh37 (Past release): http://grch37.rest.ensembl.org/documentation/info/variation_post