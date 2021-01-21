import os
import pandas as pd
from cyvcf2 import VCF
from pgs_harmonizer.harmonize import reversecomplement, chromosomes


class VCFResult:
    """Class to parse the results of a VCF Lookup and compare alleles"""

    def __init__(self, chr, pos, build, vcf_result=None):
        self.vcf_query = [chr, pos, build]
        self.vcf_result = vcf_result # List of variants
        if self.vcf_result is None:
            if chr in chromosomes:
                self.vcf_result = vcf_lookup(*self.vcf_query)

    def check_alleles(self, eff, ref=None):
        """Check if this variant exists in the ENSEMBL VCF files"""
        if self.vcf_result is None:
            return None
        else:
            # Collect variants that match or mismatch the VCF
            v_consistent = []
            v_flipped = []
            v_palindromic = []

            # Check all overlapping variants
            for v in self.vcf_result:
                alleles = [v.REF] + v.ALT
                alleles_rc = [reversecomplement(x) for x in alleles]
                #print('Alleles: {} | RC: {}'.format(alleles, alleles_rc))

                if eff in alleles:
                    if ref is not None:
                        if ref in alleles:
                            v_consistent.append(v) # Both alleles match the VCF
                    else:
                        v_consistent.append(v)  # Effect allele matches the VCF

                if eff in alleles_rc:
                    if ref is not None:
                        if ref in alleles_rc:
                            v_flipped.append(v) # Both alleles match the reverse complement VCF alleles
                    else:
                        v_flipped.append(v)  # Both alleles match the reverse complement VCF alleles

                if (v in v_consistent) and (v in v_flipped):
                    v_palindromic.append(v)
            #print(v_consistent, v_flipped, v_palindromic)

            # Decide on the harmonization code
            hm_result = -5

            hmcodes = {
                'mapped': 5,
                'mapped_palindromic': 4,
                'flippedstrand': -4
            }

            for v in self.vcf_result:
                if v in v_consistent:
                    if (v in v_palindromic) and (hmcodes['mapped_palindromic'] > hm_result):
                        hm_result = hmcodes['mapped_palindromic']
                    elif hmcodes['mapped'] > hm_result:
                        hm_result = hmcodes['mapped']
                elif (v in v_flipped) and (hmcodes['flippedstrand'] > hm_result):
                    hm_result = hmcodes['flippedstrand']

            # return hm_matchesVCF, hm_isPalindromic, hm_isFlipped
            if hm_result == 5:
                return True, False, False
            elif hm_result == 4:
                return True, True, False
            elif hm_result == -4:
                return True, False, True
            else:
                return False, False, False


def vcf_lookup(chromosome, position, build, loc_vcfref='map/vcf_ref/'):
    """Retrieve variants overlapping with a position in a VCF File"""
    if build not in ['GRCh37', 'GRCh38']:
        raise ValueError("Invalid Genome Build. Expected one of: GRCh37, GRCh38")

    if chromosome not in chromosomes:
        raise ValueError("Invalid Chromosome. Expected one of: {}".format(chromosomes))

    loc_vcf = loc_vcfref + '{}/homo_sapiens-chr{}.vcf.gz'.format(build, chromosome)
    vcf = VCF(loc_vcf)
    if (type(position) is str) and ('-' in position):
        return list(vcf('{}:{}'.format(chromosome, position)))
    else:
        return list(vcf('{}:{}-{}'.format(chromosome, position, position)))


class VCFs:
    """Class to open and hold all VCF files for a genome build"""
    def __init__(self, build, loc_vcfref='map/vcf_ref/'):
        self.by_chr = {}
        self.build = build
        for chr in chromosomes:
            loc_vcf = loc_vcfref + '{}/homo_sapiens-chr{}.vcf.gz'.format(self.build, chr)
            self.by_chr[chr] = VCF(loc_vcf)

    def vcf_lookup(self,chromosome, position):
        """Lookup a variant in a specific genome build"""
        if chromosome not in chromosomes:
            raise ValueError("Invalid Chromosome. Expected one of: {}".format(chromosomes))

        if (type(position) is str) and ('-' in position):
            r_lookup = list(self.by_chr[chromosome]('{}:{}'.format(chromosome, position)))
        else:
            r_lookup = list(self.by_chr[chromosome]('{}:{}-{}'.format(chromosome, position, position)))

        return VCFResult(chromosome, position, self.build, r_lookup)


class CohortVariants:
    """Class to load cohort-specific variants information and compare alleles"""
    def __init__(self, cohortname, variants_table):
        self.cohort = cohortname
        self.variants_table = variants_table

    def check_variant(self):
        """Checks whether the variant (chr, pos, alleles) is genotyped/imputed in this cohort"""

    def infer_reference_allele(self, chr, pos, eff, rsID = None):
        """Try to infer the reference_allele based on genotyped/imputed variants in this cohort"""

def load_cohortvariants(cohortname):
    loc_cohortvariants = 'map/cohort_ref/{}_variants.sorted.h5'.format(cohortname)
    if os.path.isfile(loc_cohortvariants):
        return CohortVariants(cohortname, pd.read_hdf(loc_cohortvariants, 'variants'))
    else:
        return None

#def guess_build(loc_file, vcf_root):