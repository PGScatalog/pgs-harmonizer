import os
import pandas as pd
from cyvcf2 import VCF
from pgs_harmonizer.harmonize import reversecomplement, chromosomes, DetermineHarmonizationCode


class VCFResult:
    """Class to parse the results of a VCF Lookup and compare alleles"""

    def __init__(self, chr, pos, build, vcf_result=None):
        self.vcf_query = [chr, pos, build]
        self.vcf_result = vcf_result # List of variants
        if self.vcf_result is None:
            if chr in chromosomes:
                self.vcf_result = vcf_lookup(*self.vcf_query)

    def check_alleles(self, eff, ref=None, rsID = None):
        """Check if this variant exists in the ENSEMBL VCF files"""
        if (self.vcf_result is None) or (len(self.vcf_result) == 0):
            return False, False, False
        else:
            # Collect variants that match or mismatch the VCF
            v_hm = []
            v_hm_code = []
            # v_rsIDmatch = []
            # if rsID is None:
            #     v_rsIDmatch = None

            # Loop through all overlapping variants
            for v in self.vcf_result:
                hm_matchesVCF = False
                hm_isPalindromic = False
                hm_isFlipped = False

                alleles = [v.REF] + v.ALT
                alleles_rc = [reversecomplement(x) for x in alleles]
                # print('Alleles: {} | RC: {}'.format(alleles, alleles_rc))

                # Check allele(s) against VCF
                if ref is not None:
                    if (eff in alleles) and (ref in alleles):
                        hm_matchesVCF = True
                        if (eff in alleles_rc) and (ref in alleles_rc):
                            hm_isPalindromic = True
                    elif (eff in alleles_rc) and (ref in alleles_rc):
                        hm_isFlipped = True
                else:
                    if eff in alleles:
                        hm_matchesVCF = True
                        if eff in alleles_rc:
                            hm_isPalindromic = True
                    elif eff in alleles_rc:
                        hm_isFlipped = True

                hm_tuple = (hm_matchesVCF, hm_isPalindromic, hm_isFlipped)
                v_hm.append(hm_tuple)
                v_hm_code.append(DetermineHarmonizationCode(hm_matchesVCF, hm_isPalindromic, hm_isFlipped))
        return v_hm[v_hm_code.index(max(v_hm_code))]

    def infer_reference_allele(self, eff, oa_ensembl=None):
        """Try to infer the reference_allele from a vcf lookup (assumes the vcf is sorted by variant quality)"""
        if oa_ensembl is None:
            oa_ensembl = []
        else:
            oa_ensembl = oa_ensembl.split('/')  # [list of potential alleles, e.g. from ENSEMBL]
        oa_mapped = []  # list of potential reference alleles (if there are more variants)
        oa_mapped_hmcodes = []  # list of potential reference alleles (if there are more variants)

        # Loop through cohort variants to look for acceptable other_alleles
        for v in self.vcf_result:
            alleles = [v.REF] + v.ALT
            if len(alleles) < 2:
                continue
            alleles_rc = [reversecomplement(x) for x in alleles]

            hm_matchesVCF = False
            hm_isPalindromic = False
            hm_isFlipped = False

            # Check if
            oa = None
            if v.REF == eff:
                oa = v.ALT[0]
            elif eff in v.ALT:
                oa = v.REF

            if oa is not None:
                hm_matchesVCF = True
                if eff in alleles_rc:
                    hm_isPalindromic = True
            else:
                if eff in alleles_rc:
                    hm_matchesVCF = True
                    hm_isFlipped = True
                    alleles_rc.remove(eff)
                    oa = alleles_rc[0]

            hm = (hm_matchesVCF, hm_isPalindromic, hm_isFlipped)
            oa_mapped.append([oa, hm])

            if hm == (True, False, False):
                oa_mapped_hmcodes.append(5)
            elif hm == (True, True, False):
                oa_mapped_hmcodes.append(4)
            elif hm == (True, False, True):
                oa_mapped_hmcodes.append(-4)
            else:
                oa_mapped_hmcodes.append(-5)

        # Select the best OA
        #print(list(zip(oa_mapped, oa_mapped_hmcodes)))
        if len(oa_mapped_hmcodes) >= 1:
            i_bestmap = oa_mapped_hmcodes.index(max(oa_mapped_hmcodes))
            return oa_mapped[i_bestmap][0], oa_mapped[i_bestmap][1], oa_mapped_hmcodes[i_bestmap]
        else:
            return None, (False, False, False), -5


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
    def __init__(self, build, loc_vcfref='map/vcf_ref/', cohort_name=None, loc_cohortref='map/cohort_ref/'):
        self.VCF = None
        self.by_chr = {}
        self.build = build
        if cohort_name is None:
            for chr in chromosomes:
                loc_vcf = loc_vcfref + '{}/homo_sapiens-chr{}.vcf.gz'.format(self.build, chr)
                self.by_chr[chr] = VCF(loc_vcf)
        else:
            loc_vcf = '{}/{}.vcf.gz'.format(loc_cohortref, cohort_name)
            self.VCF = VCF(loc_vcf)

    def vcf_lookup(self,chromosome, position):
        """Lookup a variant in a specific genome build"""
        if chromosome not in chromosomes:
            return VCFResult(chromosome, position, self.build, [])

        if (type(position) is str) and ('-' in position):
            query = '{}:{}'.format(chromosome, position)
        else:
            query = '{}:{}-{}'.format(chromosome, position, position)

        if self.VCF is None:
            r_lookup = list(self.by_chr[chromosome](query))
        else:
            r_lookup = list(self.VCF(query))

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

#def guess_build(loc_file, vcf_root):