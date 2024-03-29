from pyliftover import LiftOver

map_release = {
    'NCBI34' : 'hg16',
    'NCBI35' : 'hg17',
    'NCBI36' : 'hg18',
    'GRCh37' : 'hg19',
    'GRCh38' : 'hg38'
} # ENSEMBL -> UCSC

chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8',
               '9', '10', '11', '12', '13', '14', '15', '16',
               '17', '18', '19', '20', '21', '22',
               'X', 'x', 'Y', 'y', 'XY', 'xy', 'MT', 'Mt', 'mt']


class liftover:
    def __init__(self, build_from, build_to):
        # Source Genome Build
        if build_from in map_release.values():
            self.build_from = build_from
        else:
            build_mapped = map_release.get(build_from)
            if build_mapped is None:
                raise Exception('Unknown SOURCE genome build. The value was: {}'.format(build_from))
            else:
                self.build_from = build_mapped

        # Destination Genome Build
        if build_to in map_release.values():
            self.build_to = build_to
        else:
            build_mapped = map_release.get(build_to)
            if build_mapped is None:
                raise Exception('Unknown DESTINATION genome build. The value was: {}'.format(build_from))
            else:
                self.build_to = build_mapped

        # Download/Source the Chain from UCSC
        if self.build_from != self.build_to:
            self.GetChain()
        else:
            self.chain = None

    def GetChain(self):
        '''Downloads the chain from UCSC '''
        self.chain_name = 'UCSC: {} to {}'.format(self.build_from, self.build_to)
        self.chain = LiftOver(self.build_from, self.build_to)

    def lift(self, chr, pos):
        lifted = self.chain.convert_coordinate('chr{}'.format(str(chr)), int(pos)) # ToDo figure out whether this step should be adjusted for 0/1 indexing?
        if lifted is not None:
            if len(lifted) == 1:
                return lifted[0][0][3:], int(lifted[0][1]), False  # Only 1 position
            if len(lifted) > 1:
                for i in lifted:
                    if lifted[i][0][3:] in chromosomes:
                        return lifted[i][0][3:], int(lifted[0][1]), True # Multiple positions (take the first with standard chromosome)
                return lifted[0][0][3:], int(lifted[0][1]), True  # Multiple positions (take first)
            else:
                return None, None, None
        else:
            return None, None, None



