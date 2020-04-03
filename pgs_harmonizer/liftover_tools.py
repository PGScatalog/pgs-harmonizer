from pyliftover import LiftOver

map_release = {
    'NCBI34' : 'hg16',
    'NCBI35' : 'hg17',
    'NCBI36' : 'hg18',
    'GRCh37' : 'hg19',
    'GRCh38' : 'hg38'
} # ENSEMBL -> UCSC

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
        self.GetChain()

    def GetChain(self):
        '''Downloads the chain from UCSC '''
        self.chain_name = '{} to {}'.format(self.build_from, self.build_to)
        self.chain = LiftOver(self.build_from, self.build_to)

    def lift(self, chr, pos):
        lifted = self.chain.convert_coordinate('chr{}'.format(str(chr)), int(pos))
        if lifted is not None:
            if len(lifted) == 1:
                return lifted[0][0], int(lifted[0][1]), 'Mapped (liftover)'
            if len(lifted) > 1:
                return lifted[0][0], int(lifted[0][1]), 'Mapped (liftover !first result!)'
        else:
            return None, int(), None



