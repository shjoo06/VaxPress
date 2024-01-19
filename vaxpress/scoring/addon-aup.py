from vaxpress.scoring import ScoringFunction
import re
from vaxpress import linearpartition

class PairingProbFitness(ScoringFunction):

    name = 'aup'
    description = 'average unpaired probability'
    priority = 101
    uses_folding = False

    arguments = [
        ('weight', dict(metavar='WEIGHT',
            type=float, default=-10.0,
            help='weight for unpaired probability (default: -10.0)')),
    ]

    def __init__(self, weight, _length_cds):
        self.weight = weight
        self.length = _length_cds

    def score(self, seqs):
        aups = []
        #unpaired_ucounts = []    # choose 1 of 2
        scores = []
        for seq in seqs:
            aup, unpaired_ucount = linearpartition.get_pairingprob(seq)
            aups.append(aup)
            #unpaired_ucounts.append(unpaired_ucount)
            scores.append(aup * self.weight)

        return {self.name: scores}, {self.name: aups}

    def annotate_sequence(self, seq):
        aup, unpaired_ucount = linearpartition.get_pairingprob(seq)
        return {self.name: aup}
    