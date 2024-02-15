from vaxpress.scoring import ScoringFunction
import numpy as np

class UnpairedUridineFitness(ScoringFunction):

    name = 'unpaired_u'
    description = 'average unpaired uridine count'
    priority = 102
    uses_folding = False
    uses_basepairing_prob = True

    arguments = [
        ('weight', dict(metavar='WEIGHT',
            type=float, default=-5.0,
            help='weight for unpaired ucount (default: -5.0)')),
    ]


    def __init__(self, weight, _length_cds):
        self.weight = weight
        self.length = _length_cds

    def score(self, seqs, pairingprobs):
        ucounts = [seq.count('U') for seq in seqs]
        unpaired_ucounts = [(pairingprob['sup']/self.length) * ucount for pairingprob, ucount in zip(pairingprobs, ucounts)]
        scores = [unpaired_ucount * self.weight for unpaired_ucount in unpaired_ucounts]

        return {self.name: scores}, {self.name: unpaired_ucounts}

    def annotate_sequence(self, seq, pairingprob):
        ucount = seq.count('U')
        unpaired_ucount = (pairingprob['sup']/self.length) * ucount
        return {self.name: unpaired_ucount}

