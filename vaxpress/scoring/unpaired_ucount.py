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
            type=float, default=-3.0,
            help='weight for unpaired ucount (default: -3.0)')),
    ]


    def __init__(self, weight, _length_cds):
        self.weight = weight
        self.length = _length_cds

    def score(self, seqs, pairingprobs):
        U_sups = []
        for seq, pairingprob in zip(seqs, pairingprobs):
            U_idx = [i for i, base in enumerate(seq) if base == 'U']
            Pi_cooarray = pairingprob['Pi_array']
            Pi_array = Pi_cooarray.sum(axis=0)
            # sum only the Pi of the U index
            total_unpairedu_probs = sum(1 - Pi_array[i] for i in U_idx) # to be minimized
            U_sups.append(total_unpairedu_probs)
        scores = [u_sup * self.weight for u_sup in U_sups]

        return {self.name: scores}, {self.name: U_sups}

    def annotate_sequence(self, seq, pairingprob):
        U_idx = [i for i, base in enumerate(seq) if base == 'U']
        Pi_cooarray = pairingprob['Pi_array']
        Pi_array = Pi_cooarray.sum(axis=0)
        U_sup = sum(1 - Pi_array[i] for i in U_idx)
        return {self.name: U_sup}

