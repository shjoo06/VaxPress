from vaxpress.scoring import ScoringFunction
import numpy as np
from scipy import sparse

class PairingProbFitness(ScoringFunction):

    name = 'aup'
    description = 'average unpaired probability'
    priority = 101
    uses_folding = False
    uses_basepairing_prob = True

    arguments = [
        ('weight', dict(metavar='WEIGHT',
            type=float, default=-5.0,
            help='scoring weight for AUP (default: -5.0)')),
        ('u-weight', dict(metavar='U_WEIGHT',
            type=float, default=3.0,
            help='weight for U (default: 3.0)')),
        ('a-weight', dict(metavar='A_WEIGHT',
            type=float, default=1.5,
            help='weight for A (default: 1.5)')),
    ]

    def __init__(self, weight, u_weight, a_weight, _length_cds):
        self.weight = weight
        self.u_weight = u_weight
        self.a_weight = a_weight
        self.length = _length_cds
        self.trans_table = bytes.maketrans(b'ACGU', b'\x00\x01\x02\x03')

    def score(self, seqs, pairingprobs):
        weighted_aups = []
        weights = np.array([self.a_weight, 1, 1, self.u_weight]) # ACGU
        for seq, pairingprob in zip(seqs, pairingprobs):
            Pi_cooarray = pairingprob['Pi_array']
            Pi_array = Pi_cooarray.sum(axis=0)
            seqindex = np.frombuffer(seq.encode().translate(self.trans_table), dtype=np.uint8)
            # map weights to index
            wi = np.choose(seqindex, weights)
            xi = 1.0 - Pi_array
            weighted_aup = np.average(xi, weights=wi)
            weighted_aups.append(weighted_aup)
        scores = [weighted_aup * self.weight for weighted_aup in weighted_aups]
        return {self.name: scores}, {self.name: weighted_aups}

    def annotate_sequence(self, seq, pairingprob):
        weights = np.array([self.a_weight, 1, 1, self.u_weight])
        Pi_cooarray = pairingprob['Pi_array']
        Pi_array = Pi_cooarray.sum(axis=0)
        seqindex = np.frombuffer(seq.encode().translate(self.trans_table), dtype=np.uint8)
        wi = np.choose(seqindex, weights)
        xi = 1 - Pi_array
        weighted_aup = np.average(xi, weights=wi)
        return {'weighted_aup': weighted_aup, 'aup': np.average(xi)}