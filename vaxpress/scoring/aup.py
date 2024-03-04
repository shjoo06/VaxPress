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
            help='weight for AUP (default: -5.0)')),
        ('base_weights', dict(metavar='BASE_WEIGHTS',
            type=dict, default={'A': 1.5, 'C': 1, 'G': 1, 'U': 3},
            help='weights for each base (default: {"A": 1.5, "C": 1, "G": 1, "U": 3})')),
    ]

    def __init__(self, weight, base_weights, _length_cds):
        self.weight = weight
        self.base_weights = base_weights
        self.length = _length_cds

    def score(self, seqs, pairingprobs):
        weighted_aups = []
        weights = np.array([self.base_weights[base] for base in 'ACGU']) # ACGU 순서로
        for seq, pairingprob in zip(seqs, pairingprobs):
            Pi_cooarray = pairingprob['Pi_array']
            Pi_array = Pi_cooarray.sum(axis=0)
            seqindex = np.frombuffer(seq.encode('utf-8').translate(bytes.maketrans(b'ACGU', b'\x00\x01\x02\x03')), dtype=np.uint8)
            # map weights to index
            wi = np.choose(seqindex, weights)
            xi = 1.0 - Pi_array
            weighted_aup = np.average(xi, weights=wi)
            weighted_aups.append(weighted_aup)
        scores = [weighted_aup * self.weight for weighted_aup in weighted_aups]
        return {self.name: scores}, {self.name: weighted_aups}

    def annotate_sequence(self, seq, pairingprob):
        normalized_weights = np.array([self.base_weights[base]/sum(self.base_weights.values()) for base in 'ACGU'])
        Pi_cooarray = pairingprob['Pi_array']
        Pi_array = Pi_cooarray.sum(axis=0)
        seqindex = np.frombuffer(seq.encode('utf-8').translate(bytes.maketrans(b'ACGU', b'\x00\x01\x02\x03')), dtype=np.uint8)
        wi = np.choose(seqindex, normalized_weights)
        xi = 1 - Pi_array
        weighted_aup = np.average(xi, weights=wi)
        return {self.name: weighted_aup}