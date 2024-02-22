from vaxpress.scoring import ScoringFunction
import numpy as np

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
    ]

    def __init__(self, weight, _length_cds):
        self.weight = weight
        self.length = _length_cds

    def score(self, seqs, pairingprobs):
        # pairingprob['Pi_array']를 사용하도록 수정.
        # 각각의 Pi_array: index is i, value is Pi

        aups = []

        for pairingprob in pairingprobs:
            Pi_array = pairingprob['Pi_array']
            sup = sum(1 - Pi for Pi in Pi_array)
            aup = sup/self.length
            aups.append(aup)
        scores = [aup * self.weight for aup in aups]

        return {self.name: scores}, {self.name: aups}

    def annotate_sequence(self, seq, pairingprob):
        Pi_array = pairingprob['Pi_array']
        sup = sum(1 - Pi for Pi in Pi_array)
        aup = sup/self.length
        return {self.name: aup}