from vaxpress.scoring import ScoringFunction
import re
import numpy as np

class PairingProbFitness(ScoringFunction):

    name = 'aup'
    description = 'average unpaired probability'
    priority = 101
    uses_folding = False
    uses_basepairing_prob = True

    arguments = [
        ('weight', dict(metavar='WEIGHT',
            type=float, default=-8.0,
            help='weight for AUP (default: -8.0)')),
    ]

    def __init__(self, weight, _length_cds):
        self.weight = weight
        self.length = _length_cds

    def score(self, seqs, pairingprobs):
        aups = [pairingprob['sup']/self.length for pairingprob in pairingprobs]
        scores = [aup * self.weight for aup in aups]

        return {self.name: scores}, {self.name: aups}

    def annotate_sequence(self, seq, pairingprob):
        aup = pairingprob['sup']/self.length
        return {self.name: aup}
