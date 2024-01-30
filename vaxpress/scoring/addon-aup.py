from vaxpress.scoring import ScoringFunction
import re
import linearpartition as lp

class PairingProbFitness(ScoringFunction):

    name = 'aup'
    description = 'average unpaired probability'
    priority = 101
    uses_folding = False

    arguments = [
        ('weight', dict(metavar='WEIGHT',
            type=float, default=-8.0,
            help='weight for AUP (default: -8.0)')),
    ]

    def __init__(self, weight, _length_cds):
        self.weight = weight
        self.length = _length_cds

    def score(self, seqs):
        aups = []
        scores = []
        for seq in seqs:
            bpmtx, fe = lp.partition(seq)
            Pi = {i: 0 for i in range(0, self.length)}
            for ijprob in bpmtx:   # bpmtx: list of tuples (i, j, Pij)
                Pi[ijprob[0]] += ijprob[2]
                Pi[ijprob[1]] += ijprob[2]
            # calculate aup
            sup = sum(1 - Pi[i] for i in Pi)
            aup = sup / self.length

            aups.append(aup)
            scores.append(aup * self.weight)

        return {self.name: scores}, {self.name: aups}

    def annotate_sequence(self, seq):
        bpmtx, fe = lp.partition(seq)
        Pi = {i: 0 for i in range(0, self.length)}
        for ijprob in bpmtx:
            Pi[ijprob[0]] += ijprob[2]
            Pi[ijprob[1]] += ijprob[2]

        sup = sum(1 - Pi[i] for i in Pi)
        aup = sup / self.length
        return {self.name: aup}