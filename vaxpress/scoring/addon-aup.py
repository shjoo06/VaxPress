from vaxpress.scoring import ScoringFunction
import re
import linearpartition as lp
import numpy as np

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
            bpmtx, _ = lp.partition(seq)

            bpmtx_square = np.concatenate((bpmtx, np.array([(x[1], x[0], x[2]) for x in bpmtx], dtype=bpmtx.dtype)))
            _, inverse_indices = np.unique(bpmtx_square['i'], return_inverse=True)
            sup = np.sum(1 - np.bincount(inverse_indices, weights=bpmtx_square['prob']))
            aup = sup / self.length

            aups.append(aup)
            scores.append(aup * self.weight)
        return {self.name: scores}, {self.name: aups}

    def annotate_sequence(self, seq):
        bpmtx, fe = lp.partition(seq)
        bpmtx_transposed = np.array([(x[1], x[0], x[2]) for x in bpmtx], dtype=[('i', '<i4'), ('j', '<i4'), ('prob', '<f8')])
        bpmtx_square = np.concatenate((bpmtx, bpmtx_transposed))
        _, inverse_indices = np.unique(bpmtx_square['i'], return_inverse=True)
        OneMinusPi = 1 - np.bincount(inverse_indices, weights=bpmtx_square['prob'])

        sup = np.sum(OneMinusPi)
        aup = sup / self.length
        return {self.name: aup}
