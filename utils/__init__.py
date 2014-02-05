#!/usr/bin/env python

import numpy as np

def flatten_list(list_of_lists):
    """ This is faster than the one-liner version:-

    def(flatten): return list(itertools.chain(*list_of_lists)) """

    flat_list = []
    x = flat_list.extend
    for sublist in list_of_lists:
        x(sublist)
    return flat_list

def lognormal_parameters(mean, stdev):
    """
    Given the mean and standard deviation wanted from a lognormal distribution,
    returns the mean (mu) and standard deviation (sigma) of the underlying
    normal distribution.
    This is useful for generating lognormal samples -
        params = lognormal_parameters(m, s)
        sample = numpy.random.lognormal(*params, size=N)
        [sample.mean()=m; sample.std()=s]

    Alternatively, sample = [random.lognormvariate(*params) for _ in range(N)]
    """
    mean=float(mean)
    stdev=float(stdev)
    variance = stdev**2
    sigma_sq = np.log( 1 + (variance/mean**2) )
    mu = np.log(mean) - sigma_sq/2
    return mu, np.sqrt(sigma_sq)
