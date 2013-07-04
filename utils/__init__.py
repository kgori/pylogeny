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
    return mu, sqrt(sigma_sq)

def logN_correlated_rate(parent_rate, branch_length, autocorrel_param, size=1):
    """
    The log of the descendent rate, ln(Rd), is ~ N(mu, bl*ac), where
    the variance = bl*ac = branch_length * autocorrel_param, and mu is set
    so that E[Rd] = Rp:
    E[X] where ln(X) ~ N(mu, sigma^2) = exp(mu+(1/2)*sigma_sq)
    so Rp = exp(mu+(1/2)*bl*ac),
    ln(Rp) = mu + (1/2)*bl*ac,
    ln(Rp) - (1/2)*bl*ac = mu,
    so ln(Rd) ~ N(ln(Rp) - (1/2)*bl*ac, bl*ac)
    (NB: Var[Rd] = Rp^2 * (exp(bl*ac)-1),
         Std[Rd] = Rp * sqrt(exp(bl*ac)-1)

    See: H Kishino, J L Thorne, and W J Bruno (2001)
    """
    if autocorrel_param <= 0:
        raise Exception('Autocorrelation parameter must be greater than 0')
    
    variance = branch_length * autocorrel_param
    stdev = np.sqrt(variance)
    lnRd = np.random.normal(np.log(parent_rate) - 0.5*variance, scale=stdev, size=size)
    Rd = np.exp(lnRd)
    return float(Rd) if size==1 else Rd
