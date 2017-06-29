from math import log, exp
from numpy import matrix, diag
from scipy.stats import chi2 as chi2_distr


# 
# Functions used to perform the logistic fit
#


def logit(betas, xss):
    """
    logistic function
    """
    N = len(betas)
    xs = [1.] + xss
    bla = sum([betas[ii]*xs[ii] for ii in xrange(N)])
    return 1./(1.+exp(-bla))

def lnL(pops, Oi, betas):
    """
    likelihood for a given population and observations
    """
    Ei = [logit(betas, [ele]) for ele in pops]
    G2 = sum([Oi[ii]*log(Oi[ii]/ele) if Oi[ii] else 0. for ii, ele in enumerate(Ei)])
    p_value = (1 - chi2_distr.cdf(G2, len(Oi)-1))
    return G2, p_value

def Deviance(pops, ys, betas):
    """
    Deviance and its assocaited p-value
    """
    Ei = [logit(betas, [ele]) for ele in pops]
    G2 = sum([2.*ys[ii]*log(ys[ii]/ele) if ys[ii] else 2.*log(1./(1-ele)) for ii, ele in enumerate(Ei)])
    p_value = (1 - chi2_distr.cdf(G2, len(ys)-2))
    return G2, p_value



def dlnLdb(betas, Xs, ys, xs):
    """
    Derivative of the likelihood
    """
    Y = matrix(ys).T
    mu = matrix([logit(betas, [ele]) for ele in xs]).T
    return Xs.T * (Y-mu)

def dlnL2db2(betas, Xs, xs):
    """
    second derivative of the likelihood
    """
    pis = [logit(betas, [ele]) for ele in xs]
    W = matrix(diag([ele*(1.-ele) for ele in pis]))
    return -Xs.T * W * Xs

def get_betas(pops, ys, start=.000):
    """
    Obtains the parameters that maximize the likelihood
    """
    N = len(pops)
    K = 1
    Xs = matrix([[1. for ii in xrange(N)], pops]).T
    pi = sum(ys)*1./N
    oddi = log(pi/(1-pi))
    betas = matrix([oddi, start]).T
    for ii in xrange(200):
        betas = betas - .2*dlnL2db2(betas.T.tolist()[0], Xs, pops).I * dlnLdb(betas.T.tolist()[0], Xs, ys, pops)
    return betas.T.tolist()[0]


def chi2R(pops, Oi, betas):
    """
    One way to avaluate qualitity of the fit (chi2 test for each observation)
    """
    Ei = [logit(betas, [ele]) for ele in pops]
    chi2 = sum([(Oi[ii] - ele)**2 / ele for ii, ele in enumerate(Ei)])
    p_value = (1 - chi2_distr.cdf(chi2, len(Oi)-2))
    return p_value


def pvals(pops, betas):
    """
    Wald p-values and uncertainties for each parameter
    """
    kk = len(betas)
    N = len(pops)
    Xs = matrix([[1. for ii in xrange(N)], pops]).T
    H = dlnL2db2(betas, Xs, pops).I.tolist()
    incs = [-H[ii][ii] for ii in xrange(kk)]
    pvals = []
    for ii in xrange(kk):
        z = betas[ii]**2/incs[ii]
        pvals.append( (1 - chi2_distr.cdf(z, 1)) )
    return pvals, incs





