import numpy as np
import scipy as sp
import scipy.stats as st
from numpy import *
from random import *
from pylab import *
norminv = st.distributions.norm.ppf
norm = st.distributions.norm.cdf
from numpy.random import rand
import numpy.linalg


def defaultTimes(h, c, s = 1):
    if type(s) != int:
        N = len(s[0])
        y = s
    else:
        M = s
        N = len(h)
        y = rand(M,N)
    w = norminv(y)
    c = np.eye(N) * (1-c) + c * np.ones((N,N))
    ch = numpy.linalg.cholesky(c).T
    z = np.dot(w,ch)
    x = norm(z)
    tau = -np.log(1-x) / h
    return tau
    
def copulaLosses(r, T, notionals, recoveryRates, hazardRates, correlation, samples = 1):
    if type(samples) == int:
        M = samples
        N = len(hazardRates)
        samples = rand(M,N)
    tau = defaultTimes(hazardRates, correlation, samples)
    M = len(samples)
    N = len(samples[0])
    V = notionals * np.exp(-r * tau) * np.ones((M,N)) *(1-recoveryRates)
    losses = V
    losses[tau > T] = 0
    lossSamples = losses.sum(axis = 1)
    return lossSamples
    
    
def tranche(K1, K2, r, T, notionals, recoveryRates, hazardRates, correlation, samples = 1):
    lossSamples = copulaLosses(r, T, notionals, recoveryRates, hazardRates, correlation, samples)
    payoff = lossSamples - K1
    payoff[lossSamples<K1] = 0
    payoff[lossSamples > K2] = abs(K2 - K1)
    return mean(payoff)
    

    
