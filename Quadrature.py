from __future__ import division
from numpy import *

def trap(f, a, b, N):
    h = (b-a) / N
    res = 0.5 * f(a) + 0.5 * f(b)
    for i in range (1, N):
        res += f(a + i * h)
    return res * h
    
def phibs(T, u, sigma):
    res = exp(-0.5 * u * (u + 1j) * sigma ** 2 * T)
    return res
    
def quadbscall(T, S, K, sigma, N, uMin=0, uMax=20):
    def func(x):
        return (exp(-1j * x * log(K/S) )* phibs(T, x-0.5j, sigma)).real/(x**2 + 0.25)
    sum = trap(func, uMin, uMax, N)
    return S- sqrt(S*K)/math.pi * sum
    
def phiheston(T, u, sigma, sigmaBar, lamda, eta, rho):
    alpha = -0.5 * (u**2+ u*1j)
    beta = lamda - rho * eta * 1j * u
    gamma = 0.5 * eta **2
    omega = sqrt(beta **2 - 4 * alpha * gamma)
    rplus = (beta + omega) / eta**2
    rminus = (beta - omega) / eta**2
    g = rminus / rplus
    C = lamda * (rminus * T - 2 / eta ** 2 * log((1- g * exp(-omega * T))/(1- g)))
    D = rminus * (1-exp(-omega * T)) / (1 - g * exp(-omega * T))
    return exp(C * sigmaBar + D * sigma)
    
def quadhestoncall(T, S, K, sigma, sigmaBar, lamda, eta, rho,N, uMin=0, uMax=20):
    def func(x):
        return (exp(-1j * x * log(K/S) )* phiheston(T, x - 0.5j, sigma, sigmaBar, lamda, eta, rho)).real/(x**2 + 0.25)
    sum = trap(func, uMin, uMax, N)
    return S - sqrt(S*K)/math.pi * sum