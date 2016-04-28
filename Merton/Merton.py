from __future__ import division
from numpy import *
import math
import scipy
from scipy.optimize import newton
from scipy.stats import norm
from numpy.random import multivariate_normal as M
N = norm.cdf 

def BScall(S, K, sigma , r, T):
    q = 0
    t = 0
    if S == 0:
        return 0
    else:
        d1 =(log((S) / K) +(r-q+((sigma)*sigma)/2.)*(T-t))/(sigma * sqrt(T-t))
        d2 =(log((S)/ K) +(r-q-((sigma)*sigma)/2.)*(T-t))/(sigma * sqrt(T-t))
        call_optionValue = S*exp(-float(q)*(T-t))*norm.cdf (d1) - K*exp(-float(r)*(T-t))*norm.cdf (d2)
        return call_optionValue

def VStar(K,t,r, VB, sigmaV, T):
    def Call(S):
        return BScall(S, VB, sigmaV, r, T-t) - K
    if sigmaV != 0:
        sol = newton(Call, VB, tol = 0.00000001)
    else:
        sol = max(K * math.exp(r * (T-t)) + VB, 0)
    return sol
    
def Stock(V, VB, sigmaV, r, T):
    return BScall(V, VB, sigmaV, r, T)
    
def MertonConstants(K,t,r, V, VB, sigmaV, T):
    vstar = VStar(K, t, r, VB, sigmaV,T)
    a1 = (log(V/vstar) + (r+sigmaV**2 / 2) * t)/sigmaV/sqrt(t)
    a2 = (log(V/vstar) + (r- sigmaV**2 / 2) * t)/sigmaV/sqrt(t)
    b1 = (log(V/VB) + (r+sigmaV**2 / 2) * T)/sigmaV/sqrt(T)
    b2 = (log(V/VB) + (r- sigmaV**2 / 2) * T)/sigmaV/sqrt(T)
    return(a1, a2, b1, b2)
    
def M(a,b,rho):
	def bivariate(x,y):
		return 1/(2*math.pi*sqrt(1-rho**2))*math.exp(-(x**2+y**2-2*rho*x*y)/(2*(1-rho**2)))
	inf=float("inf")
	return scipy.integrate.dblquad (bivariate,-inf,a,lambda x:-inf,lambda x:b)[0]
    
def MertonFormula(callput, K, t, r, V, VB, sigmaV, T):
    if sigmaV == 0:
        if callput == 1 : return max(((V - VB) * exp(-r*(T-t)) -K) * exp(-r*t), 0)
        if callput == -1: return max((K - (V - VB) * exp(-r*(T-t)) ) * exp(-r*t), 0)
    else:
        (a1, a2, b1, b2) = MertonConstants(K,t,r, V, VB, sigmaV, T)
        square = sqrt(t/T)
        if callput == 1:
            return V * M(a1, b1, square) - VB * exp(-r*T) *M(a2, b2, square) - K * exp(-r * t) * N(a2)
        if callput == -1:
            return -V * M(-a1, b1, -square) + VB * exp(-r*T) *M(-a2, b2, -square) + K * exp(-r * t) * N(-a2)
            
def BSVol( callput,K,t,r,q,S, price ):
    def BSfunc(sigma):
        d1 =(log(float(S) / K) +(r-q+(float(sigma)*sigma)/2.)*(t))/(sigma * sqrt(t))
        d2 =(log(float(S)/ K) +(r-q-(float(sigma)*sigma)/2.)*(t))/(sigma * sqrt(t))
        if callput == 1:
            return S*exp(-float(q)*(t))*norm.cdf (d1) - K*exp(-float(r)*(t))*norm.cdf (d2) - price
        if callput == -1:
            return float(K)*exp(-float(r)*(t))*norm.cdf(-d2)-float(S)*(exp(-float(q)*(t)))*norm.cdf(-d1) - price
    sol = newton(BSfunc, 0.2)
    return sol
    
def MertonObjective( callputs,Ks,ts,rs,qs,optPrices,S, V, VB, sigmaV, T ):
    num = len(callputs)
    sig_diff = 0
    S_diff = 0
    for i in range (0, num):
        sigma_i = BSVol( callputs[i],Ks[i],ts[i],rs[i],qs[i],S, optPrices[i])
        price = MertonFormula(callputs[i], Ks[i], ts[i], rs[i], V, VB, sigmaV, T)
        sigmahat_i = BSVol( callputs[i],Ks[i],ts[i],rs[i],qs[i],S, price)
        S_i = Stock(V, VB, sigmaV, rs[i], T)
        sig_diff += (sigma_i - sigmahat_i)**2
        S_diff += (1- S_i / S) **2
    return 1./num * (sig_diff + S_diff)
        