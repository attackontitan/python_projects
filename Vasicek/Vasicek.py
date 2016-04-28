from __future__ import division
import math
import scipy.stats as s
import numpy
from numpy import *
from copy import copy, deepcopy
from  TridiagonalSolver import  TridiagonalSolve

def VasicekLimits(r0, sigma, kappa, theta, T, prob = 1e-6):
    max = s.distributions.norm.ppf(1- prob)
    min = s.distributions.norm.ppf(prob)
    
    if kappa != 0.:
        rmin = (r0 + theta*(math.exp(kappa*T)-1) + sigma * math.sqrt((math.exp(2 * kappa * T)-1)/2/kappa)*min)*math.exp(-kappa*T)
        rmax = (r0 + theta*(math.exp(kappa*T)-1) + sigma * math.sqrt((math.exp(2 * kappa * T)-1)/2/kappa)*max)*math.exp(-kappa*T) 
        
    elif kappa == 0.:
        rmin = r0 + sigma * min * math.sqrt(T)
        rmax = r0 + sigma * max * math.sqrt(T)
        
    return (rmin, rmax)
    
def VasicekParams(r0, M, sigma, kappa, theta, T, prob = 1e-6):
    dtau = float(T) / M
    strench = 0.5
    dr = (dtau / strench * 2 * sigma ** 2)**0.5
    (rmin, rmax) = VasicekLimits(r0, sigma, kappa, theta, T, prob)
    N = int((rmax - rmin)/dr) +1
    dr = (rmax - rmin) / (N-1)
    return (rmin, dr, N, dtau)
    
def VasicekDiagonals(sigma, kappa, theta, rmin, dr, N, dtau):
    diag = [0 ] * N
    sup = [0 ] * N
    sub = [0 ] * N
    for i in range (0, N):
        diag[i] = 1 + (rmin + i * dr) * dtau + sigma **2 * dtau / dr**2
        sup[i] = -dtau / dr / 2 * kappa * (theta - rmin - i * dr) - 0.5 * sigma**2 *dtau /dr **2
        sub[i] = dtau / dr * kappa * (theta - rmin - i * dr) / 2- 0.5 * sigma**2 *dtau /dr **2
    diag[0] = 1 + rmin * dtau + kappa * (theta - rmin) * dtau / dr
    diag[ N-1] = 1 - kappa * (theta - rmin - dr * (N-1)) * dtau / dr + (rmin + (N-1) * dr) * dtau
    sup[N-1] = 0
    sup[0] = -kappa * (theta - rmin)* dtau / dr
    sub[0] = 0
    sub[N-1] = kappa * (theta - rmin - (N-1) * dr) * dtau / dr
    res = numpy.array([sub, diag, sup])
    return res
    
def CheckExercise(v,eex):
    if type(eex) == float:
        eex =  [eex] * len(v) 
    dif = v-eex
    res = [False] * len(v)
    for i in range(0, len(v)):
        if dif[i] >= -1e-13 : res[i] = True
    return res
    
def CallExercise(R, ratio, tau):
    return ratio * math.exp(-R * tau)

def VasicekPolicyDiagonals(sub, dia, sup, old, new, eex):
    early = CheckExercise(new, eex)
    for i in range (0, len(dia)):
        if early[i] == True:
            sup[i] = 0
            sub[i] = 0
            dia[i] = float(old[i]) / eex[i]
    return numpy.array([sub, dia, sup])
    
#~ def TridiagonalSolve(a, b, c, d):
    #~ n=len(d)
    #~ res=zeros(n)
    #~ c[0]=c[0] / b[0]
    #~ d[0]=d[0] / b[0]
    #~ for i in range (1,n-1):
        #~ c[i]=c[i] / (b[i]-c[i-1]*a[i])
        #~ d[i]=(d[i]-d[i-1]*a[i])/(b[i]-c[i-1]*a[i])
    #~ c[n-1]=0
    #~ s= d[n-1]-d[n-2]*a[n-1]
    #~ t = b[n-1]-c[n-2]*a[n-1]
    #~ d[n-1]=s / t
    #~ res[n-1]=d[n-1]
    #~ for i in range(-2, -n-1,-1):
        #~ res[i] = d[i] - c[i]* res[i+1]
    #~ return res
    
def test(a,b):
    flag = 1
    for i in range (0, len(a)):
        if a[i] != b[i]: 
            flag = 0
            break
    return flag
    
def Iterate(sub, dia,sup, old, eex, maxPolicyIterations = 10):
    csup = copy(sup)
    csub = copy(sub)
    cdia = copy(dia)
    cold = copy(old)
    val = TridiagonalSolve(csub,cdia, csup, cold)
    ear = CheckExercise(val, eex)
    count = 0
    for i in range (0, maxPolicyIterations):
        nsup = copy(sup); nsub = copy(sub); ndia = copy(dia); nold = copy(old)
        (nsub, ndia, nsup) = VasicekPolicyDiagonals(nsub, ndia, nsup, old, val, eex)
        val_new = TridiagonalSolve(nsub, ndia, nsup, nold)
        checknew = CheckExercise(val_new, eex)
        if (checknew == ear): 
            ear = checknew
            if test(val, val_new) != 1: 
                count += 1
            val = val_new
            break
        else:
            count += 1
    for m in range(0, len(ear)):
        if ear[m] == True: 
            val[m] = eex[m]
    return (val, min(count, maxPolicyIterations))
    
    
def VasicekCallableZCBVals( r0, R, ratio, T,  sigma, kappa, theta, M, prob = 1e-6,maxPolicyIterations = 10 ):
    (rmin, dr, N, dtau) = VasicekParams(r0, M, sigma, kappa, theta, T, prob)
    r = arange(0, N) *dr + rmin
    taus = arange(0, M+1) * dtau
    (sub, dia, sup) = VasicekDiagonals(sigma, kappa, theta, rmin, dr, N, dtau)
    old = [1] * N
    for m in range (0, M):
        eex = N * [CallExercise(R, ratio, taus[m+1])]
        new = Iterate(sub, dia,sup, old, eex, maxPolicyIterations)[0]
        old = new
    return(r, new)
    
def VasicekCallableZCB( r0, R, ratio, T,  sigma, kappa, theta, M, prob=1e-6, maxPolicyIterations = 10):
    (r, v) = VasicekCallableZCBVals( r0, R, ratio, T,  sigma, kappa, theta, M, prob, maxPolicyIterations)

    for  i in range (0, len(r)-1):
        if (r0 - r[i])*(r0 - r[i+1]) <=  0: 
            v = ((r0 - r[i]) * v[i+1] + (r[i+1] - r0)* v[i]) / (r[i+1]-r[i])
            break
    return v
        
r0 = 0.04
sigma = 0.02
kappa = 0.0
theta = 0.07
T = 5
prob = 1e-2
M = 5
ratio = 1.0
R = 0.02
print VasicekCallableZCBVals( r0, R, ratio, T,  sigma, kappa, theta, M, prob=prob)
        