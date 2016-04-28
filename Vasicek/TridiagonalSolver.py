from numpy import *

def TridiagonalSolve(a, b, c, d):
    n=len(d)
    res=zeros(n)
    c[0]=c[0] / b[0]
    d[0]=d[0] / b[0]
    for i in range (1,n-1):
        c[i]=c[i] / (b[i]-c[i-1]*a[i])
        d[i]=(d[i]-d[i-1]*a[i])/(b[i]-c[i-1]*a[i])
    c[n-1]=0
    s= d[n-1]-d[n-2]*a[n-1]
    t = b[n-1]-c[n-2]*a[n-1]
    d[n-1]=s / t
    res[n-1]=d[n-1]
    for i in range(-2, -n-1,-1):
        res[i] = d[i] - c[i]* res[i+1]
    return res