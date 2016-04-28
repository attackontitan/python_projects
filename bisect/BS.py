import math
import Bisect
import datetime
from dateutil import parser
import csv
from numpy import array

k = 120.0
S = 100.0
t = 1.0
sigma = 0.5
r = 0.0
callput = "Call"
q = 0

def erf(x):#""" CDF distributed between 0 and 1 """
  
  gamma = 0.2316419
  k = 1.0 / (1.0 + x * gamma)

  a1 = 0.319381530
  a2 = -0.356563782
  a3 = 1.781477937
  a4 = -1.821255978
  a5 = 1.330274429

  q = 1.0 / math.sqrt(2 * math.pi)
  N = q * math.exp(-(x*x)/2.0)

  if x >= 0:
    return 1 - (N) * (a1 * k + a2 * math.pow(k,2) + a3 * math.pow(k,3) + a4 * math.pow(k,4) + a5 * math.pow(\
k,5))
  return 1 - erf(-x)

def d1(S,K,sigma,r,t):
  v1 = math.log(float(S)/K)
  v2 = (r + (sigma*sigma)/2.0) * t
  v3 = sigma * float(math.sqrt(t))
  return (v1+v2)/float(v3)
  

def d2(S,K,sigma,r,t):
  return d1(S,K,sigma,r,t) - sigma * math.sqrt(t)
  
def call(S,K,sigma,r,t):
  v1 = d1(S,K,sigma,r,t)
  v2 = d2(S,K,sigma,r,t)
  return S * math.exp(-q * t)*erf(v1) - K * math.exp(-r*t) * erf(v2)

  
def put(S,K,sigma,r,t):
  v1 = d1(S,K,sigma,r,t)
  v2 = d2(S,K,sigma,r,t)
  return K * math.exp(-r*t) * erf(-v2) - S * erf(-v1) * math. exp(-q * t)
  


def bsformula(callput, S, K, r, t, sigma, q = 0):
   vega = S * (math.sqrt(t)) *( math.exp( -0.5 * d1(S,K,sigma,r,t) **2)) /((2 * math.pi)**0.5)
   if callput == "Call":
        p = call(S, K, sigma, r, t)
        delta = erf(d1(S, K, sigma, r, t)) 
   elif callput == "Put":
        p = put(S,K, sigma, r,t)
        delta = erf(d1(S, K, sigma, r, t)) -1
   return (p, delta, vega)
   

   
def Newtown(tar, x, tol, callput, S , k, r, t,sigma):
    sigma = x
    a = bsformula(callput, S, k,r,t,sigma)
    y = float(a[0])
    xs = [x]
    count = 0
    while abs(y - tar) >= tol and count <= 100:
        dx = S * (math.sqrt(t)) *( math.exp( -0.5 * d1(S,k,sigma,r,t) **2)) /((2 * math.pi)**0.5)
        x = x + (tar - y ) / dx
        y = float(bsformula(callput,S, k, r, t, x)[0])
        sigma = x
        xs.append(x)
        count += 1
    return xs

def bsimpvol( callput, S0, K, r, t, price, q = 0.0, priceTolerance = 0.01, method = 'bisect', reportcalls = False):
    def bs_1(sigma):
        a  = bsformula(callput, S0, K, r, t, sigma, q = 0)
        b = float(a[0])
        return b
    sigmaTry = 0.5
    test = S0 + K + r + t + price + q + priceTolerance 
    reportcalls = not(math.isnan(test))
    reportcalls = True

    if (callput == "Call") and price < (max(S0 - K, 0)):
        reportcalls = False

    if (callput == "Put") and price < max (K - S0, 0 ):
        reportcalls = False
    if reportcalls == True :
        if method == 'newtown':
            y = Newtown( price, sigmaTry, priceTolerance,callput, S0,K,r,t,sigma)
        if method == 'bisect':
            y = Bisect.bisect(price, bs_1, sigmaTry, None, [ 0.01, priceTolerance ],100)
        count = len(y)
        vol = y[-1]
        return ( vol, count)
    else:
        return "NaN"
        
ins = open( "TSLAOptions.csv", "rb")
reader = csv.reader( ins, dialect = "excel")

rn = 0
today = datetime.date(2013,8,15)
for row in reader:
    if rn == 0:
        pass
    if 0 < rn < 105:
        #read in
        pb = float(row[4])
        pa = float(row[5])
        comp = float(row [6])
        yie = float(row[7])
        Stock = float(row[8])
        rate = float(row[9]) /100.
        strike = float(row[10])
        callPut = row[11]
        pr = 0.5 * ( pa + pb )
        #for date
        de = row [12]
        end = parser.parse(de)
        a = datetime.datetime.date(end)
        dt = a - today
        time = float(dt.days/365.0)
        #for date end
        res_NT = bsimpvol(callPut, Stock, strike, rate, time, pr, yie, 0.01, 'newtown', False)
        res_BS = bsimpvol( callPut, Stock, strike, rate, time, pr, yie, 0.01, 'bisect', False)
        print "bisect:", res_BS
        print "newton:" ,res_NT

    rn += 1 
