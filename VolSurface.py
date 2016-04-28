import math
import numpy as np
import csv
import datetime
from dateutil import parser
from scipy.optimize import newton
from scipy.optimize import fmin_bfgs
from scipy.stats import norm

def d1(S,K,sigma,r,t,q):
  v1 = math.log(float(S)/K)
  v2 = (r -q + (sigma*sigma)/2.0) * t
  v3 = sigma * float(math.sqrt(t))
  return (v1+v2)/float(v3)
  

def d2(S,K,sigma,r,t,q):
  return d1(S,K,sigma,r,t,q) - sigma * math.sqrt(t)
  
def call(S,K,sigma,r,t,q):
  v1 = d1(S,K,sigma,r,t,q)
  v2 = d2(S,K,sigma,r,t,q)
  return S * math.exp(-q * t)*norm.cdf(v1) - K * math.exp(-r*t) * norm.cdf(v2)

  
def put(S,K,sigma,r,t,q):
  v1 = d1(S,K,sigma,r,t,q)
  v2 = d2(S,K,sigma,r,t,q)
  return K * math.exp(-r*t) * norm.cdf(-v2) - S * norm.cdf(-v1) * math. exp(-q * t)
  
def bsPrice(cp, S, K, sigma, r, t, q):
    if cp == 1: return call(S, K, sigma, r, t, q)
    if cp == -1: return put(S, K, sigma, r, t, q)

#print call(25, 27.5, 0.05167,0.06, 0.0390, .4)
#print call(30.5, 30, 5, 0.04, 1.1, q=0.0040)
  ## Problem 1
def bsDelta(callput, S, K, r, T, sigma, q = 0., t = 0):
    if callput == 1:
        delta = math.exp(-q *(T-t))*norm.cdf(d1(S, K, sigma, r, T-t, q)) 
    elif callput == -1:
        delta = (norm.cdf(d1(S, K, sigma, r, T-t, q)) -1)* math.exp(-q * (T-t))
    return delta
   
 ## Problem 2  
def bsGamma(callput, S, K, r, T, sigma, q = 0.,t = 0):
    x = d1(S, K, sigma, r, T-t, q)
    Nd1 = 1.0/math.sqrt(2 * math.pi) * math.exp(-x*x / 2.0)
    gamma = Nd1 / (S * sigma * math.sqrt(T))
    return math.exp(-q * (T-t))*gamma
 
#print bsDelta(-1, 30, 30, 0, 4, 1, q=0.0040, t=0)
## Problem 3
def impliedBS( V, callput, S, K, R, T, q = 0., t = 0):
    tau = T - t
    def BSsigma(sigma):
        if callput == 1.0:
            #print 7
            p = call(S, K, sigma, R, tau, q)
        if callput == -1.0:
            #print 8
            p = put(S, K, sigma, R, tau, q)
        return p - V
    def vega(sigma):
         return math.exp(-q * tau) * S * (math.sqrt(tau)) *( math.exp( -0.5 * d1(S,K,sigma,R,tau,q) **2)) /((2 * math.pi)**0.5)
    if callput == 1 and (S - K < V):
        #print "hehe"
        #print BSsigma(0.0003), BSsigma(4)
        Vol = newton(BSsigma, 1, vega, (),1.e-08, 100, None)
    elif callput == -1 and (K - S < V) :Vol = newton(BSsigma, 1, vega, (), 1.e-08, 100,None)
    else: 
        return np.nan
    #rint "impliedBS",V, callput, S, K, R, T, q, t,"ig",Vol/math.sqrt(365)
    return Vol
    

    
#print impliedBS(4.4, 1, 30.5, 30, 0.04, T=1.1, q=0.0040)

## Problem 4
def sign ( x ):
    if x < 0:
        return -1
    elif x >= 0:
        return 1
        
#print sign(-3)
def PeizerPratt(z, n):
    return CumuDisAprox (z,n)

def CumuDisAprox ( z, n ):
    f = 0.5 + sign(z) * math.sqrt(0.25-0.25 * math.exp(-((z/(n+1/3.0))**2)*(n+1/6.0)))
    return f
    
#print CumuDisAprox(0.85, 5)

def lrtree(callput, S, K, r, T, sigma, q = 0., t = 0, params = {'stepCount':199}):
    n = params['stepCount']
    #print n
    tau = T - t

    if n % 2 == 0:
        n += 1
    dt = (float(tau)) / n
    dPlus = d1(S,K,sigma,r,tau,q)
    dMinus = d2(S, K, sigma, r, tau,q)
    #print "biubiu", dPlus, dMinus
    pBar = CumuDisAprox( dPlus, n)
    p = CumuDisAprox(dMinus, n)
    #print "momo", pBar, p
    u = math.exp( (r-q) * dt) * pBar / p
    d = (math.exp((r-q) * dt) - p * u) / (1- p)
    #print u,d
    
    price_tree = np.zeros([n+3, n+3])
    price_tree[0, 0] = S / u / d
    
    for i in range (1, n+3):
        price_tree[ 0 : i , i ] = price_tree[ 0 : i , i - 1] * u
        price_tree[i, i] = price_tree[ i-1,i-1 ] * d
    #print price_tree
    
    value_tree = np.zeros([n+3, n+3])
    if callput == 1:
        for j in range (0, n+2) : #indicce from 0 to n-1
            value_tree[j, n +1] = call (price_tree[ j,n+1], K, sigma, r, dt, q)
    if callput == -1:
        for j in range (0, n+2) : #indicce from 0 to n-1
            value_tree[ j, n+1] = put (price_tree[ j, n+1], K, sigma, r, dt, q)
    #print value_tree
    ##discoutn
    for i in range ( n +1, -1, -1):
        if i == n+1:
            if callput == 1:
                for j in range (0, i+1):
                    value_tree[j,i] = max(max(price_tree[j,i] -K, 0), value_tree[j,i])
                    #print i, value_tree[j,i], price_tree[j, i]-K
            elif callput == -1:
                #print "go"
                for j in range(0,i+1):
                    #nt value_tree[j, n+1], K-price_tree[j,i] 
                    value_tree[j,i] = max(max(K- price_tree[j, i], 0), value_tree[j,i])
                    #print value_tree[j,i], K-price_tree[j,i]
        if i <= n:
            value_tree [0 : i +1, i] = math.exp(-r * dt)* (p * value_tree[0: i+1, i +1] + (1-p) * value_tree[ 1:i+2, i+1])            
            if callput == 1:
                for j in range (0, i+1):
                    value_tree[ j, i ] = max(max(price_tree[ j, i ] - K,0), value_tree[ j, i])
            elif callput == -1:
                for j in range ( 0, i +1 ):
                    #print value_tree[j, i],K- price_tree[j,i] 
                    value_tree[ j, i ] = max(max(K - price_tree[ j, i ],0), value_tree[ j, i])
    #print value_tree
    p = value_tree[ 1, 2 ]
    delta = (value_tree[0,2] - value_tree[2,2]) / (S * ( u / d- d / u ))
    vu = value_tree[0,2]
    vm = value_tree[1,2]
    vd = value_tree[2,2]
    Su = S / d * u
    Sd = S / u * d
    gamma= ((vu - vm) / (Su - S) - (vm - vd) / (S - Sd)) / (Su - Sd) * 2
    tuple = (p, delta, gamma)
    #print tuple[0]
        
    return tuple
#print lrtree(-1, 25, 27, 0.06, .4, .65, q=0.0390, t=0,params={'stepCount':3})
    
##Problem 5,6
def lrtreeGamma(callput, S, K, r, T, sigma, q = 0., t = 0, params = {'stepCount':200}):
    a = lrtree(callput, S, K, r, T, sigma, q , t , params)
    return a[2]
    
def lrtreeDelta(callput, S, K, r, T, sigma, q =0., t = 0, params = {'stepCount':200}):
    a = lrtree(callput, S, K, r, T, sigma, q , t, params)
    return a[1]
    
#print lrtree(1, 25, 27.5, 0.06,.4,  0.0390, 0, params={'stepCount':3})

def lrtreePrice(callput, S, K, r, T, sigma, q =0., t = 0, params = {'stepCount':200}):
    return float(lrtree(callput, S, K, r, T, sigma, q, t, params)[0])

#~ ## Problem 7
def impliedTree( V, callput, S, K, r, T, q = 0., t = 0, params = {'stepCount' : 200}):
    #print params 
    #print V, callput, S, K, r, T, q
    def Tree( sigma ):
        a = lrtree( callput, S, K, r*365, T/365., sigma, q*365, t, params)
        return float(a[0]) - V
    if callput == 1 and S - K < V: 
        #print "line 188"
        vol = newton(Tree, 0.5,None,() , 1.e-08, 100, None)
    elif callput == -1 and K - S < V: 
        vol = newton(Tree, 0.5,None,() , 1.e-08, 100, None)
    else: 
        return np.nan
    #vol = brentq(Tree, 0.003, 4, args = (), xtol = 1.e-08)
    return vol / math.sqrt(365)
    
## Problem 8
def isnumber(x):
    try :
        float (x)
        return True
    except ValueError:
        return False
        
def isCP(x):
    if x == "Call":
        return 1
    elif x == "Put":
        return -1
        
def gettime(s):
        formatDate = parser.parse(s)
        d = datetime.datetime.date(formatDate)
        today = datetime.date(2013, 9,15)
        t = d - today
        dt = t.days
        return dt

def volParse( file ):
    
    input = open( file, "rb")
    reader = csv.reader( input, dialect = 'excel')
    
    nrow = sum(1 for line in input)
    #print nrow
    input.seek(0)
    reader = csv.reader(input, dialect = 'excel')
    chart = zeros ([ 84 - 2, 8])
    today = datetime.date(2013,9,15)
    rn = 0
    for row in reader:
        i = rn - 2
        if rn <= 1:
            pass
        if rn > 1 and  rn < 80+2:
            if isnumber(row[2]): #bid
                chart[i, 6] = float(row[2])
            
            if isnumber(row[3]): #ask
                chart[i, 7] = float(row[3])
            chart[i, 4] = float(row[4]) #q
            chart[i, 1] = float(row[5]) #S
            chart[i, 3] = float(row[6]) /100#r
            chart[i, 2] = float(row[7]) #K
            chart[i, 0] = isCP(row[8])
            chart[i, 5] = gettime(row[9])
        rn += 1
    #print chart
    return chart
#print volParse("BA.csv")

##Problem 9
def volImply( parsed, impliedVolFcn=impliedBS, params = None) :
            
    chart = parsed
    s = np.size(chart)
    line = s / 8
    #print line
    vols = np.zeros([line,3])
    d = math.sqrt(365.)
    if impliedVolFcn == impliedBS:
        for i in range (0, line):
            #print i
             # bid vol
            if chart[i, 6] == None:
                vols[i, 0] = None
            else:
                #print i
                vols[i,0] = impliedVolFcn( chart[i,6], #V
                chart[i,0], #callput 
                chart[i,1],#S
                chart[i,2],  #K, 
                chart[i,3] *365, #R, 
                float(chart[i, 5]) / 365. ,# T, 
                chart[i, 4] * 365 , #q = 0., 
                t = 0) /d
                #print "test2",chart[i,6],  chart[i,0], chart[i,1],chart[i,2],  chart[i,3] , float(chart[i, 5]),chart[i, 4] 
            # ask vol
            if chart [ i, 7] == None: vols[i, 2] = None
            else:
                vols[i, 2] = impliedVolFcn( chart[i,7], #V
                chart[i,0],#callput 
                chart[i,1],#S
                chart[i,2],  #K, 
                chart[i,3] * 365, #R, 
                float(chart[i, 5]) / 365. ,# T, 
                chart[i, 4] * 365, #q = 0., 
                t = 0) / d
                #~ print chart[i,7]
                #~ print "cp", chart[i,0]
                #~ print "s", chart[i,1]
                #print impliedBS(5.25,-1, 25, 27.5, 0.06/365., .4*365,  0.0390/365.)
            # mid vol
            if chart[i,6] != None and chart[i,7] != None:
                mid = (chart[i,6] + chart[i,7]) / 2
                vols[i,1] = impliedVolFcn( mid, #V
                chart[i,0],#callput 
                chart[i,1],#S
                chart[i,2],  #K, 
                chart[i,3] * 365, #R, 
                float(chart[i, 5])/365.,# T, 
                chart[i, 4] * 365, #q = 0., 
                t = 0) / d
            else: vols[i,1] = None

    elif impliedVolFcn == impliedTree:
        for i in range (0, line):
            print i
             # bid vol
            if chart[i, 6] == None:
                vols[i, 0] = None
            else:
                #print i
                vols[i,0] = impliedVolFcn( chart[i,6], #V
                chart[i,0], #callput 
                chart[i,1],#S
                chart[i,2],  #K, 
                chart[i,3] , #R, 
                float(chart[i, 5])  ,# T, 
                chart[i, 4] , #q = 0., 
                0 , params) 
                #print "test2",chart[i,6],  chart[i,0], chart[i,1],chart[i,2],  chart[i,3] , float(chart[i, 5]),chart[i, 4] 
            # ask vol
            if chart [ i, 7] == None: vols[i, 1] = None
            else:
                vols[i, 2] = impliedVolFcn( chart[i,7], #V
                chart[i,0],#callput 
                chart[i,1],#S
                chart[i,2],  #K, 
                chart[i,3] , #R, 
                float(chart[i, 5])  ,# T, 
                chart[i, 4] , #q = 0., 
                0, params) 
                #~ print chart[i,7]
                #~ print "cp", chart[i,0]
                #~ print "s", chart[i,1]
                #print impliedBS(5.25,-1, 25, 27.5, 0.06/365., .4*365,  0.0390/365.)
            # mid vol
            if chart[i,7] != None and chart[i,6] != None and vols[i, 0] > 0 and vols[i, 2] > 0:
                print vols[i, 0], vols[i,2]
                mid = (chart[i,6] + chart[i,7]) / 2
                vols[i,1] = impliedVolFcn( mid, #V
                chart[i,0],#callput 
                chart[i,1],#S
                chart[i,2],  #K, 
                chart[i,3] , #R, 
                float(chart[i, 5]),# T, 
                chart[i, 4] , #q = 0., 
                0, params) 
            else: vols[i, 1] = np.nan
    big_chart = np.zeros([line, 11])
    for i in range (0, line):
        for j in range (0, 6):
            big_chart[i,j] = chart[i,j]
        for j in range (6, 9):
            if vols[i, j-6] != 0:big_chart[i,j] = vols[i, j-6]
            else: big_chart[i, j] = np.nan
        for j in range (9, 11):
            if chart[i, j-3] != None: big_chart[i, j] = chart[i, j-3]
            else:big_chart[i, j] = np.nan

    return big_chart
    
parsed = np.array(
                  [
                    [-1, 25, 27.5, 0.06/365., 0.0390/365., .4*365, 5.25, 5.35],
                    [ 1, 25, 24, 0.06/365., 0.0390/365., .4*365, 2.25, 0.1],
                    [ 1, 25, 20, 0.06/365., 0.0390/365., .4*365, 0.1, 5.35],
                  ]
                  )
print volImply(parsed, impliedVolFcn=impliedTree,params={'stepCount':3})
#print volImply(parsed, impliedVolFcn=impliedTree,params={'stepCount':3})

def forwardDiff(S, t ,r, q, K):
    return S * math.exp((r-q) * t) -K

#Problem 10
def atmfVols(implied, atmfType = 'mid'):
    a = implied
    last = 0
    for row in a:
        last += 1
    #print "fa",last
    r = {}
    forward = {}
    volatility = {}
    #print r
    # end sorting
    tenor = a[0,5]
    begin = 0
    
    for i in range ( 0, last ):
        #print i
        DiffLeft = 1000
        DiffRight = -1000
        if a[i, 5] == tenor and i != last- 1:
            pass
        else:
            fw = forwardDiff (a[begin, 1], a[begin,5], a[begin,3], a[begin, 4], a[begin,2]) + a[begin, 2]
            #print"sou", i
            left = -1
            right = -1
            if begin == i:
                dis = forwardDiff (a[i, 1], a[i,5], a[i,3], a[i, 4], a[i,2])
                if dis> 0: 
                    left = i
                    DiffLeft = dis
                if dis < 0: 
                    right = i
                    DiffRight = dis
                    #print "right"
            if i == last - 1: end = i+1
            else: end = i
            for m in range(begin, end):
                #print "go", last, i, begin, end, m
                diff = forwardDiff(a[m, 1], a[m,5], a[m,3], a[m, 4], a[m,2]) 
                #if m == 5: print "??",diff
                if diff > 0 and diff < DiffLeft :
                    left = m
                    DiffLeft = diff
                    #print "left", diff
                elif diff < 0 and diff > DiffRight:
                    right = m
                    DiffRight = diff
                    #print "right",diff
            #print left, right,"puio"
            #print "lr",DiffLeft, DiffRight
            Strike1=[-1,-1]
            Strike2 = [-1,-1]
            
            #print DiffLeft,DiffRight,"jijoij"
            # pl and pr        
            if DiffLeft == 1000: 
                #print "ko"
                pl = 0
                Strike2[0] = right
            elif DiffRight == -1000: 
                pl = 1
                Strike1[0] = left
            else:
                pl = -DiffRight / (DiffLeft - DiffRight)
                Strike1[0] = left
                Strike2[0] = right
            pr = 1 - pl
            
            print pl, pr, "jiji"
            
            #print "miao",left, right
            
            # if there are 2 options for 1 K
            
            left2 = -1
            right2 = -1
            for m in range (begin, i):
                if a[m, 0] == -a[left, 0] and a[m, 2] == a[left, 2]:
                    left2 = m
                if a[m, 0] == -a[right, 0] and a[m, 2] == a[right, 2]:
                    right2 = m
            if left2 != -1: Strike1[1] = left2
            if right2 != -1: Strike2[1] = right2
            # usable strikes found  

            if atmfType == 'mid':
                sumleft = 0
                numberL = 0
                for j in range(0,2):
                    c = a[Strike1[j], 7]
                    if c > 0 and Strike1[j] != -1: sumleft += c
                    numberL += float(Strike1[j] != -1) * (float(c > 0))
                sumright = 0
                numberR = 0
                for j in range (0,2):
                    c = a[Strike2[j], 7]
                    if c > 0 and Strike2[j] != -1: sumright += c
                    numberR += float(Strike2[j] != -1) * (float(c > 0))
            
            if atmfType == 'spread':
                sumleft = 0
                numberL = 0
                #print "str",Strike1,Strike2
                for j in range (0,2):
                    #print "loop", j
                    #print "noodle",Strike1[j]
                    a1 = a[Strike1[j], 8]
                    b1 = a[Strike1[j], 6]
                    if a1 > 0 and Strike1[j] != -1: sumleft += a1
                    if b1 > 0 and Strike1[j] != -1: sumleft += b1
                    #print "woo",a1,b1
                    #sumleft += (Strike1[i] == -1) *( (a1== None) * a1+(b1 == None) * b1)
                    #print "nono",float(Strike1[j] == -1)
                    numberL += float(Strike1[j] != -1) * (float(a1 > 0) + float(b1 > 0))
                    #print "numl",numberL
                    #print "sum", sumleft
                sumright = 0
                numberR = 0
                for j in range (0,2):
                    a1 = a[Strike2[j], 8]
                    b1 = a[Strike2[j], 6]
                    if a1 > 0 and Strike2[j] != -1: sumright += a1
                    if b1 > 0 and Strike2[j] != -1: sumright += b1
                    #sumright += (Strike2[i] == -1) *( (a1== None) * a1+(b1 == None) * b1)
                    numberR += float(Strike2[j] != -1) * (float(a1 > 0) + float(b1 > 0))
                    #print "number", numberR, sumright
            #print "final", pl, sumleft, pr, sumright, numberL, numberR
            sigma_atm =   (pl*sumleft + pr * sumright) / (numberL * pl + numberR * pr)
            volatility[tenor] = sigma_atm
            forward[tenor] = fw
            r['Forwards'] = forward
            r['ATMF Volatilities'] = volatility
            tenor = a[i,5]
            #print "change", i
            begin = i
            #print "begin", i
        
    return r
    
implied = np.array(
                  [
                    [ 1, 25, 20, 0.06/365., 0.0390/365., .2*365,
                            None, None,   1.69823669e-02, 0.1, 5.35],
                    [-1, 25, 27.5, 0.06/365., 0.0390/365., .2*365,
                            .015,   .032,   .035, 5.25, 5.35],
                    [ 1, 25, 24, 0.06/365., 0.0390/365., .2*365,
                            0.01,   .011, .09, 2.25, 0.1],
                    [ 1, 25, 20, 0.06/365., 0.0390/365., .4*365,
                            None, None,   1.69823669e-02, 0.1, 5.35],
                    [-1, 25, 27.5, 0.06/365., 0.0390/365., .4*365,
                            .015,   .032,   .035, 5.25, 5.35],
                    [ 1, 25, 24, 0.06/365., 0.0390/365., .4*365,
                            0.01,   .011, .09, 2.25, 0.1],
                    [-1, 25, 27.5, 0.06/365., 0.0390/365., 290.1,
                            .015,   .032,   .035, 5.25, 5.35],
                    [ 1, 25, 24, 0.06/365., 0.0390/365., 290.1,
                            0.01,   .011, .09, 2.25, 0.1],
                    [ 1, 25, 20, 0.02/365., 0.01/365., 290.1,
                            None, None,   1.69823669e-02, 0.1, 5.35],
                    [ -1, 25, 20, 0.02/365., 0.01/365., 290.1,
                            0.01, .011,   1.69823669e-02, 0.1, 5.35],
                    [ 1, 25, 30, 0.02/365., 0.01/365., 290.1,
                            0.02, .031,   1.69823669e-02, 0.1, 5.35],
                    [ -1, 25, 30, 0.02/365., 0.01/365., 290.1,
                            0.04, .041,   1.69823669e-02, 0.1, 5.35],
                  ]
                  )
#print atmfVols(implied,'spread')

#print type(atmfVols(implied,'mid'))
#print atmfVols(implied,'spread')    
#print atmfVols(volImply("BA.csv"))

##Promblem 11
def weightedData(implied, atmfvols, forwardDic, weightType = 'delta', singleSidePenalty = 0.25, deltaFunc = bsDelta, GammaFunc = bsGamma, params = {None}):
    a = implied
    bigger_chart = []
    def judge(x):
        #print "ewer", x
        if weightType == 'delta':
            #print "go"
            if x[0] == 1:
                #print "gogo"
                weight = 1 - deltaFunc(x[0], x[1], x[2], x[3], x[5], atmfvols[x[5]], q = x[4], t = 0)
                #print weight
            if x[0] == -1:
                weight = 1 + deltaFunc(x[0], x[1], x[2], x[3], x[5], atmfvols[x[5]], q = x[4], t = 0)
        if weightType == 'gamma':
            weight = GammaFunc(x[0], x[1], x[2], x[3], x[5], atmfvols[x[5]], x[4], t = 0)
        if weightType == 'constant': weight = 1.
        #if weightType == 'moneyness': weight = g(x[1], x[2], x[3], x[5], x[4], atmfvols[x[5]])
        return weight
    for row in a:
        #print row
        if row[6] != None: 
            if row[8] != None:
                w = judge(row) 
            if row[8] == None:
                w = judge(row)  * singleSidePenalty
            line = []
            line[0:6] = row[0:6]
            line.extend([row[9], forwardDic[row[5]], atmfvols[row[5]],row[6],w])
            bigger_chart. append(line)
        if row[8] != None : 
            if row[6] != None:  w = judge(row) 
            if row[6] == None:  w = judge(row) * singleSidePenalty
            line = []
            #print w
            line = []
            line[0:6] = row[0:6]
            line.extend([row[10], forwardDic[row[5]], atmfvols[row[5]],row[8],w])
            bigger_chart. append(line)
    #~ for i in range (1, line):
            #~ if a[i,6] ==0: 
                #~ S = singleSidePenalty * a[i,2]
                #~ vols[i] = impliedBS(S, a[i,0], a[i,1], a[i,2], a[i,3], a[i,5], a[i,4])
            #~ else:vols[i] = atmfvols[a[i,5]]
    return np.array(bigger_chart)

#print weightedData(implied, atmfV)

#input = volImply("BA.csv")
#print weightedData(input, atmfVols(input)) 

##Problem 12
def fit(F, K, s_atm, Xi, Kappa, w, T):
    K_hat = 10. * w * np.arctan(math.log(F/ K)/ (w * math.sqrt(T)))
    #print K_hat
    #print "praras", Xi, Kappa, s_atm
    sigma_fitted = s_atm * ( 1 + Xi * K_hat + Kappa * K_hat**2)
    return sigma_fitted
    
def fitErrors(wdata, Xi, Kappa, w, fitType, pricingFunc = lrtree, errorResponseFunc = np.abs, params = None):
    sum = 0
    for row in wdata :
        F = float(row[7])
        K = float(row[2])
        T = float(row[5] )
        true_sigma = row[9]
        atm = row[8]
        s_fit = fit(F, K, atm, Xi, Kappa, w, T)
        print s_fit
        #print F,K,T,atm,w,s_fit
        if fitType == 'price':
            #print pricingFunc, "noiouoi"
            if pricingFunc == lrtreePrice: 
                price_fit = lrtreePrice(row[0], row[1], K, row[3], T, s_fit, row[4], 0, params)
            if pricingFunc == bsPrice: price_fit = pricingFunc(row[0], row[1], K, s_fit, row[3], T, row[4])
            error = errorResponseFunc(row[6] - price_fit) * row[10]
            sum += error
        if fitType == 'vol':
            error = errorResponseFunc(s_fit - true_sigma) * row[10]
            sum += error
    return sum
    
parsed = np.array(
                  [
                    [-1, 25, 27.5, 0.06/365., 0.0390/365., .4*365, 5.25, 5.35],
                    [ 1, 25, 26, 0.06/365., 0.0390/365., .4*365, 2.25, None],
                    [ 1, 25, 20, 0.06/365., 0.0390/365., .4*365, None, 5.35],
                  ]
                  )
#print volImply(parsed, impliedVolFcn=impliedTree,params={'stepCount':3})

##Problem 13
def fitWeightData(wdata, fitType, pricingFunc = lrtree, errorResponseFunc = np.abs, params = None):
    def targetFunc(x):
        return fitErrors(wdata, x[0], x[1], x[2],fitType, pricingFunc, errorResponseFunc, params)
    return fmin_bfgs(targetFunc, [0,0.8,0.8])
wdata = np.array([[1, 25, 20, 0.00016438356164383562, 0.00010684931506849315, 73.0,
                            5.35, 108, 0.03, 0.0169823669, 0.036901167529658427],
                           [-1, 25, 27.5, 0.00016438356164383562, 0.00010684931506849315, 73.0,
                            5.25, 108, 0.03, 0.015, 0.060197912638701861],
                           [-1, 25, 27.5, 0.00016438356164383562, 0.00010684931506849315, 73.0,
                            5.35, 108, 0.03, 0.035, 0.060197912638701861]])
fitType = 'price'
pricingFunction = lrtreePrice
#print  fitWeightData(wdata, fitType, pricingFunction, params={'stepCount':3})
