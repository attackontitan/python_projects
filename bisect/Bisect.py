
import math
from numpy import array


    
def NBseekbound(s,x,t,tarFun):#s = start ,x = tol = 0.01
    left = s-x
    right = s+x
    yLeft = tarFun (left)
    yRight = tarFun( right)
    while (yLeft-t) * (yRight-t) >= 0:
        left = 2 * (left-s) + s
        right = 2 * (right-s) + s
        yLeft =tarFun (left)
        yRight = tarFun (right)
    bound = [left, right]
    return bound

def test(u,l,t,tarFun):
    yu = tarFun(u)
    yl = tarFun(l)
    if (yu-t)*(yl-t) <= 0:
        return 1
    else:
        return 0
        
def find(l,r,t,tarFun):#left and right
    halfwaypt = 1/2.0* (l+r)
    yl = tarFun(l)
    yr = tarFun(r)
    yh = tarFun(halfwaypt)
    productr = ( yl - t ) * (yh - t)
    productl = ( yr - t ) * ( yh - t)
    if yh == t:
        bracket = [l,r,halfwaypt]
        return bracket
    if productr < 0:
        r = halfwaypt
    if productl < 0:
        l = halfwaypt
    halfway = 1/2.0* (l+r)
    bracket = [l , r , halfwaypt]
    return bracket
    
def bisect(target, targetfunction, start = None, bounds = None, tols = [0.01, 0.01], maxiter = 100):
    flag = 0
    if (bounds == None)and (start != None):
        flag = 1
        b = NBseekbound(start,tols[0],target,targetfunction)
        upper = b[1] 
        lower = b[0]
        left = lower
        right = upper
    elif (start == None)and (bounds != None):
        upper = bounds[1] 
        lower = bounds[0]
        start = 1/2.0*(upper+lower)
    elif(start != None) and (bounds != None):
        upper = bounds[1]
        lower = bounds[0]
        if test( start, lower, target, targetfunction ) == 1:
            upper = start
        if test( upper, start, target, targetfunction ) == 1:
            lower = start
    else:
        raise Exception ('cannot calculate')
    # find bound and start end
    left = lower
    right = upper
    yleft = targetfunction(left) 
    yright = targetfunction(right)
    if test(upper, lower, target, targetfunction) ==0:
        raise Exception('no solution in bounds')
    else:
        xs= [ left, right, 1/2.0 *(left + right)]
        tryvalue = xs[2]
        if abs (yleft - target) <  tols[1] and ( flag == 1 ):
            xs.append(left)
            return xs
        for m in range (2, 100) :
            y = targetfunction(tryvalue)
            if (abs(y - target)) > tols[1]:
                next = find (left,right, target, targetfunction)
                xs.append(next[2])
                left = next[0]
                right = next[1]
                tryvalue = xs[m]
                if m == 97:
                    raise Exception('max iteration')
            elif abs(tryvalue-target) < tols[1]:
                break
            
        del xs[2]
        del xs[-1]
        return xs

