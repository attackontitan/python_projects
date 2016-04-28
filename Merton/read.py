from __future__ import division
import math
from numpy import *
import numpy as np
import csv
from datetime import datetime

def UseCsv( filename):
    file = open(filename,'rb')
    table = csv.reader(file)
    datatable = [row for row in table]
    starttime = datetime.strptime('11/23/2011','%m/%d/%Y')
    header = ['Option', 'Stock', 'Strike', 'Days to expiry', 'interest rate','optionprice']
    output = []
    for i  in range(1,len(datatable)):
        call=1
        put=-1
        S=datatable[i][12]
        K=datatable[i][11]
        tau=datatable[i][9]
        r=datatable[i][8]
        cprice=datatable[i][4]
        pprice=datatable[i][5]        
        if cprice=='#N/A N/A':
            cprice=np.nan
        else:
            if S<K: #only out-of-money call
                copt=[float(call),float(S),float(K),float(tau)/365,round(float(r)/100,4),float(cprice)]
                output.append(copt)
                #print copt
        if pprice=='#N/A N/A':
            pprice=np.nan
        else:
            if S>K: #only out-of-money put
                popt=[float(put),float(S),float(K),float(tau)/365,round(float(r)/100,4),float(pprice)]
                output.append(popt)
                #print popt
    return array(output)
    
#print    UseCsv('9days.csv')

