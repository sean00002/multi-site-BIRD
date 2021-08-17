#====================================================================
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
# This is OPEN SOURxCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
#====================================================================
from __future__ import (absolute_import, division, print_function,     
   unicode_literals, generators, nested_scopes, with_statement)       
from builtins import (bytes, dict, int, list, object, range, str, ascii, 
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
import math

######################################################################
# Attributes:
#    
# Methods:
#    [mean,SD,min,max]=SummaryStats.summaryStats(array)
#    [mean,SD,min,max]=SummaryStats.roundedSummaryStats(array)
#    sum=SummaryStats.sum(array)
#    r=SummaryStats.correlation(array1,array2)
######################################################################

class SummaryStats:
    """SummaryStats computes simple mean, variance, and correlation
       statistics
    """

    @classmethod
    def summaryStats(self,array):
        n=len(array)
        minX=None
        maxX=None
        sumX=0
        sumXX=0
        for i in range(0,n):
            x=array[i]
            sumX+=x
            sumXX+=x*x
            if(i==0): minX=maxX=x
            if(x<minX): minX=x
            if(x>maxX): maxX=x
        meanX=sumX/n
        varX=None if n<2 else (sumXX-sumX*sumX/n)/(n-1)
        if(varX is not None and varX<0): varX=0
        stddevX=math.sqrt(varX) if varX is not None else None
        return [meanX,stddevX,minX,maxX]

    @classmethod
    def roundedSummaryStats(self,array):
        [mean,stddev,min,max]=SummaryStats.summaryStats(array)
        mean=int(100.0*mean+5.0/9.0)/100.0
        stddev=int(100.0*stddev+5.0/9.0)/100.0
        min=int(100.0*min+5.0/9.0)/100.0
        max=int(100.0*max+5.0/9.0)/100.0
        return [mean,stddev,min,max]

    @classmethod
    def sum(self,array):
        s=0.0
        n=len(array)
        for i in range(0,n):
            s+=array[i]
        return s;

    @classmethod
    def correlation(self,Xs,Ys):
        sumX=0.0
        sumY=0.0
        sumXY=0.0
        sumXX=0.0
        sumYY=0.0
        n=len(Xs)
        for i in range(0,n):
            x=Xs[i]
            y=Ys[i]
            sumX+=x
            sumY+=y
            sumXY+=x*y
            sumXX+=x*x
            sumYY+=y*y
        r=(sumXY-sumX*sumY/n)/math.sqrt((sumXX-sumX*sumX/n)*(sumYY-sumY*sumY/n))
        return r

