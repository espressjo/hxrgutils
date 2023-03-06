#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 11:53:34 2020

@author: espressjo
"""
from numpy import fliplr,shape
from hxrgutils.stats.stats import stats

class hxrg_amp(stats):
    def __init__(self,data,namp=-1):
        stats.__init__(self,data)
        self.namp = namp
        self.data = data
        self.shape = shape(data)
        self.w = self.shape[1]
        self.h = self.shape[0]
    def ravel(self):
        if self.namp<0:
            print('[Warning] amp number is undefined')
            return self.data.ravel()
        if self.namp%2==0:
            return fliplr(self.data).ravel()
        return self.data.ravel()
    def __itruediv__(self,other):
        self.data = self.data/other
        return self
if '__main__' in __name__:
    from numpy import zeros
    data = zeros((24,24))+1
    print(data)
    amp = hxrg_amp(data)
    print(amp.stats('mean'))