#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:08:57 2020

@author: espressjo
"""
from astropy.stats import sigma_clipped_stats as sc
from numpy import percentile,isnan

class stats():
    def __init__(self,data):
        self.data = data
    def stats(self,stats=''):
        if self.data.size<524300:
            sigma=5
        elif self.data.size<16777216:
            sigma = 5.7
        else:
            sigma = 6
        if 'std' in stats:
            return sc(self.data[~isnan(self.data)],sigma=sigma)[2]
        elif 'median' in stats:
            return sc(self.data[~isnan(self.data)],sigma=sigma)[1]
        elif 'mean' in stats:
            return sc(self.data[~isnan(self.data)],sigma=sigma)[0]
        else:
            return sc(self.data[~isnan(self.data)],sigma=sigma)
    def vline(self,ax,x,label='',color=''):
        lim = ax.get_ylim()
        ax.set_ylim(lim)
        ax.plot((x,x),lim,label=label,color=color)
    def dist(self,show=True):
      '''
      Plot the distribution of the data, will return fig
      and show==False will not show the graph
      '''
      from matplotlib import pyplot as plt
      plt.rcParams['axes.facecolor'] = '#242925'
      f,ax = plt.subplots(facecolor=(.31,.31,.31))
      ma = percentile(self.data,99.8)
      mi = percentile(self.data,0.2)
      ax.hist(self.data.ravel(),bins=150,range=(mi,ma),color='#87CEEB')
      ax.set_ylabel('Number of counts')
      ax.set_xlabel('Value')
      mn,md,st = self.stats()
      self.vline(ax,mn,label='mean',color='springgreen')
      self.vline(ax,md,label='median',color='peachpuff')
      self.vline(ax,mn-st,label='$-1 \sigma$',color='tomato')
      self.vline(ax,mn+st,label='$+1 \sigma$',color='tomato')
      print('Mean:\t %2.2f'%mn)
      print('Median:\t %2.2f'%md)
      print('Std.Dev: %2.2f'%st)
      
      plt.legend()
      plt.tight_layout()
      if show:
          plt.show()
      return f
if '__main__' in __name__:
    from astropy.io import fits
    data = fits.getdata('/home/noboru/NIRPS_R01_R01.fits')
    data = data[0:128,:]
    stats = stats(data)
    stats.dist()
    
    