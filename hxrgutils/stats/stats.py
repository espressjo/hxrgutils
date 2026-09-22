#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:08:57 2020

@author: espressjo
"""
from astropy.stats import sigma_clipped_stats as sc
from numpy import percentile,nanpercentile,isfinite,asarray,nan

class stats():
    def __init__(self,data):
        self.data = data
    def stats(self,stats=''):
        '''
        Sigma clipped mean, median and standard deviation. Non-finite pixels
        (NaN, inf) are dropped first. If no finite pixel is left we return NaN
        rather than letting numpy work on an empty array.
        '''
        data = asarray(self.data)
        if data.size<524300:
            sigma=5
        elif data.size<16777216:
            sigma = 5.7
        else:
            sigma = 6
        data = data[isfinite(data)]
        if data.size<1:
            print('[Warning] no finite pixel, statistics are undefined')
            if any(k in stats for k in ('std','median','mean')):
                return nan
            return (nan,nan,nan)
        if 'std' in stats:
            return sc(data,sigma=sigma)[2]
        elif 'median' in stats:
            return sc(data,sigma=sigma)[1]
        elif 'mean' in stats:
            return sc(data,sigma=sigma)[0]
        else:
            return sc(data,sigma=sigma)
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
      data = asarray(self.data).ravel()
      data = data[isfinite(data)]
      if data.size<1:
          print('[Warning] no finite pixel, nothing to plot')
          return f
      ma = nanpercentile(data,99.8)
      mi = nanpercentile(data,0.2)
      ax.hist(data,bins=150,range=(mi,ma),color='#87CEEB')
      ax.set_ylabel('Number of counts')
      ax.set_xlabel('Value')
      mn,md,st = self.stats()
      self.vline(ax,mn,label='mean',color='springgreen')
      self.vline(ax,md,label='median',color='peachpuff')
      self.vline(ax,mn-st,label=r'$-1 \sigma$',color='tomato')
      self.vline(ax,mn+st,label=r'$+1 \sigma$',color='tomato')
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
    
    