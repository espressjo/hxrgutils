#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 11:53:18 2020

@author: espressjo
"""

from numpy import shape,asarray,zeros,mean,std
from scipy.ndimage import zoom
from astropy.io import fits
from astropy.stats import sigma_clipped_stats as sc
from os.path import basename,join
from os import getcwd
from hxrgutils.type.amp import hxrg_amp
from hxrgutils.fits2ramp.refpixcorr import refpxcorr
from hxrgutils.stats.stats import stats
from numpy import ndarray,isnan,isfinite,interp

class hxread(refpxcorr,stats):
    def __init__(self,x=4):
        self.x = x
        self.im = zeros((1024*x,1024*x),dtype=float)
        self.H = fits.Header
        self.shape = shape(self.im)
        self.w = self.shape[1]
        self.h = self.shape[0]
        self.fname = ''
        self.low_pixels = 0
        self.hi_pixels = 0
    def __call__(self,data,header=None):
        if isinstance(data,str):
            if '/' not in data:
                data = join(getcwd(),data)
            hdu = fits.open(data)        
            self.im = asarray(hdu[0].data,dtype=float)
            self.H = hdu[0].header
            self.fname = basename(data)
            hdu.close()
            self.shape = (self.H['NAXIS1'],self.H['NAXIS2'])
            self.w = self.shape[1]
            self.h = self.shape[0]
            self.x = int(self.w/1024)
        else:
            self.im = asarray(data,dtype=float)
            self.shape = data.shape
            self.w = self.shape[1]
            self.h = self.shape[0]
            self.x = int(self.w/1024)
            if header==None:
                self.H = fits.Header()
            else:
                self.H = header
            self.fname = "tmp.fits"
        stats.__init__(self,self.im)
        refpxcorr.__init__(self,self.im)
        #implement for data
    def get_bias(self,sigma=5):
        return sc(self.im,sigma=sigma)
    def refpx(self,**kwargs):
        
        if 'toponly' in kwargs and kwargs.get('toponly'):
            self.refpxcorrtop()
            self.H.add_comment('refpixcorr performed [top only]')
        elif 'sideonly' in kwargs and kwargs.get('sideonly'):
            self.refpxcorrside()
            self.H.add_comment('refpixcorr performed [side only]')
        else:
            self.refpxcorr()
            self.H.add_comment('refpixcorr performed')
    def writeto(self,fname):
        fits.PrimaryHDU(data=self.im,header=self.H).writeto(fname,overwrite=True)
    def ds9(self):
        fits.PrimaryHDU(data=self.im,header=self.H).writeto("/var/tmp/read.fits",overwrite=True)
        from os import popen 
        popen("ds9 -zscale /var/tmp/read.fits").read()
    def get_deviant_pixels(self,pixel_bin=32,low_threshold=0.5,hi_threshold=1.3):
        #pixel size of the binning box
        import numpy as np
        pix_bin = pixel_bin
        nbin = (self.x*1024)//pix_bin
        percentile = 50
        im = np.copy(self.im)
        box = np.zeros([nbin,nbin])
        for i in range(nbin):
            for j in range(nbin):
                # -1 sigma of distrubution. Should be good to remove illuminated pixels
                box[i,j] = np.nanpercentile(im[i*pix_bin:i*pix_bin+pix_bin,j*pix_bin:j*pix_bin+pix_bin],percentile)
        lowf = zoom(box,pix_bin)
        im = im/lowf        
        mask = np.zeros((np.shape(self.im)))
        #mask[im<0.5] = 1
        mask[im>hi_threshold] = 1
        mask[0:4,:] = 0
        mask[:,0:4] = 0
        mask[-4:,:] = 0
        mask[:,-4:] = 0
        rav = mask.ravel()
        self.hi_pixels = len(rav[rav==1])
        mask = np.zeros((np.shape(self.im)))
        mask[im<low_threshold] = 1
       
        mask[0:4,:] = 0
        mask[:,0:4] = 0
        mask[-4:,:] = 0
        mask[:,-4:] = 0
        rav = mask.ravel()
        self.low_pixels = len(rav[rav==1])
        return self.low_pixels,self.hi_pixels
    def get_amp(self,n):
        x_start = (n-1)*int(self.w/32)
        x_stop = (n)*int(self.w/32)
        return hxrg_amp(self.im[:,x_start:x_stop],namp=n)
    def writeto(self,fname):
        hdul = fits.PrimaryHDU(data=self.im,header=self.H)
        hdul.writeto(fname,overwrite=True)
    def __sub__(self,other):
        if isinstance(other,hxread):
            F = hxread()
            F(self.im-other.im)
            F.H = self.H
            newFile = "%s-%s"%(self.fname,other.fname)
            F.fname = newFile
            F.H.add_comment(newFile)
            return F
        else:
            temp = self.copy()
            temp.im = self.im-other
        return temp
    def __add__(self,other):
        if isinstance(other,hxread):
            F = hxread()
            F(self.im+other.im)
            F.H = self.H
            newFile = "%s+%s"%(self.fname,other.fname)
            F.fname = newFile
            F.H.add_comment(newFile)
            return F
        else:
            temp = self.copy()
            temp.im = self.im+other
        return temp
    def __imul__(self,other):
        if isinstance(other,hxread):
            self.im = self.im*other.im
        else:
            self.im = self.im*other  
        return self
    def __itruediv__(self,other):
        if isinstance(other,hxread):
            self.im = self.im/other.im
        else:
            self.im = self.im/other 
        return self
    def __iadd__(self,other):
        if isinstance(other,hxread):
            self.im = self.im+other.im
        else:
            self.im = self.im+other 
        return self
    def __isub__(self,other):
        if isinstance(other,hxread):
            self.im = self.im-other.im
        else:
            self.im = self.im-other 
        return self
    def __truediv__(self,other):
        if isinstance(other,hxread):
            F = hxread()
            F(self.im/other.im)
            F.H = self.H
            newFile = "%s/%s"%(self.fname,other.fname)
            F.fname = newFile
            F.H.add_comment(newFile)
            return F
        else:
            temp = self.copy()
            temp.im = self.im/other
        return temp
    def __mul__(self,other):
        if isinstance(other,hxread):
            F = hxread()
            F(self.im*other.im)
            F.H = self.H
            newFile = "%s*%s"%(self.fname,other.fname)
            F.fname = newFile
            F.H.add_comment(newFile)
            return F
        else:
            temp = self.copy()
            temp.im = self.im*other
        return temp
    def __pow__(self,other):
        temp = self.copy()
        for i in range(other-1):
            temp.im = temp.im*self.im
        return temp
    def show_fft(self,show=True,max_freq=200):
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_theme()

        f,ax = plt.subplots(figsize=(10,5))
        x,y = self.fft(max_freq)
        ax.plot(x,y)
        xx,yy = self.get_freq()
        if 'SEQID' in self.H:
            ax.set_title('Unique ID: %s'%self.H['SEQID'])
        ax.plot(xx,yy,'o',markersize=5)
        for i in range(len(xx)):
            if i>7:
                break
            ax.text(xx[i],yy[i],"%.1fHz"%xx[i],fontsize=9)
        
        ax.set_ylabel('Amplitude')
        ax.set_xlabel('Frequency (Hz)')
        if show:
            plt.show()
        return f,ax
    def fft(self,max_freq=200):
        from scipy.fft import fft
        from numpy import linspace,abs
        
        t = asarray(range(int( (self.w/32)*self.w) ))*0.00001
        y = zeros(( int( (self.w/32)*self.w )  ))
        for i in range(32):
            amp = self.get_amp(i+1)
            amp/=amp.stats('median')
            y+= amp.ravel()
        
        if any(isnan(y)):
            print("\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            print(":::   [Warning] NaN detected, fix solution not yet fully tested.   :::")
            print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n")
            mask = isfinite(y)
            y = interp(t,t[mask],y[mask])    
            t = t[mask]
        
        tt = t[-1]-t[0]
        N = len(t)
        T = tt / N
        
        yf = fft(y)
        yf = asarray(2.0/N * abs(yf[0:N//2]))
        xf = linspace(0.0, 1.0/(2.0*T), N//2)
        xf = xf[xf<=max_freq]
        yf = yf[0:len(xf)]
        
        return xf[1:],yf[1:]
    def get_freq(self):
        from scipy.signal import argrelextrema
        from numpy import greater,zeros
        from astropy.stats import sigma_clipped_stats as sc
        
        x,y = self.fft()        
        mn,_,st = sc(y)
        peaks_id = argrelextrema(y,greater,order=10)
        mask1 = zeros((len(y)),dtype=bool)
        mask2 = zeros((len(y)),dtype=bool)
        mask1[peaks_id]=True
        mask2[y>(mn+20*st)]=True
        peaksX = x[mask1&mask2]
        peaksY = y[mask1&mask2]
        if len(peaksY)<1:
            return [],[]
        Y,X = zip(*sorted(zip(peaksY,peaksX),reverse=True))
        return X,Y
    def copy(self):
        h = hxread()
        h.im = self.im
        h.H = self.H
        h.fname = self.fname
        h.shape = self.shape
        h.w = self.w
        h.h = self.h
        h.x = self.x
        return h
    def get_ro(self):
        st = [self.get_amp(i+1).stats('std') for i in range(32)]
        return mean(st),std(st)
    def show(self,show=True,**kwargs):
        from matplotlib import pyplot as plt
        import seaborn as sns
        from astropy.visualization import (SqrtStretch,ImageNormalize,ZScaleInterval)
        #sns.set_theme()#sns.set_style("dark")
        sns.set_style("white")
        f,ax = plt.subplots()
        if 'title' in kwargs:
            title = kwargs.get('title')
            ax.set_title(title)
        else:
            ax.set_title(basename(self.fname))
        norm = ImageNormalize(self.im,interval=ZScaleInterval(),stretch=SqrtStretch()) 
        ax.imshow(self.im,cmap='Greys_r',norm=norm)     
        if show:
            plt.show()
        return f,ax
    def __str__(self):
        return self.H.tostring(sep='\n')
    def __eq__(self):
        return self.H.tostring(sep='\n')
if '__main__' in __name__:
    read1 = hxread()
    read1('/home/noboru/NIRPS_R01_R01.fits')
    
    
