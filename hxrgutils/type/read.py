#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 11:53:18 2020

@author: espressjo
"""

from numpy import shape,asarray,zeros,mean,std,fliplr
from scipy.ndimage import zoom
from astropy.io import fits
from astropy.stats import sigma_clipped_stats as sc
from os.path import basename,join
from os import getcwd
from hxrgutils.type.amp import hxrg_amp
from hxrgutils.fits2ramp.refpixcorr import refpxcorr
from hxrgutils.stats.stats import stats
from hxrgutils.type.header import as_header
from numpy import isfinite,interp,arange

class hxread(refpxcorr,stats):
    """
    One detector read (a single 2D FITS frame) from an H2RG/H4RG.

    The object is built empty and then *called* to load data:

        r = hxread()
        r('read.fits')                 # from a FITS file
        r(array, header=some_header)   # or from a 2D array (header optional)

    Inherits the reference-pixel corrections from refpxcorr (in place on .im)
    and the sigma-clipped statistics/histogram helpers from stats.

    Attributes
    ----------
    im : ndarray
        The image (float). Also reachable as .data.
    H : fits.Header
        Header of the read. Always a real header, empty when the data came
        from an array without a header.
    fname : str
        Base name of the source file ("tmp.fits" for an array).
    shape, h, w, x : tuple/int
        Geometry, always re-derived from the data loaded (x = 2 for an H2RG,
        4 for an H4RG). The x argument of __init__ only sizes the empty image.
    low_pixels, hi_pixels : int
        Filled by get_deviant_pixels().
    """
    def __init__(self,x=4):
        """
        Parameters
        ----------
        x : INT, optional
            Detector size in units of 1024 pixels (2: H2RG, 4: H4RG). Only
            used to size the empty placeholder image; loading data overrides
            it. The default is 4.
        """
        self.x = x
        self.im = zeros((1024*x,1024*x),dtype=float)
        self.H = None#the H property turns this into an empty fits.Header
        self.fname = ''
        self.low_pixels = 0
        self.hi_pixels = 0
        self._update_geometry()
        stats.__init__(self,self.im)
        refpxcorr.__init__(self,self.im)
    @property
    def H(self):
        """
        FITS header of the read. Always a fits.Header instance, even when the
        read was initialized from a numpy array without a header; in that case
        it is simply empty. Anything assigned here (None, a dict, an old-style
        fits.Header class, ...) is coerced into a header.
        """
        return self._H
    @H.setter
    def H(self,value):
        self._H = as_header(value)
    @property
    def data(self):
        """
        stats() and refpxcorr() operate on .data; we alias it to .im so both
        always see the current image, even after .im has been reassigned
        (e.g., by the in-place operators).
        """
        return self.im
    @data.setter
    def data(self,value):
        self.im = asarray(value,dtype=float)
    def _update_geometry(self):
        """
        Refresh shape,h,w and x from the image currently held in .im. Called
        every time a new image is loaded, whatever the source.
        """
        self.shape = shape(self.im)
        self.h = self.shape[0]
        self.w = self.shape[1]
        if self.w%1024!=0 or self.h%1024!=0:
            print("[Warning] %dx%d is not a standard HxRG format"%(self.w,self.h))
        self.x = int(self.w/1024)
    def __call__(self,data,header=None):
        """
        Load a read, either from a FITS file or from a numpy array.

        Parameters
        ----------
        data : STR/ndarray
            FITS file name, or a 2D array of the read.
        header : fits.Header, optional
            Header to attach when data is an array. Without it the read simply
            carries an empty header.

        Returns
        -------
        None.

        """
        if isinstance(data,str):
            if '/' not in data:
                data = join(getcwd(),data)
            hdu = fits.open(data)
            if hdu[0].data is None:
                hdu.close()
                raise ValueError("no data in the primary HDU of %s"%data)
            im = asarray(hdu[0].data,dtype=float)
            self.H = hdu[0].header if header is None else header
            self.fname = basename(data)
            hdu.close()
        else:
            im = asarray(data,dtype=float)
            self.H = header#None gives an empty header
            self.fname = "tmp.fits"
        if im.ndim!=2:
            raise ValueError("hxread handles a single 2D read, got a %dD array. Use make_ramp() for a cube of reads."%im.ndim)
        self.im = im
        self._update_geometry()
        stats.__init__(self,self.im)
        refpxcorr.__init__(self,self.im)
    def unfold(self):
        '''
        Put the image in readout-time order, in place. Adjacent amplifiers
        are read in opposite directions, so every other ribbon (amp. 1, 3,
        5, ...) is flipped left/right. Useful before an fft or an IPC
        analysis. Calling it twice restores the original orientation.

        Returns
        -------
        None.

        '''
        w = int(self.x*1024/32.)
        for i in range(32):
            if (i+1)%2!=0:
                amp = self.im[:,i*w:(i*w)+w]
                self.im[:,i*w:(i*w)+w] = fliplr(amp)                
    def get_bias(self,sigma=5):
        """
        Sigma-clipped (mean, median, std) of the whole image, e.g., to
        characterize a bias read.

        Parameters
        ----------
        sigma : FLOAT, optional
            Clipping threshold. The default is 5.

        Returns
        -------
        tuple
            (mean, median, std) as returned by astropy's sigma_clipped_stats.
        """
        return sc(self.im,sigma=sigma)
    def refpx(self,**kwargs):
        """
        Reference-pixel correction of the read, in place on self.im, and a
        comment added to the header.

        Parameters
        ----------
        toponly : BOOL, optional
            Side + top reference pixels (refpxcorrtop).
        sideonly : BOOL, optional
            Side (left/right) reference pixels only (refpxcorrside).
            toponly has priority if both are given. With neither, the full
            correction (side + top + bottom, refpxcorr) is applied.

        NOTE
        ----
            oddeven and medfiltersize are not forwarded; call refpxcorr(),
            refpxcorrtop() or refpxcorrside() directly to use them.

        Returns
        -------
        None.
        """
        if 'toponly' in kwargs and kwargs.get('toponly'):
            self.refpxcorrtop()
            self.H.add_comment('refpixcorr performed [top only]')
        elif 'sideonly' in kwargs and kwargs.get('sideonly'):
            self.refpxcorrside()
            self.H.add_comment('refpixcorr performed [side only]')
        else:
            self.refpxcorr()
            self.H.add_comment('refpixcorr performed')
    def ds9(self):
        """Write the read to /var/tmp/read.fits and open it in ds9 (zscale). Needs ds9 on the PATH."""
        fits.PrimaryHDU(data=self.im,header=self.H).writeto("/var/tmp/read.fits",overwrite=True)
        from os import popen 
        popen("ds9 -zscale /var/tmp/read.fits").read()
    def get_deviant_pixels(self,pixel_bin=32,low_threshold=0.5,hi_threshold=1.3):
        """
        Count pixels far from the local level of the image. The image is
        divided by a low-frequency map (median of each pixel_bin x pixel_bin
        box, zoomed back to full size); pixels whose ratio is above
        hi_threshold or below low_threshold are counted. The 4-pixel
        reference border is ignored.

        Parameters
        ----------
        pixel_bin : INT, optional
            Size of the binning box in pixels. The default is 32.
        low_threshold : FLOAT, optional
            Ratio under which a pixel is "low". The default is 0.5.
        hi_threshold : FLOAT, optional
            Ratio over which a pixel is "high". The default is 1.3.

        Returns
        -------
        low_pixels, hi_pixels : INT
            Counts, also stored in self.low_pixels and self.hi_pixels.
        """
        #pixel size of the binning box
        import numpy as np
        pix_bin = pixel_bin
        nbin = self.w//pix_bin
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
        """
        Return the nth amplifier (1 to 32) as a hxrg_amp instance.
        """
        if n<1 or n>32:
            raise ValueError("amp. number must be between 1 and 32, got %s"%n)
        x_start = (n-1)*int(self.w/32)
        x_stop = (n)*int(self.w/32)
        return hxrg_amp(self.im[:,x_start:x_stop],namp=n)
    def writeto(self,fname):
        """Write the image and its header to a single-HDU FITS file (overwrites)."""
        hdul =fits.PrimaryHDU(data=self.im,header=self.H)
        hdul.writeto(fname,overwrite=True)
    # Arithmetic. The binary operators (+ - * / **) return a NEW hxread and
    # leave the operands untouched; with another hxread the header is copied
    # from the left operand and a comment names the operation. The in-place
    # operators (+= -= *= /=) modify self.im. other is an hxread or a scalar/array.
    def __sub__(self,other):
        """self - other (CDS when both are reads). Returns a new hxread."""
        if isinstance(other,hxread):
            F = hxread()
            F(self.im-other.im)
            F.H = self.H.copy()
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
            F.H = self.H.copy()
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
            F.H = self.H.copy()
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
            F.H = self.H.copy()
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
        """
        Plot the amplitude spectrum from fft() and label (up to 8) detected
        peaks from get_freq().

        Parameters
        ----------
        show : BOOL, optional
            Call plt.show(). The default is True.
        max_freq : FLOAT, optional
            Highest frequency plotted, in Hz. The default is 200.

        Returns
        -------
        f, ax : matplotlib figure and axes.
        """
        from matplotlib import pyplot as plt
        import seaborn as sns
        sns.set_theme()

        f,ax = plt.subplots(figsize=(10,5))
        x,y = self.fft(max_freq)
        if len(x)<1:
            print("[Warning] nothing to plot")
            return f,ax
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
        '''
        Amplitude spectrum of the read in the readout time domain. Each of the
        32 amp. is normalized by its median and co-added (in readout order) to
        boost the SNR of the noise pattern common to all amp.

        NaN are handled as follow; a NaN pixel is simply not co-added and the
        sample is rescaled by the number of amp. that did contribute. A sample
        where every amp. is NaN is linearly interpolated from its finite
        neighbours, and an amp. with no valid normalization is skipped. If
        nothing is finite, two empty arrays are returned.

        Parameters
        ----------
        max_freq : FLOAT, optional
            Highest frequency returned, in Hz. The default is 200.

        Returns
        -------
        xf : ndarray
            Frequencies (Hz).
        yf : ndarray
            Amplitudes.

        '''
        from scipy.fft import fft
        from numpy import linspace,abs

        nsample = int(self.h*int(self.w/32))
        t = arange(nsample)*0.00001
        y = zeros(nsample)
        nco = zeros(nsample)#number of amp. co-added on each sample
        for i in range(32):
            amp = self.get_amp(i+1)
            norm = amp.stats('median')
            if not isfinite(norm) or norm==0:
                print("[Warning] amp. %d has no valid median, skipped in the fft"%(i+1))
                continue
            amp/=norm
            a = amp.ravel()
            good = isfinite(a)
            y[good]+=a[good]
            nco[good]+=1

        good = nco>0
        if good.sum()<2:
            print("[Warning] no valid pixel left, no fft computed")
            return asarray([]),asarray([])
        #keep the amplitude scale of a full 32 amp. co-add
        y[good]*= 32./nco[good]
        if not good.all():
            print("\n::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
            print(":::   [Warning] NaN detected, interpolated before the fft.        :::")
            print("::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n")
            #interpolate on the full time grid, t must keep its sampling
            y = interp(t,t[good],y[good])

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
        """
        Find the noise peaks of the fft(): local maxima that stand more than
        20 sigma (sigma-clipped) above the mean amplitude.

        Returns
        -------
        X, Y : tuple/list
            Frequencies (Hz) and amplitudes of the peaks, strongest first.
            Two empty lists when there is no peak or no valid spectrum.
        """
        from scipy.signal import argrelextrema
        from numpy import greater,zeros
        from astropy.stats import sigma_clipped_stats as sc
        
        x,y = self.fft()
        if len(y)<3:
            return [],[]
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
        """
        Deep copy of the read. The new instance shares neither its image nor
        its header with the original.
        """
        h = hxread()
        h(self.im.copy(),header=self.H.copy())
        h.fname = self.fname
        return h
    def get_ro(self):
        """
        Readout noise estimate: sigma-clipped std of each of the 32
        amplifiers; returns (mean, std) of those 32 values, in ADU. Most
        meaningful on a bias-like read or a CDS.
        """
        st = [self.get_amp(i+1).stats('std') for i in range(32)]
        return mean(st),std(st)
    def show(self,show=True,**kwargs):
        """
        Display the read (zscale, sqrt stretch).

        Parameters
        ----------
        show : BOOL, optional
            Call plt.show(). The default is True.
        title : STR, optional (kwarg)
            Plot title; defaults to the file name.

        Returns
        -------
        f, ax : matplotlib figure and axes.
        """
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
        """The header, one card per line."""
        return self.H.tostring(sep='\n')
if '__main__' in __name__:
    read1 = hxread()
    read1('/home/noboru/NIRPS_R01_R01.fits')
    
    
