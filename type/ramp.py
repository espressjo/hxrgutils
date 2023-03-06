#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:03:14 2020

@author: espressjo
"""
from astropy.io import fits
from hxrgutils.type.read import hxread
from os.path import join,isfile
from os import getcwd
from os import stderr
import numpy as np

def show(im):
    fits.PrimaryHDU(data=im).writeto("/var/tmp/tmp.fits",overwrite=True)
    from os import popen
    popen("ds9 -zscale /var/tmp/tmp.fits").read()
def make_ramp(fname,errors=False):
    from tqdm import tqdm
    R = hxramp()
    r = hxread()
    print("Performing fits2ramp...")
    if isinstance(fname, list):
        for f in tqdm(fname):
            r(f)
            R<<r
    elif isinstance(fname, str):
        cube = fits.getdata(fname)
        for im in tqdm(cube):
            r(im)
            R<<r
    R.fit()
    if errors:
        print("Calculating errors...")
        R.create_error(fname)
    return R
class hxramp():
    def __init__(self,**kwargs):
        self.i = 0
        self.H = fits.Header()
        self.timestamp = [];
        self.inttime = 5.24288
        self.bias = np.zeros((4096,4096))
        self.err_calulated = False#flag for the error calculation
        self.errslope = np.zeros([4096,4096],dtype=float)+np.nan
        if 'saturation' in kwargs:
            self.saturation = kwargs.get('saturation')
        else:
            self.saturation = 45000
        
        self.refpxeachread = False
        self.selfbias = True
        self.toponly = False
        self.notopnobottom = False
        self.oddeven = False
        self.nl = False
        if 'toponly' in kwargs:
            self.toponly = kwargs.get('toponly')
        if 'notopnobottom' in kwargs:
            self.notopnobottom = kwargs.get('notopnobottom')
        if 'oddeven' in kwargs:
            self.oddeven = kwargs.get('oddeven')
        if 'refpxeachread' in kwargs:
            self.refpxeachread = kwargs.get('refpxeachread')

    def nlcorr(self,array):
        """
        Apply the non-lineariy correction. If no NL file
        is provided, this function will return its inputs

        Parameters
        ----------
        array : 
            4096X4096 array.

        Returns
        -------
        TYPE
            NL corrected array.

        """
        if self.nl:
            array = np.asarray(array,dtype=float)        
            return array+self.c3*(array**3)+self.c2*(array**2)
        else:
            return array
    def upload_nonlin(self,f:str):
        """
        Upload a non-linearity calibration file. The NL file
        should have 2 extension where ext=0 is the C3 coefficient
        and ext=1 is the C2 coefficient.

        Parameters
        ----------
        f : str
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if isfile(f):
            hdu = fits.open(f)            
            self.c3 = hdu[0].data
            self.c2 = hdu[1].data
            self.nl = True
            if 'UID' in hdu[0].header:
                self.H["NLUID"] = (hdu[0].header["UID"],"Unique ID number for the NL cal. file.")
        else:
            print('%s not found, no nl-correction will be apply'%f)
    def upload_bias(self,fname:str):
        '''
        Upload a super bias
        '''
        if isfile(fname):
            self.bias = fits.getdata(fname)
            if self.bias.shape!=(4096,4096):
                print("Bias file has wrong format")
                exit(1)
            self.selfbias = False
            if "UID" in self.getheader(fname):
                self.H["BIASUID"] = (self.getheader(fname)["UID"],"Unique ID number for the bias cal. file.")
        else:
            print("[Warning] bias file not found, will perform selfbias")            
    def unfold(self):
        '''
        Flip [l/r] every odd amp. This can be usefull to perform fft, 
        or IPC analysis. fit() method should be called before this.

        Returns
        -------
        None.

        '''
        for i in range(32):
            if (i+1)%2!=0:
                amp = self.a[:,i*128:(i*128)+128]
                self.a[:,i*128:(i*128)+128] = np.fliplr(amp)                
                amp = self.b[:,i*128:(i*128)+128]
                self.b[:,i*128:(i*128)+128] = np.fliplr(amp)
                amp = self.n[:,i*128:(i*128)+128]
                self.n[:,i*128:(i*128)+128] = np.fliplr(amp)
    def __call__(self,fname:str):
        '''
        Open a ramp file for NIRPS, SPIP or SPIRou.
        '''
        
        if '/' not in fname:
            fname = join(getcwd(),fname)
        hdu = fits.open(fname)
        self.H = hdu[0].header
        #NIRPS Case
        if 'HIERARCH ESO DET SEQ1 DIT' in self.H:
            self.total_integration = self.H['HIERARCH ESO DET SEQ1 DIT']
        if len(hdu)==5:
            #SPIP or SPIRou
            self.a = np.asarray(hdu[1].data,dtype=float)
            self.b = np.asarray(hdu[2].data,dtype=float)
            self.errslope = np.asarray(hdu[3].data,dtype=float)
            self.n = hdu[4].data
        elif len(hdu)==4:
            #NIRPS
            self.a = np.asarray(hdu[1].data,dtype=float)
            self.b = np.asarray(hdu[2].data,dtype=float)
            self.n = hdu[3].data
        else:
            print("Shape of ramp file not recongnized.",file=stderr)
    def init_from_read(self,read):
        '''
        Initialize a new ramp using a 1st read

        Parameters
        ----------
        read : hxread
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        self.shape = (read.x*1024,read.x*1024)
        self.sx = np.zeros(self.shape)
        self.sx2 = np.zeros(self.shape)
        self.n = np.zeros(self.shape,dtype=np.int16)#number of read used to perform fits2ramp
        self.sy = np.zeros(self.shape)
        self.sxy = np.zeros(self.shape)
        self.a = np.zeros(self.shape)#slope
        self.b = np.zeros(self.shape)#intercept
    def applyrefpxcorr(self,array):
        '''
        Apply ref pix. correction given the options

        Parameters
        ----------
        array : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        from hxrgutils.fits2ramp.refpixcorr import refpxcorr as rpc
        if self.toponly:
            return rpc(np.asarray(array,dtype=float)).refpxcorrtop(oddeven=self.oddeven)
        elif self.notopnobottom:
            return rpc(np.asarray(array,dtype=float)).refpxcorrside(oddeven=self.oddeven)
        else:
            return rpc(np.asarray(array,dtype=float)).refpxcorr(oddeven=self.oddeven)
    def get_time(self):
        return self.timestamp[:-1]
    def get_coeff(self):
        '''
        Return slope and intercept as hxread datatype

        Returns
        -------
        a : hxread
            slope.
        b : hxread
            intercept.

        '''
        a = hxread()
        a(self.a)
        a.H = self.H
        a.fname = 'coeff_a.fits'
        b = hxread()
        b(self.b)
        b.H = self.H        
        b.fname = 'coeff_b.fits'
        return a,b
    def writeto(self,fname):
        """
        Write the fits2ramp to file.

        Parameters
        ----------
        fname : STR
            filename.

        Returns
        -------
        None.

        """
        phdu = fits.PrimaryHDU(header=self.H)
        if self.err_calulated:
            hdul = fits.HDUList([phdu,fits.ImageHDU(data=self.a,name="slope"),fits.ImageHDU(data=self.b,name="intercept"),fits.ImageHDU(data=self.errslope,name="errors"),fits.ImageHDU(data=self.n,name="n")])
        else:
            hdul = fits.HDUList([phdu,fits.ImageHDU(data=self.a,name="slope"),fits.ImageHDU(data=self.b,name="intercept"),fits.ImageHDU(data=self.n,name="n")])
        hdul.writeto(fname,overwrite=True)
    def __lshift__(self, o):  
        """
        Calculate the fits2ramp as a streaming algorithm. The nth read is 
        feed as a hxread to the hxramp. e.g.,
        
        for read in [read1,read2,...]:
            Ramp << read
        ramp.fit()
        
        Parameters
        ----------
        o : hxread
            Ramp nth read as a hxread instance.

        Returns
        -------
        None.

        """
        if not isinstance(o,hxread):
            _o = hxread()
            _o(o)
            o = _o
        if self.i==0:
            self.init_from_read(o)
            if self.selfbias:
                self.bias = o.im
            #self.timestamp.append(self.inttime)
            #self.i+=1
            self.total_integration = 0
            self.goodmask = np.full((o.w,o.h),True,dtype=bool)
            #self.n+=self.goodmask
            self.H = o.H
        self.i+=1
        self.total_integration +=self.inttime
        self.timestamp.append(self.inttime*(self.i))
        self.goodmask = (o.im <= self.saturation)*self.goodmask
        self.n+=self.goodmask    
        if self.nl:
            
            imc = self.nlcorr(o.im-self.bias)
        else:
            imc = (o.im-self.bias)
        if self.refpxeachread:
             imc = self.applyrefpxcorr(imc)
        
        self.sy[self.goodmask]+=imc[self.goodmask]
        self.sxy[self.goodmask]+=(imc[self.goodmask]*self.inttime*(self.i))
    def create_error(self,fname):
        """
        Once the fits2ramp is calculated either by calling  << operator
        or by uploading a SPIP, SPIRou or NIRPS ramp file.

        Parameters
        ----------
        fname : STR/LIST
            fname can be a list of individual read or a cube of reads.
            If fname is a list of individual file, we us natsorted to make
            sure the list is in order.
        Returns
        -------
        None.

        """
        if isinstance(fname, list):
            from natsort import natsorted
            fname = natsorted(fname)
            Cube = np.asarray([np.asarray(fits.getdata(f),dtype=float) for f in fname])
        elif isinstance(fname, str):
            Cube = np.asarray(fits.getdata(fname),dtype=float)
        else:
            print("create_error only accept list of file or stream cube")
            exit(1)
        dim3 = len(Cube)
        sx = self.create_sx()
        
        varx2 = np.zeros([4096,4096],dtype=float)
        vary2 = np.zeros([4096,4096],dtype=float)
        xp = np.zeros([4096,4096],dtype=float)
        goodmask = np.full((4096,4096),True,dtype=bool)
        valid = (self.n>2)
        xp[valid]=sx[valid]/self.n[valid] # used in the determination of error below
        print('we now compute the standard error on the slope')


        #bias = self.nlcorr(np.copy(np.copy(Cube[0])-self.bias))
        bias = np.copy(Cube[0])-self.bias 
        bias = self.nlcorr(np.copy(bias))
        print(f"bias {np.median(bias)}")
        #show(bias)
        for i in range(1,dim3-1):

            # we read the npz as this file has been linearized (if the -linearize keyword has been set)
            # and we subtracted the reference regions on the array
            im = np.copy(Cube[i])
            im = im - self.bias
            print(f"im {np.median(im)}")
            im = self.nlcorr(np.copy(im))
            #show(im)
            im =    self.applyrefpxcorr(np.copy(im-bias))
            
            
            print(f"Image {i+1} (corr): {np.nanmedian(im)} ADU")
            #show(im)
            goodmask = (self.n > i)
            #print(self.timestamp[i])
            yp = self.a*self.timestamp[i] - self.a*self.timestamp[0]
            print(f"Reconstruction: {np.nanmedian(yp)} ADU")
            #show(im-yp)
            print(i+1,'/',dim3,' ~~~> Computing slope error')

            varx2+= ((self.timestamp[i-1]-xp)**2)*goodmask # we multiply by goodmask so that only
            vary2+= ((im-yp)**2)*goodmask

        valid*=(varx2!=0) # avoid diving by zero
        self.err_calulated = True
        self.errslope[valid] = np.sqrt(vary2[valid]/(self.n[valid]-2))/np.sqrt(varx2[valid])
    def create_sx(self):
        """
        Utility function in case a ramp file is uploaded (sx not in this case)

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if len(self.timestamp)!=0:
            
            return np.where(self.n>0,(np.cumsum(self.timestamp))[self.n-1],0)
        else:
            max_n = np.amax(self.n)
            for i in range(max_n):
                self.timestamp.append((i+1)*self.inttime)
            return np.where(self.n>0,(np.cumsum(self.timestamp))[self.n-1],0)
    def fit(self):
        """
        One of the main function. Once all the FITS file (hxread) are
        streamed, you can call fit to compute the slope and intercept.

        Returns
        -------
        None.

        """
        self.sx=np.where(self.n>0,(np.cumsum(self.timestamp))[self.n-1],0)
        self.sx2=np.where(self.n>0,(np.cumsum(np.asarray(self.timestamp)**2))[self.n-1],0)
        valid = self.n>1
        self.b[valid] = (self.sx*self.sxy-self.sx2*self.sy)[valid]/(self.sx**2-self.n*self.sx2)[valid] # algebra of the linear fit
        self.a[valid] = (self.sy-self.n*self.b)[valid]/self.sx[valid]        
        if not self.refpxeachread:
            self.a = self.applyrefpxcorr(self.a)
            
        if self.refpxeachread:
            print("Ref. px. corrected after each read")
        else:
            print("Ref. px. corrected on the slope")
        if self.selfbias:
            print("Self biased")
        else:
            print("super bias used")
        
        if self.toponly:
            print("Only using top and side ref. px.")
        elif self.notopnobottom:
            print("Only side ref. px. used")
        else:
            print("All ref. px. used")
        
        if self.oddeven:
            print("odd and even column processed separetly")
   
        if self.nl:
            print("nl correction applied.")
        self.b+=self.bias            
    def show(self):
        """
        Show slope and intercept in a plot.

        Returns
        -------
        None.

        """
        from matplotlib import pyplot as plt
        from astropy.visualization import (SqrtStretch,ImageNormalize,ZScaleInterval)
        f,ax = plt.subplots()
        ax.set_title('Coeff. a')
        norm = ImageNormalize(self.a,interval=ZScaleInterval(),stretch=SqrtStretch()) 
        ax.imshow(self.a,cmap='Greys_r',norm=norm)   
        
        f,ax = plt.subplots()
        ax.set_title('Coeff. b')
        norm = ImageNormalize(self.b,interval=ZScaleInterval(),stretch=SqrtStretch()) 
        ax.imshow(self.b,cmap='Greys_r',norm=norm) 
        plt.show()
    def preprocess(self):
        """
        Apply a quick preprocessing to the slope.

        Returns
        -------
        None.

        """
        from hxrgutils.preprocess import preprocess
        self.a = preprocess.preprocess(self.a)
    