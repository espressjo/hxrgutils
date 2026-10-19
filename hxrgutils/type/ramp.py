#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:03:14 2020

@author: espressjo
"""
from astropy.io import fits
from hxrgutils.type.read import hxread
from hxrgutils.type.header import as_header
from os.path import join,isfile
from os import getcwd
import numpy as np

def show(im):
    """
    Given an ndarray, this function will show 
    the image in ds9.

    Parameters
    ----------
    im : ndarray
        Image.

    Returns
    -------
    None.

    """
    fits.PrimaryHDU(data=im).writeto("/var/tmp/tmp.fits",overwrite=True)
    from os import popen
    popen("ds9 -zscale /var/tmp/tmp.fits").read()
def make_ramp(fname,errors=False,nonlin="",bias="",**kwargs):
    """
    Generate a ramp object from a list of filename or 
    a data cube of individual read.

    Parameters
    ----------
    fname : STR/LIST/ndarray
        Either a list of individual read (file names or arrays), the file
        name of a cube of individual read, or a 3D array of individual read.
    errors : BOOL, optional
        Weither we compute the error on the slope or not. 
        The default is false.
    nonlin : STR, optional
        filename of the nonlinearity reference file
    bias : STR, optional
        filename of the bias reference file
    NOTE
    -----
        **kwargs is pass directly to hxramp()
    Returns
    -------
    R : hxramp
        Ramp object.

    """
    from tqdm import tqdm
    R = hxramp(**kwargs)
    if nonlin!="" and isfile(nonlin):
        R.upload_nonlin(nonlin)
    if bias!="" and isfile(bias):
        R.upload_bias(bias)
    r = hxread()
    print("Performing fits2ramp...")
    if isinstance(fname, str):
        cube = np.asarray(fits.getdata(fname),dtype=float)
        if cube.ndim!=3:
            raise ValueError("%s does not hold a cube of reads"%fname)
    elif isinstance(fname,(list,tuple)):
        cube = fname
    else:
        cube = np.asarray(fname,dtype=float)
        if cube.ndim!=3:
            raise ValueError("make_ramp accepts a list of reads, a cube file name or a 3D array")
    for im in tqdm(cube):
        r(im)#hxread handles both a file name and an array
        R<<r
    R.fit()
    if errors:
        print("Calculating errors...")
        R.create_error(fname)
    return R
class hxramp():
    """
    Up-the-ramp (slope) fit of a series of reads, or a ramp file loaded from disk.

    Two ways to use it:

    1. Build a ramp from reads: stream them one at a time with <<, then fit().

           R = hxramp(h2rg=True, saturation=45000)
           for read in reads:      # hxread, 2D array or file name, in time order
               R << read
           R.fit()                 # fills R.a (slope, ADU/s) and R.b (intercept)
           R.writeto('ramp.fits')

    2. Load an existing ramp product: R('ramp.fits'). The layout is told by
       the number of HDUs: 4 = NIRPS (slope, intercept, n), 5 = SPIRou/SPIP
       (slope, intercept, errslope, n). A single HDU or a 3D array is a cube
       of reads and is streamed instead (fit() is then left to the caller).

    Keyword arguments (all optional)
    --------------------------------
    h2rg : BOOL
        Hint for a 2048x2048 array (else 4096x4096). The data loaded always
        overrides it.
    saturation : FLOAT
        ADU above which a read is excluded, for that pixel and all later
        reads. Default 45000.
    norefpxcorr : BOOL
        Skip the reference-pixel correction.
    refpxeachread : BOOL
        Correct each (read - bias) while streaming instead of correcting the
        slope once in fit().
    toponly : BOOL
        Side + top reference pixels only (no bottom).
    notopnobottom : BOOL
        Side reference pixels only. toponly has priority if both are set.
    oddeven : BOOL
        Treat odd and even columns separately in the correction.

    Main attributes
    ---------------
    a, b : ndarray
        Slope (ADU/s) and intercept (ADU, bias added back).
    n : ndarray (int16)
        Number of non-saturated reads used per pixel.
    errslope : ndarray
        Formal slope error; NaN until create_error() has run (or the file
        loaded had one).
    H : fits.Header
        Header, copied from the first read (empty for array input).
    inttime : float
        Time between reads, (w*h)*1e-5/32 s.
    nl : BOOL
        True once a non-linearity file has been uploaded.
    """
    def __init__(self,**kwargs):
        """See the class docstring for the accepted keyword arguments."""
        if 'h2rg' in kwargs and kwargs.get('h2rg'):
            self.x = 2
        else:
            self.x = 4
        self.i = 0#number of read streamed so far
        self.loaded = False#True when the ramp comes from an already fitted file
        self._geometry_set = False#True once we know the real format of the data
        self.H = None#the H property turns this into an empty fits.Header
        self.timestamp = [];
        self.shape = (self.x*1024,self.x*1024)
        self.inttime = (self.x*1024)**2*1e-5/32.
        self.bias = np.zeros(self.shape)
        self.total_integration = 0.
        self.err_calulated = False#flag for the error calculation
        self.errslope = np.zeros(self.shape,dtype=float)+np.nan
        if 'saturation' in kwargs:
            self.saturation = kwargs.get('saturation')
        else:
            self.saturation = 45000
        self.norefpxcorr = False
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
        if 'saturation' in kwargs:
            self.saturation = kwargs.get('saturation')
        if 'norefpxcorr' in kwargs:
            self.norefpxcorr = kwargs.get('norefpxcorr')
    @property
    def H(self):
        """
        FITS header of the ramp. Always a fits.Header instance; a ramp built
        from numpy arrays instead of FITS files simply carries an empty one,
        so every method that writes a keyword or a comment keeps working.
        """
        return self._H
    @H.setter
    def H(self,value):
        self._H = as_header(value)
    def _set_geometry(self,shape):
        """
        Adapt every size dependent attribute to the format of the data we are
        actually given (H2RG, H4RG, ...), no matter if it comes from a file or
        from a numpy array. The h2rg keyword is therefore only a hint; the
        data always wins.

        Parameters
        ----------
        shape : tuple
            Shape of the detector array we are about to handle.

        Returns
        -------
        None.

        """
        shape = tuple(shape)
        if len(shape)!=2:
            raise ValueError("expecting a 2D detector array, got a %dD one"%len(shape))
        if not self.selfbias and self.bias.shape!=shape:
            raise ValueError("bias is %s but the data is %s"%(str(self.bias.shape),str(shape)))
        self.shape = shape
        h,w = shape
        if w%1024!=0 or h%1024!=0:
            print("[Warning] %dx%d is not a standard HxRG format"%(w,h))
        self.x = int(w/1024)
        self.inttime = (w*h)*1e-5/32.
        if self.bias.shape!=shape:
            self.bias = np.zeros(shape)
        if self.errslope.shape!=shape:
            self.errslope = np.zeros(shape,dtype=float)+np.nan
        if self.nl and np.shape(self.c3)!=shape:
            print("[Warning] nl-coeff. are %s while the data is %s, no nl-correction will be applied"%(str(np.shape(self.c3)),str(shape)))
            self.nl = False
        self._geometry_set = True
    def _stream_cube(self,cube,header=None):
        """
        Stream every read of a cube through the << operator. fit() is left to
        the caller.

        Parameters
        ----------
        cube : ndarray/list
            3D array, or list of 2D arrays, of individual reads.
        header : fits.Header, optional
            Header to attach to the ramp. Without it the ramp keeps an empty
            header.

        Returns
        -------
        None.

        """
        if isinstance(cube,(list,tuple)) and all([isinstance(c,str) for c in cube]):
            reads = cube#list of file name
        else:
            reads = np.asarray(cube,dtype=float)
            if reads.ndim!=3:
                raise ValueError("expecting a cube of individual read (3D) or a list of read, got a %dD array. A single read is handled by hxread."%reads.ndim)
        r = hxread()
        for im in reads:
            r(im,header=header)
            self<<r
    def nlcorr(self,array):
        """
        Apply the non-linearity correction, array + C3*array**3 + C2*array**2.
        The array should already be bias-subtracted. If no NL file was
        uploaded, the input is returned unchanged.

        Parameters
        ----------
        array : ndarray
            2D array with the same shape as the NL coefficients (e.g.,
            4096x4096 for an H4RG).

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
        and ext=1 is the C2 coefficient (2D images). Call it before
        streaming reads. A missing file only prints a message (no NL
        correction). If the ramp geometry is already known and differs from
        the coefficients, a ValueError is raised. A "UID" keyword in the file
        is copied to the ramp header as NLUID.

        Parameters
        ----------
        f : str
            Path of the non-linearity file.

        Returns
        -------
        None.

        """
        if isfile(f):
            hdu = fits.open(f)            
            self.c3 = np.asarray(hdu[0].data,dtype=float)
            self.c2 = np.asarray(hdu[1].data,dtype=float)
            self.nl = True
            if self._geometry_set and self.c3.shape!=self.shape:
                self.nl = False
                raise ValueError("nl-coeff. are %s but the ramp is %s"%(str(self.c3.shape),str(self.shape)))
            if 'UID' in hdu[0].header:
                self.H["NLUID"] = (hdu[0].header["UID"],"Unique ID number for the NL cal. file.")
        else:
            print('%s not found, no nl-correction will be apply'%f)
    def upload_bias(self,fname:str):
        '''
        Upload a super bias (2D FITS image in the primary HDU), used instead
        of the default self-bias (first read). Call it before streaming
        reads; it also fixes the ramp geometry. A missing file prints a
        warning and self-bias is kept; a wrong shape or ndim raises
        ValueError. A "UID" keyword in the file is copied to the header as
        BIASUID.

        Parameters
        ----------
        fname : str
            Path of the bias file.
        '''
        if isfile(fname):
            bias = np.asarray(fits.getdata(fname),dtype=float)
            if bias.ndim!=2:
                raise ValueError("bias file must hold a 2D image, got a %dD one"%bias.ndim)
            if self._geometry_set and bias.shape!=self.shape:
                raise ValueError("bias file is %s but the ramp is %s"%(str(bias.shape),str(self.shape)))
            self.bias = bias
            self.selfbias = False
            self._set_geometry(bias.shape)#no read streamed yet: the bias sets the format
            if "UID" in fits.getheader(fname):
                self.H["BIASUID"] = (fits.getheader(fname)["UID"],"Unique ID number for the bias cal. file.")
        else:
            print("[Warning] bias file not found, will perform selfbias")            
    def unfold(self):
        '''
        Put slope, intercept and n in readout-time order, in place: every
        other amplifier ribbon (amp. 1, 3, 5, ...) is flipped left/right.
        Useful for an fft or an IPC analysis. fit() (or loading a ramp) must
        be done first; calling it twice restores the original orientation.

        Returns
        -------
        None.

        '''
        if not hasattr(self,'a'):
            raise RuntimeError("nothing to unfold, compute or load a ramp first")
        w = int(self.shape[1]/32.)
        for i in range(32):
            if (i+1)%2!=0:
                amp = self.a[:,i*w:(i*w)+w]
                self.a[:,i*w:(i*w)+w] = np.fliplr(amp)                
                amp = self.b[:,i*w:(i*w)+w]
                self.b[:,i*w:(i*w)+w] = np.fliplr(amp)
                amp = self.n[:,i*w:(i*w)+w]
                self.n[:,i*w:(i*w)+w] = np.fliplr(amp)
    def __call__(self,data,header=None):
        '''
        Load a ramp. Data can come from a file or from memory:

            * file name of a ramp already fitted (NIRPS, SPIP or SPIRou);
            * file name of a cube of individual read;
            * a 3D array, or a list of 2D arrays, of individual read. Those
              are streamed exactly like "ramp << read", fit() is left to the
              caller.

        Parameters
        ----------
        data : STR/ndarray/list
            See above.
        header : fits.Header, optional
            Header to attach. Mostly usefull when data is an array; without
            it the ramp simply carries an empty header.

        Returns
        -------
        None.

        '''
        if not isinstance(data,str):
            self._stream_cube(data,header=header)
            return
        fname = data
        if '/' not in fname:
            fname = join(getcwd(),fname)
        hdu = fits.open(fname)
        nhdu = len(hdu)
        if nhdu==1 or np.ndim(hdu[0].data)==3:
            #cube of individual read rather than an already fitted ramp
            cube = np.asarray(hdu[0].data,dtype=float)
            H = hdu[0].header if header is None else header
            hdu.close()
            self._stream_cube(cube,header=H)
            return
        self.H = hdu[0].header if header is None else header
        if nhdu==5:
            #SPIP or SPIRou
            print("SPIRou/SPIP")
            self.a = np.asarray(hdu[1].data,dtype=float)
            self.b = np.asarray(hdu[2].data,dtype=float)
            self.errslope = np.asarray(hdu[3].data,dtype=float)
            self.n = np.asarray(hdu[4].data,dtype=np.int16)
            self.err_calulated = True
        elif nhdu==4:
            #NIRPS
            self.a = np.asarray(hdu[1].data,dtype=float)
            self.b = np.asarray(hdu[2].data,dtype=float)
            self.n = np.asarray(hdu[3].data,dtype=np.int16)
        else:
            hdu.close()
            raise ValueError("Shape of ramp file not recongnized (%d HDU in %s)."%(nhdu,fname))
        hdu.close()
        self.loaded = True
        self.timestamp = []
        #the format of the file wins over the h2rg keyword
        self._set_geometry(np.shape(self.a))
        #NIRPS Case
        if 'HIERARCH ESO DET SEQ1 DIT' in self.H:
            self.total_integration = self.H['HIERARCH ESO DET SEQ1 DIT']
        else:
            self.total_integration = self.inttime*float(np.amax(self.n))
    def init_from_read(self,read):
        '''
        Initialize the accumulators (sx, sx2, sy, sxy, n, a, b), all zeros, for
        a new ramp using a 1st read to set the geometry. Called by << on the
        first read; there is normally no need to call it yourself.

        Parameters
        ----------
        read : hxread
            First read; only its shape is used.

        Returns
        -------
        None.

        '''
        self._set_geometry(np.shape(read.im))
        self.sx = np.zeros(self.shape)
        self.sx2 = np.zeros(self.shape)
        self.n = np.zeros(self.shape,dtype=np.int16)#number of read used to perform fits2ramp
        self.sy = np.zeros(self.shape)
        self.sxy = np.zeros(self.shape)
        self.a = np.zeros(self.shape)#slope
        self.b = np.zeros(self.shape)#intercept
    def applyrefpxcorr(self,array):
        '''
        Apply the reference-pixel correction selected by the ramp options:
        toponly -> refpxcorrtop (side + top), else notopnobottom ->
        refpxcorrside (side only), else refpxcorr (side + top + bottom);
        oddeven is forwarded to all of them. norefpxcorr is *not* checked
        here (fit() does that).

        Parameters
        ----------
        array : ndarray
            2D array to correct. It is wrapped, not copied, so a float array
            is modified in place; pass np.copy(array) to keep the original.

        Returns
        -------
        ndarray
            The corrected array.

        '''
        from hxrgutils.fits2ramp.refpixcorr import refpxcorr as rpc
        if self.toponly:
            return rpc(np.asarray(array,dtype=float)).refpxcorrtop(oddeven=self.oddeven)
        elif self.notopnobottom:
            return rpc(np.asarray(array,dtype=float)).refpxcorrside(oddeven=self.oddeven)
        else:
            return rpc(np.asarray(array,dtype=float)).refpxcorr(oddeven=self.oddeven)
    def get_time(self):
        """Time stamps (s) of the reads streamed so far, without the last one."""
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
        a(self.a,header=self.H.copy())
        a.fname = 'coeff_a.fits'
        b = hxread()
        b(self.b,header=self.H.copy())
        b.fname = 'coeff_b.fits'
        return a,b
    def writeto(self,fname):
        """
        Write the ramp to a FITS file (overwrites). Layout: primary header,
        then slope, intercept, [errors,] n. The errors extension (5 HDUs,
        SPIRou/SPIP layout) is present only if the errors were computed or
        loaded; otherwise 4 HDUs (NIRPS layout).

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

        Per read: update the saturation mask and n, subtract the bias (the
        first read unless a bias file was uploaded), apply the NL correction
        and (if refpxeachread) the reference-pixel correction, then add to the
        running sums sy and sxy. Reads must come in time order; read i is
        stamped i*inttime. The header of the ramp is a copy of the first
        read's. Streaming into a ramp loaded from a file raises RuntimeError.

        Parameters
        ----------
        o : hxread/ndarray/str
            Ramp nth read as a hxread instance. A 2D array or a file name is
            wrapped into an hxread automatically.

        Returns
        -------
        None.

        """
        if self.loaded:
            raise RuntimeError("this ramp was loaded from a file, create a new hxramp() to stream reads")
        if not isinstance(o,hxread):
            _o = hxread()
            _o(o)#hxread accepts a file name or a 2D array
            o = _o
        if self.i==0:
            self.init_from_read(o)
            if self.selfbias:
                self.bias = np.copy(o.im)
            #self.timestamp.append(self.inttime)
            #self.i+=1
            self.total_integration = 0
            self.goodmask = np.full(self.shape,True,dtype=bool)
            #self.n+=self.goodmask
            self.H = o.H.copy()#empty header when the read came from an array
        elif np.shape(o.im)!=self.shape:
            raise ValueError("read #%d is %s while the ramp is %s"%(self.i+1,str(np.shape(o.im)),str(self.shape)))
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
        Compute the formal error on the slope (self.errslope) from the reads,
        once the fits2ramp is calculated either by calling  << operator
        or by uploading a SPIP, SPIRou or NIRPS ramp file. Each read is
        bias-subtracted, NL-corrected and reference-pixel corrected with the
        same options as the fit, and compared with the fitted line. Sets
        err_calulated so writeto() then emits the 5-HDU layout. Needs at
        least 3 reads (otherwise only a warning is printed) and reads with
        the same shape as the ramp.

        Parameters
        ----------
        fname : STR/LIST/ndarray
            fname can be a list of individual read (file names or arrays),
            the file name of a cube of reads, or a 3D array of reads. If
            fname is a list of individual file, we us natsorted to make
            sure the list is in order.
        Returns
        -------
        None.

        """
        if isinstance(fname,(list,tuple)) and all([isinstance(f,str) for f in fname]):
            from natsort import natsorted
            fname = natsorted(fname)
            Cube = np.asarray([np.asarray(fits.getdata(f),dtype=float) for f in fname])
        elif isinstance(fname, str):
            Cube = np.asarray(fits.getdata(fname),dtype=float)
        else:
            Cube = np.asarray(fname,dtype=float)
        if Cube.ndim!=3:
            raise ValueError("create_error only accept a list of read, a cube file name or a 3D array")
        if not hasattr(self,'n'):
            raise RuntimeError("no ramp yet, compute or load a ramp before calling create_error()")
        if Cube.shape[1:]!=self.shape:
            raise ValueError("the reads are %s while the ramp is %s"%(str(Cube.shape[1:]),str(self.shape)))
        dim3 = len(Cube)
        if dim3<3:
            print("[Warning] %d read(s) is not enough to compute an error on the slope"%dim3)
            return
        sx = self.create_sx()

        varx2 = np.zeros(self.shape,dtype=float)
        vary2 = np.zeros(self.shape,dtype=float)
        xp = np.zeros(self.shape,dtype=float)
        goodmask = np.full(self.shape,True,dtype=bool)
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
        Per-pixel sum of the time stamps of the n valid reads, rebuilt from n
        (and the time stamps, generated from inttime if missing). Needed by
        create_error() when a ramp file was loaded, since sx is not stored.

        Returns
        -------
        ndarray
            sx, 0 where n == 0.

        """
        if len(self.timestamp)!=0:
            
            return np.where(self.n>0,(np.cumsum(self.timestamp))[self.n-1],0)
        else:
            max_n = np.amax(self.n)
            print(max_n)
            for i in range(max_n):
                self.timestamp.append((i+1)*self.inttime)
            
            print(self.timestamp)
            return np.where(self.n>0,(np.cumsum(self.timestamp))[self.n-1],0)
    def fit(self):
        """
        One of the main function. Once all the FITS file (hxread) are
        streamed, you can call fit to compute the slope and intercept.

        Solves the linear fit from the running sums for pixels with n > 1
        (others keep slope and intercept 0), then, unless norefpxcorr or
        refpxeachread is set, applies the reference-pixel correction to the
        slope (not to the intercept). The bias is added back to the
        intercept. Prints which options were applied. Raises RuntimeError if
        no read was streamed.

        Returns
        -------
        None.

        """
        if not hasattr(self,'sy'):
            raise RuntimeError("no read streamed yet, feed the ramp with 'ramp << read' before calling fit()")
        self.sx=np.where(self.n>0,(np.cumsum(self.timestamp))[self.n-1],0)
        self.sx2=np.where(self.n>0,(np.cumsum(np.asarray(self.timestamp)**2))[self.n-1],0)
        valid = self.n>1
        self.b[valid] = (self.sx*self.sxy-self.sx2*self.sy)[valid]/(self.sx**2-self.n*self.sx2)[valid] # algebra of the linear fit
        self.a[valid] = (self.sy-self.n*self.b)[valid]/self.sx[valid]        
        if not self.norefpxcorr:
            if not self.refpxeachread:
                self.a = self.applyrefpxcorr(self.a)
        
        if self.refpxeachread and not self.norefpxcorr:
            print("Ref. px. corrected after each read")
        else:
            print("Ref. px. corrected on the slope")
        if self.selfbias:
            print("Self biased")
        else:
            print("super bias used")
        
        if self.toponly and not self.norefpxcorr:
            print("Only using top and side ref. px.")
        elif self.notopnobottom and not self.norefpxcorr:
            print("Only side ref. px. used")
        else:
            print("All ref. px. used")
        
        if self.oddeven and not self.norefpxcorr:
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
        Apply a quick preprocessing to the slope (horizontal-line and
        "butterfly" pattern removal, see hxrgutils.preprocess.preprocess).
        Replaces self.a.

        Returns
        -------
        None.

        """
        from hxrgutils.preprocess import preprocess
        if not hasattr(self,'a'):
            raise RuntimeError("nothing to preprocess, compute or load a ramp first")
        self.a = preprocess.preprocess(self.a)
    
