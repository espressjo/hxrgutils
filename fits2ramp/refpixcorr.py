import numpy as np
import scipy.ndimage.filters


class refpxcorr():
    def __init__(self,data):
        self.im = data
    def nanmedian(self,data):
        # this function returns the median of finite values within
        # a vector. This is for python2 only and we will replace by
        # the python3 version np.nanmedian that does exactly the same
        # thing. When swithing to python3, we will simply add :
        #
        # import np.nanmedian as nanmedian 
        # 
        # and it should be completely transparent for the rest of the code.
        #
        data2=np.asarray(data)
        g=np.isfinite(data2)
        if np.max(g)==False:
            return(np.nan)
        return(np.median(data2[g]))
    
    def refpxcorr(self,oddeven=False):
        # function that corrects with reference pixels on the sides of the H2RG and H4RG.
        #
        # On the periphery of the arrays, there are 4 pixels that are not light-sensitive 
        # and that track drifts in the amplifiers. These are reference pixels and they can
        # reduce the effective readout noise by a factor of at least 2 if properly used.
        #
        # The top and bottom pixels of each output (one of 32 vertical "ribbons") see the
        # start and end of each readout. To filter noise on a readout timescale, we measure
        # the median of the top and bottom reference pixels. We then define a "slope" 
        # (matrix y_frac) that forces interpolates the gradient through the light-sensitive
        # pixels.
        #
        # For some arrays (e.g., the H2RG used for the AT4), the odd and even pixels within
        # each amplifier differ in behaviour. We therefore measure and correct this "slope"
        # independently for odd and even pixels. This is done by setting oddeven=True in the
        # function call. The default is oddeven=False
        #
        # The side (x=0-3 and x=N-3:N) of the HxRG arrays see the "faster" 1/f noise that 
        # affects all amplifier. We therefore need to subtract the mean of the side reference
        # pixels to remove (most of) the 1/f noise. As the reference pixels are themselves 
        # noisy, we apply a median filter to these pixels before subtracting.
        # The size of this running median filter is set with the "medfilterwidth" 
        # variable.
        #
    
        imdim=(np.shape(self.im))[0]
    
    
        # x position of the side deference pixels
        ref_sides = [0, 1, 2, 3,imdim - 4, imdim - 3, imdim - 2, imdim - 1]
        
        # filtering with ref pixels on either side of image
        medfilterwidth = 15 # value used for JWST H2RGs. Could be modified
    
        ref=np.zeros(imdim) # contains the median-filter, mean value of the vertical ref pixels
        for xpix in ref_sides:
            ref+=scipy.ndimage.filters.median_filter(self.im[:,xpix], medfilterwidth)/np.size(ref_sides)
        
        # pad the ref pixel value into a imdim x imdim square and subtract from image
        self.im-=np.repeat(ref,imdim).reshape(imdim,imdim) # correct an error, used to be "tile" instead of "repeat", which pads in the wrong direction
    
        # we filter independently the odd and even pixels in the bottom and top reference regions
        if oddeven:
            odd_bottom=np.zeros([imdim,imdim//32],dtype=float) # contains a range from 0 to 1 on odd pixels, 1 at bottom, 0 at top
            even_bottom=np.zeros([imdim,imdim//32],dtype=float)
            odd_top=np.zeros([imdim,imdim//32],dtype=float)
            even_top=np.zeros([imdim,imdim//32],dtype=float)
            
            g_odd_bottom=np.zeros([imdim,imdim//32],dtype=bool) # contains a range from 0 to 1 on odd pixels, 1 at bottom, 0 at top
            g_even_bottom=np.zeros([imdim,imdim//32],dtype=bool)
            g_odd_top=np.zeros([imdim,imdim//32],dtype=bool)
            g_even_top=np.zeros([imdim,imdim//32],dtype=bool)
        else:
            bottom=np.zeros([imdim,imdim//32],dtype=float) # contains a range from 0 to 1 on odd pixels, 1 at bottom, 0 at top
            top=np.zeros([imdim,imdim//32],dtype=float)
            
            g_bottom=np.zeros([imdim,imdim//32],dtype=bool) # contains a range from 0 to 1 on odd pixels, 1 at bottom, 0 at top
            g_top=np.zeros([imdim,imdim//32],dtype=bool)
        
        frac=np.asarray(range(imdim))/(imdim-1.0)
        for j in range(imdim//64):
            if oddeven:
                odd_bottom[:,j*2+1]=1-frac
                even_bottom[:,j*2]=1-frac
                odd_top[:,j*2+1]=frac
                even_top[:,j*2]=frac
            else:
                bottom[:,j*2+1]=1-frac
                bottom[:,j*2]=1-frac
                top[:,j*2+1]=frac
                top[:,j*2]=frac
            if oddeven:
                g_odd_bottom[0:4,j*2+1]=True # contains a range from 0 to 1 on odd pixels, 1 at bottom, 0 at top
                g_even_bottom[0:4,j*2]=True
                g_odd_top[imdim-4:imdim,j*2+1]=True
                g_even_top[imdim-4:imdim,j*2]=True
            else:
                g_bottom[0:4,j*2+1]=True
                g_bottom[0:4,j*2]=True
                g_top[imdim-4:imdim,j*2+1]=True
                g_top[imdim-4:imdim,j*2]=True
                    
        for j in range(32): # looping through the 32 outputs
            # subtract median value of ref unilluminated pixels
            ribbon = self.im[:,j*imdim//32:(j+1)*imdim//32]
            if oddeven:
                y_even_bottom = self.nanmedian(  ribbon[g_even_bottom])
                y_odd_bottom = self.nanmedian(  ribbon[g_odd_bottom])
                y_even_top = self.nanmedian(  ribbon[g_even_top])
                y_odd_top = self.nanmedian(  ribbon[g_odd_top])
                self.im[:,j*imdim//32:(j+1)*imdim//32]-=( y_even_bottom*even_bottom + y_odd_bottom*odd_bottom + y_odd_top*odd_top + y_even_top*even_top)

            else:
                y_bottom = self.nanmedian(  ribbon[g_bottom])
                y_top = self.nanmedian(  ribbon[g_top])
                self.im[:,j*imdim//32:(j+1)*imdim//32]-=( y_bottom*bottom  + y_top*top )
    
        return(self.im)
    def refpxcorrtop(self,oddeven=False):
        # function that corrects with reference pixels on the sides of the H2RG and H4RG.
        #
        # On the periphery of the arrays, there are 4 pixels that are not light-sensitive 
        # and that track drifts in the amplifiers. These are reference pixels and they can
        # reduce the effective readout noise by a factor of at least 2 if properly used.
        #
        # The top and bottom pixels of each output (one of 32 vertical "ribbons") see the
        # start and end of each readout. To filter noise on a readout timescale, we measure
        # the median of the top and bottom reference pixels. We then define a "slope" 
        # (matrix y_frac) that forces interpolates the gradient through the light-sensitive
        # pixels.
        #
        # For some arrays (e.g., the H2RG used for the AT4), the odd and even pixels within
        # each amplifier differ in behaviour. We therefore measure and correct this "slope"
        # independently for odd and even pixels. This is done by setting oddeven=True in the
        # function call. The default is oddeven=False
        #
        # The side (x=0-3 and x=N-3:N) of the HxRG arrays see the "faster" 1/f noise that 
        # affects all amplifier. We therefore need to subtract the mean of the side reference
        # pixels to remove (most of) the 1/f noise. As the reference pixels are themselves 
        # noisy, we apply a median filter to these pixels before subtracting.
        # The size of this running median filter is set with the "medfilterwidth" 
        # variable.
        #
    
        imdim=(np.shape(self.im))[0]
    
    
        # x position of the side deference pixels
        ref_sides = [0, 1, 2, 3,imdim - 4, imdim - 3, imdim - 2, imdim - 1]
        
        # filtering with ref pixels on either side of image
        medfilterwidth = 15 # value used for JWST H2RGs. Could be modified
    
        ref=np.zeros(imdim) # contains the median-filter, mean value of the vertical ref pixels
        for xpix in ref_sides:
            ref+=scipy.ndimage.filters.median_filter(self.im[:,xpix], medfilterwidth)/np.size(ref_sides)
        
        # pad the ref pixel value into a imdim x imdim square and subtract from image
        self.im-=np.repeat(ref,imdim).reshape(imdim,imdim) # correct an error, used to be "tile" instead of "repeat", which pads in the wrong direction
    
        # we filter independently the odd and even pixels in the bottom and top reference regions
        #odd_bottom=np.zeros([imdim,imdim//32],dtype=float) # contains a range from 0 to 1 on odd pixels, 1 at bottom, 0 at top
        #even_bottom=np.zeros([imdim,imdim//32],dtype=float)
        if oddeven:
            odd_top=np.zeros([imdim,imdim//32],dtype=float)
            even_top=np.zeros([imdim,imdim//32],dtype=float)
            
            #g_odd_bottom=np.zeros([imdim,imdim//32],dtype=bool) # contains a range from 0 to 1 on odd pixels, 1 at bottom, 0 at top
            #g_even_bottom=np.zeros([imdim,imdim//32],dtype=bool)
            g_odd_top=np.zeros([imdim,imdim//32],dtype=bool)
            g_even_top=np.zeros([imdim,imdim//32],dtype=bool)
        else:
            top=np.ones([imdim,imdim//32],dtype=float)
            g_top=np.zeros([imdim,imdim//32],dtype=bool)
        
        #frac=np.asarray(range(imdim))/(imdim-1.0)
        for j in range(imdim//64):
            #odd_bottom[:,j*2+1]=1-frac
            #even_bottom[:,j*2]=1-frac
            if oddeven:
                odd_top[:,j*2+1]=1.0
                even_top[:,j*2]=1.0
                g_odd_top[imdim-4:imdim,j*2+1]=True
                g_even_top[imdim-4:imdim,j*2]=True
            else:
                g_top[imdim-4:imdim,j*2+1]=True
                g_top[imdim-4:imdim,j*2]=True
            #g_odd_bottom[0:4,j*2+1]=True # contains a range from 0 to 1 on odd pixels, 1 at bottom, 0 at top
            #g_even_bottom[0:4,j*2]=True
            
            
                    
        for j in range(32): # looping through the 32 outputs
            # subtract median value of ref unilluminated pixels
            ribbon = self.im[:,j*imdim//32:(j+1)*imdim//32]
            #y_even_bottom = self.nanmedian(  ribbon[g_even_bottom])
            #y_odd_bottom = self.nanmedian(  ribbon[g_odd_bottom])
            if oddeven:
                y_even_top = self.nanmedian(  ribbon[g_even_top])
                y_odd_top = self.nanmedian(  ribbon[g_odd_top])
                self.im[:,j*imdim//32:(j+1)*imdim//32]-=( y_odd_top*odd_top + y_even_top*even_top)
            else:
                y_top = self.nanmedian(  ribbon[g_top])

                self.im[:,j*imdim//32:(j+1)*imdim//32]-=( y_top*top )
        return(self.im)
    def refpxcorrside(self,oddeven=False):
        # function that corrects with reference pixels on the sides of the H2RG and H4RG.
        #
        # On the periphery of the arrays, there are 4 pixels that are not light-sensitive 
        # and that track drifts in the amplifiers. These are reference pixels and they can
        # reduce the effective readout noise by a factor of at least 2 if properly used.
        #
        # The top and bottom pixels of each output (one of 32 vertical "ribbons") see the
        # start and end of each readout. To filter noise on a readout timescale, we measure
        # the median of the top and bottom reference pixels. We then define a "slope" 
        # (matrix y_frac) that forces interpolates the gradient through the light-sensitive
        # pixels.
        #
        # For some arrays (e.g., the H2RG used for the AT4), the odd and even pixels within
        # each amplifier differ in behaviour. We therefore measure and correct this "slope"
        # independently for odd and even pixels. This is done by setting oddeven=True in the
        # function call. The default is oddeven=False
        #
        # The side (x=0-3 and x=N-3:N) of the HxRG arrays see the "faster" 1/f noise that 
        # affects all amplifier. We therefore need to subtract the mean of the side reference
        # pixels to remove (most of) the 1/f noise. As the reference pixels are themselves 
        # noisy, we apply a median filter to these pixels before subtracting.
        # The size of this running median filter is set with the "medfilterwidth" 
        # variable.
        #
    
        imdim=(np.shape(self.im))[0]
    
    
        # x position of the side deference pixels
        ref_sides = [0, 1, 2, 3,imdim - 4, imdim - 3, imdim - 2, imdim - 1]
        
        # filtering with ref pixels on either side of image
        medfilterwidth = 15 # value used for JWST H2RGs. Could be modified
    
        ref=np.zeros(imdim) # contains the median-filter, mean value of the vertical ref pixels
        for xpix in ref_sides:
            ref+=scipy.ndimage.filters.median_filter(self.im[:,xpix], medfilterwidth)/np.size(ref_sides)
        
        # pad the ref pixel value into a imdim x imdim square and subtract from image
        self.im-=np.repeat(ref,imdim).reshape(imdim,imdim) # correct an error, used to be "tile" instead of "repeat", which pads in the wrong direction
        return(self.im)
    def patch_shift(self,im,bias):

        # this bit of code 
        index=np.asarray(range(4,60,2))
        cut1 = 0.2 # max CC for shifts that are invalid
        cut2 = 0.9 # min CC for shifts that is valid
        ccs = np.zeros(3)
    
        print(np.shape(im))
        i=0
        for off in range(-1,2):
            ccs[i]= (pearsonr(im[0,index],bias[0,off+index]))[0]
            i+=1
    
        message = 'Ambiguous Pearson correlation with bias... suspicious data!'
    
        if (ccs[2] >= cut2) and (ccs[1]<=cut1) and (ccs[0]<=cut1):
            message='We have a pixel shift problem... we correct it!'
    
            xpix2=np.asarray(range(2048))
            xpix=np.asarray(range(2048))
            x64=np.asarray(range(64))
            for i in range(32):
                xpix[i*64:i*64+64]=(i*32)+((x64+(2*(i % 2)-1) ) % 64)
            im[:,xpix2]=im[:,xpix]
    
    
        if (ccs[1] >= cut2) and (ccs[2]<=cut1) and (ccs[0]<=cut1):
            message = 'all good, there is no mischievous pixel shift in your data!'
    
        print(message)    
        return(im)
if '__main__' in __name__:
    f2r = fits2ramp()
    from fits_images import loadIm
    f = '/home/noboru/NIRPS_R01_R01.fits' 
    print(type(loadIm(f)[0,0]))
    f2r.refpixcorr(loadIm(f))
    
    
    
#from scipy import stats
#import numpy as np
#from array import *
#import glob
#import os
## import pyfits --> rendered obsolete by the use of the more recent astropy.io.fits
#import time
#import sys
#import scipy.ndimage.filters
#from astropy.io import fits as pyfits
#from scipy.stats.stats import pearsonr
#
#
#def nanmedian(data):
#    # this function returns the median of finite values within
#    # a vector. This is for python2 only and we will replace by
#    # the python3 version np.nanmedian that does exactly the same
#    # thing. When swithing to python3, we will simply add :
#    #
#    # import np.nanmedian as nanmedian 
#    # 
#    # and it should be completely transparent for the rest of the code.
#    #
#    data2=np.asarray(data)
#    g=np.isfinite(data2)
#    if np.max(g)==False:
#        return(np.nan)
#    return(np.median(data2[g]))
#
#def refpixcorr(im,oddeven=False):
#    # function that corrects with reference pixels on the sides of the H2RG and H4RG.
#    #
#    # On the periphery of the arrays, there are 4 pixels that are not light-sensitive 
#    # and that track drifts in the amplifiers. These are reference pixels and they can
#    # reduce the effective readout noise by a factor of at least 2 if properly used.
#    #
#    # The top and bottom pixels of each output (one of 32 vertical "ribbons") see the
#    # start and end of each readout. To filter noise on a readout timescale, we measure
#    # the median of the top and bottom reference pixels. We then define a "slope" 
#    # (matrix y_frac) that forces interpolates the gradient through the light-sensitive
#    # pixels.
#    #
#    # For some arrays (e.g., the H2RG used for the AT4), the odd and even pixels within
#    # each amplifier differ in behaviour. We therefore measure and correct this "slope"
#    # independently for odd and even pixels. This is done by setting oddeven=True in the
#    # function call. The default is oddeven=False
#    #
#    # The side (x=0-3 and x=N-3:N) of the HxRG arrays see the "faster" 1/f noise that 
#    # affects all amplifier. We therefore need to subtract the mean of the side reference
#    # pixels to remove (most of) the 1/f noise. As the reference pixels are themselves 
#    # noisy, we apply a median filter to these pixels before subtracting.
#    # The size of this running median filter is set with the "medfilterwidth" 
#    # variable.
#    #
#
#    imdim=(np.shape(im))[0]
#
#
#    # x position of the side deference pixels
#    ref_sides = [0, 1, 2, 3,imdim - 4, imdim - 3, imdim - 2, imdim - 1]
#    
#    # filtering with ref pixels on either side of image
#    medfilterwidth = 15 # value used for JWST H2RGs. Could be modified
#
#    ref=np.zeros(imdim) # contains the median-filter, mean value of the vertical ref pixels
#    for xpix in ref_sides:
#        ref+=scipy.ndimage.filters.median_filter(im[:,xpix], medfilterwidth)/np.size(ref_sides)
#    
#    # pad the ref pixel value into a imdim x imdim square and subtract from image
#    im-=np.repeat(ref,imdim).reshape(imdim,imdim) # correct an error, used to be "tile" instead of "repeat", which pads in the wrong direction
#
#    # we filter independently the odd and even pixels in the bottom and top reference regions
#    odd_bottom=np.zeros([imdim,imdim//32],dtype=float) # contains a range from 0 to 1 on odd pixels, 1 at bottom, 0 at top
#    even_bottom=np.zeros([imdim,imdim//32],dtype=float)
#    odd_top=np.zeros([imdim,imdim//32],dtype=float)
#    even_top=np.zeros([imdim,imdim//32],dtype=float)
#    
#    g_odd_bottom=np.zeros([imdim,imdim//32],dtype=bool) # contains a range from 0 to 1 on odd pixels, 1 at bottom, 0 at top
#    g_even_bottom=np.zeros([imdim,imdim//32],dtype=bool)
#    g_odd_top=np.zeros([imdim,imdim//32],dtype=bool)
#    g_even_top=np.zeros([imdim,imdim//32],dtype=bool)
#    
#    frac=np.asarray(range(imdim))/(imdim-1.0)
#    for j in range(imdim//64):
#        odd_bottom[:,j*2+1]=1-frac
#        even_bottom[:,j*2]=1-frac
#        odd_top[:,j*2+1]=frac
#        even_top[:,j*2]=frac
#    
#        g_odd_bottom[0:4,j*2+1]=True # contains a range from 0 to 1 on odd pixels, 1 at bottom, 0 at top
#        g_even_bottom[0:4,j*2]=True
#        g_odd_top[imdim-4:imdim,j*2+1]=True
#        g_even_top[imdim-4:imdim,j*2]=True
#        
#                
#    for j in range(32): # looping through the 32 outputs
#        # subtract median value of ref unilluminated pixels
#        ribbon = im[:,j*imdim//32:(j+1)*imdim//32]
#        y_even_bottom = nanmedian(  ribbon[g_even_bottom])
#        y_odd_bottom = nanmedian(  ribbon[g_odd_bottom])
#        y_even_top = nanmedian(  ribbon[g_even_top])
#        y_odd_top = nanmedian(  ribbon[g_odd_top])
#        im[:,j*imdim//32:(j+1)*imdim//32]-=( y_even_bottom*even_bottom+y_odd_bottom*odd_bottom+y_odd_top*odd_top+y_even_top*even_top)
#
#    return(im)
#
#
#def patch_shift(im,bias):
#
#    # this bit of code 
#    index=np.asarray(range(4,60,2))
#    cut1 = 0.2 # max CC for shifts that are invalid
#    cut2 = 0.9 # min CC for shifts that is valid
#    ccs = np.zeros(3)
#
#    print(np.shape(im))
#    i=0
#    for off in range(-1,2):
#        ccs[i]= (pearsonr(im[0,index],bias[0,off+index]))[0]
#        i+=1
#
#    message = 'Ambiguous Pearson correlation with bias... suspicious data!'
#
#    if (ccs[2] >= cut2) and (ccs[1]<=cut1) and (ccs[0]<=cut1):
#        message='We have a pixel shift problem... we correct it!'
#
#        xpix2=np.asarray(range(2048))
#        xpix=np.asarray(range(2048))
#        x64=np.asarray(range(64))
#        for i in range(32):
#            xpix[i*64:i*64+64]=(i*32)+((x64+(2*(i % 2)-1) ) % 64)
#        im[:,xpix2]=im[:,xpix]
#
#
#    if (ccs[1] >= cut2) and (ccs[2]<=cut1) and (ccs[0]<=cut1):
#        message = 'all good, there is no mischievous pixel shift in your data!'
#
#    print(message)    
#    return(im)
#
## will be set to True if selfbias=True. If we use a file for bias (later update?) then this will also
## change the dobias to True
#
#dobias = False
#
#
#arg=np.asarray(sys.argv)
#arg=arg[1:] # first argument is simply the name of the program and needs to be removed
#
#write_cube = sum(arg=='-cube') ==1. # if set, then we will write cube, if not, then we skip this step that may be long
#skip_error = sum(arg=='-noerror') ==1. # if set, we skip slope error
#skip_ref   = sum(arg=='-noref') ==1. # if set, we skip reference pixel corrections
#linearize = sum(arg=='-linearize') ==1. # if set, we correct for non-linearity
#selfbias = sum(arg=='-selfbias') ==1. # if set, we correct ref pixels on a frame-to-frame basis
#
#nmax_set=False
#for argn in arg:
#    if (argn)[0:3] == '-n=':
#        nmax_set=True
#        dim3=np.int( (argn)[3:] )
#
## here we remove arguments with a "-"
#keep=np.zeros(len(arg))
#for i in range(len(arg)):
#    keep[i] = (arg[i])[0] != '-'
#
#arg=arg[keep ==1]  # keep only params not beginning with a "-"
#
#if len(arg)>=1:
#    odometer = arg[0] # first argument after program and flags is the output name
#fic = arg[1:]
#
#if len(fic)>=1:
#    h = pyfits.getheader(fic[0])
#    h2=h
#
#mef_flag=0 # file is a MEF flag
#cubefits_flag=0 # file is a CUBE flag
#
#if len(fic) ==1:
#    naxis =h['naxis']
#    if naxis ==0:
#        mef_flag=1# we have a flag to know that the input file is a MEF and that extensions need to be read from there
#    if naxis==3:
#        cubefits_flag=1#this is a cuube
#
#
#exists = np.zeros(len(fic),dtype=bool)
#for i in range(len(fic)):
#    exists[i] = os.path.isfile(fic[i])
#
#if np.sum(exists ==0) !=0:
#    print('some files given as inputs do not exist')
#    print('missing file(s) --')
#    print('')
#    missing=fic[exists !=1]
#    for i in range(len(missing)):
#        print(missing[i])
#    print('')
#    print('... you way also have given some erroneous input, double check your inputs dude!')
#    sys.exit()
#
#if len(sys.argv) <=2:
#    print('***** !!! warning, something went wrong !!! *****')
#    print('')
#    print('            ----- you can provide a list of files as an input -----')
#    print('')
#    print('syntax : python fits2ramp.py outname directory/file*.fits -cube -noerror -linearize')
#    print('')
#    print('')
#    print(' the argument after the "outname" must be the files to combine')
#    print(' with the ramp-fitting algorithm. ex: 20170322140210/H2RG_R01_M01_N08*.fits ')
#    print(' should also accept *.fits.gz files')
#    print(' you need at least two files in the wildcard. You can also expliclty')
#    print(' name the files you combine.')
#    print(' The syntax would be :')
#    print('  python fits2ramp.py outname file1.fits file2.fits ... fileN.fits')
#    print('')
#    print('    ----- you can also provide a single file that has a MEF format -----')
#    print('')
#    print('syntax : python fits2ramp.py outname mef_file*.fits -cube -noerror -linearize')
#    print('')
#    print(' if you provide an outname and a single fits file, then we know its a MEF')
#    print('')
#    print(' if you provide a -n=XXXX then only the first XXXX readouts within the MEF')
#    print('')
#    print(' will be used for slope fitting')
#    print('                      ---- some more options ----'  )
#    print('')
#    print('  -cube saves all slices in a cube. This is slower and takes disk space')
#    print('  -noerror does not compute the slope error. This is faster.' )
#    print('  -linearize corrects for non-linearity. This is slower but more accurate.')
#    print('')
#    print(' If all goes well, the programs outputs 2 files: ')
#    print(' outnameo.fits ')
#    print('   ... ext=1, ramp frame' )
#    print('   ... ext=2, ramp intercept')
#    print('   ... ext=3, ramp error' )
#    print('   ... ext=4, ramp # valid frames')
#    print('   ... every where, NaN values trace saturated pixel')
#    print(' outnamer.fits.gz')
#    print('   ... cube with as many slices as there are files in the wildcard above')
#    print('   ... outnamer.fits.gz contains the same info as the files' )
#    print('   ... this is only done if we pass the "-cube" argument')
#    print('')
#    sys.exit()
#
#
#
#
##################################################################
##################################################################
#
## We need the size of the image. Should be 2048 or 4096 (H2RG/H4RG)
#imdim=(np.shape(pyfits.getdata(fic[0])))[1]
#
#if (imdim!=2048) and (imdim!=4096):
#    print('')
#    print('')
#    print('     something is really wrong with the size of the input image')
#    print('     the image    '+fic[0]+'    has a width of :',imdim,' pixel(s)')
#    print('     and we should only have values of 2048 or 4096 pixels')
#    print('')
#    print('')
#    sys.exit()
#    
#
#
## reading the relevant calibrations
##mask = getdata(calibdir+'/mask.fits') # 0/1 mask defining the area of the science array used as pose-meter
#mask=np.zeros([imdim,imdim],dtype=float) # dummy ~~~>>> will need to be changed for the H4RG
#
##  this is the region used for the posemeter
## For SPIRou, we will have a binary mask selecting the H-band orders (science and not ref channel)
#mask[1912:1938,572:777]=1
#
#mask=np.where(mask ==1)
#
## non-linearity cube with 4 slices. The linearized flux will be derived from the measured flux with the
## following relation :
##  F_lin = a0 + a1*(F_mea - bias) + a2*(F_mea - bias)**2 + a3*(F_mea - bias)**3
## where aN is the Nth slice of the linearity cube
## ... bias is the super-bias
## ... F_lin is the linearised flux
## ... F_mea is the measured flux
##linearity = getdata(calibdir+'/non_lin.fits')  # we will use files with non-linearity correction here
#
## This is an operation that may be done if we do not have a bias in hand and want to
## correct non-linearity. Lets consider this under development and set it to False for now
##
#
#
#linearity_saturation = pyfits.getdata('nonlin.fits')
#
## Slice 1 - 2nd ordre term of non-linearity correction
## Slice 2 - 3rd ordre term of non-linearity correction
#linearity = linearity_saturation[0:2,:,:]
#
## Slice 3 - dynamical range for <20% non-linearity
#saturation = linearity_saturation[2,:,:]
#
#
#if mef_flag==0 and cubefits_flag==0:
#    if nmax_set == False:
#        dim3 = len(fic)
#    else:
#        if len(fic) < dim3:
#            print('You requested a ramp of ',dim3,' readouts... ')
#            print(' ... but you have only ',len(fic),' files')
#            sys.exit()
#
#if mef_flag==1:
#    hdulist = pyfits.open(fic[0],memmap=False) ## We will use memmap when CFHT gets rid of  BZERO/BSCALE/BLANK header keywords
#    dims=np.shape(hdulist[1])
#    if nmax_set == False:
#        dim3= len(hdulist)-1
#    else:
#        if (len(hdulist)-1) < dim3:
#            print('You requested a ramp of ',dim3,' readouts... ')
#            print(' ... but you have only ',len(hdulist)-1,' slices in your MEF')
#            sys.exit()
#
#if cubefits_flag==1:
#    if nmax_set == False:
#        dim3 = h['naxis3']
#    else:
#        if (h['naxis3']) < dim3:
#            print('You requested a ramp of ',dim3,' readouts... ')
#            print(' ... but you have only ',len(hdulist)-1,' slices in your cube')
#            sys.exit()
#
#
## delete all keywords  from the reference file
#del_keywords=['DATLEVEL', 'ASICGAIN', 'NOMGAIN', 'AMPRESET', 'KTCREMOV', 'SRCCUR',\
#    'AMPINPUT', 'V4V3V2V1', 'PDDECTOR', 'CLKOFF', 'NADCS', 'INTTIME',\
#    'TSTATION', 'SEQNUM_N', 'SEQNUM_M', 'CLOCKING', 'NEXTRAP','NEXTRAL', 'SEQNNAME']
#
#
#for key in del_keywords:
#    if key in h: # as keywords may change from version to version, we check if the keyword we want to delete is present
#        del h[key]
#del h['bias*']
#
#timestamp=np.zeros(dim3,dtype=float)
#
#
## loop to check image size and populate header with time stamps
#for i in range(dim3):
#    if mef_flag==0 and cubefits_flag==0: # we have a mef file, info is in the ith extension
#        h_tmp = pyfits.getheader(fic[i])
#
#        if 'frmtime' not in h_tmp:
#            h_tmp['frmtime'] = 5.24288, 'assumed integration time (s)'
#
#        if 'inttime' not in h_tmp:
#            h_tmp['inttime'] = 5.24288*(i+1), 'assumed frame time (s)'
#
#        timestamp[i]=h_tmp['inttime']
#
#    if cubefits_flag==1: # we have a cube, calculate from FRMTIME
#        timestamp[i]= (i+1)*h['frmtime'] # sets zero time at the time of reset
#
#    if mef_flag==1: # we read the ith extension
#        h_tmp = hdulist[i+1].header
#        timestamp[i]=h_tmp['inttime']
#
#if mef_flag==0 and cubefits_flag==0:
#    order = np.argsort(timestamp) # who knows, the files may not be in the right order! Lets sort them according to their timestamps
#    fic=fic[order]
#    timestamp=timestamp[order]
#
#for i in range(dim3):
#    tag0 = str(i+1)
#    if len(tag0) < 4:
#        tag = '0'*(4-len(tag0))+tag0
#    tag = 'INTT'+tag
#    h[tag] = (timestamp[i],'Timestamp, '+tag0+'/'+str(dim3))
#
#
#if mef_flag==1:
#    write_cube=False
#
#if write_cube:
#    cube=np.zeros([dim3,dim2,dim1],dtype=float)
#    print('loading all files in cube')
#    for i in range(dim3):
#        print(i+1,'/',len(fic),fic[i])
#        im=pyfits.getdata(fic[i])
#        cube[i,:,:] = im
#
#    print('writing the cube file --> '+odometer+'r.fits ')
#    t1 = time.time()
#
#    hcube=h2
#    hcube['NAXIS'] = 3
#    hcube['NAXIS3'] = dim3
#    pyfits.writeto(odometer+'r.fits', cube,header=hcube)
#
#    # This operation is somewhat long and could lead to back-log of files on a slow machine
#    # ... for the code development, we time it. This may be remove at a later point.
#    print('Duration of file writting : '+str(float(time.time()-t1))+' s')
#
#    # zipping the .fits file. Normally this could be done within pyfits.writeto, but its much, much slower
#    os.system('gzip -f '+odometer+'r.fits  &')
#
#    print('done writing the cube file --> '+odometer+'r.fits')
#    print(' compressing file in background ... ')
#    del cube # removing cube from memory to make things lighter... unclear in necessary
#
#else:
#    print('we do not write the cube file for this ramp')
#
## place htimestampolders for some arithmetics for the linear fit
##sx = 0#np.zeros([dim2,dim1])
##sx2 = 0#np.zeros([dim2,dim1])
#sy = np.zeros([imdim,imdim],dtype=float)
#n = np.zeros([imdim,imdim],dtype=np.int16)
#
#sxy = np.zeros([imdim,imdim],dtype=float)
#fmask = np.zeros(dim3,dtype=float)
#
## mask for pixels that are valid
#goodmask = np.full((imdim,imdim),True,dtype=bool)
## when a pixels goes above saturation, it remains invalid for the rest of the ramp
#
#
#if skip_error == False:
#    savname=['']*dim3
#
#
#print(mef_flag,cubefits_flag,linearize)
#
#
#
#t_start=time.time()
#for i in range(dim3):
#    t0=time.time()
#    print(i+1,'/',dim3,' ~~~> Computing slope')
#
#    if mef_flag==0 and cubefits_flag==0: # this is a set with N files
#        im = pyfits.getdata(fic[i])
#
#    if mef_flag==1:
#        im=hdulist[i+1].data # reading the Nth extension
#
#    if cubefits_flag==1:
#        if i ==0:
#            bigcube=pyfits.getdata(fic[0]) # that's dangerous as it may overfill memory
#        im=bigcube[i,:,:]
#
#    im = np.array(im,dtype='float')
#
#    
#    if selfbias and (i ==0):
#        bias = np.array(im)
#        print('setting 1st extension as a bias file')
#        dobias=True
#
#    goodmask = (im <= saturation)*goodmask
#
#    if dobias:
#        if selfbias:
#            print('bias subtraction with 1st readout')
#        else:
#            print('bias subtraction with provided bias file')
#        im-=bias
#        
#    if linearize:
#        print('applying non-lin correction')
#        # first we linearize the data by applying the non-linearity coefficients and bias correction
#        for j  in range(2):
#            im += linearity[j,:,:]*(im)**(j+2)
#
#    if selfbias and (skip_ref == False):
#        print('as we applied self-bias, we correct ref pixels')
#        im=refpixcorr(im)
#
#    
#    n+= goodmask
#    fmask[i]=np.nanmean( im[mask])
#    # m*=goodmask # starting now, only the product of the two is needed. saves one multipltication
#    # Actually, best not fill what used to be saturated elements in the array with
#    # 0, which is what this did.  Then, if the errslope calculation wants to check
#    # im <= saturation as it used to do, it will come up with the wrong answer.
#    # Since the first check for im <= saturation (about 20 lines above) does so
#    # before linearity correction and this check would be after, they could also
#    # come up with different answers though, unless the linearity function is
#    # is guaranteed to apply a correction that keeps saturation values at the same
#    # ADU.  Since we already have n[], when the errslope calculation happens, it
#    # uses that, now with a simple "goodmask = (n > i)" for each i on that pass.
#    sy[goodmask]+= im[goodmask]#*goodmask
#    sxy[goodmask]+=(im[goodmask]*timestamp[i])
#
#    # here we save the non-linearity corrected images as python npz files
#    # we could just dump everything into a big cube to be used in the slope
#    # error determination. We opt to write these files to disk to avoid overfilling
#    # the memory. This should be safer for very large number of reads.
#    #
#    # We cannot simply re-read the fits files are the "im" variable saved in the npz has been corrected for
#    # non-linearity, which is NOT the case for the .fits.gz. We save the NPZ only if the data is linearized
#    #
#    # We also corrected for the bias regions of the detector, so a temporary file is necessary if we want to properly compute slope error
#    # and cannot afford to keep everything in memory. Keeping everything in memory may be fine for small datasets, but we want
#    # to avoid having a code that crashes for long sequences or on machines with less memory!
#
#
#    if skip_error == False:
#        savname[i]='.tmp'+str(i)+'.npz'
#        np.savez(savname[i],im=im)  # this file is temporary and will be deleted after computing the slope error
#        
#    dt=(time.time()-t_start)/(i+1.0)
#    print('dt[last image] ','{:5.2f}'.format(time.time()-t0),'s; dt[mean/image] ','{:5.2f}'.format(dt),'s; estimated time left '+'{:3.0f}'.format(np.floor((dim3-i)*dt/60))+'m'+'{:2.0f}'.format(np.floor((dim3-i)*dt % 60))+'s')
#
#
#
## we now have these variables outside the loop. We keep n that contains the 
## number of valid reads, and directely interpolate the vector with the cumulative
## sum of timestamp and timestamp**2. Previously, we added these values to the sx and sx2
## matrices for each frame. This operation is much, much faster and equivalent.
#sx=np.where(n>0,(np.cumsum(timestamp))[n-1],0)
#sx2=np.where(n>0,(np.cumsum(timestamp**2))[n-1],0)
#
#if mef_flag==1:
#    hdulist.close()
#
#
#fmask-=fmask[0]
#for i in range(dim3):
#    tag0 = str(i+1)
#    if len(tag0) < 4:
#        tag = '0'*(4-len(tag))+tag0
#    tag = 'POSE'+tag
#    h[tag] =  (fmask[i],'Posemeter, '+tag0+'/'+str(len(fic)))
#
#
#a = np.zeros([imdim,imdim],dtype=float)+np.nan # slope, NaN if not enough valid readouts
#b = np.zeros([imdim,imdim],dtype=float)+np.nan # intercept
#
#valid=n>1 # only valid where there's more than one good readout(s)
#b[valid] = (sx*sxy-sx2*sy)[valid]/(sx**2-n*sx2)[valid] # algebra of the linear fit
#a[valid] = (sy-n*b)[valid]/sx[valid]
#
## For the sake of consistency, we fix the slope, error and intercept to NaN for
## pixels that have 0 or 1 valid (i.e., not saturated) values and for which
## one cannot determine a valid slope
#
#
#
#errslope = np.zeros([imdim,imdim],dtype=float)+np.nan
#
#goodmask = np.full((imdim,imdim),True,dtype=bool)
#
#if skip_error == False:
#    varx2 = np.zeros([imdim,imdim],dtype=float)
#    vary2 = np.zeros([imdim,imdim],dtype=float)
#    xp = np.zeros([imdim,imdim],dtype=float)
#
#    valid = (n>2)
#    xp[valid]=sx[valid]/n[valid] # used in the determination of error below
#
#    print('we now compute the standard error on the slope')
#    for i in range(dim3):
#
#        # we read the npz as this file has been linearized (if the -linearize keyword has been set)
#        # and we subtracted the reference regions on the array
#        data=np.load(savname[i])
#        os.system('rm '+savname[i])
#
#        im=data['im']
#
#        goodmask = (n > i)
#
#        yp = b+a*timestamp[i]
#        print(i+1,'/',dim3,' ~~~> Computing slope error')
#
#        varx2+= ((timestamp[i]-xp)**2)*goodmask # we multiply by goodmask so that only
#        vary2+= ((im-yp)**2)*goodmask
#
#    valid*=(varx2!=0) # avoid diving by zero
#    errslope[valid] = np.sqrt(vary2[valid]/(n[valid]-2))/np.sqrt(varx2[valid])
#    # deleting the temporary npz
#else:
#    print(' We do not calculate the error on slope.')
#    print(' This is faster and intended for debugging but ')
#    print(' ultimately we will want to compute slope error ')
#    print(' for all files')
#
#
#h['satur1']=(nanmedian(saturation),'median saturation limit in ADU')
#h['satur2']=(nanmedian(saturation)/max(timestamp),'median saturation limit in ADU/s')
#
#dfmask = fmask[1:]-fmask[0:-1] # flux received between readouts
#dtimestamp = timestamp[1:]+0.5*(timestamp[-1]-timestamp[0])/(len(timestamp)-1) # mid-time of Nth readout
#
#### we estimate the RON by checking the slope error in pixels receiving little flux
#### as the orders cover ~50% of the science array, we take the median slope error of
#### pixels that are below the median slope. We assume that these pixels have an RMS that is
#### dominated by readout noise (TO BE CONFIRMED).
#### we also clip pixels that are above 3x the median RMS
#
#pseudodark = 0.0 # (a < np.median(a))*(errslope < 3*np.median(errslope))
#ron_estimate = 0.0 #np.median(errslope[pseudodark])*(max(timestamp)-min(timestamp)) # converted into ADU instead of ADU/s
#
#####   Standard FITS Keywords BITPIX = 16 / 16bit
#h['BSCALE']=(1.0 , 'Scale factor')
#####  FITS keyword related to the detector
#h['RON_EST']=(ron_estimate , '[ADU] read noise estimate')
#h['NSUBEXPS']=(len(fic) , 'Total number of sub-exposures of 5.5s ')
##h['TMID']= (np.sum(dtimestamp*dfmask)/np.sum(dfmask) , '[s] Flux-weighted mid-exposure time ' )
##h['CMEAN']= ( np.mean(dfmask)/(timestamp[1]-timestamp[0]), '[ADU/s] Average count posemeter' )
#
#if skip_ref == False:
#    a=refpixcorr(a,oddeven=True)
#
#a=np.float32(a)
#
#if dobias:
#    # we subtracted the bias from all frames, we need to add it to the intercept
#    b+=bias
#
#b=np.float32(b)
#errslope=np.float32(errslope)
#
#hdu1 = pyfits.PrimaryHDU()
#hdu1.header = h
#hdu1.header['NEXTEND'] = 4
#hdu2 = pyfits.ImageHDU(a)
#hdu2.header['UNITS'] = ('ADU/S','Slope of fit, flux vs time')
#hdu2.header['EXTNAME'] = ('slope','Slope of fit, flux vs time')
#
#hdu3 = pyfits.ImageHDU(b)
#hdu3.header['UNITS'] =  ('ADU','Intercept of the pixel/time fit.')
#hdu3.header['EXTNAME'] = ('intercept','Intercept of the pixel/time fit.')
#
#hdu4 = pyfits.ImageHDU(errslope)
#hdu4.header['UNITS'] = ('ADU/S','Formal error on slope fit')
#hdu4.header['EXTNAME'] = ('errslope','Formal error on slope fit')
#
#hdu5 = pyfits.ImageHDU(n)
#hdu5.header['UNITS'] = ('Nimages','N readouts below saturation')
#hdu5.header['EXTNAME'] = ('count','N readouts below saturation')
#
#new_hdul = pyfits.HDUList([hdu1, hdu2, hdu3, hdu4, hdu5])
#
## just to avoid an error message with writeto
#if os.path.isfile(odometer+'.fits'):
#    print('file : '+odometer+'.fits exists, we are overwriting it')
#    os.system('rm '+odometer+'.fits')
#
#new_hdul.writeto(odometer +'.fits', clobber=True)
#
#print('Elapsed time for entire fits2ramp : '+str(float(time.time()-t0))+' s')