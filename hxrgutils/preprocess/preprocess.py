#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:20:47 2023

@author: espressjo
"""

import warnings

import numpy as np

def f_h_lines(im:np.ndarray,start=0,stop=None):
    """
    Remove horizontal line commonly found after 1/f noise is removed due
    to capacitive coupling.

    Parameters
    ----------
    im : np.ndarray
        HXRG ndarray, any format (H2RG, H4RG).
    start,stop : INT, optional
        Column range used to measure the level of each row. By default the
        full width of the array is used.

    Returns
    -------
    im : ndarray
        Image without its horizontal lines.

    """
    im = np.asarray(im,dtype=np.float64)
    if stop is None:
        stop = im.shape[1]
    with warnings.catch_warnings():
        #an all-NaN row simply gets no correction
        warnings.simplefilter("ignore",category=RuntimeWarning)
        level = np.nanmedian(im[:,start:stop],axis=1)
    level = np.where(np.isfinite(level),level,0.0)
    return im-level[:,None]
def butterfly(im:np.ndarray,amps=32):
    """
    Calculate and subtract the butterfly effect we see on raw ramps.

    Parameters
    ----------
    im : np.ndarray
        Ramps, any format (H2RG, H4RG).
    amps : INT, optional
        Number of amplifiers. The default is 32.

    Returns
    -------
    ndarray
        Corrected image.
    butterfly : ndarray
        Calculated butterfly pattern.

    """
    im = np.asarray(im,dtype=np.float64)
    if im.shape[1]%amps!=0:
        raise ValueError("width %d is not a multiple of the %d amp."%(im.shape[1],amps))
    w = im.shape[1]//amps#128 px on an H4RG, 64 px on an H2RG
    butterfly = np.zeros(im.shape)
    strip = [];
    for i in range(amps):
        if i%2==0:
            strip.append(np.fliplr(im[:,i*w:(i+1)*w]))
        else:
            strip.append(im[:,i*w:(i+1)*w])
    with warnings.catch_warnings():
        #an all-NaN position across every amp. gets no correction
        warnings.simplefilter("ignore",category=RuntimeWarning)
        _strip = np.nanmedian(strip,axis=0)
    _strip = np.where(np.isfinite(_strip),_strip,0.0)
    for i in range(amps):
        if i%2==0:
            butterfly[:,i*w:(i+1)*w] = np.fliplr(_strip)
        else:
            butterfly[:,i*w:(i+1)*w] = _strip
    return im-butterfly,butterfly
def preprocess(im:np.ndarray):
    """
    Remove the horizontal lines and the butterfly pattern on the given image

    Parameters
    ----------
    im : np.ndarray
        Ramps.

    Returns
    -------
    im : ndarray
        Preproceesed image.

    """
    im = f_h_lines(im)
    im,_ = butterfly(im)
    return im