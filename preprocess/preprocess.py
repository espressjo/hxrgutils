#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:20:47 2023

@author: espressjo
"""

import numpy as np

def f_h_lines(im:np.ndarray):
    """
    Remove horizontal line commonly found after 1/f noise is removed due
    to capacitive coupling.

    Parameters
    ----------
    im : np.ndarray
        HXRG ndarray.

    Returns
    -------
    im : TYPE
        DESCRIPTION.

    """
    if im.dtype!=np.float64:
        im = im.astype(np.float64)
    h_l = np.median(im,axis=1)
    for i in range(im.shape[0]):
        im[i,:]-=h_l[i]
    return im
def butterfly(im:np.ndarray):
    """
    Calculate and subtract the butterfly effect we see on raw ramps.

    Parameters
    ----------
    im : np.ndarray
        Ramps.

    Returns
    -------
    ndarray
        Corrected image.
    butterfly : ndarray
        Calculated butterfly pattern.

    """
    if im.dtype!=np.float64:
        im = im.astype(np.float64)
    butterfly = np.zeros(im.shape)
    strip = [];
    for i in range(32):
        if i%2==0:
            strip.append(np.fliplr(im[:,i*128:(i+1)*128]))
        else:
            strip.append(im[:,i*128:(i+1)*128])
    _strip = np.median(strip,axis=0)
    for i in range(32):
        if i%2==0:
            butterfly[:,i*128:(i+1)*128] = np.fliplr(_strip)
        else:
            butterfly[:,i*128:(i+1)*128] = _strip
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