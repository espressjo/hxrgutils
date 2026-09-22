#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 22 2026

Small helper shared by hxread/hxramp. Both can be initialized either from a
FITS file (we get a header for free) or from a numpy array (we have no header
at all). Every method that touches a header goes through as_header() so that
the no-header case behaves exactly like an empty header instead of crashing.

@author: espressjo
"""
from astropy.io import fits


def as_header(value=None):
    """
    Coerce anything into a usable astropy Header *instance*.

    Parameters
    ----------
    value : fits.Header, dict, list of cards, None, ...
        Whatever the user (or an older version of this code) stored in the
        .H attribute. None, an unusable value or the fits.Header class itself
        all give an empty header.

    Returns
    -------
    fits.Header
        Always a header instance, never a class, never None.

    """
    if isinstance(value, fits.Header):
        return value
    if value is None:
        return fits.Header()
    try:
        return fits.Header(value)
    except Exception:
        return fits.Header()
