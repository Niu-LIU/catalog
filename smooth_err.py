#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: smoothed_error.py
"""
Created on Fri Aug 31 18:59:49 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np


__all__ = ["smooth_by_dec", "smooth_by_time"]


# -----------------------------  FUNCTIONS -----------------------------
def smooth_by_dec(dec, sig_pos_max, interv_size=50):
    """Calculate the smoothed standard error as a function of declination

    Parameters
    ----------
    dec : array_like of float type
        declination in degree
    sig_pos_max : array_like of float type
        positional uncertainty (ellpise major axis) in mas
    interv_size : int
        number of source in a subset, default is 50.

    Returns
    ----------
    med_dec : array of float
        median declination of every subset
    med_poserr : array of int
        median positional uncertainty of every subset
    """

    # Sort the data according to the declination
    ind = np.argsort(dec)
    dec_sort = np.take(dec, ind)
    poserr_sort = np.take(sig_pos_max, ind)

    # Then we calculate the median ellpise major axis for subset of 50 sources.
    if dec_sort.size <= interv_size:
        print("# Too smaller sample!")
        exit()
    else:
        interv_num = dec_sort.size - interv_size + 1

    med_dec = np.zeros(interv_num)
    med_poserr = np.zeros(interv_num)

    for i in range(interv_num):
        ind_b, ind_e = i, i + interv_size
        deci = dec_sort[ind_b: ind_e]
        pos_erri = poserr_sort[ind_b: ind_e]

        med_dec[i] = np.median(deci)
        med_poserr[i] = np.median(pos_erri)

    return med_dec, med_poserr


def smooth_by_time(t, x, tb=None, te=None, ts=0.1):
    """Smooth the series by time coordinate.

    First bin the series by the time coordinates with a bin size of 0.1,
    then calculate the median of each bin.

    Parameters
    ----------
    t : array
        time coordinate of series
    x : array
        value of series
    tb : float
        beginning time of the first bin
    te : float
        end time of the last bin
    ts : float
        bin size. The default value is 0.1

    Returns
    -------
    t_med : array
        median time coordinate of the bin
    med : array
        smoothed median of original series
    """

    if tb is None:
        tb = np.min(t)  # Beginning
    if te is None:
        te = np.max(t)  # End

    bin_num = int((te - tb) // ts) + 1
    med = np.zeros(bin_num)
    t_med = np.zeros(bin_num)

    for i in range(bin_num-1):
        t1i = tb + i * ts
        t2i = t1i + ts

        mask = (t1i <= t) & (t < t2i)
        med[i] = np.median(x[mask])
        t_med[i] = np.median(t[mask])

    # Last element
    t1e = tb + (bin_num-1) * ts
    t2e = t1e + ts
    mask = (t1e <= t) & (t <= t2e)
    med[-1] = np.median(x[mask])
    t_med[-1] = np.median(t[mask])

    return t_med, med


# --------------------------------- END --------------------------------


if __name__ == "__main__":
    dec_smoothed_err(dec, sig_pos_max, interv_size=50)
