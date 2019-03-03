#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: pos_diff.py
"""
Created on Fri Sep 21 15:39:02 2018

@author: Neo(liuniu@smail.nju.edu.cn)
"""

import numpy as np
from numpy import cos, sqrt
from functools import reduce


__all__ = ["nor_sep", "pos_diff_calc", "pa_calc"]


# -----------------------------  FUNCTIONS -----------------------------
def nor_sep(dRA, dRA_err, dDC, dDC_err, C):
    '''Calculate the normalized seperation.

    Parameters
    ----------
    dRA/dDC : Right Ascension / Declination differences in micro-as
    dRA_err/dDC_err : formal uncertainty of dRA*cos(Dec)/dDC in micro-as
    C : correlation coeffient between dRA*cos(Dec) and dDC.

    Returns
    ----------
    ang_sep : angular seperation, in micro-as
    X_a / X_d : normalized coordinate differences in RA / DC, unit-less
    X : Normalized separations, unit-less.
    '''

    # Angular seperations
    ang_sep = sqrt(dRA**2 + dDC**2)

    # Normalised coordinate difference
    X_a = dRA / dRA_err
    X_d = dDC / dDC_err

    # Normalised separation - Mignard's statistics (considering covariance)
    X = np.zeros_like(X_a)

    for i, (X_ai, X_di, Ci) in enumerate(zip(X_a, X_d, C)):
        if Ci == -1.:
            Ci = -0.99
        if Ci == 1.0:
            Ci = 0.99
        # print(Ci)

        wgt = np.linalg.inv(np.mat([[1, Ci],
                                    [Ci, 1]]))
        Xmat = np.mat([X_ai, X_di])
        X[i] = sqrt(reduce(np.dot, (Xmat, wgt, Xmat.T)))

    # # Normalised separation (normal way)
    X2 = sqrt(X_a**2 + X_d**2)

    # return ang_sep, X_a, X_d, X
    return ang_sep, X_a, X_d, X, X2


def pos_diff_calc(RA1, RA1_err, DC1, DC1_err, Cor1,
                  RA2, RA2_err, DC2, DC2_err, Cor2,
                  arccof=None):
    '''Calculate the normalized seperation between VLBI and Gaia positions.


    Parameters
    ----------
    RA / DC : Right Ascension / Declination, degress
    e_RA / e_DC : formal uncertainty of RA * cos(Dec) / DC, mas
    Cor : correlation coeffient between RA and DC.
    arccof : cos(Dec.)

    Note: suffix 'G' stands for GaiaDR1 and I for VLBI catalog.

    Returns
    ----------
    ang_sep : angular seperation in micro-as
    X_a / X_d : normalized seperation in RA / DC, unit-less
    X : Normalized separations, unit-less.
    '''

    if arccof is None:
        arccof = cos(np.deg2rad(DC1))

    # # deg -> uas
    # dRA = (RA1 - RA2) * 3.6e9 * arccof
    # dRA_err = sqrt(RA1_err**2 + RA2_err**2)
    # dDC = (DC1 - DC2) * 3.6e9
    # dDC_err = sqrt(DC1_err**2 + DC2_err**2)

    # deg -> mas
    dRA = (RA1 - RA2) * 3.6e6 * arccof
    dRA_err = sqrt(RA1_err**2 + RA2_err**2)
    dDC = (DC1 - DC2) * 3.6e6
    dDC_err = sqrt(DC1_err**2 + DC2_err**2)

    # Correlation coefficient of combined errors
    cov = RA1_err * DC1_err * Cor1 + RA2_err * DC2_err * Cor2
    corf = cov / (dRA_err * dDC_err)

    # Normalised separation
    ang_sep, X_a, X_d, X, X2 = nor_sep(dRA, dRA_err,
                                       dDC, dDC_err, corf)

    # return ang_sep, X_a, X_d, X
    return dRA, dDC, dRA_err, dDC_err, cov, ang_sep, X_a, X_d, X, X2


def pa_calc(dra, ddec):
    """Calcaulte the position angle (East of North) of positional offsets.

    Parameters
    ----------
    dra/ddec : array
        positional offset in RA/decl.

    Return
    ------
    pa : array
        position angles in degree.
    """

    pa = np.zeros_like(dra)

    for i, (drai, ddeci) in enumerate(zip(dra, ddec)):
        ang_sep = np.sqrt(drai**2 + ddeci**2)
        pai = np.rad2deg(np.arccos(ddeci/ang_sep))

        if drai < 0:
            pai = 360 - pai

        pa[i] = pai

    return pa


def main():
    """Main function.
    """

    print("This code is just a module of functions.")


if __name__ == "__main__":
    main()
# --------------------------------- END --------------------------------
