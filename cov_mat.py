#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# File name: cov_mat.py
"""
Created on Mon Oct 29 11:12:29 2018

@author: Neo(liuniu@smail.nju.edu.cn)

Calculate the covariance matrix and weighted matrix
"""

import numpy as np
from numpy import concatenate
import sys


__all__ = ["calc_cov_mat", "calc_wgt_mat", "read_cov_mat"]


# -----------------------------  FUNCTIONS -----------------------------
def calc_cov_mat(x_err, y_err, cov=None, corr=None):
    '''Generate the covriance matrix.

    Parameters
    ----------
    x_err/y_err : array of float
        formal uncertainty of two components, usually dRA(*cos(Decl.))/dDecl.
    cov : array of float
        covariance between two components, default is None
    corr : array of float
        correlation coefficient between two components, default is None

    Returns
    ----------
    cov_mat : matrix
        covriance matrix.
    '''

    # check if cov and corr are both given used
    if corr is not None and cov is not None:
        print("corr and cov couldn't be used simultaneously!")
        sys.exit()

    err = concatenate((x_err, y_err), axis=0)

    # Covariance matrix.
    cov_mat = np.diag(err**2)

    # Take the correlation into consideration.
    if corr is not None:
        cov = corr * x_err * y_err

    if cov is not None:
        num = x_err.size
        for i, covi in enumerate(cov):
            cov_mat[i, i + num] = covi
            cov_mat[i + num, i] = covi

    return cov_mat


def calc_wgt_mat(x_err, y_err, cov=None, corr=None):
    '''Generate the weighted matrix.

    Parameters
    ----------
    x_err/y_err : array of float
        formal uncertainty of two components, usually dRA(*cos(DE))/dDE
    cov : array of float
        covariance between two components, default is None
    corr : array of float
        correlation coefficient between two components, default is None

    Returns
    ----------
    wgt_mat : matrix
        weighted matrix used in the least squares fitting.
    '''

    # check if cov and corr are both given used
    if corr is not None and cov is not None:
        print("corr and cov couldn't be used simultaneously!")
        sys.exit()

    # Calculate the covariance matrix
    cov_mat = calc_cov_mat(x_err, y_err, cov, corr)

    # Inverse the covariance matrix to obtain the weighted matrix.
    wgt_mat = np.linalg.inv(cov_mat)

    return wgt_mat


def read_cov_mat(cov_mat):
    """Obtain formal uncertainties and correlations from the covariance matrix.

    Parameter
    ---------
    cov_mat : matrix of (N,N)
            covariance matrix

    Returns
    -------
    sig : 1-d array
        formal uncertainties
    corr_mat : matrix of N * N
        matrix of correlation coefficients
    """

    sig = np.sqrt(cov.diagonal())
    num = sig.size()

    # Correlation coefficient.
    corr_mat = np.array([cov[i, j] / sig[i] / sig[j]
                         for j in range(num) for i in range(num)])

    corr_mat.resize((num, num))

    return sig, corr_mat


def main():
    pass


if __name__ == '__main__':
    main()
# --------------------------------- END --------------------------------
