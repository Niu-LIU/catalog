# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 14:12:08 2017

Write the output information (estimates associated with their formal errors
and covariance) into a text file or on screen.

@author: Neo
"""


import numpy as np
import sys

__all__ = ["htable", "vtable", "print_vsh1_corr", "print_vsh2_corr",
           "print_vsh1_result", "print_vsh2_result"]


# There are some codes written in the past.
# I will check them later.
# ----------- Horizontal Table -------------------
def htable(names, pmts, sigs, opt=sys.stdout, fmt="%8.2f"):
    '''Print estmates and corresponding formal errors into a horizontal table.

    Parameters
    ----------
    names : array of string
        names or labels of the parameters
    pmts : array of float
        estimates of parameters
    sigs : array of float
        formal errors of estimates
    opt : file handling
        Default to print all the output on the screen.
    fmt : string
        specifier of the output format
    '''

    pmt_nb = pmts.size
    tmp = np.vstack((pmts, sigs)).flatten("F")

    # Table header
    head_line = ["  $%10s$" % names[0]]

    data_fmt1 = "  ${0:s} \\pm {0:s}$".format(fmt)
    data_fmt2 = "  &${0:s} \\pm {0:s}$".format(fmt)

    # Estimate with formal errors
    data_line = [data_fmt1 % (pmts[0], sigs[0])]

    for i in range(1, pmt_nb-1):
        head_line.append("  &$%10s$" % names[i])
        data_line.append(data_fmt2 % (pmts[i], sigs[i]))

    head_line.append(" \\\\")
    data_line.append(" \\\\")

    print("".join(head_line), file=opt)
    print("\\hline", file=opt)
    print("".join(data_line), file=opt)


# ----------- Vertical Table -------------------
def vtable(names, pmts, sigs, opt=sys.stdout, fmt="%8.2f"):
    '''Print estmates and corresponding formal errors into a horizontal table.

    Parameters
    ----------
    names : array of string
        names or labels of the parameters
    pmts : array of float
        estimates of parameters
    sigs : array of float
        formal errors of estimates
    opt : file handling
        Default to print all the output on the screen.
    fmt : string
        specifier of the output format
    '''

    line_fmt = "$%10s$  &${0:s} \\pm {0:s}$ \\\\".format(fmt)

    for pmt_namei, pmti, pmt_erri in (names, pmts, sigs):
        print(line_fmt % (pmt_namei, pmti, pmt_erri), file=opt)


#######################################################################
# New codes
# ----------- Correlation coefficients -------------------
def print_corr(names, corrs, deci_digit=2, included_one=True,
               opt=sys.stdout):
    """Print the correlation coefficient in the screen.

    Parameters
    ----------
    names : array
        names or labels of the parameters
    corrs : array of N * N
        matrix or correlation coefficient
    deci_digit : int
        decimal digits. Only 1 and 2 are supported. Default is 2.
    included_one : boolean
        to include the correlation between the same parameters (of course the
        value equals to 1) or not. True for yes.
    opt : file handling
        Default to print all the output on the screen.
    """

    pmt_nb = names.size

    # If the correlation matrix is an object of np.mat type,
    # convert it to np.ndarray
    if type(corrs) is np.mat or type(corrs) is np.matrix:
        corrs = np.array(corrs)

    if deci_digit == 1:
        # Now begin to print the correlation coefficient
        if included_one:
            # The first line
            print(("  %4s" * (pmt_nb+1)) % ("    ", *names), file=opt)

            # Include the correlation coefficient of one
            for i, pmt_namei in enumerate(names):
                line_fmt = "  %4s" + "  %+4.1f" * (i + 1)
                print(line_fmt % (pmt_namei, *corrs[i, :i+1]), file=opt)
        else:
            # The first line
            print(("  %4s" * pmt_nb) % ("    ", *names[:-1]), file=opt)

            for i in range(1, pmt_nb):
                line_fmt = "  %4s" + "  %+4.1f" * i
                print(line_fmt % (names[i], *corrs[i, :i]), file=opt)

    elif deci_digit == 2:

        # Now begin to print the correlation coefficient
        if included_one:
            # The first line
            print(("  %5s" * (pmt_nb+1)) % ("    ", *names), file=opt)

            for i, pmt_namei in enumerate(names):
                line_fmt = "  %5s" + "  %+5.2f" * (i + 1)
                print(line_fmt % (pmt_namei, *corrs[i, :i+1]), file=opt)
        else:
            # The first line
            print(("  %5s" * pmt_nb) % ("    ", *names[:-1]), file=opt)

            for i in range(1, pmt_nb):
                line_fmt = "  %5s" + "  %+5.2f" * i
                print(line_fmt % (names[i], *corrs[i, :i]), file=opt)
    else:
        print("The decimal digit of only 1 or 2 is supported!")
        sys.exit()


def print_vsh1_corr(corrs, names=None,
                    deci_digit=2, included_one=True, opt=sys.stdout):
    """Print the correlation coefficient of VSH01 parameters in the screen.

    Parameters
    ----------
    corrs : array of N * N
        matrix or correlation coefficient
    deci_digit : int
        decimal digits. Only 1 and 2 are supported. Default is 2.
    included_one : boolean
        to include the correlation between the same parameters (of course the
        value equals to 1) or not. True for yes.
    """

    if names is None:
        names = np.array(["R1", "R2", "R3", "D1", "D2", "D3"])

    pmt_nb = names.size

    # Check the shape of the matrix
    a, b = corrs.shape
    if a != b or a != pmt_nb:
        print("The shape of the correlation matrix should be (N, N)(N=6)!")
        sys.exit()

    print_corr(names, corrs, deci_digit, included_one, opt)


def print_vsh2_corr(corrs, names=None,
                    deci_digit=2, included_one=True, opt=sys.stdout):
    """Print the correlation coefficient of VSH02 parameters in the screen.

    Parameters
    ----------
    corrs : array of N * N
        matrix or correlation coefficient
    deci_digit : int
        decimal digits. Only 1 and 2 are supported. Default is 2.
    included_one : boolean
        to include the correlation between the same parameters (of course the
        value equals to 1) or not. True for yes.
    """

    if names is None:
        names = np.array(["R1", "R2", "R3", "D1", "D2", "D3",
                          "E22R", "E22I", "E21R", "E21I", "E20",
                          "M22R", "M22I", "M21R", "M21I", "M20"])

    pmt_nb = names.size

    # Check the shape of the matrix
    a, b = corrs.shape
    if a != b or a != pmt_nb:
        print("The shape of the correlation matrix should be (N, N)(N=16)!")
        sys.exit()

    print_corr(names, corrs, deci_digit, included_one, opt)


def print_vsh1_result(pmts, sigs, corrs,
                      names=None, opt=sys.stdout, fmt="%8.2f"):
    '''Print estmates and corresponding formal errors of vsh01 parameters.

    Parameters
    ----------
    names : array of string
        names or labels of the parameters
    pmts : array of float
        estimates of parameters
    sigs : array of float
        formal errors of estimates
    opt : file handling
        Default to print all the output on the screen.
    fmt : string
        specifier of the output format
    '''

    if names is None:
        names = np.array(["R1", "R2", "R3", "D1", "D2", "D3"])

    htable(names, pmts, sigs, opt, fmt)
    print_vsh1_corr(corrs, names,
                    deci_digit=2, included_one=True, opt=opt)


def print_vsh2_result(pmts, sigs, corrs,
                      names=None, opt=sys.stdout, fmt="%8.2f"):
    '''Print estmates and corresponding formal errors of vsh02 parameters.

    Parameters
    ----------
    names : array of string
        names or labels of the parameters
    pmts : array of float
        estimates of parameters
    sigs : array of float
        formal errors of estimates
    opt : file handle
        Default to print all the output on the screen.
    fmt : string
        specifier of the output format
    '''

    if names is None:
        names = np.array(["R1", "R2", "R3", "D1", "D2", "D3",
                          "E22R", "E22I", "E21R", "E21I", "E20",
                          "M22R", "M22I", "M21R", "M21I", "M20"])

    htable(names, pmts, sigs, opt, fmt)
    print_vsh2_corr(corrs, names,
                    deci_digit=2, included_one=True, opt=opt)


def save_vsh1_result(pmts, sigs, ofile, fmt="%5.0f"):
    """Save the VSH results into a text file.

    Parameters
    ----------
    pmts : array of float
        estimates of parameters
    sigs : array of float
        formal errors of estimates
    opt : file handle
        Default to print all the output on the screen.
    fmt : string
        specifier of the output format
    """

    names = Column(["D1", "D2", "D3",
                    "R1", "R2", "R3"])
    tvsh1 = Table([names, pmts, sig2], names=[
        "Names", "Estimate", "Error"])
    tvsh1["Estimate"].format = fmt
    tvsh1["Error"].format = fmt

    tvsh1.write(ofile, format="ascii", overwrite=True)


def save_vsh2_result(pmts, sigs, ofile, fmt="%5.0f"):
    """Save the VSH results into a text file.

    Parameters
    ----------
    pmts : array of float
        estimates of parameters
    sigs : array of float
        formal errors of estimates
    opt : file handle
        Default to print all the output on the screen.
    fmt : string
        specifier of the output format
    """

    names = Column(["D1", "D2", "D3",
                    "R1", "R2", "R3",
                    "ER22", "EI22", "ER21", "EI21", "E20",
                    "MR22", "MI22", "MR21", "MI21", "M20"])
    tvsh2 = Table([names, pmts, sig2], names=[
        "Names", "Estimate", "Error"])
    tvsh2["Estimate"].format = fmt
    tvsh2["Error"].format = fmt

    # tvsh2["Estimate"].unit = u.uas
    # tvsh2["Error"].unit = u.uas

    tvsh2.write(ofile, format="ascii", overwrite=True)


# ------------------------- MAIN ---------------------------------------
if __name__ == "__main__":
    print("Sorry but there is nothing to do.")
# ------------------------- END ----------------------------------------
