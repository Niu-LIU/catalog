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

__all__ = [""]


# There are some codes written in the past.
# I will check them later.
# ----------- Horizontal Table -------------------
def htable(names, estimates, errors, fout):
    '''save estmates and corresponding formal errors into a horizontal table.
    '''
    nameline = '$%10s$' % names[0]
    # dataline = '%+8.3f $\\pm$ %8.3f' % (estimates[0], errors[0])
    dataline = '$%+8.1f \\pm %8.1f$' % (estimates[0], errors[0])
    for i in range(1, len(names)):
        nameline += '  &$%10s$' % names[i]
        # dataline += '  &%+8.3f $\\pm$ %8.3f' % (estimates[i], errors[i])
        dataline += '  &$%+8.1f \\pm %8.1f$' % (estimates[i], errors[i])
    nameline += ' \\\\'
    dataline += ' \\\\'
    print(nameline + '\n \\hline \n' + dataline, file=fout)


def htable_int(names, estimates, errors, fout):
    '''save estmates and corresponding formal errors into a horizontal table.
    '''
    nameline = '$%10s$' % names[0]
    # dataline = '%+8.3f $\\pm$ %8.3f' % (estimates[0], errors[0])
    # dataline = '$%+8.1f \\pm %8.1f$' % (estimates[0], errors[0])
    dataline = '$%+7.0f \\pm %7.0f$' % (estimates[0], errors[0])
    for i in range(1, len(names)):
        nameline += '  &$%10s$' % names[i]
        dataline += '  &$%+7.0f \\pm %7.0f$' % (estimates[i], errors[i])
        # dataline += '  &$%+8.1f \\pm %8.1f$' % (estimates[i], errors[i])
    nameline += ' \\\\'
    dataline += ' \\\\'
    print(nameline + '\n \\hline \n' + dataline, file=fout)


# ----------- Vertical Table -------------------
def vtable(names, estimates, errors, fout):
    '''save estmates and corresponding formal errors into a vertical table.
    '''
    for i in range(len(names)):
        # print('$%10s$  &%+8.3f $\\pm$ %8.3f \\\\'\
        print('$%10s$  &$%+8.1f \\pm %8.1f$ \\\\'
              % (names[i], estimates[i], errors[i]), file=fout)


# ----------- Vertical Table -------------------
def vtable_int(names, estimates, errors, fout):
    '''save estmates and corresponding formal errors into a vertical table.
    '''
    for i in range(len(names)):
        # print('$%10s$  &%+8.3f $\\pm$ %8.3f \\\\'\
        # print('$%10s$  &$%+8.1f \\pm %8.1f$ \\\\'
        print('$%10s$  &$%+7.0f \\pm %7.0f$ \\\\'
              % (names[i], estimates[i], errors[i]), file=fout)


# ----------- Correlation coefficients -------------------
def cor_table(names, cor, fout):
    print('## correlation coefficients.', file=fout)
    num = len(names)
# Heading
    headline = ' ' * 10
    for name in names:
        headline += '  &$%10s$' % name
    headline += '\\\\ \\hline'
    print(headline, file=fout)
# Coefficient
    for i, name in enumerate(names):
        headpart = '  $%10s$' % name
        blankpart = '  &      ' * i
        datapart = ''
        for j in range(i, num):
            datapart += '  &$%+5.2f$' % cor[i, j]
        print(headpart + blankpart + datapart + '\\\\', file=fout)


def write_result_deg1(x1name, x1, sig1, corr1, flog):
    htable(x1name, x1, sig1, flog)
    cor_table(x1name, corr1, flog)


def write_result_deg2(x1name, x2name, x2, sig2, corr2, flog):
    x21, sig21 = x2[: 6], sig2[: 6]
    x22, sig22 = x2[6:], sig2[6:]
    # htable(x1name, x21, sig21, flog)
    htable(x2name, x22, sig22, flog)
    # htable(x1name, x21, sig21, flog)
    htable_int(x2name, x22, sig22, flog)
    cor_table(x1name + x2name, corr2, flog)


def write_result_deg2_int(x1name, x2name, x2, sig2, corr2, flog):
    x21, sig21 = x2[: 6], sig2[: 6]
    x22, sig22 = x2[6:], sig2[6:]
    htable_int(x1name, x21, sig21, flog)
    htable_int(x2name, x22, sig22, flog)
    cor_table(x1name + x2name, corr2, flog)


#######################################################################
# New codes

def print_corr(pmt_name, pmt_corr, deci_digit=2, included_one=True):
    """Print the correlation coefficient in the screen.

    Parameters
    ----------
    pmt_name : array
        names or labels of the parameters
    pmt_corr : array of N * N
        matrix or correlation coefficient
    deci_digit : int
        decimal digits. Only 1 and 2 are supported. Default is 2.
    included_one : boolean
        to include the correlation between the same parameters (of course the
        value equals to 1) or not. True for yes.
    """

    pmt_nb = pmt_name.size

    # If the correlation matrix is an object of np.mat type,
    # convert it to np.ndarray
    if type(pmt_corr) is np.mat or type(pmt_corr) is np.matrix:
        pmt_corr = np.array(pmt_corr)

    if deci_digit == 1:
        # Now begin to print the correlation coefficient
        if included_one:
            # The first line
            print(("  %4s" * (pmt_nb+1)) % ("    ", *pmt_name))

            # Include the correlation coefficient of one
            for i, pmt_namei in enumerate(pmt_name):
                line_fmt = "  %4s" + "  %+4.1f" * (i + 1)
                print(line_fmt % (pmt_namei, *pmt_corr[i, :i+1]))
        else:
            # The first line
            print(("  %4s" * pmt_nb) % ("    ", *pmt_name[:-1]))

            for i in range(1, pmt_nb):
                line_fmt = "  %4s" + "  %+4.1f" * i
                print(line_fmt % (pmt_name[i], *pmt_corr[i, :i]))

    elif deci_digit == 2:

        # Now begin to print the correlation coefficient
        if included_one:
            # The first line
            print(("  %5s" * (pmt_nb+1)) % ("    ", *pmt_name))

            for i, pmt_namei in enumerate(pmt_name):
                line_fmt = "  %5s" + "  %+5.2f" * (i + 1)
                print(line_fmt % (pmt_namei, *pmt_corr[i, :i+1]))
        else:
            # The first line
            print(("  %5s" * pmt_nb) % ("    ", *pmt_name[:-1]))

            for i in range(1, pmt_nb):
                line_fmt = "  %5s" + "  %+5.2f" * i
                print(line_fmt % (pmt_name[i], *pmt_corr[i, :i]))
    else:
        print("The decimal digit of only 1 or 2 is supported!")
        sys.exit()


def print_vsh1_corr(pmt_corr, deci_digit=2, included_one=True):
    """Print the correlation coefficient of VSH01 parameters in the screen.

    Parameters
    ----------
    pmt_corr : array of N * N
        matrix or correlation coefficient
    deci_digit : int
        decimal digits. Only 1 and 2 are supported. Default is 2.
    included_one : boolean
        to include the correlation between the same parameters (of course the
        value equals to 1) or not. True for yes.
    """

    pmt_name = np.array(["R1", "R2", "R3", "D1", "D2", "D3"])

    pmt_nb = pmt_name.size

    # Check the shape of the matrix
    a, b = pmt_corr.shape
    if a != b or a != pmt_nb:
        print("The shape of the correlation matrix should be (N, N)(N=6)!")
        sys.exit()

    print_corr(pmt_name, pmt_corr, deci_digit, included_one)


def print_vsh2_corr(pmt_corr, deci_digit=2, included_one=True):
    """Print the correlation coefficient of VSH02 parameters in the screen.

    Parameters
    ----------
    pmt_corr : array of N * N
        matrix or correlation coefficient
    deci_digit : int
        decimal digits. Only 1 and 2 are supported. Default is 2.
    included_one : boolean
        to include the correlation between the same parameters (of course the
        value equals to 1) or not. True for yes.
    """

    pmt_name = np.array(["R1", "R2", "R3", "D1", "D2", "D3",
                         "E22R", "E22I", "E21R", "E21I", "E20",
                         "M22R", "M22I", "M21R", "M21I", "M20"])

    pmt_nb = pmt_name.size

    # Check the shape of the matrix
    a, b = pmt_corr.shape
    if a != b or a != pmt_nb:
        print("The shape of the correlation matrix should be (N, N)(N=16)!")
        sys.exit()

    print_corr(pmt_name, pmt_corr, deci_digit, included_one)


# ------------------------- MAIN ---------------------------------------
if __name__ == "__main__":
    print("Sorry but there is nothing to do.")
# ------------------------- END ----------------------------------------
