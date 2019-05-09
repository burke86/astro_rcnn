#!/usr/bin/env python
##
# @package cwfs
# @file cwfsTools.py
##
# @authors: Bo Xin & Chuck Claver
# @       Large Synoptic Survey Telescope

# getCenterAndR() is partly based on the EF wavefront sensing software
# by Laplacian Optics

##

import sys
import numpy as np


def padArray(inArray, dim):
    m, n = inArray.shape
    if (m != n):
        raise('padArray: array is not square')

    if m > dim:
        raise('padArray: array is larger than dimension')

    out = np.zeros((dim, dim))
    i = np.floor((dim - m) / 2)
    j = i + m
    out[i:j, i:j] = inArray

    return out


def ZernikeAnnularGrad(Z, x, y, e, type):
    '''Gradient of the Annular Zernicke'''
    m1, n1 = x.shape
    m2, n2 = y.shape

    if(m1 != m2 or n1 != n2):
        print('x & y are not the same size')
        exit()

    if(len(Z) > 22):
        print('ZernikeAnnularEval() is not implemented with >22 terms')
        return
    elif len(Z) < 22:
        Z[21] = 0

    x2 = x * x
    y2 = y * y
    x4 = x2 * x2
    y4 = y2 * y2
    xy = x * y
    r2 = x2 + y2
    r4 = r2 * r2
    e2 = e * e
    e4 = e2 * e2
    e6 = e4 * e2
    e8 = e6 * e2
    e10 = e8 * e2
    e12 = e10 * e2

    if (type == 'dx'):
        d = Z[0] * 0 * x  # to make d an array with the same size as x
        den = np.sqrt(1 + e2)
        d = d + Z[1] * 2 * 1 / den
        d = d + Z[2] * 2 * 0
        den = 1 - e**2
        d = d + Z[3] * np.sqrt(3) * 4 * x / den
        den = np.sqrt(1 + e2 + e4)
        d = d + Z[4] * np.sqrt(6) * 2 * y / den
        d = d + Z[5] * np.sqrt(6) * 2 * x / den
        den = np.sqrt((1 - e2)**2 * (1 + e2) * (1 + 4 * e2 + e4))
        d = d + Z[6] * np.sqrt(8) * 6 * xy * (1 + e2) / den
        d = d + Z[7] * np.sqrt(8) * ((9 * x2 + 3 * y2 - 2) *
                                     (1 + e2) - 2 * e4) / den
        den = np.sqrt(1 + e2 + e4 + e6)
        d = d + Z[8] * np.sqrt(8) * 6 * xy / den
        d = d + Z[9] * np.sqrt(8) * (3 * x2 - 3 * y2) / den
        den = (1 - e2)**2
        d = d + Z[10] * np.sqrt(5) * 12 * x * (2 * r2 - 1 - e2) / den
        den = (1 - e2)**3 * (1 + e2 + e4)
        num = np.sqrt((1 - e2)**4 * (1 + e2 + e4) /
                      (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8))
        d = d + Z[11] * np.sqrt(10) * (x * (16 * x2 - 6) *
                                       (1 + e2 + e4) - 6 * x * e6) * num / den
        d = d + Z[12] * np.sqrt(10) * (y * (24 * x2 + 8 * y2 - 6) *
                                       (1 + e2 + e4) - 6 * y * e6) * num / den
        den = np.sqrt(1 + e2 + e4 + e6 + e8)
        d = d + Z[13] * np.sqrt(10) * 4 * x * (x2 - 3 * y2) / den
        d = d + Z[14] * np.sqrt(10) * 4 * y * (3 * x2 - y2) / den
        den = (1 - e2)**3 * (1 + 4 * e2 + e4)
        num = np.sqrt((1 - e2)**2 * (1 + 4 * e2 + e4) /
                      (1 + 9 * e2 + 9 * e4 + e6))
        d = d + Z[15] * np.sqrt(12) * (
            3 * e8 - 36 * e6 * x2 - 12 * e6 * y2 + 12 * e6 +
            50 * e4 * x4 + 60 * e4 * x2 * y2 - 144 * e4 * x2 +
            10 * e4 * y4 - 48 * e4 * y2 + 30 * e4 + 200 * e2 * x4 + 240 *
            e2 * x2 * y2 - 144 * e2 * x2 + 40 * e2 * y4 - 48 * e2 * y2 +
            12 * e2 + 50 * x4 + 60 * x2 * y2 - 36 * x2 +
            10 * y4 - 12 * y2 + 3) * num / den
        d = d + Z[16] * np.sqrt(12) * (
            8 * xy * (5 * r2 * (1 + 4 * e2 + e4) -
                      (3 + 12 * e2 + 12 * e4 + 3 * e6))) * num / den
        den = (1 - e2)**4 * (1 + e2) * (1 + e4)
        num = np.sqrt((1 - e2)**6 * (1 + e2) * (1 + e4) /
                      (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 +
                       4 * e10 + e12))
        d = d + Z[17] * np.sqrt(12) * (
            25 * (e6 + e4 + e2 + 1) * x4 +
            (- 12 * e8 - 30 * e6 * y2 - 12 * e6 - 30 * e4 * y2 - 12 * e4 -
             30 * e2 * y2 - 12 * e2 - 30 * y2 - 12) * x2 + 12 * e8 * y2 -
            15 * e6 * y4 + 12 * e6 * y2 - 15 * e4 * y4 + 12 * e4 * y2 -
            15 * e2 * y4 + 12 * e2 * y2 - 15 * y4 + 12 * y2) * num / den
        d = d + Z[18] * np.sqrt(12) * (
            4.0 * xy * (15 * (e6 + e4 + e2 + 1) * x2 - 6 * e8 + 5 * e6 * y2 -
                        6 * e6 + 5 * e4 * y2 - 6 * e4 + 5 * e2 * y2 -
                        6 * e2 + 5 * y2 - 6)) * num / den
        den = np.sqrt(1 + e2 + e4 + e6 + e8 + e10)
        d = d + Z[19] * np.sqrt(12) * 5 * (x2 * (x2 - 6 * y2) + y4) / den
        d = d + Z[20] * np.sqrt(12) * 20 * xy * (x2 - y2) / den
        den = (1 - e2)**3
        d = d + Z[21] * np.sqrt(7) * 24 * x * (
            e4 - e2 * (5 * y2 - 3) + 5 * x4 - 5 * y2 + 5 * y4 -
            x2 * (5 * e2 - 10 * y2 + 5) + 1) / den
    elif (type == 'dy'):
        d = Z[0] * 0 * x
        den = np.sqrt(1 + e2)
        d = d + Z[1] * 2 * 0
        d = d + Z[2] * 2 * 1 / den
        den = 1 - e**2
        d = d + Z[3] * np.sqrt(3) * 4 * y / den
        den = np.sqrt(1 + e2 + e4)
        d = d + Z[4] * np.sqrt(6) * 2 * x / den
        d = d + Z[5] * np.sqrt(6) * (-2) * y / den
        den = np.sqrt((1 - e2)**2 * (1 + e2) * (1 + 4 * e2 + e4))
        d = d + Z[6] * np.sqrt(8) * ((1 + e2) *
                                     (3 * x2 + 9 * y2 - 2) - 2 * e4) / den
        d = d + Z[7] * np.sqrt(8) * 6 * xy * (1 + e2) / den
        den = np.sqrt(1 + e2 + e4 + e6)
        d = d + Z[8] * np.sqrt(8) * (3 * x2 - 3 * y2) / den
        d = d + Z[9] * np.sqrt(8) * (-6) * xy / den
        den = (1 - e2)**2
        d = d + Z[10] * np.sqrt(5) * 12 * y * (2 * r2 - 1 - e2) / den
        den = (1 - e2)**3 * (1 + e2 + e4)
        num = np.sqrt((1 - e2)**4 * (1 + e2 + e4) /
                      (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8))
        d = d + Z[11] * np.sqrt(10) * (y * (6 - 16 * y2) *
                                       (1 + e2 + e4) + 6 * y * e6) * num / den
        d = d + Z[12] * np.sqrt(10) * (x * (8 * x2 + 24 * y2 - 6) *
                                       (1 + e2 + e4) - 6 * x * e6) * num / den
        den = np.sqrt(1 + e2 + e4 + e6 + e8)
        d = d + Z[13] * np.sqrt(10) * 4 * y * (y2 - 3 * x2) / den
        d = d + Z[14] * np.sqrt(10) * 4 * x * (x2 - 3 * y2) / den
        den = (1 - e2)**3 * (1 + 4 * e2 + e4)
        num = np.sqrt((1 - e2)**2 * (1 + 4 * e2 + e4) /
                      (1 + 9 * e2 + 9 * e4 + e6))
        d = d + Z[15] * np.sqrt(12) * (
            -x * (24 * y + 4 * e2 * (24 * y - 40 * y * r2) +
                  2 * e4 * (48 * y - 20 * y * r2) + 24 * e6 * y -
                  40 * y * r2)) * num / den
        d = d + Z[16] * np.sqrt(12) * (
            3 * e8 - 12 * e6 * x2 - 36 * e6 * y2 + 12 * e6 + 10 * e4 * x4 +
            60 * e4 * x2 * y2 - 48 * e4 * x2 +
            50 * e4 * y4 - 144 * e4 * y2 + 30 * e4 + 40 * e2 * x4 + 240 *
            e2 * x2 * y2 - 48 * e2 * x2 + 200 * e2 * y4 - 144 * e2 * y2 +
            12 * e2 + 10 * x4 + 60 * x2 * y2 - 12 * x2 +
            50 * y4 - 36 * y2 + 3) * num / den
        den = (1 - e2)**4 * (1 + e2) * (1 + e4)
        num = np.sqrt((1 - e2)**6 * (1 + e2) * (1 + e4) /
                      (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 +
                       4 * e10 + e12))
        d = d + Z[17] * np.sqrt(12) * (
            4.0 * xy * ((- 5) * (e6 + e4 + e2 + 1) * x2 + 6 * e8 -
                        15 * e6 * y2 + 6 * e6 - 15 * e4 * y2 +
                        6 * e4 - 15 * e2 * y2 + 6 * e2 -
                        15 * y2 + 6)) * num / den
        d = d + Z[18] * np.sqrt(12) * (
            - 12 * e8 * x2 + 12 * e8 * y2 + 15 * e6 * x4 +
            30 * e6 * x2 * y2 - 12 * e6 * x2 - 25 * e6 * y4 +
            12 * e6 * y2 + 15 * e4 * x4 + 30 * e4 * x2 * y2 - 12 * e4 * x2 -
            25 * e4 * y4 + 12 * e4 * y2 + 15 * e2 * x4 + 30 * e2 * x2 * y2 -
            12 * e2 * x2 - 25 * e2 * y4 + 12 * e2 * y2 + 15 * x4 +
            30 * x2 * y2 - 12 * x2 - 25 * y4 + 12 * y2) * num / den
        den = np.sqrt(1 + e2 + e4 + e6 + e8 + e10)
        d = d + Z[19] * np.sqrt(12) * 20 * xy * (y2 - x2) / den
        d = d + Z[20] * np.sqrt(12) * 5 * (x2 * (x2 - 6 * y2) + y4) / den
        den = (1 - e2)**3
        d = d + Z[21] * np.sqrt(7) * 24 * y * (
            e4 - e2 * (5 * x2 - 3) - 5 * x2 + 5 * x4 + 5 * y4 -
            y2 * (5 * e2 - 10 * x2 + 5) + 1) / den
    elif (type == 'dx2'):
        d = Z[0] * 0 * x  # to make d an array with the same size as x
        d = d + Z[1] * 0
        d = d + Z[2] * 0
        den = 1 - e**2
        d = d + Z[3] * np.sqrt(3) * 4 / den
        d = d + Z[4] * 0
        den = np.sqrt(1 + e2 + e4)
        d = d + Z[5] * np.sqrt(6) * 2 / den
        den = np.sqrt((1 - e2)**2 * (1 + e2) * (1 + 4 * e2 + e4))
        d = d + Z[6] * np.sqrt(8) * 6 * y * (1 + e2) / den
        d = d + Z[7] * np.sqrt(8) * 18 * x * (1 + e2) / den
        den = np.sqrt(1 + e2 + e4 + e6)
        d = d + Z[8] * np.sqrt(8) * 6 * y / den
        d = d + Z[9] * np.sqrt(8) * 6 * x / den
        den = (1 - e2)**2
        d = d + Z[10] * np.sqrt(5) * 12 * (6 * x2 + 2 * y2 - e2 - 1) / den
        den = (1 - e2)**3 * (1 + e2 + e4)
        num = np.sqrt((1 - e2)**4 * (1 + e2 + e4) /
                      (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8))
        d = d + Z[11] * np.sqrt(10) * ((48 * x2 - 6) *
                                       (1 + e2 + e4) - 6 * e6) * num / den
        d = d + Z[12] * np.sqrt(10) * 48 * xy * (1 + e2 + e4) * num / den
        den = np.sqrt(1 + e2 + e4 + e6 + e8)
        d = d + Z[13] * np.sqrt(10) * 12 * (x2 - y2) / den
        d = d + Z[14] * np.sqrt(10) * 24 * xy / den
        den = (1 - e2)**3 * (1 + 4 * e2 + e4)
        num = np.sqrt((1 - e2)**2 * (1 + 4 * e2 + e4) /
                      (1 + 9 * e2 + 9 * e4 + e6))
        d = d + Z[15] * np.sqrt(12) * (
            -8 * x * (9 * e6 - 25 * e4 * x2 - 15 * e4 * y2 + 36 * e4 -
                      100 * e2 * x2 - 60 * e2 * y2 + 36 * e2 - 25 * x2 -
                      15 * y2 + 9)) * num / den
        d = d + Z[16] * np.sqrt(12) * (
            -8 * y * (3 * e6 - 15 * e4 * x2 - 5 * e4 * y2 + 12 * e4 -
                      60 * e2 * x2 - 20 * e2 * y2 + 12 * e2 - 15 * x2 -
                      5 * y2 + 3)) * num / den
        den = (1 - e2)**4 * (1 + e2) * (1 + e4)
        num = np.sqrt((1 - e2)**6 * (1 + e2) * (1 + e4) /
                      (1 + 4 * e2 + 10 * e4 + 20 * e6 +
                       10 * e8 + 4 * e10 + e12))
        d = d + Z[17] * np.sqrt(12) * (
            -4 * x * (6 * e8 - 25 * e6 * x2 + 15 * e6 * y2 + 6 * e6 -
                      25 * e4 * x2 + 15 * e4 * y2 + 6 * e4 - 25 * e2 * x2 +
                      15 * e2 * y2 + 6 * e2 - 25 * x2 +
                      15 * y2 + 6)) * num / den
        d = d + Z[18] * np.sqrt(12) * (
            -4 * y * (6 * e8 - 45 * e6 * x2 - 5 * e6 * y2 + 6 * e6 -
                      45 * e4 * x2 - 5 * e4 * y2 + 6 * e4 - 45 * e2 * x2 -
                      5 * e2 * y2 + 6 * e2 - 45 * x2 - 5 * y2 + 6)) * num / den
        den = np.sqrt(1 + e2 + e4 + e6 + e8 + e10)
        d = d + Z[19] * np.sqrt(12) * 20 * x * (x2 - 3 * y2) / den
        d = d + Z[20] * np.sqrt(12) * 20 * y * (3 * x2 - y2) / den
        den = (1 - e2)**3
        d = d + Z[21] * np.sqrt(7) * (
            480 * x2 * r2 + 120 * r4 + 24 * e4 - 360 * x2 - 120 * y2 -
            3 * e2 * (120 * x2 + 40 * y2 - 24) + 24) / den

    elif (type == 'dy2'):
        d = Z[0] * 0 * x  # to make d an array with the same size as x
        d = d + Z[1] * 0
        d = d + Z[2] * 0
        den = 1 - e**2
        d = d + Z[3] * np.sqrt(3) * 4 / den
        d = d + Z[4] * 0
        den = np.sqrt(1 + e2 + e4)
        d = d + Z[5] * np.sqrt(6) * (-2) / den
        den = np.sqrt((1 - e2)**2 * (1 + e2) * (1 + 4 * e2 + e4))
        d = d + Z[6] * np.sqrt(8) * (1 + e2) * 18 * y / den
        d = d + Z[7] * np.sqrt(8) * 6 * x * (1 + e2) / den
        den = np.sqrt(1 + e2 + e4 + e6)
        d = d + Z[8] * np.sqrt(8) * (-6) * y / den
        d = d + Z[9] * np.sqrt(8) * (-6) * x / den
        den = (1 - e2)**2
        d = d + Z[10] * np.sqrt(5) * 12 * (2 * x2 + 6 * y2 - e2 - 1) / den
        den = (1 - e2)**3 * (1 + e2 + e4)
        num = np.sqrt((1 - e2)**4 * (1 + e2 + e4) /
                      (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8))
        d = d + Z[11] * np.sqrt(10) * ((6 - 48 * y2) *
                                       (1 + e2 + e4) + 6 * e6) * num / den
        d = d + Z[12] * np.sqrt(10) * 48 * xy * (1 + e2 + e4) * num / den
        den = np.sqrt(1 + e2 + e4 + e6 + e8)
        d = d + Z[13] * np.sqrt(10) * 12 * (y2 - x2) / den
        d = d + Z[14] * np.sqrt(10) * (-24) * xy / den
        den = (1 - e2)**3 * (1 + 4 * e2 + e4)
        num = np.sqrt((1 - e2)**2 * (1 + 4 * e2 + e4) /
                      (1 + 9 * e2 + 9 * e4 + e6))
        d = d + Z[15] * np.sqrt(12) * (
            -8 * x * (3 * e6 - 5 * e4 * x2 - 15 * e4 * y2 + 12 * e4 -
                      20 * e2 * x2 - 60 * e2 * y2 + 12 * e2 - 5 * x2 -
                      15 * y2 + 3)) * num / den
        d = d + Z[16] * np.sqrt(12) * (
            -8 * y * (9 * e6 - 15 * e4 * x2 - 25 * e4 * y2 + 36 * e4 -
                      60 * e2 * x2 - 100 * e2 * y2 + 36 * e2 - 15 * x2 -
                      25 * y2 + 9)) * num / den
        den = (1 - e2)**4 * (1 + e2) * (1 + e4)
        num = np.sqrt((1 - e2)**6 * (1 + e2) * (1 + e4) /
                      (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 +
                       4 * e10 + e12))
        d = d + Z[17] * np.sqrt(12) * (
            4 * x * (6 * e8 - 5 * e6 * x2 - 45 * e6 * y2 + 6 * e6 -
                     5 * e4 * x2 - 45 * e4 * y2 + 6 * e4 - 5 * e2 * x2 -
                     45 * e2 * y2 + 6 * e2 - 5 * x2 - 45 * y2 +
                     6)) * num / den
        d = d + Z[18] * np.sqrt(12) * (
            4 * y * (6 * e8 + 15 * e6 * x2 - 25 * e6 * y2 + 6 * e6 +
                     15 * e4 * x2 - 25 * e4 * y2 + 6 * e4 + 15 * e2 * x2 -
                     25 * e2 * y2 + 6 * e2 + 15 * x2 - 25 * y2 +
                     6)) * num / den
        den = np.sqrt(1 + e2 + e4 + e6 + e8 + e10)
        d = d + Z[19] * np.sqrt(12) * 20 * x * (3 * y2 - x2) / den
        d = d + Z[20] * np.sqrt(12) * 20 * y * (y2 - 3 * x2) / den
        den = (1 - e2)**3
        d = d + Z[21] * np.sqrt(7) * (
            480 * y2 * r2 + 120 * r4 + 24 * e4 - 120 * x2 - 360 * y2 -
            3 * e2 * (40 * x2 + 120 * y2 - 24) + 24) / den

    elif (type == 'dxy'):
        d = Z[0] * 0 * x  # to make d an array with the same size as x
        d = d + Z[1] * 0
        d = d + Z[2] * 0
        d = d + Z[3] * 0
        den = np.sqrt(1 + e2 + e4)
        d = d + Z[4] * np.sqrt(6) * 2 / den
        d = d + Z[5] * 0
        den = np.sqrt((1 - e2)**2 * (1 + e2) * (1 + 4 * e2 + e4))
        d = d + Z[6] * np.sqrt(8) * (1 + e2) * (6 * x) / den
        d = d + Z[7] * np.sqrt(8) * 6 * y * (1 + e2) / den
        den = np.sqrt(1 + e2 + e4 + e6)
        d = d + Z[8] * np.sqrt(8) * 6 * x / den
        d = d + Z[9] * np.sqrt(8) * (-6) * y / den
        den = (1 - e2)**2
        d = d + Z[10] * np.sqrt(5) * 48 * xy / den
        den = (1 - e2)**3 * (1 + e2 + e4)
        num = np.sqrt((1 - e2)**4 * (1 + e2 + e4) /
                      (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8))
        d = d + Z[11] * np.sqrt(10) * 0
        d = d + Z[12] * np.sqrt(10) * ((24 * x2 + 24 * y2 - 6) *
                                       (1 + e2 + e4) - 6 * e6) * num / den
        den = np.sqrt(1 + e2 + e4 + e6 + e8)
        d = d + Z[13] * np.sqrt(10) * (-24) * xy / den
        d = d + Z[14] * np.sqrt(10) * 12 * (x2 - y2) / den
        den = (1 - e2)**3 * (1 + 4 * e2 + e4)
        num = np.sqrt((1 - e2)**2 * (1 + 4 * e2 + e4) /
                      (1 + 9 * e2 + 9 * e4 + e6))
        d = d + Z[15] * np.sqrt(12) * (
            -8 * y * (3 * e6 - 15 * e4 * x2 - 5 * e4 * y2 + 12 * e4 -
                      60 * e2 * x2 - 20 * e2 * y2 + 12 * e2 - 15 * x2 -
                      5 * y2 + 3)) * num / den
        d = d + Z[16] * np.sqrt(12) * (
            -8 * x * (3 * e6 - 5 * e4 * x2 - 15 * e4 * y2 + 12 * e4 -
                      20 * e2 * x2 - 60 * e2 * y2 + 12 * e2 - 5 * x2 -
                      15 * y2 + 3)) * num / den
        den = (1 - e2)**4 * (1 + e2) * (1 + e4)
        num = np.sqrt((1 - e2)**6 * (1 + e2) * (1 + e4) /
                      (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 +
                       4 * e10 + e12))
        d = d + Z[17] * np.sqrt(12) * (
            12 * y * (2 * e8 - 5 * e6 * r2 + 2 * e6 - 5 * e4 * r2 + 2 * e4 -
                      5 * e2 * r2 + 2 * e2 - 5 * r2 + 2)) * num / den
        d = d + Z[18] * np.sqrt(12) * (
            -12 * x * (2 * e8 - 5 * e6 * r2 + 2 * e6 -
                       5 * e4 * r2 + 2 * e4 - 5 * e2 * r2 + 2 * e2 -
                       5 * r2 + 2)) * num / den
        den = np.sqrt(1 + e2 + e4 + e6 + e8 + e10)
        d = d + Z[19] * np.sqrt(12) * 20 * y * (y2 - 3 * x2) / den
        d = d + Z[20] * np.sqrt(12) * 20 * x * (x2 - 3 * y2) / den
        den = (1 - e2)**3
        d = d + Z[21] * np.sqrt(7) * 240 * xy * (2 * r2 - 1 - e2) / den

    return d


def ZernikeGrad(Z, x, y, atype):

    m1, n1 = x.shape
    m2, n2 = y.shape
    if(m1 != m2 or n1 != n2):
        print('x & y are not the same size')

    if(len(Z) > 22):
        print('ZernikeGrad() is not implemented with >22 terms')
        return
    elif len(Z) < 22:
        Z[21] = 0

    x2 = x * x
    y2 = y * y
    xy = x * y
    r2 = x2 + y2

    if (atype == 'dx'):
        d = Z[0] * 0 * x  # to make d an array with the same size as x
        d = d + Z[1] * 2 * 1
        d = d + Z[2] * 2 * 0
        d = d + Z[3] * np.sqrt(3) * 4 * x
        d = d + Z[4] * np.sqrt(6) * 2 * y
        d = d + Z[5] * np.sqrt(6) * 2 * x
        d = d + Z[6] * np.sqrt(8) * 6 * xy
        d = d + Z[7] * np.sqrt(8) * (9 * x2 + 3 * y2 - 2)
        d = d + Z[8] * np.sqrt(8) * 6 * xy
        d = d + Z[9] * np.sqrt(8) * (3 * x2 - 3 * y2)
        d = d + Z[10] * np.sqrt(5) * 12 * x * (2 * (x2 + y2) - 1)
        d = d + Z[11] * np.sqrt(10) * x * (16 * x2 - 6)
        d = d + Z[12] * np.sqrt(10) * y * (24 * x2 + 8 * y2 - 6)
        d = d + Z[13] * np.sqrt(10) * 4 * x * (x2 - 3 * y2)
        d = d + Z[14] * np.sqrt(10) * 4 * y * (3 * x2 - y2)
        d = d + Z[15] * np.sqrt(12) * (
            x2 * (50.0 * x2 + 60.0 * y2 - 36.0) + y2 * (10.0 * y2 - 12.0) + 3)
        d = d + Z[16] * np.sqrt(12) * (xy * (40.0 * r2 - 24.0))
        d = d + Z[17] * np.sqrt(12) * (
            x2 * (25.0 * x2 - 12.0 - 30.0 * y2) + y2 * (12.0 - 15.0 * y2))
        d = d + Z[18] * np.sqrt(12) * (4.0 * xy *
                                       (-6.0 + 15.0 * x2 + 5.0 * y2))
        d = d + Z[19] * np.sqrt(12) * 5 * (x2 * (x2 - 6 * y2) + y2 * y2)
        d = d + Z[20] * np.sqrt(12) * 20 * xy * (x2 - y2)
        d = d + Z[21] * np.sqrt(7) * 24 * x * (
            1 + x2 * (10 * y2 - 5 + 5 * x2) + y2 * (5 * y2 - 5))

    elif (atype, 'dy'):

        d = Z[0] * 0 * x
        d = d + Z[1] * 2 * 0
        d = d + Z[2] * 2 * 1
        d = d + Z[3] * np.sqrt(3) * 4 * y
        d = d + Z[4] * np.sqrt(6) * 2 * x
        d = d + Z[5] * np.sqrt(6) * (-2) * y
        d = d + Z[6] * np.sqrt(8) * (3 * x2 + 9 * y2 - 2)
        d = d + Z[7] * np.sqrt(8) * 6 * xy
        d = d + Z[8] * np.sqrt(8) * (3 * x2 - 3 * y2)
        d = d + Z[9] * np.sqrt(8) * (-6) * xy
        d = d + Z[10] * np.sqrt(5) * 12 * y * (2 * (x2 + y2) - 1)
        d = d + Z[11] * np.sqrt(10) * y * (6 - 16 * y2)
        d = d + Z[12] * np.sqrt(10) * x * (8 * x2 + 24 * y2 - 6)
        d = d + Z[13] * np.sqrt(10) * 4 * y * (y2 - 3 * x2)
        d = d + Z[14] * np.sqrt(10) * 4 * x * (x2 - 3 * y2)
        d = d + Z[15] * np.sqrt(12) * (xy * (40.0 * r2 - 24.0))
        d = d + Z[16] * np.sqrt(12) * (
            x2 * (10.0 * x2 + 60.0 * y2 - 12.0) + y2 * (50.0 * y2 - 36.0) + 3)
        d = d + Z[17] * np.sqrt(12) * (4.0 * xy * (6.0 - 5.0 * x2 - 15.0 * y2))
        d = d + Z[18] * np.sqrt(12) * (
            y2 * (-25.0 * y2 + 12.0 + 30.0 * x2) + x2 * (-12.0 + 15.0 * x2))
        d = d + Z[19] * np.sqrt(12) * 20 * xy * (y2 - x2)
        d = d + Z[20] * np.sqrt(12) * 5 * (x2 * (x2 - 6 * y2) + y2 * y2)
        d = d + Z[21] * np.sqrt(7) * 24 * y * (
            1 + y2 * (10 * x2 - 5 + 5 * y2) + x2 * (5 * x2 - 5))

    return d


def ZernikeJacobian(Z, x, y, atype):

    m1, n1 = x.shape
    m2, n2 = y.shape
    if(m1 != m2 or n1 != n2):
        print('x & y are not the same size')

    if(len(Z) > 22):
        print('ZernikeGrad() is not implemented with >22 terms')
        return
    elif len(Z) < 22:
        Z[21] = 0

    x2 = x * x
    y2 = y * y
    xy = x * y
    r2 = x2 + y2

    if (atype == '1st'):
        j = Z[0] * 0 * x  # to make d an array with the same size as x
        j = j + Z[1] * 0
        j = j + Z[2] * 0
        j = j + Z[3] * np.sqrt(3) * 8
        j = j + Z[4] * np.sqrt(6) * 0
        j = j + Z[5] * np.sqrt(6) * 0
        j = j + Z[6] * np.sqrt(8) * 24 * y  # W8 in Roddier's 1993 table
        j = j + Z[7] * np.sqrt(8) * 24 * x  # W7 in Roddier's 1993 table
        j = j + Z[8] * np.sqrt(8) * 0
        j = j + Z[9] * np.sqrt(8) * 0
        j = j + Z[10] * np.sqrt(5) * (96 * r2 - 24)
        j = j + Z[11] * np.sqrt(10) * 48 * (x2 - y2)
        j = j + Z[12] * np.sqrt(10) * 96 * xy
        j = j + Z[13] * np.sqrt(10) * 0
        j = j + Z[14] * np.sqrt(10) * 0
        j = j + Z[15] * np.sqrt(12) * x * (240.0 * r2 - 96.0)
        j = j + Z[16] * np.sqrt(12) * y * (240.0 * (x2 + y2) - 96.0)
        j = j + Z[17] * np.sqrt(12) * 80.0 * x * (x2 - 3.0 * y2)
        j = j + Z[18] * np.sqrt(12) * 80.0 * y * (3 * x2 - y2)
        j = j + Z[19] * np.sqrt(12) * 0
        j = j + Z[20] * np.sqrt(12) * 0
        j = j + Z[21] * np.sqrt(7) * 48 * (
            1 + x2 * (30 * y2 + 15 * x2 - 10) + y2 * (15 * y2 - 10))

    elif (atype == '2nd'):

        j = Z[0]**2 * 0 * x  # to make d an array with the same size as x
        j = j + Z[1]**2 * 0
        j = j + Z[2]**2 * 0
        j = j + Z[3]**2 * (3) * 16
        j = j + Z[4]**2 * (6) * (-4)
        j = j + Z[5]**2 * (6) * (-4)
        # W8 in Roddier's 1993 table
        j = j + Z[6]**2 * (8) * (108 * y2 - 36 * x2)
        # W7 in Roddier's 1993 table
        j = j + Z[7]**2 * (8) * (108 * x2 - 36 * y2)
        j = j + Z[8]**2 * (8) * (-36 * r2)
        j = j + Z[9]**2 * (8) * (-36 * r2)
        j = j + Z[10]**2 * (5) * 144 * (12 * r2**2 - 8 * r2 + 1)
        j = j + Z[11]**2 * (10) * 36 * (8 * x2 - 1) * (1 - 8 * y2)
        j = j + Z[12]**2 * (10) * 36 * (8 * r2 - 16 * (x2 - y2)**2 - 1)
        j = j + Z[13]**2 * (10) * (-144) * r2**2
        j = j + Z[14]**2 * (10) * (-144) * r2**2
        j = j + Z[15]**2 * (12) * 64 * (5.0 * (x2 + y2) - 3) * \
            (x2 * (25.0 * x2 + 20.0 * y2 - 9) - y2 * (5.0 * y2 - 3.0))
        j = j + Z[16]**2 * (12) * 64 * (5.0 * (x2 + y2) - 3) * \
            (y2 * (25.0 * y2 + 20.0 * x2 - 9) - x2 * (5.0 * x2 - 3.0))
        j = j + Z[17]**2 * (12) * 16.0 * (
            x2 * (-36.0 + 360 * y2 + x2 * (
                180 - 1275 * y2 - 125 * x2)) + y2 * (
                -36 + y2 * (180 + 225 * (x2 - y2))))
        j = j + Z[18]**2 * (12) * 16.0 * (
            y2 * (-36.0 + 360 * x2 + y2 * (
                180 - 1275 * x2 - 125 * y2)) + x2 * (
                -36 + x2 * (180 + 225 * (y2 - x2))))
        j = j + Z[19]**2 * (12) * (-400) * r2**3
        j = j + Z[20]**2 * (12) * (-400) * r2**3
        j = j + Z[21]**2 * (7) * 576 * (
            1 + x2 * (10 * y2 + 5 * x2 - 5) + y2 * (
                5 * y2 - 5)) * (1 + x2 * (25 * x2 + 50 * y2 - 15) +
                                y2 * (25 * y2 - 15))

    return j


def ZernikeAnnularJacobian(Z, x, y, e, atype):

    m1, n1 = x.shape
    m2, n2 = y.shape

    if(m1 != m2 or n1 != n2):
        print('x & y are not the same size')
        exit()

    if(len(Z) > 22):
        print('ZernikeAnnularEval() is not implemented with >22 terms')
        return
    elif len(Z) < 22:
        Z[21] = 0

    x2 = x * x
    y2 = y * y
    xy = x * y
    r2 = x2 + y2
    x4 = x2 * x2
    x6 = x4 * x2
    y4 = y2 * y2
    y6 = y4 * y2
    e2 = e * e
    e4 = e2 * e2
    e6 = e4 * e2
    e8 = e6 * e2
    e10 = e8 * e2
    e12 = e10 * e2
    e14 = e12 * e2
    e16 = e14 * e2

    if (atype == '1st'):
        j = Z[0] * 0 * x  # to make d an array with the same size as x
        j = j + Z[1] * 0
        j = j + Z[2] * 0
        den = 1 - e**2
        j = j + Z[3] * np.sqrt(3) * 8 / den
        j = j + Z[4] * np.sqrt(6) * 0
        j = j + Z[5] * np.sqrt(6) * 0
        den = np.sqrt((1 - e2)**2 * (1 + e2) * (1 + 4 * e2 + e4))
        j = j + Z[6] * np.sqrt(8) * 24 * y * (1 + e2) / den
        j = j + Z[7] * np.sqrt(8) * 24 * x * (1 + e2) / den
        j = j + Z[8] * np.sqrt(8) * 0
        j = j + Z[9] * np.sqrt(8) * 0
        den = (1 - e2)**2
        j = j + Z[10] * np.sqrt(5) * (96 * r2 - 24 * (1 + e2)) / den
        den = (1 - e2)**3 * (1 + e2 + e4)
        num = np.sqrt((1 - e2)**4 * (1 + e2 + e4) /
                      (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8))
        j = j + Z[11] * np.sqrt(10) * 48 * (x2 - y2) * \
            (1 + e2 + e4) * num / den
        j = j + Z[12] * np.sqrt(10) * 96 * xy * (1 + e2 + e4) * num / den
        j = j + Z[13] * np.sqrt(10) * 0
        j = j + Z[14] * np.sqrt(10) * 0
        den = (1 - e2)**3 * (1 + 4 * e2 + e4)
        num = np.sqrt((1 - e2)**2 * (1 + 4 * e2 + e4) /
                      (1 + 9 * e2 + 9 * e4 + e6))
        j = j + Z[15] * np.sqrt(12) * 48 * x * (
            5 * r2 * (1 + 4 * e2 + e4) - 2 *
            (1 + 4 * e2 + 4 * e4 + e6)) * num / den
        j = j + Z[16] * np.sqrt(12) * 48 * y * (
            5 * r2 * (1 + 4 * e2 + e4) - 2 *
            (1 + 4 * e2 + 4 * e4 + e6)) * num / den
        den = (1 - e2)**4 * (1 + e2) * (1 + e4)
        num = np.sqrt((1 - e2)**6 * (1 + e2) * (1 + e4) /
                      (1 + 4 * e2 + 10 * e4 + 20 * e6 +
                       10 * e8 + 4 * e10 + e12))
        j = j + Z[17] * np.sqrt(12) * 80.0 * x * \
            (x2 - 3.0 * y2) * (1 + e2) * (1 + e4) * num / den
        j = j + Z[18] * np.sqrt(12) * 80.0 * y * \
            (3 * x2 - y2) * (1 + e2) * (1 + e4) * num / den
        j = j + Z[19] * np.sqrt(12) * 0
        j = j + Z[20] * np.sqrt(12) * 0
        den = (1 - e2)**3
        j = j + Z[21] * np.sqrt(7) * 48 * (
            e4 - 10 * e2 * x2 - 10 * e2 * y2 +
            3 * e2 + 15 * x4 + 30 * x2 * y2 - 10 * x2 +
            15 * y4 - 10 * y2 + 1) / den
    elif (atype == '2nd'):

        j = Z[0]**2 * 0 * x  # to make d an array with the same size as x
        j = j + Z[1]**2 * 0
        j = j + Z[2]**2 * 0
        den = 1 - e**2
        j = j + Z[3]**2 * (3) * 16 / den / den
        den = (1 + e2 + e4)
        j = j + Z[4]**2 * (6) * (-4) / den
        j = j + Z[5]**2 * (6) * (-4) / den
        den = (1 - e2)**2 * (1 + e2) * (1 + 4 * e2 + e4)
        j = j + Z[6]**2 * (8) * (108 * y2 - 36 * x2) * (1 + e2) / den
        j = j + Z[7]**2 * (8) * (108 * x2 - 36 * y2) * (1 + e2) / den
        den = (1 + e2 + e4 + e6)
        j = j + Z[8]**2 * (8) * (-36 * r2) / den
        j = j + Z[9]**2 * (8) * (-36 * r2) / den
        den = (1 - e2)**4
        j = j + Z[10]**2 * (5) * 144 * (1 + e2 - 2 * r2) * \
            (1 + e2 - 6 * r2) / den
        den = (1 - e2)**6 * (1 + e2 + e4)**2
        num = ((1 - e2)**4 * (1 + e2 + e4) /
               (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8))
        j = j + Z[11]**2 * (10) * 36 * (
            8 * (1 + e2 + e4) * x2 - 1 - e2 - e4 - e6) * \
            (1 + e2 + e4 + e6 - 8 * (1 + e2 + e4) * y2) * num / den
        j = j + Z[12]**2 * (10) * 36 * (
            -4 * (x - y)**2 * (e4 + e2 + 1) + 1 + e2 + e4 + e6) * \
            (4 * (x + y)**2 * (e4 + e2 + 1) - 1 - e2 - e4 - e6) * num / den
        den = (1 + e2 + e4 + e6 + e8)
        j = j + Z[13]**2 * (10) * (-144) * r2**2 / den
        j = j + Z[14]**2 * (10) * (-144) * r2**2 / den
        den = (1 - e2)**6 * (1 + 4 * e2 + e4)**2
        num = (1 - e2)**2 * (1 + 4 * e2 + e4) / (1 + 9 * e2 + 9 * e4 + e6)
        j = j + Z[15]**2 * (12) * 64 * (
            (3 * e6 - 5 * e4 * r2 + 12 * e4 - 20 * e2 * r2 +
             12 * e2 - 5 * r2 + 3) *
            (9 * e6 * x2 - 3 * e6 * y2 - 25 * e4 * x4 - 20 * e4 * x2 * y2 +
             36 * e4 * x2 + 5 * e4 * y4 - 12 * e4 * y2 - 100 * e2 * x4 -
             80 * e2 * x2 * y2 + 36 * e2 * x2 + 20 * e2 * y4 -
             12 * e2 * y2 - 25 * x4 - 20 * x2 * y2 +
             9 * x2 + 5 * y4 - 3 * y2)) * num / den
        j = j + Z[16]**2 * (12) * 64 * (
            -(3 * e6 - 5 * e4 * r2 + 12 * e4 - 20 * e2 * r2 + 12 * e2 -
              5 * r2 + 3) * (3 * e6 * x2 - 9 * e6 * y2 - 5 * e4 * x4 +
                             20 * e4 * x2 * y2 + 12 * e4 * x2 + 25 * e4 * y4 -
                             36 * e4 * y2 - 20 * e2 * x4 + 80 * e2 * x2 * y2 +
                             12 * e2 * x2 + 100 * e2 * y4 -
                             36 * e2 * y2 - 5 * x4 +
                             20 * x2 * y2 + 3 * x2 + 25 * y4 -
                             9 * y2)) * num / den
        den = (1 - e2)**8 * (1 + e2)**2 * (1 + e4)**2
        num = (1 - e2)**6 * (1 + e2) * (1 + e4) / \
            (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 + 4 * e10 + e12)
        j = j + Z[17]**2 * (12) * 16.0 * (
            - 36 * e16 * x2 - 36 * e16 * y2 + 180 * e14 * x4 +
            360 * e14 * x2 * y2 - 72 * e14 * x2 +
            180 * e14 * y4 - 72 * e14 * y2 - 125 * e12 * x6 -
            1275 * e12 * x4 * y2 + 360 * e12 * x4 + 225 * e12 * x2 * y4 +
            720 * e12 * x2 * y2 - 108 * e12 * x2 - 225 * e12 *
            y6 + 360 * e12 * y4 - 108 * e12 * y2 - 250 * e10 * x6 -
            2550 * e10 * x4 * y2 + 540 * e10 * x4 + 450 * e10 * x2 *
            y4 + 1080 * e10 * x2 * y2 - 144 * e10 * x2 - 450 * e10 * y6 +
            540 * e10 * y4 - 144 * e10 * y2 - 375 * e8 * x6 - 3825 *
            e8 * x4 * y2 + 720 * e8 * x4 + 675 * e8 * x2 * y4 +
            1440 * e8 * x2 * y2 - 180 * e8 * x2 - 675 * e8 * y6 + 720 *
            e8 * y4 - 180 * e8 * y2 - 500 * e6 * x6 - 5100 * e6 * x4 * y2 +
            720 * e6 * x4 + 900 * e6 * x2 * y4 + 1440 * e6 * x2 * y2 -
            144 * e6 * x2 - 900 * e6 * y6 + 720 * e6 * y4 - 144 * e6 * y2 -
            375 * e4 * x6 - 3825 * e4 * x4 * y2 + 540 * e4 * x4 + 675 * e4 *
            x2 * y4 + 1080 * e4 * x2 * y2 - 108 * e4 * x2 - 675 * e4 * y6 +
            540 * e4 * y4 - 108 * e4 * y2 - 250 * e2 * x6 - 2550 * e2 * x4 *
            y2 + 360 * e2 * x4 + 450 * e2 * x2 * y4 + 720 * e2 * x2 * y2 -
            72 * e2 * x2 - 450 * e2 * y6 + 360 * e2 * y4 - 72 * e2 *
            y2 - 125 * x6 - 1275 * x4 * y2 + 180 * x4 + 225 * x2 * y4 +
            360 * x2 * y2 - 36 * x2 - 225 * y6 + 180 * y4 -
            36 * y2) * num / den
        j = j + Z[18]**2 * (12) * 16.0 * ((
            - 225 * e12 - 450 * e10 - 675 * e8 - 900 * e6 - 675 * e4 -
            450 * e2 - 225) * x6 +
            (180 * e14 + 225 * e12 * y2 + 360 * e12 + 450 * e10 * y2 +
             540 * e10 + 675 * e8 * y2 + 720 * e8 + 900 * e6 * y2 +
             720 * e6 + 675 * e4 * y2 + 540 * e4 + 450 * e2 * y2 +
             360 * e2 + 225 * y2 + 180) * x4 +
            (- 36 * e16 + 360 * e14 * y2 - 72 * e14 - 1275 * e12 * y4 +
             720 * e12 * y2 - 108 * e12 - 2550 * e10 * y4 +
             1080 * e10 * y2 - 144 * e10 - 3825 * e8 * y4 + 1440 *
             e8 * y2 - 180 * e8 - 5100 * e6 * y4 + 1440 * e6 * y2 -
             144 * e6 - 3825 * e4 * y4 + 1080 * e4 * y2 - 108 * e4 -
             2550 * e2 * y4 + 720 * e2 * y2 - 72 * e2 - 1275 * y4 +
             360 * y2 - 36) * x2 - 36 * e16 * y2 + 180 * e14 * y4 -
            72 * e14 * y2 - 125 * e12 * y6 + 360 * e12 * y4 - 108 * e12 * y2 -
            250 * e10 * y6 + 540 * e10 * y4 - 144 * e10 * y2 - 375 *
            e8 * y6 + 720 * e8 * y4 - 180 * e8 * y2 - 500 * e6 * y6 +
            720 * e6 * y4 - 144 * e6 * y2 - 375 * e4 * y6 + 540 * e4 * y4 -
            108 * e4 * y2 - 250 * e2 * y6 + 360 * e2 * y4 - 72 * e2 * y2 -
            125 * y6 + 180 * y4 - 36 * y2) * num / den
        den = (1 + e2 + e4 + e6 + e8 + e10)
        j = j + Z[19]**2 * (12) * (-400) * r2**3 / den
        j = j + Z[20]**2 * (12) * (-400) * r2**3 / den
        den = (1 - e2)**6
        j = j + Z[21]**2 * (7) * 576 * ((
            e4 - 5 * e2 * x2 - 5 * e2 * y2 + 3 * e2 + 5 * x4 +
            10 * x2 * y2 - 5 * x2 + 5 * y4 - 5 * y2 + 1) *
            (e4 - 15 * e2 * x2 - 15 * e2 * y2 + 3 * e2 +
             25 * x4 + 50 * x2 * y2 - 15 * x2 + 25 * y4 - 15 * y2 + 1)) / den

    return j


def poly10_2D(c, data, y=None):

    if (y is None):
        x = data[0, :]
        y = data[1, :]
    else:
        x = data

    F = c[0] + c[1] * x + c[2] * y + c[3] * x * x + \
        c[4] * x * y + c[5] * y * y + c[6] * x**3 + \
        c[7] * x**2 * y + c[8] * x * y**2 + c[9] * y**3 + \
        c[10] * x**4 + c[11] * x**3 * y + c[12] * x**2 * y**2 + \
        c[13] * x * y**3 + c[14] * y**4 + c[15] * x**5 + \
        c[16] * x**4 * y + c[17] * x**3 * y**2 + c[18] * x**2 * y**3 + \
        c[19] * x * y**4 + c[20] * y**5 + c[21] * x**6 + \
        c[22] * x**5 * y + c[23] * x**4 * y**2 + c[24] * x**3 * y**3 + \
        c[25] * x**2 * y**4 + c[26] * x * y**5 + c[27] * y**6 + \
        c[28] * x**7 + c[29] * x**6 * y + c[30] * x**5 * y**2 + \
        c[31] * x**4 * y**3 + c[32] * x**3 * y**4 + c[33] * x**2 * y**5 + \
        c[34] * x * y**6 + c[35] * y**7 + c[36] * x**8 + \
        c[37] * x**7 * y + c[38] * x**6 * y**2 + c[39] * x**5 * y**3 + \
        c[40] * x**4 * y**4 + c[41] * x**3 * y**5 + c[42] * x**2 * y**6 + \
        c[43] * x * y**7 + c[44] * y**8 + c[45] * x**9 + \
        c[46] * x**8 * y + c[47] * x**7 * y**2 + c[48] * x**6 * y**3 + \
        c[49] * x**5 * y**4 + c[50] * x**4 * y**5 + c[51] * x**3 * y**6 + \
        c[52] * x**2 * y**7 + c[53] * x * y**8 + c[54] * y**9 + \
        c[55] * x**10 + c[56] * x**9 * y + c[57] * x**8 * y**2 + \
        c[58] * x**7 * y**3 + c[59] * x**6 * y**4 + c[60] * x**5 * y**5 + \
        c[61] * x**4 * y**6 + c[62] * x**3 * y**7 + c[63] * x**2 * y**8 + \
        c[64] * x * y**9 + c[65] * y**10

    return F


def poly10Grad(c, x, y, atype):

    if (atype == 'dx'):
        out = c[1] + c[3] * 2 * x + c[4] * y + c[6] * 3 * x**2 + \
            c[7] * 2 * x * y + c[8] * y**2 + c[10] * 4 * x**3 + \
            c[11] * 3 * x**2 * y + c[12] * 2 * x * y**2 + c[13] * y**3 + \
            c[15] * 5 * x**4 + c[16] * 4 * x**3 * y + \
            c[17] * 3 * x**2 * y**2 + \
            c[18] * 2 * x * y**3 + c[19] * y**4 + c[21] * 6 * x**5 + \
            c[22] * 5 * x**4 * y + c[23] * 4 * x**3 * y**2 + \
            c[24] * 3 * x**2 * y**3 + c[25] * 2 * x * y**4 + c[26] * y**5 + \
            c[28] * 7 * x**6 + c[29] * 6 * x**5 * y + \
            c[30] * 5 * x**4 * y**2 + \
            c[31] * 4 * x**3 * y**3 + c[32] * 3 * x**2 * y**4 + \
            c[33] * 2 * x * y**5 + c[34] * y**6 + \
            c[36] * 8 * x**7 + c[37] * 7 * x**6 * y + \
            c[38] * 6 * x**5 * y**2 + \
            c[39] * 5 * x**4 * y**3 + c[40] * 4 * x**3 * y**4 + \
            c[41] * 3 * x**2 * y**5 + c[42] * 2 * x * y**6 + c[43] * y**7 + \
            c[45] * 9 * x**8 + c[46] * 8 * x**7 * y + \
            c[47] * 7 * x**6 * y**2 + \
            c[48] * 6 * x**5 * y**3 + c[49] * 5 * x**4 * y**4 + \
            c[50] * 4 * x**3 * y**5 + c[51] * 3 * x**2 * y**6 + \
            c[52] * 2 * x * y**7 + c[53] * y**8 + c[55] * 10 * x**9 + \
            c[56] * 9 * x**8 * y + c[57] * 8 * x**7 * y**2 + \
            c[58] * 7 * x**6 * y**3 + c[59] * 6 * x**5 * y**4 + \
            c[60] * 5 * x**4 * y**5 + c[61] * 4 * x**3 * y**6 + \
            c[62] * 3 * x**2 * y**7 + c[63] * 2 * x * y**8 + c[64] * y**9

    elif (atype == 'dy'):
        out = c[2] + c[4] * x + c[5] * 2 * y + c[7] * x**2 + \
            c[8] * x * 2 * y + c[9] * 3 * y**2 + c[11] * x**3 + \
            c[12] * x**2 * 2 * y + c[13] * x * 3 * y**2 + c[14] * 4 * y**3 + \
            c[16] * x**4 + c[17] * x**3 * 2 * y + c[18] * x**2 * 3 * y**2 + \
            c[19] * x * 4 * y**3 + c[20] * 5 * y**4 + c[22] * x**5 + \
            c[23] * x**4 * 2 * y + c[24] * x**3 * 3 * y**2 + \
            c[25] * x**2 * 4 * y**3 + c[26] * x * 5 * y**4 + \
            c[27] * 6 * y**5 + c[29] * x**6 + c[30] * x**5 * 2 * y + \
            c[31] * x**4 * 3 * y**2 + c[32] * x**3 * 4 * y**3 + \
            c[33] * x**2 * 5 * y**4 + c[34] * x * 6 * y**5 + \
            c[35] * 7 * y**6 + \
            c[37] * x**7 + c[38] * x**6 * 2 * y + c[39] * x**5 * 3 * y**2 + \
            c[40] * x**4 * 4 * y**3 + c[41] * x**3 * 5 * y**4 + \
            c[42] * x**2 * 6 * y**5 + c[43] * x * 7 * y**6 + \
            c[44] * 8 * y**7 + \
            c[46] * x**8 + c[47] * x**7 * 2 * y + c[48] * x**6 * 3 * y**2 + \
            c[49] * x**5 * 4 * y**3 + c[50] * x**4 * 5 * y**4 + \
            c[51] * x**3 * 6 * y**5 + c[52] * x**2 * 7 * y**6 + \
            c[53] * x * 8 * y**7 + c[54] * 9 * y**8 + c[56] * x**9 + \
            c[57] * x**8 * 2 * y + c[58] * x**7 * 3 * y**2 + \
            c[59] * x**6 * 4 * y**3 + c[60] * x**5 * 5 * y**4 + \
            c[61] * x**4 * 6 * y**5 + c[62] * x**3 * 7 * y**6 + \
            c[63] * x**2 * 8 * y**7 + c[64] * x * 9 * y**8 + c[65] * 10 * y**9

    return out


def extractArray(inArray, dim):

    m, n = inArray.shape
    if m != n:
        print('extractArray: array is not square')

    if m < dim:
        print('extractArray: array is smaller than dimension')

    # print "DIMEN", dim
    i = np.floor((m - dim) / 2)
    j = i + dim
    out = inArray[i:j, i:j]

    return out


def ZernikeMaskedFit(S, x, y, numTerms, mask, e):

    j, i = np.nonzero(mask[:])
    S = S[i, j]
    x = x[i, j]
    y = y[i, j]
    if (e > 0):
        Z = ZernikeAnnularFit(S, x, y, numTerms, e)
    else:
        Z = ZernikeFit(S, x, y, numTerms)
    return Z


def ZernikeFit(S, x, y, numTerms):
    # print x.shape
    # if x,y are 2D, m1,m2 are lists, still (m1!=m2) below works
    m1 = x.shape
    m2 = y.shape
    if((m1 != m2)):
        print('x & y are not the same size')

    S = S[:].copy()
    x = x[:].copy()
    y = y[:].copy()

    i = np.isfinite(S + x + y)
    S = S[i]
    x = x[i]
    y = y[i]

    H = np.zeros((len(S), int(numTerms)))

    for i in range(int(numTerms)):
        Z = np.zeros(int(numTerms))
        Z[i] = 1
        H[:, i] = ZernikeEval(Z, x, y)

    Z = np.dot(np.linalg.pinv(H), S)

    return Z


def ZernikeEval(Z, x, y):
    '''Evaluate Zernicke'''

    # if x,y are 2D, m1,m2 are lists, still (m1!=m2) below works
    m1 = x.shape
    m2 = y.shape

    # print Z.shape
    # print Z

    if((m1 != m2)):
        print('x & y are not the same size')
        exit()

    if(len(Z) > 22):
        print('ZernikeAnnularEval() is not implemented with >22 terms')
        return
    elif len(Z) < 22:
        Z[21] = 0

    r2 = x * x + y * y
    r = np.sqrt(r2)
    r3 = r2 * r
    r4 = r2 * r2
    r5 = r3 * r2
    r6 = r3 * r3

    t = np.arctan2(y, x)
    s = np.sin(t)
    c = np.cos(t)
    s2 = np.sin(2 * t)
    c2 = np.cos(2 * t)
    s3 = np.sin(3 * t)
    c3 = np.cos(3 * t)
    s4 = np.sin(4 * t)
    c4 = np.cos(4 * t)
    s5 = np.sin(5 * t)
    c5 = np.cos(5 * t)

    S = Z[0] * (1 + 0 * x)  # 0*x to set NaNs properly
    S = S + Z[1] * 2 * r * c
    S = S + Z[2] * 2 * r * s
    S = S + Z[3] * np.sqrt(3) * (2 * r2 - 1)
    S = S + Z[4] * np.sqrt(6) * r2 * s2
    S = S + Z[5] * np.sqrt(6) * r2 * c2
    S = S + Z[6] * np.sqrt(8) * (3 * r3 - 2 * r) * s
    S = S + Z[7] * np.sqrt(8) * (3 * r3 - 2 * r) * c
    S = S + Z[8] * np.sqrt(8) * r3 * s3
    S = S + Z[9] * np.sqrt(8) * r3 * c3
    S = S + Z[10] * np.sqrt(5) * (6 * r4 - 6 * r2 + 1)
    S = S + Z[11] * np.sqrt(10) * (4 * r4 - 3 * r2) * c2
    S = S + Z[12] * np.sqrt(10) * (4 * r4 - 3 * r2) * s2
    S = S + Z[13] * np.sqrt(10) * r4 * c4
    S = S + Z[14] * np.sqrt(10) * r4 * s4
    S = S + Z[15] * np.sqrt(12) * (10 * r5 - 12 * r3 + 3 * r) * c
    S = S + Z[16] * np.sqrt(12) * (10 * r5 - 12 * r3 + 3 * r) * s
    S = S + Z[17] * np.sqrt(12) * (5 * r5 - 4 * r3) * c3
    S = S + Z[18] * np.sqrt(12) * (5 * r5 - 4 * r3) * s3
    S = S + Z[19] * np.sqrt(12) * r5 * c5
    S = S + Z[20] * np.sqrt(12) * r5 * s5
    S = S + Z[21] * np.sqrt(7) * (20 * r6 - 30 * r4 + 12 * r2 - 1)

    return S


def ZernikeAnnularFit(S, x, y, numTerms, e):

    m1 = x.shape
    m2 = y.shape
    if (m1 != m2):
        print('x & y are not the same size')

    S = S[:].copy()
    x = x[:].copy()
    y = y[:].copy()

    i = np.isfinite(S + x + y)
    S = S[i]
    x = x[i]
    y = y[i]

    H = np.zeros((len(S), numTerms))

    for i in range(numTerms):
        Z = np.zeros((numTerms))
        Z[i] = 1
        H[:, i] = ZernikeAnnularEval(Z, x, y, e)

    Z = np.dot(np.linalg.pinv(H), S)

    return Z


def ZernikeAnnularEval(Z, x, y, e):
    '''Evaluate the Annular Zernickes'''

    m1 = x.shape
    m2 = y.shape

    if (m1 != m2):
        print('x & y are not the same size')
        exit()

    if(len(Z) > 22):
        print('ZernikeAnnularEval() is not implemented with >22 terms')
        return
    elif len(Z) < 22:
        Z[21] = 0

    r2 = x * x + y * y
    r = np.sqrt(r2)
    r3 = r2 * r
    r4 = r2 * r2
    r5 = r3 * r2
    r6 = r3 * r3

    t = np.arctan2(y, x)
    s = np.sin(t)
    c = np.cos(t)
    s2 = np.sin(2 * t)
    c2 = np.cos(2 * t)
    s3 = np.sin(3 * t)
    c3 = np.cos(3 * t)
    s4 = np.sin(4 * t)
    c4 = np.cos(4 * t)
    s5 = np.sin(5 * t)
    c5 = np.cos(5 * t)

    e2 = e * e
    e4 = e2 * e2
    e6 = e4 * e2
    e8 = e6 * e2
    e10 = e8 * e2
    e12 = e10 * e2

    S = Z[0] * (1 + 0 * x)  # 0*x to set NaNs properly

    den = np.sqrt(1 + e2)
    S = S + Z[1] * 2 * r * c / den
    S = S + Z[2] * 2 * r * s / den

    den = 1 - e**2
    S = S + Z[3] * np.sqrt(3) * (2 * r2 - 1 - e2) / den

    den = np.sqrt(1 + e2 + e4)
    S = S + Z[4] * np.sqrt(6) * r2 * s2 / den
    S = S + Z[5] * np.sqrt(6) * r2 * c2 / den

    den = np.sqrt((1 - e2)**2 * (1 + e2) * (1 + 4 * e2 + e4))
    S = S + Z[6] * np.sqrt(8) * (3 * r3 - 2 * r - 2 *
                                 e4 * r + e2 * r * (3 * r2 - 2)) * s / den
    S = S + Z[7] * np.sqrt(8) * (3 * r3 - 2 * r - 2 *
                                 e4 * r + e2 * r * (3 * r2 - 2)) * c / den

    den = np.sqrt(1 + e2 + e4 + e6)
    S = S + Z[8] * np.sqrt(8) * r3 * s3 / den
    S = S + Z[9] * np.sqrt(8) * r3 * c3 / den

    den = (1 - e2)**2
    S = S + Z[10] * np.sqrt(5) * (6 * r4 - 6 * r2 + 1 +
                                  e4 + e2 * (4 - 6 * r2)) / den

    den = (1 - e2)**3 * (1 + e2 + e4)
    num = np.sqrt((1 - e2)**4 * (1 + e2 + e4) /
                  (1 + 4 * e2 + 10 * e4 + 4 * e6 + e8))
    S = S + Z[11] * np.sqrt(10) * (
        4 * r4 - 3 * r2 - 3 * e6 * r2 - e2 * r2 * (3 - 4 * r2) -
        e4 * r2 * (3 - 4 * r2)) * c2 * num / den
    S = S + Z[12] * np.sqrt(10) * (
        4 * r4 - 3 * r2 - 3 * e6 * r2 - e2 * r2 * (3 - 4 * r2) -
        e4 * r2 * (3 - 4 * r2)) * s2 * num / den

    den = np.sqrt(1 + e2 + e4 + e6 + e8)
    S = S + Z[13] * np.sqrt(10) * r4 * c4 / den
    S = S + Z[14] * np.sqrt(10) * r4 * s4 / den

    den = (1 - e2)**3 * (1 + 4 * e2 + e4)
    numE = np.sqrt((1 - e2)**2 * (1 + 4 * e2 + e4) /
                   (1 + 9 * e2 + 9 * e4 + e6))
    numQ = 10 * r5 - 12 * r3 + 3 * r + 3 * e8 * r - 12 * e6 * r * (r2 - 1) + \
        2 * e4 * r * (15 - 24 * r2 + 5 * r4) + \
        4 * e2 * r * (3 - 12 * r2 + 10 * r4)
    S = S + Z[15] * np.sqrt(12) * numE * numQ * c / den
    S = S + Z[16] * np.sqrt(12) * numE * numQ * s / den

    den = (1 - e2)**4 * (1 + e2) * (1 + e4)
    numE = np.sqrt((1 - e2)**6 * (1 + e2) * (1 + e4) /
                   (1 + 4 * e2 + 10 * e4 + 20 * e6 + 10 * e8 + 4 * e10 + e12))
    numQ = r3 * (5 * r2 - 4 - 4 * e8 - e2 * (4 - 5 * r2) -
                 e4 * (4 - 5 * r2) - e6 * (4 - 5 * r2))
    S = S + Z[17] * np.sqrt(12) * numE * numQ * c3 / den
    S = S + Z[18] * np.sqrt(12) * numE * numQ * s3 / den

    den = np.sqrt(1 + e2 + e4 + e6 + e8 + e10)
    S = S + Z[19] * np.sqrt(12) * r5 * c5 / den
    S = S + Z[20] * np.sqrt(12) * r5 * s5 / den

    den = (1 - e2)**3
    S = S + Z[21] * np.sqrt(7) * (
        20 * r6 - 30 * r4 + 12 * r2 - 1 - e6 +
        3 * e4 * (-3 + 4 * r2) - 3 * e2 * (3 - 12 * r2 + 10 * r4)) / den

    return S


def outParam(filename, algo, inst, I1, I2, model):
    if (filename == ''):
        fout = sys.stdout
    else:
        fout = open(filename, 'w')

    fout.write('intra image: \t %s \t field in deg =(%6.3f, %6.3f)\n' %
               (I1.filename, I1.fieldX, I1.fieldY))
    fout.write('extra image: \t %s \t field in deg =(%6.3f, %6.3f)\n' %
               (I2.filename, I2.fieldX, I2.fieldY))
    fout.write('Using optical model:\t %s\n' % model)
    fout.write('\n')
    finst = open(inst.filename)
    fout.write('---instrument file: --- %s ----------\n' % inst.filename)
    iscomment = False
    for line in finst:
        line = line.strip()
        if (line.startswith('###')):
            iscomment = ~iscomment
        if (not(line.startswith('#')) and (not iscomment) and len(line) > 0):
            fout.write(line + '\n')
    finst.close()

    fout.write('\n')
    falgo = open(algo.filename)
    fout.write('---algorithm file: --- %s ----------\n' % inst.filename)
    iscomment = False
    for line in falgo:
        line = line.strip()
        if (line.startswith('###')):
            iscomment = ~iscomment
        if (not(line.startswith('#')) and (not iscomment) and len(line) > 0):
            fout.write(line + '\n')
    falgo.close()

    if not (filename == ''):
        fout.close()

def outZer4Up(z, unit, filename=''):
    try:
        if unit == 'm':
            z = z * 1e-9
        elif unit == 'nm':
            pass
        elif unit == 'um':
            z = z * 1e-3
        else:
            raise(unknownUnitError)
    except unknownUnitError:
        print('Unknown unit: %s' % unit)
        print('Known options are: m, nm, um')
        sys.exit()

    if (filename == ''):
        f = sys.stdout
    else:
        f = open(filename, 'w')

    for i in range(4, len(z) + 4):
        f.write('%d\t %8.0f\n' % (i, z[i - 4]))
    if not (filename == ''):
        f.close()
         
