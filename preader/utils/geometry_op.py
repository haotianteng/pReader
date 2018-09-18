#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 01:01:35 2018

@author: Haotian.Teng
"""
import numpy as np
def get_angle(P1,P2,P3,out_format = 'deg'):
    """
    Given three points, calculate their angle by:
        P12 = P1-P2
        P32 = P3-P2
        cos(theta) = P12 . P32 / (|P12||P32|)
    """
    P1 = np.asarray(P1)
    P2 = np.asarray(P2)
    P3 = np.asarray(P3)
    P21 = P1-P2
    P32 = P3-P2
    angle = np.arccos(np.dot(P21,P32)/np.linalg.norm(P21)/np.linalg.norm(P32))
    if out_format == 'deg':
        return np.rad2deg(angle)
    elif out_format == 'rad':
        return angle

def get_len(P1,P2):
    """
    Get the length of vector of P1-P2
    """
    P1 = np.asarray(P1)
    P2 = np.asarray(P2)
    return np.linalg.norm(P1-P2)

def _norm_v(Vector):
    return Vector/np.linalg.norm(Vector)

def get_diherdral(Points,out_format = 'deg'):
    """
    Calculate the diherdral angle given the four atoms
    out_format can be 'deg' or 'rad'
    """
    Ps = [np.asarray(point) for point in Points]
    V1 = Ps[1] - Ps[0]
    V2 = Ps[2] - Ps[1]
    V3 = Ps[3] - Ps[2]
    a = np.cross(V1,V2)
    b = np.cross(V2,V3)
    c = np.cross(a,b)
    if out_format ==  'deg':
        return np.rad2deg(np.arctan2(np.dot(c,_norm_v(V2)),np.dot(a,b)))
    elif out_format == 'rad':
        return np.arctan2(np.dot(c,_norm_v(V2)),np.dot(a,b))