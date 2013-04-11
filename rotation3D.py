import math as ma
from math import sin, cos
import numpy as np
import scipy as sc

'''Short script of rotation matrices
Matthew Newby, November 15,2010'''


def rotx (psi):
    psi = psi*(ma.pi / 180.0)
    R = sc.matrix([
        [1.0,  0.0,       0.0],
        [0.0,  cos(psi), -sin(psi)],
        [0.0,  sin(psi),  cos(psi)]])
    return R

def roty (theta):
    theta = theta*(ma.pi / 180.0)
    R = sc.matrix([
        [cos(theta),  0.0,  sin(theta)],
        [0.0,         1.0,  0.0],
        [-sin(theta), 0.0,  cos(theta)]])
    return R

def rotz (phi):
    phi = phi*(ma.pi / 180.0)
    R = sc.matrix([
        [cos(phi), -sin(phi), 0.0],
        [sin(phi),  cos(phi), 0.0],
        [0.0,       0.0,      1.0]])
    return R

def rot3D (xang, yang, zang):
    """Tait-Bryan angles Convention"""
    xang = xang*(ma.pi / 180.0)
    yang = yang*(ma.pi / 180.0)
    zang = zang*(ma.pi / 180.0)
    xrot = rotx(xang)
    yrot = roty(yang)
    zrot = rotz(zang)
    return (zrot*xrot*yrot)