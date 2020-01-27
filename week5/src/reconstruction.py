import cv2
import numpy as np

import utils as h
import maths as mth
from fundamental import make_homogeneous

def compute_proj_camera(F, i):
    # Result 9.15 of MVG (v = 0, lambda = 1). It assumes P1 = [I|0]
    # your code here
    v = np.zeros(3)
    lamb = 1
    u,d,v = np.linalg.svd(F.T)
    ep = v[:,-1]
    ep = ep/ep[2]
    ep_x = np.array([[0, -ep[2], ep[1]],
                     [ep[2], 0, -ep[0]],
                     [-ep[1], ep[0], 0]])
    P = np.c_[np.dot(ep_x,F), lamb*ep]

    return P

def estimate_3d_points(P1, P2, xr1, xr2):
    # Triangulate 3D points from camera matrices
    Xprj = cv2.triangulatePoints(P1, P2, xr1, xr2)

    # Divide by the last column
    Xprj = Xprj / Xprj[3, :]

    if h.debug >2:
        print("  X estimated:\n", Xprj)

    return Xprj

def compute_reproj_error(xproj, cam1, cam2, x1, x2):

    x_hat1 = cam1 @ xproj
    x_hat2 = cam2 @ xproj

    x_hat1_h = (x_hat1 / x_hat1[2, :])[:-1, :]
    x_hat2_h = (x_hat2 / x_hat2[2, :])[:-1, :]

    error = np.sqrt(((x_hat1_h - x1)**2 + (x_hat2_h - x2)**2).sum())

    return error

def transform(aff_hom, Xprj, cams_pr):
    # your code here
    cams_aff = np.dot(cams_pr,np.linalg.inv(aff_hom.T))

    Xaff = np.dot(aff_hom, Xprj)
    Xaff = Xaff / Xaff[3, :]

    return Xaff, cams_aff

def resection(tracks, img):
    # your code here


    return P
