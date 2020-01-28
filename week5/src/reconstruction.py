import cv2
import numpy as np

import utils as h
import maths as mth
from fundamental import make_homogeneous


def compute_proj_camera(F, i):
    # Result 9.15 of MVG (v = 0, lambda = 1). It assumes P1 = [I|0]
    # your code here
    lamb = 1
    u, d, v = np.linalg.svd(F.T)  # were looking for the epipole on the second image e'
    ep = v.T[:, -1]
    ep = ep / ep[-1]

    ep_x = np.array([[0, -ep[2], ep[1]],
                     [ep[2], 0, -ep[0]],
                     [-ep[1], ep[0], 0]])
    P = np.c_[ep_x @ F, lamb * ep]

    return P


def estimate_3d_points(P1, P2, xr1, xr2):
    # Triangulate 3D points from camera matrices
    Xprj = cv2.triangulatePoints(P1, P2, xr1, xr2)

    # Divide by the last column
    Xprj = Xprj / Xprj[3, :]

    if h.debug > 2:
        print("  X estimated:\n", Xprj)

    return Xprj


def compute_reproj_error(xproj, cam1, cam2, x1, x2):
    x_hat1 = cam1 @ xproj
    x_hat2 = cam2 @ xproj

    x_hat1_h = (x_hat1 / x_hat1[2, :])[:-1, :]
    x_hat2_h = (x_hat2 / x_hat2[2, :])[:-1, :]

    error = np.sqrt(((x_hat1_h - x1) ** 2 + (x_hat2_h - x2) ** 2)).sum()

    return error


def transform(aff_hom, Xprj, cams_pr):
    # your code here
    cams_aff = cams_pr @ np.linalg.inv(aff_hom)
    Xaff = aff_hom @ Xprj
    Xaff = Xaff / Xaff[-1, :]

    return Xaff, cams_aff


def resection(tracks, n_img):
    pts_3d = []
    pts1 = []
    pts2 = []
    for t in tracks:
        if len(t.ref_views) >= 3:
            pts_3d.append(t.pt)
            pts1.append(t.ref_views[n_img - 1])
            pts2.append(t.ref_views[n_img])

    pts1 = make_homogeneous(np.array(pts1))
    pts2 = make_homogeneous(np.array(pts2))

    sig_x = np.std(pts2[:, 0]) / np.sqrt(2)
    sig_y = np.std(pts2[:, 1]) / np.sqrt(2)
    mu_x = np.mean(pts2[:, 0])
    mu_y = np.mean(pts2[:, 1])

    T = np.array([[sig_x, 0, -1 * mu_x],
                  [0, sig_y, -1 * mu_y],
                  [0, 0, 1]])

    pts_3d = np.array(pts_3d)
    sig_x = np.std(pts_3d[:, 0]) / np.sqrt(2)
    sig_y = np.std(pts_3d[:, 1]) / np.sqrt(2)
    sig_z = np.std(pts_3d[:, 2]) / np.sqrt(2)

    mu_x = np.mean(pts_3d[:, 0])
    mu_y = np.mean(pts_3d[:, 1])
    mu_z = np.mean(pts_3d[:, 2])

    UU = np.array([[sig_x, 0, 0, -1 * mu_x],
                   [0, sig_y, 0, -1 * mu_y],
                   [0, 0, sig_z, -1 * mu_z],
                   [0, 0, 0, 1]])

    pt2_n = pts2 @ T
    pt3d_n = pts_3d @ UU

    n_points = np.shape(pts2)[0]
    A = np.zeros((2 * n_points, 12))
    c = 0
    for i in range(0, n_points):
        x2i = pt2_n[i, :]
        x3di = pt3d_n[i, :]
        A[c, :] = np.array([0, 0, 0, 0,
                            -x2i[2] * x3di[0], -x2i[2] * x3di[1], -x2i[2] * x3di[2], -x2i[2] * x3di[3],
                            x2i[1] * x3di[0], x2i[1] * x3di[1], x2i[1] * x3di[2], x2i[1] * x3di[3]])
        A[c + 1, :] = np.array([x2i[2] * x3di[0], x2i[2] * x3di[1], x2i[2] * x3di[2], x2i[2] * x3di[3],
                                0, 0, 0, 0,
                                -x2i[0] * x3di[0], -x2i[0] * x3di[1], -x2i[0] * x3di[2], -x2i[0] * x3di[3]])
        c = c + 2
    u, d, v = np.linalg.svd(A)
    P = v.T[:, -1]
    P = np.array([[P[0], P[1], P[2], P[3]],
                  [P[4], P[5], P[6], P[7]],
                  [P[8], P[9], P[10], P[11]]])

    P = np.linalg.inv(T) @ P
    P = P @ UU
    P = P / P[2, 3]

    return P
