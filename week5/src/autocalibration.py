import numpy as np

import utils as h
import reconstruction as rc
import maths as mth

def estimate_aff_hom(cams, vps):
    # your code here
    P1 = cams[0]
    P2 = cams[1]

    vps1 = vps[0]
    vps2 = vps[1]

    #find 3D coordinates by triangulation
    vp_3d = rc.estimate_3d_points(P1, P2, vps1.T, vps2.T)
    u,d,v = np.linalg.svd(np.transpose(vp_3d))
    p = 1*v[:,-1]
    p = p / p[-1]
    aff_hom = np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,1,0],
                        [p[0],p[1],p[2],1]])

    return aff_hom


def estimate_euc_hom(cams,vps):
    P1 = cams[0]
    P2 = cams[1]

    vps1 = vps[0]
    vps2 = vps[1]

    #find 3D coordinates by triangulation
    vp_3d = rc.estimate_3d_points(P1, P2, vps1.T, vps2.T)
    vp_3d = vp_3d.T
    u = vp_3d[0]
    v = vp_3d[1]
    z = vp_3d[1]
    A = np.array([[u[0]*v[0], u[0]*v[1]+u[1]*v[0], u[0]*v[2]+u[2]*v[0], u[1]*v[1], u[1]*v[2]+u[2]*v[1], u[2]*v[2]],
                  [u[0]*z[0], u[0]*z[1]+u[1]*z[0], u[0]*z[2]+u[2]*z[0], u[1]*z[1], u[1]*z[2]+u[2]*z[1], u[2]*z[2]],
                  [v[0]*z[0], v[0]*z[1]+v[1]*z[0], v[0]*z[2]+v[2]*z[0], v[1]*z[1], v[1]*z[2]+v[2]*z[1], v[2]*z[2]],
                  [0,1,0,0,0,0],
                  [1,0,0,-1,0,0]])
    u,d,v = np.linalg.svd(A)
    w = v[:,-1]
    w_arr = np.array([[w[0], w[1], w[2]],
                      [w[1], w[3], w[4]],
                      [w[2], w[4], w[5]]])
    M = P2[:3,:3]
    prod = np.linalg.inv(np.dot(np.dot(np.linalg.inv(M),w_arr),M))
    A = np.linalg.cholesky(prod)
    A_it = np.linalg.inv(A.T)
    H_metric = np.array([[A_it[0,0],A_it[0,1],A_it[0,2], 0],
                  [A_it[1,0], A_it[1,1], A_it[1,2], 0],
                  [A_it[2,0], A_it[2,1], A_it[2,2], 0],
                  [0, 0, 0, 1]])
    return H_metric
