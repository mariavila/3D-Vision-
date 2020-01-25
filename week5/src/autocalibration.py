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
    vp_3d = vp_3d.T
    u,d,v = np.linalg.svd(np.transpose(vp_3d))
    p = v[:,-1]
    aff_hom = np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,1,0],
                        [p[0],p[1],p[2],1]])
    return aff_hom 


