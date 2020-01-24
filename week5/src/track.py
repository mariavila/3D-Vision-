import numpy as np

import utils as h

class Track:
    # A class for storing a 3d point and its views from cameras
    def __init__(self, view, rview, img, hs_vs, pt = None):
        self.views = {}     # dictionary: {#cam, 2d point}
        self.ref_views = {} # dictionary: {#cam, 2d point}. Refined views according to F
        self.views[img] = view
        self.ref_views[img] = rview
        hs_vs[view] = self 

        self.pt = np.zeros(4, dtype=np.float32)
        if pt is not None:
            self.pt = pt

    def add_view(self, view, rview, img, hs_vs=None, pt=None):
        #for key, view in self.views.items():
        #    print("view[",key,"]: ",view)
        #for key, rview in self.ref_views.items():
        #    print("refined view[",key,"]: ", rview)
        if hs_vs is not None:
            # The addition to the hash table is done here
            hs_vs[view] = self 
            if h.debug > 2:
                print("    hs_vs added to Track")

        if img not in self.views: 
            # self doesn't have a view in img
            self.views[img] = view
            self.ref_views[img] = rview
            if h.debug > 2:
                print("    view and rview added to Track")
        elif view is not self.views[img]:
            #conflict of views
            self.deal_with_conflicts(view, self.views[img])
            self.deal_with_conflicts(rview, self.ref_views[img])
            if h.debug > 2:
                print("    Conflict managed")
            # deal with pt
            if pt is not None:
                if pt is None:
                    self.pt = pt
                    if h.debug > 2:
                        print("    pt assigned")
                else:
                    if (np.linalg.norm(np.array(pt) - np.array(self.pt)) < 1):
                        self.pt = (pt + self.pt)*0.5
                        if h.debug > 2:
                            print("    pt averaged")

    def is_in_views(self, img):
        if img not in self.views:
            return False
        else: 
            return True

    def is_in_views(self, view, img):
        if is_in_views(img):
            if view is self.views[img]:
                return True
            else: 
                return False
        else: 
            return False

    def deal_with_conflicts(self, v1, v2):
        new_vw = v1
        if (np.linalg.norm(np.array(v1) - np.array(v2)) < 0.5):
            new_vw = ((v1[0] + v2[0])*0.5, (v1[1] + v2[1])*0.5)
        return new_vw

    def merge (self, views):
        for key in views.views.keys():
            self.add_view(views.views[key], views.ref_views[key], key, None, views.pt)

def add_tracks(xi, xj, xri, xrj, i, j, tracks, hs_vs):
    # Add views matched to the tracks and to the hash table of tracks

    for fi, fj, rfi, rfj in zip(xi, xj, xri, xrj):
        tfi = tuple(fi)
        tfj = tuple(fj)
        fi_is_v = (tfi in hs_vs)
        fj_is_v = (tfj in hs_vs)

        #print ("tfi:",tfi)
        #print ("tfj:",tfj)
        #print ("trfi:",tuple(rfi))
        #print ("trfj:",tuple(rfj))

        if h.debug > 2:
            print ("x1 is in hs_vs?", fi_is_v)
            print ("x2 is in hs_vs?", fj_is_v)
            
        if not fi_is_v and not fj_is_v: 
            # create a new view and add the matches to it and to the hash table of tracks
            v = Track(tfi, tuple(rfi), i, hs_vs)
            v.add_view(tfj, tuple(rfj), j, hs_vs)
            tracks.append(v)
            if h.debug > 2:
                print ("Track created and view added")
                print("track", v,":")
                for key, view in v.views.items():
                    print("view[",key,"]: ",view)
                for key, rview in v.ref_views.items():
                    print("refined view[",key,"]: ", rview)
        elif not fi_is_v: 
            # retrieve the view and add the orphan to it and to the hash table
            v = hs_vs[tfj]
            v.add_view(tfi, tuple(rfi), i, hs_vs)
            if h.debug > 2:
                print ("view 1 added")
        elif not fj_is_v: 
            # retrieve the view and add the orphan to it and to the hash table
            v = hs_vs[tfi]
            v.add_view(tfj, tuple(rfj), j, hs_vs)
            if h.debug > 2:
                print ("view 2 added")
        else:
            # both are in hash table
            if hs_vs[tfi] is not hs_vs[tfj]: 
                # They are different tracks: merge their views into one and update the hash table
                v = hs_vs[tfi]
                w = hs_vs[tfj]
                if v.views.keys() > w.views.keys():
                    # The smaller is added to the bigger
                    v.merge(w)
                    hs_vs[tfj] = v
                    if h.debug > 2:
                        print ("view 2 merged into 1")
                else:
                    w.merge(v)
                    hs_vs[tfi] = w
                    if h.debug > 2:
                        print ("view 1 merged into 2")

            #v = hs_vs[tfi]
            #if not v.is_in_views(tfj, j):
            #    v.add_view(tfj, tuple(rfj), j)

def add_pts_tracks(Xaff, x1, x2, tracks, hs_vs):
    # your code here

