import numpy as np
import matplotlib.pyplot as plt
import scipy

def _s(phi):
    F = scipy.special.ellipkinc(np.arccos((0.5-phi) / (0.5+phi)), 0.5)
    return (F + np.sqrt(2*phi*(4*phi**2+1))) / np.sqrt(18)

s = np.vectorize(_s)

class TrackSequence:
    def __init__(self):
        self.tracks=[]

    def append_track(self, track):
        self.tracks.append(track)

    def total_arclength(self):
        return sum([track.arclength for track in self.tracks])

    def vertices_normals(self, arclengths):
        arclength_thresholds = np.append(np.zeros(1),np.cumsum([track.arclength for track in self.tracks]))
        arclength_conditions = np.array([(arclength_thresholds[i] <= arclengths) & (arclengths < arclength_thresholds[i+1]) \
                                for i in range(len(arclength_thresholds) - 1)])
        result = np.zeros((4,len(arclengths)))
        for i in range(len(arclengths)):
            track_idx = np.where(arclength_conditions[:,i])[0][0]
            result[:,i] = self.tracks[track_idx]._vertices_normals(arclengths[i] - arclength_thresholds[track_idx])
        return result

class Track:
    def __init__(self, arclength, next_track=None):
        self.arclength, self.next_track = arclength, next_track

    def _vertices_normals(self, tline_xs):
        pass

class FermatSpiralTrack(Track):
    #number of turns, phi_start should not be needed
    def __init__(self, r_max, n_turns, phi_start, neg_branch, next_track=None):
        self.phi_max = n_turns * 2 * np.pi
        self.a = r_max / np.sqrt(self.phi_max)
        self.phi_start=phi_start
        self.neg_branch=neg_branch
        self.phi_sampled = np.append(np.linspace(0, 1e-4, 1000000), np.linspace(1e-4, 2*np.pi*3, 1000000))
        try:
            with open('arclength.npy') as f:
                pass
        except FileNotFoundError:
            self.arclength_sampled = s(self.phi_sampled)
            np.save('arclength', self.arclength_sampled)
        self.arclength_sampled = np.load('arclength.npy') * self.a
        #plt.plot(self.phi_sampled[:100], self.arclength_sampled[:100]); plt.show()
        Track.__init__(self, np.interp(self.phi_max, self.phi_sampled, self.arclength_sampled))

    def _vertices_normals(self, tline_xs):
        if not self.neg_branch:
            tline_xs = self.arclength - tline_xs
        self.phis = np.interp(tline_xs, self.arclength_sampled, self.phi_sampled)
        self.rs = self.a * np.sqrt(self.phis)
        self.normal_angles = self.phis + self.phi_start - np.arctan(0.5 / self.phis)
        v = np.array((self.rs * np.cos(self.phis + self.phi_start),
                    self.rs * np.sin(self.phis + self.phi_start),
                      np.cos(self.normal_angles),
                      np.sin(self.normal_angles)))
        if self.neg_branch:
            v[0:2]*=-1
        return v

    def _outer_point_normal(self):
        r_max = self.a * np.sqrt(self.phi_max)
        outer_point = np.array((r_max * np.cos(self.phi_max), r_max * np.sin(self.phi_max)))
        if self.neg_branch:
            outer_point *= -1

        outer_normal_phi = self.phi_max - np.arctan(0.5/self.phi_max)
        outer_normal = np.array((np.cos(outer_normal_phi), np.sin(outer_normal_phi)))

        if self.neg_branch:
            outer_normal *= -1

        return outer_point, outer_normal, outer_normal_phi

    def construct_suitable_ArcTrack(self, Rlaunch):
        outer_point, outer_normal, outer_normal_phi = self._outer_point_normal()
        arc_center_point = outer_point + outer_normal * Rlaunch
        #final_angle = self.phi_max - np.arctan(0.5/self.phi_max)
        if not self.neg_branch:
            return ArcTrack(arc_center_point, Rlaunch, np.pi / 2, outer_normal_phi % np.pi, False)
        else:
            return ArcTrack(arc_center_point, Rlaunch, -outer_normal_phi % np.pi, np.pi / 2, True)

class ArcTrack(Track):
    def __init__(self, center, r, phi_start, phi_end, invert_phi):
        self.center, self.r, self.phi_start, self.phi_end = center, r, phi_start, phi_end
        self.invert_phi=invert_phi
        Track.__init__(self, r * (phi_end-phi_start))

    def _vertices_normals(self, tline_xs):
        self.phis = tline_xs / self.r + self.phi_start
        if self.invert_phi:
            self.phis *= -1
        #     # self.phis = np.flip(self.phis, axis=0) #added so arc sections not in reverse
        v = np.array((self.r * np.cos(self.phis) + self.center[0],
                         self.r * np.sin(self.phis) + self.center[1],
                         np.cos(self.phis),
                         np.sin(self.phis)))
        if not self.invert_phi:
            v[2:4] *= -1
        return v

if __name__ == '__main__':
    Rspiral = 20e-3
    launch_angle = 0
    turns = 1.9
    final_angle = 2 * np.pi * turns
    #final_point = R
    #final_normal_angle = final_angle % (2*np.pi) - np.arctan(0.5/final_angle)
    Rlaunch = 30e-3

    fermat1 = FermatSpiralTrack(Rspiral, turns, 0, True)
    fermat2 = FermatSpiralTrack(Rspiral, turns, 0, False)
    arc2 = fermat2.construct_suitable_ArcTrack(Rlaunch)
    arc1 = fermat1.construct_suitable_ArcTrack(Rlaunch)

    track1 = arc2
    track2 = fermat2
    track3 = fermat1
    track4 = arc1

    track1.next_track, track2.next_track, track3.next_track = track2, track3, track4

    initial_track = track1
    arclengths = np.arange(0, initial_track.total_arclength(), 1e-6)
    xs, ys, us, vs = initial_track.vertices_normals(arclengths)
    fig, axs = plt.subplots()
    axs.set_aspect('equal')
    axs.quiver(xs, ys, us, vs, scale_units='width', scale=100)

    xs,ys,_,_=track1.vertices_normals(np.array([0, track1.arclength, track1.arclength+track2.arclength, track1.arclength+track2.arclength+track3.arclength, track1.arclength+track2.arclength+track3.arclength+track4.arclength-1e-12]))
    axs.scatter(xs, ys, c='r', marker='.')
    A=40e-3
    axs.set_xlim(-A, A), axs.set_ylim(-A,A)
    plt.show()
