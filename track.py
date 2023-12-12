import numpy as np
import matplotlib.pyplot as plt
import scipy

def _s(phi):
    F = scipy.special.ellipkinc(np.arccos((0.5-phi) / (0.5+phi)), 0.5)
    return (F + np.sqrt(2*phi*(4*phi**2+1))) / np.sqrt(18)

s = np.vectorize(_s)

class Track:
    def __init__(self, arclength, next_track=None):
        self.arclength, self.next_track = arclength, next_track

    def vertices_normals(self, arclengths):
        local_vertices_normals = self._vertices_normals(arclengths[arclengths < self.arclength])
        if self.next_track is None:
            return local_vertices_normals
        else:
            return np.append(local_vertices_normals,
                         self.next_track.vertices_normals(arclengths[arclengths >= self.arclength] - self.arclength),
                         axis = 1)

    def total_arclength(self):
        if self.next_track is None:
            return self.arclength
        else:
            return self.arclength + self.next_track.total_arclength()

    def _vertices_normals(self, tline_xs):
        pass

class FermatSpiralTrack(Track):

    #number of turns, phi_start should not be needed
    def __init__(self, r_max, n_turns, phi_start, neg_branch, next_track):
        self.phi_max = n_turns * 2 * np.pi
        self.a = r_max / np.sqrt(self.phi_max)
        self.phi_start=phi_start
        self.neg_branch=neg_branch
        self.phi_sampled = np.append(np.linspace(0, 1e-4, 100000), np.linspace(1e-4, 2*np.pi*n_turns, 100000))
        try:
            with open('arclength.npy') as f:
                pass
        except FileNotFoundError:
            self.arclength_sampled = s(self.phi_sampled)
            np.save('arclength', self.arclength_sampled)
        self.arclength_sampled = np.load('arclength.npy') * self.a
        #plt.plot(self.phi_sampled[:100], self.arclength_sampled[:100]); plt.show()
        Track.__init__(self, self.arclength_sampled[-1], next_track)

    def _vertices_normals(self, tline_xs):
        if not self.neg_branch:
            tline_xs = self.arclength - tline_xs
        self.phis = np.interp(tline_xs, self.arclength_sampled, self.phi_sampled, left=0, right=self.phi_max)
        self.rs = self.a * np.sqrt(self.phis)
        self.normal_angles = self.phis + self.phi_start - np.arctan(0.5 / self.phis)
        v = np.array((self.rs * np.cos(self.phis + self.phi_start),
                    self.rs * np.sin(self.phis + self.phi_start),
                      np.cos(self.normal_angles),
                      np.sin(self.normal_angles)))
        if self.neg_branch:
            v[0:2,:]*=-1
        return v

class ArcTrack(Track):
    def __init__(self, center, r, phi_start, phi_end, invert_phi, next_track):
        self.center, self.r, self.phi_start, self.phi_end = center, r, phi_start, phi_end
        self.invert_phi=invert_phi
        Track.__init__(self, r * np.abs((phi_end-phi_start)), next_track)

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
            v[2:4,:] *= -1
        return v

if __name__ == '__main__':
    Rspiral = 20e-3
    launch_angle = 0
    turns = 1.95
    final_angle = turns * 2 * np.pi - (2 * np.pi) * (turns // 1)
    Rlaunch = 10e-3

    arc1 = ArcTrack((-Rspiral * np.cos(final_angle) - Rlaunch * np.cos(final_angle),
                           -Rspiral * np.sin(final_angle) - Rlaunch * np.sin(final_angle)),
                          Rlaunch, 3/2*np.pi, final_angle, False, next_track=None)
    fermat1 = FermatSpiralTrack(Rspiral, turns, 0, True, next_track=arc1)
    fermat2 = FermatSpiralTrack(Rspiral, turns, 0, False, next_track=fermat1)
    arc2 = ArcTrack(
        ((Rspiral + Rlaunch) * np.cos(final_angle), Rspiral * np.sin(final_angle) + Rlaunch * np.sin(final_angle)),
        Rlaunch,
        np.pi / 2, final_angle % np.pi, False, next_track=fermat2)

    track1 = arc2
    track2 = fermat2
    track3 = fermat1
    track4 = arc1

    initial_track = track1
    arclengths = np.arange(0, initial_track.total_arclength(), 1e-5)
    xs, ys, us, vs = initial_track.vertices_normals(arclengths)
    fig, axs = plt.subplots()
    axs.set_aspect('equal')
    axs.quiver(xs, ys, us, vs, scale_units='width', scale=100)


    # #Rlaunch = Rspiral * np.sin(launch_angle) / (1-np.sin(launch_angle))
    # #Rlaunch = Rspiral * np.sin(launch_angle + 2*np.pi*turns) / (1 - np.sin(launch_angle + 2*np.pi*turns))
    # #track4 = ArcTrack(((Rlaunch+Rspiral)*np.cos(launch_angle),-Rlaunch), Rlaunch, np.pi+launch_angle, 3/2*np.pi, True, next_track=None)
    # track4 = ArcTrack(((Rlaunch + Rspiral), -Rlaunch), Rlaunch,
    #                   0, -2*np.pi/3, True, next_track=None)
    # track3 = FermatSpiralTrack(Rspiral, turns, np.pi-launch_angle, True, next_track=track4)
    # track2 = FermatSpiralTrack(Rspiral, turns, np.pi-launch_angle, False, next_track=track3)
    # #track1 = ArcTrack((-(Rlaunch+Rspiral)*np.cos(launch_angle), Rlaunch), Rlaunch, 3/2 * np.pi, 2*np.pi-launch_angle, False, next_track=track2)
    # track1 = ArcTrack((-(Rlaunch + Rspiral), Rlaunch), Rlaunch, 0,
    #                   final_angle, False, next_track=track2)

    arclengths = np.arange(0, track1.total_arclength(), 1e-5) #track1.arclength+track2.arclength*0.99, track1.arclength+track2.arclength*1.01, 1e-6)
    xs, ys, us, vs = track1.vertices_normals(arclengths)
    fig, axs = plt.subplots()
    axs.set_aspect('equal')
    axs.quiver(xs,ys,us,vs,scale_units='width', scale=100)
    xs,ys,_,_=track1.vertices_normals(np.array([0, track1.arclength, track1.arclength+track2.arclength, track1.arclength+track2.arclength+track3.arclength, track1.arclength+track2.arclength+track3.arclength+track4.arclength-1e-12]))
    axs.scatter(xs, ys, c='r', marker='.')
    A=40e-3
    axs.set_xlim(-A, A), axs.set_ylim(-A,A)
    plt.show()
