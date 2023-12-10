import track, tline
import numpy as np, matplotlib.pyplot as plt

class Bender: #uh oh
    def __init__(self, track, tline):
        self.track, self.tline = track, tline

    def bend(self):
        pass

fishboneA, fishboneB = tline.FishboneUnitCell(5e-6, 10e-6, 1e-6), tline.FishboneUnitCell(7.5e-6, 20e-6, 1e-6)
floquet = tline.FloquetUnitCell()
floquet.append_fishbones(fishboneA, 5)
floquet.append_fishbones(fishboneB, 3)
floquet.append_fishbones(fishboneA, 5)
floquet_vertices = floquet.vertices()
floquet_cell_length = floquet_vertices[0,-1]

Rspiral = 20e-3
launch_angle = np.pi/6
Rlaunch = Rspiral * np.sin(launch_angle) / (1-np.sin(launch_angle))
track4 = track.ArcTrack(((Rlaunch+Rspiral)*np.cos(launch_angle),-Rlaunch), Rlaunch, np.pi+launch_angle, 3/2*np.pi, True, next_track=None)
track3 = track.FermatSpiralTrack(Rspiral, 6, np.pi-launch_angle, True, next_track=track4)
track2 = track.FermatSpiralTrack(Rspiral, 6, np.pi-launch_angle, False, next_track=track3)
track1 = track.ArcTrack((-(Rlaunch+Rspiral)*np.cos(launch_angle), Rlaunch), Rlaunch, 3/2 * np.pi, 2*np.pi-launch_angle, False, next_track=track2)

arclength = track1.total_arclength()
n_floquet_cells = int(np.floor(arclength / floquet_cell_length))

try:
    with open('tline_vertices.npy') as f:
        pass
except FileNotFoundError:
    tline_vertices = np.array([[],[]])
    offset = np.array([[0.0],[0.0]])
    for i in range(n_floquet_cells):
        offset += np.array([[floquet_cell_length], [0.0]])
        tline_vertices = np.append(tline_vertices, floquet_vertices + offset, axis=1)
    np.save('tline_vertices', tline_vertices)

tline_vertices = np.load('tline_vertices.npy')

arclengths, heights = tline_vertices
print(arclengths, heights)
heights = heights[:-1]
xs,ys,us,vs = track1.vertices_normals(arclengths)

bent_xs = xs + us * heights
bent_ys = ys + vs * heights

fig, axs = plt.subplots()
axs.set_aspect('equal')
A=40e-3
axs.set_xlim(-A, A), axs.set_ylim(-A,A)
plt.plot(bent_xs, bent_ys, c='b', marker='.')

xs,ys,_,_=track1.vertices_normals(np.array([0, track1.arclength, track1.arclength+track2.arclength, track1.arclength+track2.arclength+track3.arclength, track1.arclength+track2.arclength+track3.arclength+track4.arclength-1e-12]))
axs.scatter(xs, ys, c='r', marker='.')
plt.show()
