import track, tline
import numpy as np, matplotlib.pyplot as plt


class Bender:  # uh oh
    def __init__(self, track, tline):
        self.track, self.tline = track, tline

    def bend(self):
        pass


if __name__ == '__main__':
    fishboneA, fishboneB = tline.FishboneUnitCell(4e-6, 2e-6, 25e-6, 2e-6), tline.FishboneUnitCell(4e-6, 2e-6, 100e-6,
                                                                                                   2e-6)

    # fishboneA, fishboneB = tline.FishboneUnitCell(5e-6, 10e-6, 1e-6), tline.FishboneUnitCell(7.5e-6, 20e-6, 1e-6)
    floquet = tline.FloquetUnitCell()
    floquet.append_fishbones(fishboneA, 5)
    floquet.append_fishbones(fishboneB, 3)
    floquet.append_fishbones(fishboneA, 5)
    floquet_vertices = floquet.vertices()
    # print("floquet vertices")
    # print(floquet_vertices)
    # print("end")
    floquet_cell_length = floquet_vertices[0, -1]
    #print(floquet_cell_length)

    # arbitrary turn version

    Rspiral = 20e-3
    launch_angle = 0
    turns = 1.95
    final_angle = turns * 2 * np.pi - (2 * np.pi) * (turns //1)
    Rlaunch = 10e-3

    arc1 = track.ArcTrack((-Rspiral * np.cos(final_angle) - Rlaunch * np.cos(final_angle),
                     -Rspiral * np.sin(final_angle) - Rlaunch * np.sin(final_angle)),
                    Rlaunch, 3 / 2 * np.pi, final_angle, False, next_track=None)
    fermat1 = track.FermatSpiralTrack(Rspiral, turns, 0, True, next_track=arc1)
    fermat2 = track.FermatSpiralTrack(Rspiral, turns, 0, False, next_track=fermat1)
    arc2 = track.ArcTrack(
        ((Rspiral + Rlaunch) * np.cos(final_angle), Rspiral * np.sin(final_angle) + Rlaunch * np.sin(final_angle)),
        Rlaunch,
        np.pi / 2, final_angle % np.pi, False, next_track=fermat2)



    track1 = arc2
    track2 = fermat2
    track3 = fermat1
    track4 = arc1

    # Rspiral = 20e-3
    # launch_angle = np.pi/6
    # Rlaunch = Rspiral * np.sin(launch_angle) / (1-np.sin(launch_angle))
    # track4 = track.ArcTrack(((Rlaunch+Rspiral)*np.cos(launch_angle),-Rlaunch), Rlaunch, np.pi+launch_angle, 3/2*np.pi, True, next_track=None)
    # track3 = track.FermatSpiralTrack(Rspiral, 6, np.pi-launch_angle, True, next_track=track4)
    # track2 = track.FermatSpiralTrack(Rspiral, 6, np.pi-launch_angle, False, next_track=track3)
    # track1 = track.ArcTrack((-(Rlaunch+Rspiral)*np.cos(launch_angle), Rlaunch), Rlaunch, 3/2 * np.pi, 2*np.pi-launch_angle, False, next_track=track2)
    #


    arclength = track1.total_arclength()
    n_floquet_cells = int(np.floor(arclength / floquet_cell_length))

    try:
        with open('tline_vertices.npy') as f:
            pass
    except FileNotFoundError:
        tline_vertices = np.array([[], []])
        offset = np.array([[0.0], [0.0]])
        for i in range(n_floquet_cells):
            offset += np.array([[floquet_cell_length], [0.0]])
            tline_vertices = np.append(tline_vertices, floquet_vertices + offset, axis=1)
        np.save('tline_vertices', tline_vertices)

    tline_vertices = np.load('tline_vertices.npy')

    arclengths, heights = tline_vertices
    arclengths = arclengths
    heights = heights
    xs, ys, us, vs = track1.vertices_normals(arclengths)


    bent_xs = xs + us * heights[:len(us)] #for some reason heights is too long
    bent_xs_mirror = xs - us * heights[:len(us)]
    bent_ys = ys + vs * heights[:len(vs)]
    bent_ys_mirror = ys - vs*heights[:len(vs)]

    fig, axs = plt.subplots()
    axs.set_aspect('equal')
    A = 40e-3
    axs.set_xlim(-A, A), axs.set_ylim(-A, A)
    plt.plot(bent_xs, bent_ys, c='b', marker='.')
    plt.plot(bent_xs_mirror, bent_ys_mirror, c='r', marker='.')

    xs, ys, _, _ = track1.vertices_normals(np.array([0, track1.arclength, track1.arclength + track2.arclength,
                                                     track1.arclength + track2.arclength + track3.arclength,
                                                     track1.arclength + track2.arclength + track3.arclength + track4.arclength - 1e-12]))
    axs.scatter(xs, ys, c='g', marker='*')
    plt.show()
