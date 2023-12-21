import track, tline
import numpy as np, matplotlib.pyplot as plt
import ezdxf
from ezdxf import units
import math


class Bender:  # uh oh
    def __init__(self, track, tline):
        self.track, self.tline = track, tline

    def bend(self):
        pass


if __name__ == '__main__':
    fishboneA, fishboneB = tline.FishboneUnitCell(4e-6, 2e-6, 25e-6, 2e-6), tline.FishboneUnitCell(4e-6, 2e-6, 100e-6,
                                                                                                   2e-6)

    floquet = tline.FloquetUnitCell()
    floquet.append_fishbones(fishboneA, 5)
    floquet.append_fishbones(fishboneB, 10)
    floquet.append_fishbones(fishboneA, 5)
    xs, ys = floquet.vertices()
    unit_cell_vertices_count = len(xs)
    xgnd, ygnd = floquet.vertices_gnd()
    unit_cell_vertices_count_gnd = len(xgnd)
    floquet_cell_length = xs[-1]
    floquet_vertices = floquet.vertices()
    floquet_vertices_gnd = floquet.vertices_gnd()


    Rspiral = 10e-3
    launch_angle = 0
    turns = 1.95#1.950620200000362
    final_angle = 2 * np.pi * turns

    # arbitrary turn version


    Rlaunch = 10e-3

    # arc1 = track.ArcTrack((-Rspiral * np.cos(final_angle) - Rlaunch * np.cos(final_angle),
    #                  -Rspiral * np.sin(final_angle) - Rlaunch * np.sin(final_angle)),
    #                 Rlaunch, -final_angle % np.pi, 1 / 2 * np.pi, True, next_track=None)
    # fermat1 = track.FermatSpiralTrack(Rspiral, turns, 0, True, next_track=None)
    # fermat2 = track.FermatSpiralTrack(Rspiral, turns, 0, False, next_track=fermat1)
    # arc2 = track.ArcTrack(
    #     ((Rspiral + Rlaunch) * np.cos(final_angle), Rspiral * np.sin(final_angle) + Rlaunch * np.sin(final_angle)),
    #     Rlaunch,
    #     np.pi / 2, final_angle % np.pi, False, next_track=fermat2)

    fermat1 = track.FermatSpiralTrack(Rspiral, turns, 0, True)
    fermat2 = track.FermatSpiralTrack(Rspiral, turns, 0, False)
    arc2 = fermat2.construct_suitable_ArcTrack(Rlaunch)
    arc1 = fermat1.construct_suitable_ArcTrack(Rlaunch)
    x, y = arc1.center


    track1 = arc2
    track2 = fermat2
    track3 = fermat1
    track4 = arc1

    track1.next_track, track2.next_track, track3.next_track = track2, track3, track4
    initial_track = track1

    # Rspiral = 20e-3
    # launch_angle = np.pi/6
    # Rlaunch = Rspiral * np.sin(launch_angle) / (1-np.sin(launch_angle))
    # track4 = track.ArcTrack(((Rlaunch+Rspiral)*np.cos(launch_angle),-Rlaunch), Rlaunch, np.pi+launch_angle, 3/2*np.pi, True, next_track=None)
    # track3 = track.FermatSpiralTrack(Rspiral, 6, np.pi-launch_angle, True, next_track=track4)
    # track2 = track.FermatSpiralTrack(Rspiral, 6, np.pi-launch_angle, False, next_track=track3)
    # track1 = track.ArcTrack((-(Rlaunch+Rspiral)*np.cos(launch_angle), Rlaunch), Rlaunch, 3/2 * np.pi, 2*np.pi-launch_angle, False, next_track=track2)
    #


    arclength = initial_track.total_arclength()
    n_floquet_cells = int(np.floor(arclength / floquet_cell_length)-1)
    print(n_floquet_cells)
    print((arclength / floquet_cell_length))


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

    try:
        with open('tline_vertices_gnd.npy') as f:
            pass
    except FileNotFoundError:
        tline_vertices_gnd = np.array([[], []])
        offset = np.array([[(n_floquet_cells)*floquet_cell_length], [0.0]])
        for i in range(n_floquet_cells):
            tline_vertices_gnd = np.append(tline_vertices_gnd, floquet_vertices_gnd + offset, axis=1)
            offset -= np.array([[floquet_cell_length], [0.0]])
        np.save('tline_vertices_gnd', np.flip(tline_vertices_gnd, axis=1))

    tline_vertices = np.load('tline_vertices.npy')
    tline_vertices_gnd = np.load('tline_vertices_gnd.npy')

    arclengths, heights = tline_vertices
    arclengths_gnd, heights_gnd = tline_vertices_gnd
    xs, ys, us, vs = track1.vertices_normals(arclengths)
    xgnd, ygnd, ugnd, vgnd = track1.vertices_normals(arclengths_gnd)


    bent_xs = xs + us* heights #for some reason heights is too long
    bent_xs_mirror = xs - us* heights
    bent_ys = ys + vs* heights
    bent_ys_mirror = ys - vs*heights

    bent_xs_gnd = xgnd + ugnd * heights_gnd
    bent_xs_gnd_mirror = xgnd - ugnd * heights_gnd
    bent_ys_gnd = ygnd + vgnd * heights_gnd
    bent_ys_gnd_mirror = ygnd - vgnd * heights_gnd



    doc = ezdxf.new()
    msp = doc.modelspace()
    doc.units = units.M


    # bent_xs = list(dict.fromkeys(bent_xs))
    # bent_ys = list(dict.fromkeys(bent_ys))
    # bent_xs_gnd = list(dict.fromkeys(bent_xs_gnd))
    # bent_ys_gnd = list(dict.fromkeys(bent_ys_gnd))
    # bent_xs_mirror = list(dict.fromkeys(bent_xs_mirror))
    # bent_ys_mirror = list(dict.fromkeys(bent_ys_mirror))
    # bent_xs_gnd_mirror = list(dict.fromkeys(bent_xs_gnd_mirror))
    # bent_ys_gnd_mirror = list(dict.fromkeys(bent_ys_gnd_mirror))



    for k in range(n_floquet_cells-1):
        j = (k+1)
        merged_list = [(bent_xs[i], bent_ys[i]) for i in range((j-1) * unit_cell_vertices_count, j * unit_cell_vertices_count)]
        merged_list_gnd = [(bent_xs_gnd[i], bent_ys_gnd[i]) for i in range((j-1) * unit_cell_vertices_count_gnd, j * unit_cell_vertices_count_gnd)]

        merged_list_mir = [(bent_xs_mirror[i], bent_ys_mirror[i]) for i in range((j-1) * unit_cell_vertices_count, j * unit_cell_vertices_count)]
        merged_list_gnd_mir = [(bent_xs_gnd_mirror[i], bent_ys_gnd_mirror[i]) for i in
                           range((j - 1) * unit_cell_vertices_count_gnd, j * unit_cell_vertices_count_gnd)]

        merged_list_gnd_mir.reverse()
        merged_list_gnd.reverse()

        # if math.dist(merged_list[-1], merged_list_gnd_mir[0]) > math.dist(merged_list[-1], merged_list_gnd_mir[1]):
        #     print("Error")
        #     print(j)
        # merged_list_gnd_mir.insert(0, merged_list_mir[-1])
        # merged_list_gnd_mir.append(merged_list_mir[0])
        # print(merged_list_gnd_mir)
        msp.add_lwpolyline(merged_list + merged_list_gnd, close=True)
        # msp.add_lwpolyline(merged_list_gnd, close = False)
        msp.add_lwpolyline(merged_list_mir + merged_list_gnd_mir, close=True)
        # print(merged_list_gnd)
        # print(merged_list[-1][0])
        # print(merged_list_gnd[0][0])
        # print(merged_list_gnd[-1][0])
        # print("***")




    # Save the DXF file
    doc.saveas("justspiral.dxf")



    fig, axs = plt.subplots()
    axs.set_aspect('equal')
    A = 20e-3
    axs.set_xlim(-A, A), axs.set_ylim(-A, A)
    xx = np.append(bent_xs, bent_xs_gnd)
    yy = np.append(bent_ys, bent_ys_gnd)
    plt.plot(xx, yy, c='b', marker='.')
    plt.plot(xs, ys, c='r', marker='.')
    # plt.scatter( np.arange(len(arclengths_gnd)) , arclengths_gnd)
    plt.show()
    # axs.quiver(xs, ys, us, vs, scale_units='width', scale=100)
    # plt.plot(arclengths_gnd, heights_gnd, c='g', marker='.')
    # for i in range(len(arclengths_gnd)):
    #     # plt.text(bent_xs_gnd[i], bent_ys_gnd[i]+(0.1*bent_ys_gnd[i]*(i%2)), str(i))
    #     # plt.text(xgnd[i], ygnd[i] + (1e-5* (i % 2)), str(i))
    #     plt.text(arclengths_gnd[i], heights_gnd[i] + (1e-9 * (i % 2)), str(i))
    # axs.scatter(arclengths_gnd[len(arclengths_gnd)//2], heights_gnd[len(heights_gnd)//2], c='r', marker='*')
    # plt.show()

    # xs, ys, _, _ = track1.vertices_normals(np.array([0, track1.arclength, track1.arclength + track2.arclength,
    #                                                          track1.arclength + track2.arclength + track3.arclength,
    #                                                          track1.arclength + track2.arclength + track3.arclength + track4.arclength - 1e-12]))

    # xs, ys, _, _ = fermat1.vertices_normals(np.array([0, fermat1.arclength, fermat1.arclength + fermat1.arclength]))
    #
    # axs.scatter(xs, ys, c='g', marker='*')
    # plt.show()
