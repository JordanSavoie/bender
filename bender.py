import track, tline
import numpy as np, matplotlib.pyplot as plt
import ezdxf
from ezdxf import units


class Bender:  # uh oh
    def __init__(self, track, tline):
        self.track, self.tline = track, tline

    def bend(self):
        pass


if __name__ == '__main__':
    fishboneA, fishboneB = tline.FishboneUnitCell(40e-6, 2e-6, 25e-6, 2e-6), tline.FishboneUnitCell(4e-6, 2e-6, 100e-6,
                                                                                                   2e-6)

    # fishboneA, fishboneB = tline.FishboneUnitCell(5e-6, 10e-6, 1e-6), tline.FishboneUnitCell(7.5e-6, 20e-6, 1e-6)
    floquet = tline.FloquetUnitCell()
    floquet.append_fishbones(fishboneA, 3)
    # floquet.append_fishbones(fishboneB, 1)
    # floquet.append_fishbones(fishboneA, 10)
    xs, ys = floquet.vertices()
    unit_cell_vertices_count = len(xs)
    xgnd, ygnd = floquet.vertices_gnd()
    unit_cell_vertices_count_gnd = len(xgnd)
    floquet_cell_length = xs[-1]
    floquet_vertices = floquet.vertices()
    floquet_vertices_gnd = floquet.vertices_gnd()


    Rspiral = 10e-3
    launch_angle = 0
    turns = 0.01#1.950620200000362
    final_angle = turns * 2 * np.pi - (2 * np.pi) * (turns // 1)

    # arbitrary turn version


    Rlaunch = 10e-3

    # arc1 = track.ArcTrack((-Rspiral * np.cos(final_angle) - Rlaunch * np.cos(final_angle),
    #                  -Rspiral * np.sin(final_angle) - Rlaunch * np.sin(final_angle)),
    #                 Rlaunch, -final_angle % np.pi, 1 / 2 * np.pi, True, next_track=None)
    fermat1 = track.FermatSpiralTrack(Rspiral, turns, 0, True, next_track=None)
    fermat2 = track.FermatSpiralTrack(Rspiral, turns, 0, False, next_track=fermat1)
    # arc2 = track.ArcTrack(
    #     ((Rspiral + Rlaunch) * np.cos(final_angle), Rspiral * np.sin(final_angle) + Rlaunch * np.sin(final_angle)),
    #     Rlaunch,
    #     np.pi / 2, final_angle % np.pi, False, next_track=fermat2)



    # track1 = arc2
    # track2 = fermat2
    # track3 = fermat1
    # track4 = arc1

    track1 = fermat2


    # Rspiral = 20e-3
    # launch_angle = np.pi/6
    # Rlaunch = Rspiral * np.sin(launch_angle) / (1-np.sin(launch_angle))
    # track4 = track.ArcTrack(((Rlaunch+Rspiral)*np.cos(launch_angle),-Rlaunch), Rlaunch, np.pi+launch_angle, 3/2*np.pi, True, next_track=None)
    # track3 = track.FermatSpiralTrack(Rspiral, 6, np.pi-launch_angle, True, next_track=track4)
    # track2 = track.FermatSpiralTrack(Rspiral, 6, np.pi-launch_angle, False, next_track=track3)
    # track1 = track.ArcTrack((-(Rlaunch+Rspiral)*np.cos(launch_angle), Rlaunch), Rlaunch, 3/2 * np.pi, 2*np.pi-launch_angle, False, next_track=track2)
    #


    arclength = track1.total_arclength()
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
        offset = np.array([[0.0], [0.0]])
        for i in range(n_floquet_cells):
            offset += np.array([[floquet_cell_length], [0.0]])
            tline_vertices_gnd = np.append(tline_vertices_gnd, floquet_vertices_gnd + offset, axis=1)
        np.save('tline_vertices_gnd', tline_vertices_gnd)

    tline_vertices = np.load('tline_vertices.npy')
    tline_vertices_gnd = np.load('tline_vertices_gnd.npy')

    arclengths, heights = tline_vertices
    arclengths_gnd, heights_gnd = tline_vertices_gnd
    xs, ys, us, vs = track1.vertices_normals(arclengths)
    xgnd, ygnd, ugnd, vgnd = track1.vertices_normals(arclengths_gnd)
    print(arclengths_gnd)


    bent_xs = xs + us* heights #for some reason heights is too long
    bent_xs_mirror = xs - us* heights
    bent_ys = ys + vs* heights
    bent_ys_mirror = ys - vs*heights

    bent_xs_gnd = xgnd + ugnd * heights_gnd
    bent_xs_gnd_mirror = xgnd - ugnd * heights_gnd
    bent_ys_gnd = ygnd + vgnd * heights_gnd
    bent_ys_gnd_mirror = ygnd - vgnd * heights_gnd

    # xs, ys = arclengths, heights
    # xgnd, ygnd = arclengths_gnd, heights_gnd
    #
    # bent_xs = xs
    # bent_xs_mirror = xs
    # bent_ys = ys
    # bent_ys_mirror = ys
    #
    # bent_xs_gnd = xgnd
    # bent_xs_gnd_mirror = xgnd
    # bent_ys_gnd = ygnd
    # bent_ys_gnd_mirror = ygnd

    doc = ezdxf.new()
    msp = doc.modelspace()
    doc.units = units.M

    print(arclengths_gnd)


    for k in range(n_floquet_cells-1):
        j = k+1
        merged_list = [(bent_xs[i], bent_ys[i]) for i in range((j-1) * unit_cell_vertices_count, j * unit_cell_vertices_count)]
        merged_list_gnd = [(bent_xs_gnd[i], bent_ys_gnd[i]) for i in range((j-1) * unit_cell_vertices_count_gnd, j * unit_cell_vertices_count_gnd)]

        merged_list_mir = [(bent_xs_mirror[i], bent_ys_mirror[i]) for i in range((j-1) * unit_cell_vertices_count, j * unit_cell_vertices_count)]
        merged_list_gnd_mir = [(bent_xs_gnd_mirror[i], bent_ys_gnd_mirror[i]) for i in
                           range((j - 1) * unit_cell_vertices_count_gnd, j * unit_cell_vertices_count_gnd)]

        # merged_list_gnd[0] = (merged_list[-1][0], merged_list_gnd[0][1])
        # merged_list_gnd_mir[0] = (merged_list_mir[-1][0], merged_list_gnd_mir[0][1])
        msp.add_lwpolyline(merged_list, close=False)
        msp.add_lwpolyline(merged_list_gnd, close = False)
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
    # axs.set_xlim(-A, A), axs.set_ylim(-A, A)
    xx = np.append(bent_xs, bent_xs_gnd)
    yy = np.append(bent_ys, bent_ys_gnd)
    # plt.plot(xx, yy, c='b', marker='.')
    # plt.plot(xs, ys, c='r', marker='.')
    plt.plot( np.arange(len(arclengths_gnd)) , arclengths_gnd)
    plt.show()
    # axs.quiver(xs, ys, us, vs, scale_units='width', scale=100)
    # plt.plot(bent_xs_gnd, bent_ys_gnd, c='g', marker='.')

    # xs, ys, _, _ = track1.vertices_normals(np.array([0, track1.arclength, track1.arclength + track2.arclength,
    #                                                          track1.arclength + track2.arclength + track3.arclength,
    #                                                          track1.arclength + track2.arclength + track3.arclength + track4.arclength - 1e-12]))

    # xs, ys, _, _ = fermat1.vertices_normals(np.array([0, fermat1.arclength, fermat1.arclength + fermat1.arclength]))
    #
    # axs.scatter(xs, ys, c='g', marker='*')
    # plt.show()
