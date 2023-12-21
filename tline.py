# load .dat file with unit cell parameters
# construct polygon arrays for each relevant line
# break into fishbone unit cells

# construct spiral track polygon (Fermat spiral with diametrically opposed launch curves)
#     parametrized by arclength
# sample spiral track at arclengths given by x-values of artwork polygon fishbone cell starts
# compute spiral track tangents and normals at same sample arclengths

# offset from track sample points by fishbone unit cell rotated into (tangent, normal) coordinate system

# generate svg file 

import numpy as np
import matplotlib.pyplot as plt
import ezdxf

class FishboneUnitCell:
    def __init__(self, cell_length, fishbone_length, fishbone_height, line_width):
        self.line_width, self.cell_length, self.fishbone_height, self.fishbone_length = (
            line_width, cell_length, fishbone_height, fishbone_length)
    
    def vertices(self):
        xvals1=[0,
               self.cell_length/2 - self.fishbone_length/2, self.cell_length/2 - self.fishbone_length/2,
               self.cell_length/2 + self.fishbone_length/2, self.cell_length/2 + self.fishbone_length/2,
               self.cell_length]
        #yvals=[self.fishbone_length+self.line_width/2, self.fishbone_length+self.line_width/2, self.line_width/2, self.line_width/2]
        yvals1 = [self.line_width / 2,
                 self.line_width/2, self.fishbone_height + self.line_width / 2,
                 self.fishbone_height + self.line_width / 2, self.line_width / 2,
                 self.line_width / 2 ]
        #changed this to make fishbone unit symmetric

        return np.array((xvals1, yvals1))

    def vertices_gnd(self):
        xvals1=[0,
               self.cell_length/2 - self.fishbone_length/2, self.cell_length/2 - self.fishbone_length/2,
               self.cell_length/2 + self.fishbone_length/2, self.cell_length/2 + self.fishbone_length/2,
               self.cell_length]
        yvals1 = [self.line_width / 2,
                 self.line_width/2, self.fishbone_height + self.line_width / 2,
                 self.fishbone_height + self.line_width / 2, self.line_width / 2,
                 self.line_width / 2 ]
        xvals = [xvals1[-1], xvals1[0]]
        yvals = [self.fishbone_height + 2e-6, self.fishbone_height + 2e-6]

        return np.array((xvals, yvals))


# class FishboneUnitCellNegative:
#     def __init__(self, cell_length, fishbone_length, fishbone_height, line_width, gap_width):
#         self.line_width, self.cell_length, self.fishbone_height, self.fishbone_length, self.gap_width = (
#             line_width, cell_length, fishbone_height, fishbone_length, gap_width)
#
#     def vertices(self):
#         xvals = [0,
#                  0,
#                  0,
#                  self.cell_length / 2 - self.fishbone_length / 2,
#                  self.cell_length / 2 - self.fishbone_length / 2,
#                  self.cell_length / 2 + self.fishbone_length / 2,
#                  self.cell_length / 2 + self.fishbone_length / 2,
#                  self.cell_length,
#                  self.cell_length,
#                  0,
#                  self.cell_length,
#                  self.cell_length
#                  ]
#         # yvals=[self.fishbone_length+self.line_width/2, self.fishbone_length+self.line_width/2, self.line_width/2, self.line_width/2]
#         yvals = [self.line_width/2,
#                  self.line_width / 2 + self.fishbone_height + self.gap_width,
#                  self.line_width/2,
#                  self.line_width/2,
#                  self.line_width / 2 + self.fishbone_height,
#                  self.line_width / 2 + self.fishbone_height,
#                  self.line_width / 2,
#                  self.line_width / 2,
#                  self.line_width / 2 + self.fishbone_height + self.gap_width,
#                  self.line_width / 2 + self.fishbone_height + self.gap_width,
#                  self.line_width / 2 + self.fishbone_height + self.gap_width,
#                  self.line_width / 2
#                 ]
#
#         return np.array((xvals, yvals))





class FloquetUnitCell:
    def __init__(self):
        self.fishbones=[]
    def append_fishbones(self, fishbone_cell, n_fishbones=1):
        for n in range(n_fishbones):
            self.fishbones = np.append(self.fishbones, fishbone_cell)
            
    def vertices(self):
        v=self.fishbones[0].vertices()
        xstart=self.fishbones[0].cell_length
        for fishbone in self.fishbones[1:]:
            fishbone_vertices=fishbone.vertices()
            #print(fishbone_vertices)
            fishbone_vertices[0,:]+=xstart
            v=np.append(v, fishbone_vertices, axis = 1)
            xstart += fishbone.cell_length

        return v

    def vertices_gnd(self):
        v=self.fishbones[-1].vertices_gnd()
        xstart=(len(self.fishbones)-1)*self.fishbones[0].cell_length
        v[0,:] += xstart
        for fishbone in reversed(self.fishbones[0:-1]):
            fishbone_vertices=fishbone.vertices_gnd()
            #print(fishbone_vertices)
            # fishbone_vertices[0,:]-=xstart
            v=np.append(v, fishbone_vertices, axis = 1)
            # xstart += fishbone.cell_length
            print(xstart)
        return v


if __name__=='__main__':
    fishboneA, fishboneB = FishboneUnitCell(4e-6, 2e-6, 25e-6, 2e-6), FishboneUnitCell(4e-6, 2e-6, 100e-6, 2e-6)
    #fishboneA, fishboneB = FishboneUnitCellNegative(4e-6,2e-6, 25e-6, 2e-6, 2e-6), FishboneUnitCellNegative(4e-6, 2e-6, 100e-6, 2e-6, 2e-6)
    floquet = FloquetUnitCell()
    floquet.append_fishbones(fishboneA, 1)
    # floquet.append_fishbones(fishboneB, 1)
    # floquet.append_fishbones(fishboneA, 3)
    xs, ys = floquet.vertices()
    xgnd, ygnd = floquet.vertices_gnd()

    # #add ground boundary
    # numbones = 2
    # bone_vertices = 6
    # bones = [fishboneA, fishboneB]
    # # xs = np.append(xs, xs[-1])
    # # ys = np.append(ys, bones[-1].fishbone_height + 2e-6)
    # #
    # # for i in range(numbones):
    # #     xs = np.append(xs, [xs[-1-(i+1)*bone_vertices], xs[-1-(i+1)*bone_vertices]])#, xs[i*numbones*bone_vertices]])
    # #     ys = np.append(ys, [bones[-1].fishbone_height + 2e-6, bones[-2].fishbone_height + 2e-6])#, bones[-1-i].fishbone_height + 2e-6, ys[0]])


    # print(xs)
    mir_ys = -ys
    mir_ygnd = -ygnd

    fig,ax=plt.subplots()
    ax.scatter(xs,ys,c = np.arange(len(xs)),marker='.', cmap = 'coolwarm')
    ax.scatter(xgnd, ygnd, c = np.arange(len(xgnd)), marker='.', cmap = 'coolwarm')
    ax.set_ylim(0,None)
    plt.show()
    fig,ax = plt.subplots()
    ax.plot(xs, ys, marker = '.')
    ax.plot(xgnd, ygnd, marker = '.')
    plt.show()

    merged_list = [(xs[i], ys[i]) for i in range(0, len(xs))]
    merged_list_gnd = [(xgnd[i], ygnd[i]) for i in range(0, len(xgnd))]
    merge_list_mir = [(xs[i], mir_ys[i]) for i in range(0, len(xs))]
    merged_list_mir_gnd = [(xgnd[i], mir_ygnd[i]) for i in range(0, len(xgnd))]
    print(merged_list_gnd)


    doc = ezdxf.new()
    msp = doc.modelspace()
    #polyline = msp.add_lwpolyline(merged_list + merge_list_mir[::-1], close=True)
    polyline = msp.add_polyline2d(merged_list + merged_list_gnd, close=True)
    polyline_mir = msp.add_polyline2d(merge_list_mir + merged_list_mir_gnd, close = True)

    # Save the DXF file
    doc.saveas("unitcell.dxf")

    # print(merged_list)

