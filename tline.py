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
        xvals=[0,
               self.cell_length/2 - self.fishbone_length/2, self.cell_length/2 - self.fishbone_length/2,
               self.cell_length/2 + self.fishbone_length/2, self.cell_length/2 + self.fishbone_length/2,
               self.cell_length]
        #yvals=[self.fishbone_length+self.line_width/2, self.fishbone_length+self.line_width/2, self.line_width/2, self.line_width/2]
        yvals = [self.line_width / 2,
                 self.line_width/2, self.fishbone_height + self.line_width / 2,
                 self.fishbone_height + self.line_width / 2, self.line_width / 2,
                 self.line_width / 2 ]
        #changed this to make fishbone unit symmetric
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
            self.fishbones.append(fishbone_cell)
            
    def vertices(self):
        v=self.fishbones[0].vertices()
        xstart=self.fishbones[0].cell_length
        for fishbone in self.fishbones[1:]:
            fishbone_vertices=fishbone.vertices()
            #print(fishbone_vertices)
            fishbone_vertices[0,:]+=xstart
            v=np.append(v, fishbone_vertices, axis=1)
            xstart += fishbone.cell_length

        return v


if __name__=='__main__':
    fishboneA, fishboneB = FishboneUnitCell(4e-6, 2e-6, 25e-6, 2e-6), FishboneUnitCell(4e-6, 2e-6, 100e-6, 2e-6)
    #fishboneA, fishboneB = FishboneUnitCellNegative(4e-6,2e-6, 25e-6, 2e-6, 2e-6), FishboneUnitCellNegative(4e-6, 2e-6, 100e-6, 2e-6, 2e-6)
    floquet = FloquetUnitCell()
    floquet.append_fishbones(fishboneA, 1)
    floquet.append_fishbones(fishboneB, 1)
    floquet.append_fishbones(fishboneA, 1)
    xs, ys = floquet.vertices()

    #add ground boundary
    xs = np.append(xs, [xs[-1], xs[0], xs[0]])
    ys = np.append(ys, [fishboneB.fishbone_height + 2e-6, fishboneB.fishbone_height + 2e-6, ys[0]])
    print(xs)
    mir_ys = -ys

    fig,ax=plt.subplots()
    ax.plot(xs,ys,marker='.')
    ax.set_ylim(0,None)
    plt.show()

    merged_list = [(xs[i], ys[i]) for i in range(0, len(xs))]
    merge_list_mir = [(xs[i], mir_ys[i]) for i in range(0, len(xs))]
    print(merged_list)


    doc = ezdxf.new()
    msp = doc.modelspace()
    #polyline = msp.add_lwpolyline(merged_list + merge_list_mir[::-1], close=True)
    polyline = msp.add_polyline2d(merged_list + merge_list_mir[::-1], close=True)

    # Save the DXF file
    doc.saveas("unitcell.dxf")

    print(merged_list)

