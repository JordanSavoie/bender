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
    def __init__(self, cell_length, fishbone_length, fishbone_height, line_width, gnd_spacing=2e-6):
        self.line_width, self.cell_length, self.fishbone_height, self.fishbone_length, self.gnd_spacing = \
            line_width, cell_length, fishbone_height, fishbone_length, gnd_spacing

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
        gnd_height = self.line_width/2+self.fishbone_height+self.gnd_spacing
        return np.array(((self.cell_length, 0), (gnd_height, gnd_height)))

class FloquetUnitCell:
    def __init__(self):
        self.fishbones=[]
    
    def append_fishbones(self, fishbone_cell, n_fishbones=1):
        for n in range(n_fishbones):
            self.fishbones = np.append(self.fishbones, fishbone_cell)
            
    def vertices(self):
        v=np.array([[],[]])
        xstart=0
        for fishbone in self.fishbones:
            fishbone_vertices=fishbone.vertices()
            fishbone_vertices[0,:]+=xstart
            v=np.append(v, fishbone_vertices, axis=1)
            xstart += fishbone.cell_length

        for fishbone in reversed(self.fishbones):
            xstart -= fishbone.cell_length
            gnd_vertices = fishbone.vertices_gnd()
            gnd_vertices[0,:]+=xstart
            v=np.append(v, gnd_vertices, axis=1)

        return v

    def cell_length(self):
        return sum([fishbone.cell_length for fishbone in self.fishbones])

    def cell_min_spacing(self):
        return max(np.concatenate([fishbone.vertices()[1,:] for fishbone in self.fishbones]))

if __name__=='__main__':
    pass
    # fishboneA, fishboneB = FishboneUnitCell(4e-6, 2e-6, 25e-6, 2e-6), FishboneUnitCell(4e-6, 2e-6, 100e-6, 2e-6)
    # #fishboneA, fishboneB = FishboneUnitCellNegative(4e-6,2e-6, 25e-6, 2e-6, 2e-6), FishboneUnitCellNegative(4e-6, 2e-6, 100e-6, 2e-6, 2e-6)
    # floquet = FloquetUnitCell()
    # floquet.append_fishbones(fishboneA, 2)
    # floquet.append_fishbones(fishboneB, 1)
    # xs, ys = floquet.vertices()

    # mir_ys = -ys
    #
    # fig,ax=plt.subplots()
    # ax.scatter(xs,ys)
    # ax.scatter(xs, mir_ys)
    # for i in range(len(xs)):
    #     ax.annotate(str(i), (xs[i],ys[i]))
    # plt.show()

    # merged_list = [(xs[i], ys[i]) for i in range(0, len(xs))]
    # merged_list_gnd = [(xgnd[i], ygnd[i]) for i in range(0, len(xgnd))]
    # merge_list_mir = [(xs[i], mir_ys[i]) for i in range(0, len(xs))]
    # merged_list_mir_gnd = [(xgnd[i], mir_ygnd[i]) for i in range(0, len(xgnd))]
    # print(merged_list_gnd)


    # doc = ezdxf.new()
    # msp = doc.modelspace()
    # #polyline = msp.add_lwpolyline(merged_list + merge_list_mir[::-1], close=True)
    # polyline = msp.add_polyline2d(merged_list + merged_list_gnd, close=True)
    # polyline_mir = msp.add_polyline2d(merge_list_mir + merged_list_mir_gnd, close = True)
    #
    # # Save the DXF file
    # doc.saveas("unitcell.dxf")

    # print(merged_list)

