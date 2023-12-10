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

class FishboneUnitCell:
    def __init__(self, cell_length, fishbone_length, line_width):
        self.line_width,self.cell_length,self.fishbone_length=line_width,cell_length,fishbone_length
    
    def vertices(self):
        xvals=[0, self.line_width, self.line_width, self.cell_length]
        yvals=[self.fishbone_length+self.line_width/2, self.fishbone_length+self.line_width/2, self.line_width/2, self.line_width/2]
        return np.array((xvals, yvals))

class FloquetUnitCell:
    def __init__(self):
        self.fishbones=[]

    def append_fishbones(self,fishbone_cell, n_fishbones=1):
        for n in range(n_fishbones):
            self.fishbones.append(fishbone_cell)
            
    def vertices(self):
        v=self.fishbones[0].vertices()
        xstart=self.fishbones[0].cell_length
        for fishbone in self.fishbones[1:]:
            fishbone_vertices=fishbone.vertices()
            fishbone_vertices[0,:]+=xstart
            v=np.append(v, fishbone_vertices, axis=1)
            xstart += fishbone.cell_length
        return v

if __name__=='__main__':
    fishboneA, fishboneB = FishboneUnitCell(5e-6, 10e-6, 1e-6), FishboneUnitCell(7.5e-6, 20e-6, 1e-6)
    floquet = FloquetUnitCell()
    floquet.append_fishbones(fishboneA, 5)
    floquet.append_fishbones(fishboneB, 3)
    floquet.append_fishbones(fishboneA, 5)
    xs, ys = floquet.vertices()
    fig,ax=plt.subplots()
    ax.plot(xs,ys,marker='.')
    ax.set_ylim(0,None)
    plt.show()

