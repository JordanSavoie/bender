import ezdxf, numpy as np

points = np.array([[0,0], [1,0], [0,1], [0,0]])

doc = ezdxf.new()
msp = doc.modelspace()
msp.add_lwpolyline(points)

doc.saveas('test.dxf')
