# import ezdxf, numpy as np
#
# points = np.array([[0,0], [1,0], [0,1], [0,0]])
#
# doc = ezdxf.new()
# msp = doc.modelspace()
# hatch = msp.add_hatch(color=1)
# hatch.paths.add_polyline_path(points, is_closed=True)
#
# doc.saveas('test.dxf')

import ezdxf, numpy as np
import ezdxf

def create_filled_square(file_name):
    # Create a new DXF document
    doc = ezdxf.new()

    # Add a new model space
    msp = doc.modelspace()

    # Define the four corners of the square
    corners = [(0, 0), (5, 0), (5, 5), (0, 5), (0,0)]  # Change the coordinates as needed

    # Create a polyline for the square
    square_polyline = msp.add_lwpolyline(corners, close=True)

    # Save the DXF file
    doc.saveas(file_name)

# Specify the file name for the new DXF file
file_name = "filled_square.dxf"

# Create the filled square and save it as "filled_square.dxf"
create_filled_square(file_name)


