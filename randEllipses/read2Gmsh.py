
# ---------------------------------------------------------------------------
# The Python scripts read the txt file containing the geo infor of inclusions
# and then write to Gmsh
# 
# by Mingzhao Zhuo @ TU Delft, Oct. 2019
# ---------------------------------------------------------------------------


import numpy as np
import pygmsh
import meshio
import os

from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon

from shapely.geometry.point import Point
from shapely.geometry import LineString
from shapely import affinity

class Ellipse:
    """ Create a new Ellipse, with center, semi-axis a and b, angle """
    """ Ellipse include the five numbers and a affinity object """
    # semi major and minor axes

    def __init__(self, x=0, y=0, a=2, b=1, ang=0):
        """ Create a new center point at x, y """
        self.center = (x, y)
        self.lengths = (a, b)
        self.angle = ang
        # self.eps = create_ellipse(self.center, self.lengths, self.angle)

    def create_ellipse(self, ratio):
        """
        create a shapely ellipse. adapted from
        https://gis.stackexchange.com/a/243462
            https://gis.stackexchange.com/questions/243459/drawing-ellipse-with-shapely/243462#243462
        """
        circ = Point(self.center).buffer(1.0)
        ell = affinity.scale(circ, float(
            self.lengths[0]*ratio), float(self.lengths[1]*ratio))
        ellr = affinity.rotate(ell, self.angle)
        return ellr


allEllipses = []

#  store a and b
radius = [] 

# see this
# https://stackoverflow.com/questions/2081836/reading-specific-lines-only
# 
filename = './RVE25.txt'

with open(filename, 'r') as f:

#   Read and print the entire file line by line
    lines = f.readlines()
    lineSec = lines[1]
    spl = [x.strip() for x in lineSec.split(',')]
    Rx = float(spl[0])
    Ry = float(spl[1])

#   read ellipse data
#   cnt = 2
    for i in range( 3, len(lines) ):
        iline = lines[i]
        spl = [x.strip() for x in iline.split(',')]
        x = float(spl[0])
        y = float(spl[1])
        a = float(spl[2])
        b = float(spl[3])
        ang = float(spl[4])
        newEllip = Ellipse(x, y, a, b, ang)
        allEllipses.append(newEllip)
        if i == 3:
            radius.append(a)
            radius.append(b)
            # print(radius[0])

# plot the result
fig, ax = plt.subplots()
# these next few lines are pretty important because
# otherwise your ellipses might only be displayed partly
# or may be distorted
a = radius[0]
b = radius[1]
ax.set_xlim([-2*a, Rx+2*a])
ax.set_ylim([-2*a, Ry+2*a])
ax.set_aspect('equal')
ax.add_patch(
  patches.Rectangle(
    (0.0, 0.0),
    Rx,
    Ry,
    fill=False      # remove background
  ))
print("ellipse number:", len(allEllipses))
for test in allEllipses:
  shape = test.create_ellipse( 1.0 )
  verts = np.array(shape.exterior.coords.xy)
  patch = Polygon(verts.T, color='blue', alpha=0.5)
  ax.add_patch(patch)
plt.show()



# # output to geo file test.lengths test.angle
# geom = pygmsh.opencascade.Geometry(
#   characteristic_length_min=radius[0]/20, characteristic_length_max=radius[0]/10
# )

# rectangle = geom.add_rectangle([0.0, 0.0, 0.0], Rx, Ry)

# # to put all ellipse, called disk in opencascade
# allEllipsePS = []
# # print(len(allEllipses))

# if radius[0] == radius[1]: # output circular particles
#   for isp in allEllipses:
#       #   get the center and radius of all ellipse
#       x = isp.center[0]
#       y = isp.center[1]
#       a = isp.lengths[0]
#       disk = geom.add_disk([x, y, 0.0], a)
#       allEllipsePS.append(disk)   
# else: # output ellipse inclusions
#   for isp in allEllipses:
#       #   get the center and radius of all ellipse
#       x = isp.center[0]
#       y = isp.center[1]
#       a = isp.lengths[0]
#       b = isp.lengths[1]
#       ori = isp.angle/180*np.pi
#       # print(x, y, a, b, ori)
#       p1 = geom.add_point([x, y, 0.0])
#       p2 = geom.add_point([x+a*np.cos(ori), y+a*np.sin(ori), 0.0])
#       p3 = geom.add_point([x-b*np.sin(ori), y+b*np.cos(ori), 0.0])
#       p4 = geom.add_point([x-a*np.cos(ori), y-a*np.sin(ori), 0.0])
#       p5 = geom.add_point([x+b*np.sin(ori), y-b*np.cos(ori), 0.0])
#       arc1 = geom.add_ellipse_arc(p2, p1, p2, p3)
#       arc2 = geom.add_ellipse_arc(p4, p1, p4, p3)
#       arc3 = geom.add_ellipse_arc(p2, p1, p2, p5)
#       arc4 = geom.add_ellipse_arc(p4, p1, p4, p5)
#       arc = (arc1, arc2, arc3, arc4)
#       ellipseLL = geom.add_line_loop(arc)
#       ellipsePS = geom.add_plane_surface(ellipseLL)
#       allEllipsePS.append(ellipsePS)

# boldiff = geom.boolean_difference([rectangle], allEllipsePS, delete_first=False,delete_other=False)
# bolinte = geom.boolean_intersection([rectangle], allEllipsePS)

# geom.add_physical(boldiff, label='matrix')

# geom.add_physical(bolinte, label='inclusion')

# geom.add_raw_code("Mesh 2;")
# geom.add_raw_code("Coherence Mesh;")

# if os.path.exists('./rve.geo'):
#   os.remove("rve.geo")

# mesh = pygmsh.generate_mesh(geom, geo_filename="rve.geo") #mesh_file_type="msh"

# meshio.write("rve.msh", mesh, 'gmsh2-ascii', write_binary=False)

# # geom.add_raw_code("Periodic Curve {{{}}} = {{{}}};".format(l0.id, l2.id))

