# ---------------------------------------------------------------------------
# The Python scripts write the generated inclusions to Gmsh
# 
# features and difficulty: non-overlapping and periodic boundary conditions
# 
# by Mingzhao Zhuo @ TU Delft, Oct. 2019
# ---------------------------------------------------------------------------


import numpy as np
import pygmsh
import meshio
import os

def writeGeo (a, b, Rx, Ry, allEllipses, filename):

	geom = pygmsh.opencascade.Geometry( characteristic_length_min=a/20, characteristic_length_max=a/10 )

	rectangle = geom.add_rectangle([0.0, 0.0, 0.0], Rx, Ry)

#  to put all ellipse, called disk in opencascade
	allEllipsePS = []

	# print('test')

	if a == b: # output circular particles 
		for isp in allEllipses:
			#   get the center and radius of all ellipse
			x = isp.center[0]
			y = isp.center[1]
			a = isp.lengths[0]
			disk = geom.add_disk([x, y, 0.0], a)
			allEllipsePS.append(disk)   
	else: # output ellipse inclusions
		for isp in allEllipses:
			#   get the center and radius of all ellipse
			x = isp.center[0]
			y = isp.center[1]
			a = isp.lengths[0]
			b = isp.lengths[1]
			ori = isp.angle/180*np.pi
			# print(x, y, a, b, ori)
			p1 = geom.add_point([x, y, 0.0])
			p2 = geom.add_point([x+a*np.cos(ori), y+a*np.sin(ori), 0.0])
			p3 = geom.add_point([x-b*np.sin(ori), y+b*np.cos(ori), 0.0])
			p4 = geom.add_point([x-a*np.cos(ori), y-a*np.sin(ori), 0.0])
			p5 = geom.add_point([x+b*np.sin(ori), y-b*np.cos(ori), 0.0])
			arc1 = geom.add_ellipse_arc(p2, p1, p2, p3)
			arc2 = geom.add_ellipse_arc(p4, p1, p4, p3)
			arc3 = geom.add_ellipse_arc(p2, p1, p2, p5)
			arc4 = geom.add_ellipse_arc(p4, p1, p4, p5)
			arc = (arc1, arc2, arc3, arc4)
			ellipseLL = geom.add_line_loop(arc)
			ellipsePS = geom.add_plane_surface(ellipseLL)
			allEllipsePS.append(ellipsePS)

#   need to modify the boolean_intersection function in geometry.py to the same as boolean_difference

	boldiff = geom.boolean_difference([rectangle], allEllipsePS, delete_first=False,delete_other=False)
	bolinte = geom.boolean_intersection([rectangle], allEllipsePS)

	geom.add_physical(boldiff, label='matrix')

	geom.add_physical(bolinte, label='inclusion')

	geom.add_raw_code("Mesh 2;")
	geom.add_raw_code("Coherence Mesh;")

# -----------------------------------------------------------------
#   just get matrix
	# bol = geom.boolean_difference([rectangle], allEllipsePS)
	# geom.add_physical(bol, label='spe')
# -----------------------------------------------------------------

	filenamegeo = filename + '.geo'
	if os.path.exists( filenamegeo ):
		os.remove( filenamegeo )
#   write geo file
	mesh = pygmsh.generate_mesh(geom, geo_filename=filenamegeo) #mesh_file_type="msh"

#   write msh file
	filenamemsh = filename + '.msh'
	meshio.write(filenamemsh, mesh, 'gmsh2-ascii', write_binary=False)
