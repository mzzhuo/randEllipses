# ---------------------------------------------------------------------------
# The Python scripts use the Random  Sequential Algorithm for generation of 
# randomly packed ellipses/circles 
# 
# features and difficulty: non-overlapping and periodic boundary conditions
# 
# by Mingzhao Zhuo @ TU Delft, Oct. 2019
# ---------------------------------------------------------------------------

from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon

from shapely.geometry.point import Point
from shapely.geometry import LineString
from shapely import affinity

import numpy as np
import random
import pygmsh
import os

import toGmsh


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

# detect the edges intersect with the ellipse
def intersectEdges ( ellipShape, edges ):
	whichedges = []
	for ieg in range(4): # 0, 1, 2, 3
		iedge = edges[ieg]
		intersecEdge = ellipShape.intersects(iedge)
		if intersecEdge:
			whichedges.append(ieg)
	return whichedges


def randEllipsePack(a, b, Rx, Ry, angRang, neps, ratio, filename):

#   flag indicating whether creating one configuration success or not
	yesGene = False

#   4 corner nodes
	lbPt = Point(0.0, 0.0)		# 0
	rbPt = Point(Rx, 0.0)		# 1
	rtPt = Point(Rx, Ry)			# 2
	ltPt = Point(0.0, Ry)		# 3
	cornerNodes = (lbPt, rbPt, rtPt, ltPt)

#   4 boundary edges
	left = LineString([(0, 0), (0, Ry)])  # 0
	bottom = LineString([(0, 0), (Rx, 0)])  # 1
	right = LineString([(Rx, 0), (Rx, Ry)])  # 2
	top = LineString([(0, Ry), (Rx, Ry)])  # 3
	edges = (left, bottom, right, top)

#   count of successfully placed ellipses
	ig = 0

#   put all generated ellipses and their periodic copies
	allEllipses = []

#   for each experiment of inserting neps non-overlapping ellipses (create one config),
#   maximum trial times neps*60, return success or not

	# while ig <= neps 
	for trialTimes in range(1, neps*60): 

#		generate a random ellipse
		x = random.uniform(0.0, Rx)
		y = random.uniform(0.0, Ry)
		ang = random.uniform(angRang[0], angRang[1])
		newEllip = Ellipse(x, y, a, b, ang)
		newshape = newEllip.create_ellipse(ratio)
	
#       put the new generated ellipse and its periodic copies here 
		ellipCopies = []

	#   check whether the new ellipse contains any of the four corner nodes
		containCorner = False

	#   assume newly generated ellipse and its copies will be kept 
		keep = True

		for ipt in range(4):  # 0 1 2 3
			inode = cornerNodes[ipt]
			containCorner = newshape.contains(inode)
			if containCorner:
				cornerTrans = (([0, 0], [Rx, 0], [Rx, Ry], [0, Ry]), 		# left bottom pt
							   # right bottom pt
							   ([-Rx, 0],  [0, 0], [0, Ry], [-Rx, Ry]),
							   ([-Rx, -Ry],  [0, -Ry], [0, 0], [-Rx, 0]),  # right top pt
							   ([0, -Ry],  [Rx, -Ry], [Rx, 0], [0, 0]))  # left top pt
				# print("current ellipse contains one corner node")
				trans = np.array(cornerTrans[ipt])
				pos = np.array([x, y])
				posArray = pos + trans
				for ieps in range(4):
					newEllip_copy = Ellipse( posArray[ieps, 0], posArray[ieps, 1], a, b, ang )
					ellipCopies.append( newEllip_copy )
				break

#		not contain corner node then either intersect with edges
# 		or inside
		if not containCorner:

# 			get which edges intersect
			whichedges = intersectEdges( newshape, edges )


			if len(whichedges) == 0:

# 				no intersect = inside

				ellipCopies.append(newEllip)

# 			only intersect with one edge

			elif len(whichedges) == 1:

				edgeTrans = ( [Rx, 0], [0, Ry], [-Rx, 0], [0, -Ry] )

				trans = np.array( edgeTrans[ whichedges[0] ] )
				pos = np.array([x, y])
				pos_copy = pos + trans
				newEllip_copy = Ellipse(pos_copy[0], pos_copy[1], a, b, ang)
				newEllip_copy_shape = newEllip_copy.create_ellipse(ratio)

# 				check whether the copy also intersect with only one edge

				whichedges4copy = intersectEdges( newEllip_copy_shape, edges )

				if len(whichedges4copy) == 1:
					ellipCopies.append(newEllip_copy)
					ellipCopies.append(newEllip)
				else:
					keep = False

			else:
				keep = False

		if keep: 

#		 check whether the new ellipse overlap with existing ones
			if len(allEllipses) > 0:
				for i in allEllipses:
					itershape = i.create_ellipse(ratio)
					overlap = False
					for j in ellipCopies:
						jitershape = j.create_ellipse(ratio)
						overlap = jitershape.intersects(itershape)
						if overlap:
							keep = False
							break
					if overlap:
						keep = False
						break

# 		if all above conditions meet, store them to allEllipses and add one successful count
		if keep:
			allEllipses.extend(ellipCopies)
			ig = ig + 1

# 		if desired number achieved, stop inserting
		if ig == neps :
			yesGene = True
			break
	
	# output geometry information

	if yesGene:
		filenametxt = filename + '.txt'
	
		with open(filenametxt, 'w') as f:
			f.write( 'this is the length of box [0,0] - [Rx,Ry]\n' )
			f.write( '{:12.6f},{:12.6f}\n'.format(Rx, Ry) )
			f.write( 'center, axis a and b, angle in degree\n' )
			for test in allEllipses:
				f.write('{:15.5e},'.format(test.center[0]))
				f.write('{:15.5e},'.format(test.center[1]))
				f.write('{:15.5e},'.format(test.lengths[0]))
				f.write('{:15.5e},'.format(test.lengths[1]))
				f.write('{:15.5e}\n'.format(test.angle))
	
		#  write to geo and msh file
		
		toGmsh.writeGeo(a, b, Rx, Ry, allEllipses, filename)

	return yesGene


if __name__ == "__main__":

	# volume fraction of inclusion
	vf = 0.5
	# number of ellipses
	neps = 25
	
	# ellipse size
	r = 2 ; # a/b axis ratio
	b = 1 
	a = b * r
	
	# ratio: expand the ellipse when compute intersection so that they are not too close
	ratio = 1.05
	
	# simulation box size
	Rx = np.sqrt(np.pi * a * b * neps / vf)
	Ry = np.sqrt(np.pi * a * b * neps / vf)
	# Rx = np.sqrt(np.pi * neps / vf) * a
	# Ry = np.sqrt(np.pi * neps / vf) * b

	ang_sta = 0
	ang_end = 180
	ang = [ang_sta, ang_end]
	
	filename = './RVE' + str(neps) 
	
	okgene = False
	
	while not okgene:
		okgene = randEllipsePack(a, b, Rx, Ry, ang, neps, ratio, filename)






