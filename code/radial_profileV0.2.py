##### README ######
# Author: Gayathri Nadar, SCF MPI-CBG (nadar@mpi-cbg.de)
# Version 0.0 (14.01.2019) : radial profile test
# Version 0.1 (15.01.2019): script adapted from Alertoduina's MyOvalProfile.java on Github
# Version 0.2 (21.01.2019): additional features in script:
#								- use npoints from roi.size() i.e. all points in roi boundary
#								- split channels and get plots for each channel

# Description
# Oval_Profile takes the image region bounded by an oval region and samples the
# oval at N equal angles around the oval. The program generates a ProfilePlot
# of Pixel intensities at equiangular sample points along the oval profile.

# Usage
# Open an image and script with FIJI.
# Draw a circular roi around the region you want to explore and click run.

# Output
# a radial plot
# an image with the sampled points as overlay

# Note:
# The starting point is seen as a CIRCLE in the output image for reference!
# The circle is sampled from this starting point in the CLOCKWISE direction equi-angularly.
# Split the channels if needed!

###=============================================================================###

#@ImagePlus imp0

from ij import IJ
from ij import ImagePlus, ImageStack
from ij.process import ImageProcessor, ImageStatistics, FloatPolygon
from ij.measure import ResultsTable, Calibration
from ij.plugin import Duplicator
from ij.gui import ShapeRoi, Overlay, Roi, PolygonRoi, Toolbar, Plot, PointRoi, Overlay
from ij.gui import WaitForUserDialog, ImageWindow
from ij.plugin import ChannelSplitter
from ij.measure import Measurements
from ij.measure import ResultsTable
from ij import WindowManager
from ij.text import TextWindow
import java.awt.Color as Color
import math

global start 
start = 0

def getOvalProfile(npoints, ip, roi, start):
	profile = []
	xValues = []
	xcoords = []
	ycoords = []

	# get rectangular bounds, width, height and center
	b = roi.getBounds()
	width = b.width
	height = b.height
	w2 = width*width / 4.0
	h2 = height*height / 4.0
	cx = b.x+width/2.0 
	cy = b.y+height/2.0 

	print "center_x ", cx
	print "center_y ", cy

	# get sampling angle
	tink = 2 * math.pi / npoints
	
	theta = 0
	i2 = 0

	for i in range(npoints):

		# get starting factor
		i2 = i + start

		if i2 > npoints:
			i2 = i2 - npoints

		# multiple factor with sampling angle to get theta
		theta = i2 * tink
#		print "theta ", theta

		# get base and height
		dx = math.cos(theta)
		dy = math.sin(theta)
		x = cx
		y = cy
		hotx = 0
		hoty = 0

		while(roi.contains(int(x), int(y))):
			x += dx
			y += dy

		m = math.sqrt(w2 * h2 / (dx * dx * h2 + dy * dy * w2))

		# coordinates along the circumference to get profile
		hotx = cx + dx * m
		hoty = cy + dy * m

#		print "hotx ", hotx
#		print "hoty ", hoty
#		print "I=", ip.getPixel(int(hotx), int(hoty))

		# get intensities and theta in degrees and append to array
#		profile.append(ip.getPixel(int(hotx), int(hoty)))
		profile.append(ip.getInterpolatedPixel(hotx, hoty))
		xcoords.append(hotx)
		ycoords.append(hoty)
		xValues.append(math.degrees(theta))
		
	return profile, xValues, xcoords, ycoords


#### main ####
roi = imp0.getRoi()

if roi==None or roi.getType() != Roi.OVAL:
	IJ.error("Oval selection required!")

n_ch = imp0.getNChannels()

npts = roi.size()
print "No of points in roi boundary= ", npts

imps = ChannelSplitter.split(imp0)

for idx, imp in enumerate(imps, 1):

	imp.setRoi(roi)
#	imp.show()

	ip = imp.getProcessor()
	stats = ip.getStats()
	minI= stats.min
	maxI= stats.max
	
	p, theta, xcoords, ycoords = getOvalProfile(npts, ip, roi, start)

	print "Interpolated intensities= ", p
	print "Sampled angles= ", theta

	plot = Plot("Intensities around a circumference ch: "+ str(idx), "angle", "intensity")
	plot.setColor(Color.RED)
	plot.setLimits(-45, 380, 0, maxI)
	plot.addPoints( theta, p, Plot.CONNECTED_CIRCLES)
	plot.show()

	pp = PointRoi(xcoords, ycoords)
	pp.setPointType(1)
	starting_point = PointRoi(xcoords[0], ycoords[0])
	starting_point.setPointType(2)
	starting_point.setSize(2)
	ov = Overlay()
	ov.add(pp)
	ov.setStrokeColor(Color.MAGENTA)
	ov.add(starting_point)


imp0.killRoi()
imp_final = Duplicator().run(imp0)
imp_final.setOverlay(ov)
#imp.setRoi(roi)
imp_final.show()
imp0.setRoi(roi)

    
	
	

	
	

	
	

	


