#@ CommandService command 

from de.csbdresden.stardist import StarDist2D
from net.imglib2.img.display.imagej import ImageJFunctions;
from de.mpicbg.scf.fijiplugins.ui.roi import LabelMapToRoiManagerPlugin;

from ij import IJ
from ij import ImagePlus
from ij.plugin import Duplicator, RGBStackMerge, ChannelSplitter, ZProjector, ImageCalculator
from ij.gui import Roi, PointRoi, WaitForUserDialog, Overlay
from ij.measure import Measurements, ResultsTable
from ij.io import FileSaver 
from ij.plugin.frame import RoiManager

import loci.plugins
from loci.plugins import BF
from loci.plugins.in import ImporterOptions 
from loci.formats import ImageReader
from java.awt import Color as color 

import os, math 
from os.path import isfile, join

from inra.ijpb.binary import BinaryImages
from inra.ijpb.measure import IntensityMeasures
from inra.ijpb.label import LabelImages
from inra.ijpb.measure.region2d import Centroid, InertiaEllipse
from inra.ijpb.measure import IntensityMeasures


def filterBySize(imp_label, imp_calibration = None, diameter_in_um = 5.0):
	"""
	Filter out regions having avg diamater < 5 um. Obj greater than threshold will be kept.
	For each region compute major and minor axis of fitted ellipse and take its average as diameter 

	params: label imp, calibration of raw image, size threshold 
	returns: new label imp with obj with diameter > size threshold 
	"""

	# find all labels 
	labels = LabelImages.findAllLabels(imp_label)

	keep = []

	ellipses = InertiaEllipse.inertiaEllipses(imp_label.getProcessor(), labels, imp_calibration)
	radius1 = [e.radius1() for e in ellipses]
	radius2 = [e.radius2() for e in ellipses]

	avgdiameters = [r1 + r2 for r1, r2 in zip(radius1, radius2)]

	for l,d in zip(labels, avgdiameters):
		if d > diameter_in_um:
			keep.append(l)

	imp_keep = imp_label.duplicate() 
	imp_keep = LabelImages.keepLabels(imp_label, keep)

	return imp_keep 




imp = IJ.getImage(); 

stats = imp.getStatistics();
print stats.mean 

res = command.run(StarDist2D, False,
		"input", imp, "modelChoice", "Versatile (fluorescent nuclei)",
		"outputType", "Label Image", "nTiles", 32).get()

label = res.getOutput("label")
labelimp = ImageJFunctions.wrap(label, "labelimp");
IJ.run(labelimp, "glasbey_on_dark", "");
IJ.run(labelimp, "Enhance Contrast", "saturated=0.35");

labelimp.show()


rt = ResultsTable()
rt.reset();
im = IntensityMeasures(imp, labelimp)

rt = im.getMean()
mean = rt.getColumn(0)

print len(mean)

ids = LabelImages.findAllLabels(labelimp)

count = 0
pos = []
for m, i in zip(mean, ids):
#	print m
	if m > stats.mean:
#		print "Greater"
#		print "Mean =", m 
		count += 1
#		print m, i
		pos.append(i)



positive_antibody = [ids[i1] for i1, val1 in enumerate(mean) if val1 > stats.mean] 

print len(positive_antibody)
print len(pos)
print count 


#for p in positive_antibody:
#	print p



pos1 = LabelImages.keepLabels(labelimp, positive_antibody)
pos2 = LabelImages.keepLabels(labelimp, pos)

pos1.show()
pos2.show()



labelfilter = filterBySize(labelimp, imp.getCalibration(), diameter_in_um = 5)
labelfilter.show()











	