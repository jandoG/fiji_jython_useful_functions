#@ File (label= "Open Label1", style= "file") f1 
#@ File (label= "Open Label2", style= "file") f2 
#@ File (label= "Open Label3", style= "file") f3

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

def openImageWithBF(path, virtual= True, groupfiles = True):
	"""
	set options to open image using bio-formats- use virtual for quick loading
	"""
	options = ImporterOptions()
	options.setColorMode(ImporterOptions.COLOR_MODE_DEFAULT)
	options.setAutoscale(True)
	options.setStackFormat("Hyperstack")
	options.setVirtual(virtual)
	options.setGroupFiles(groupfiles) 
	options.setId(path)
	imp = BF.openImagePlus(options)[0]	

	return imp 
	
def getDistance(c1 , c2):
	"""
	Compute euclidean distance between 2 coordinates.

	params: c1, c2 where c1 = [x1, y1] and c2 = [x2, y2]
	retruns: distance 
	"""
	dist = ( ((c2[0] - c1[0])**2) + ((c2[1] - c1[1])**2) )

	return math.sqrt(dist)


def closest(cur_pos, positions, maxdist = 100.0):
	"""
	Find the closest point from list to a point, within the max distance threshold 

	params: cur_pos = a point (x, y)
	positions = list of positions 
	returns: closest coordinates from the list of positions
	"""
	closestpts = None
	for pos in positions:
		dist = getDistance(cur_pos , pos)
		if dist >= maxdist:
			continue
		closestpts = pos
		maxdist = dist

	return closestpts

def applyGlasbeyLUT(imp):
	"""Applies Glasbey on Dark LUT on label maps"""
	IJ.run(imp, "glasbey_on_dark", "");


l1 = openImageWithBF(f1.getPath(), False, False);
l2 = openImageWithBF(f2.getPath(), False, False);
l3 = openImageWithBF(f3.getPath(), False, False);

applyGlasbeyLUT(l1)
applyGlasbeyLUT(l2)
applyGlasbeyLUT(l3)

l1.show()
l2.show()
l3.show()

labels1 = LabelImages.findAllLabels(l1) 
labels2 = LabelImages.findAllLabels(l2)
labels3 = LabelImages.findAllLabels(l3)

centroids_1 = Centroid().centroids(l1.getProcessor(), labels1)
centroids_2 = Centroid().centroids(l2.getProcessor(), labels2)
centroids_3 = Centroid().centroids(l3.getProcessor(), labels3)

ids_c1 = []
ids_c2 = []
ids_c3 = []
usedl2 = []
usedl3 = []

distance_threshold = 5.0

# for each in src find neighbor in target
for i1, (c1, ll1) in enumerate(zip(centroids_1, labels1)):
	close = closest(c1, centroids_2, distance_threshold)	
	
	if close is not None:
		close1 = closest(close, centroids_3, distance_threshold) 
		
		if close1 is not None:
			print "l1-l2 =", close 
			print "l2-l3 =", close1
			
			idx = centroids_2.index(close)
			idx1 = centroids_3.index(close1)

			if labels2[idx] not in usedl2 and labels3[idx1] not in usedl3:
				ids_c2.append(labels2[idx])
				ids_c3.append(labels3[idx1])
				ids_c1.append(ll1)
				usedl2.append(labels2[idx])
				usedl3.append(labels3[idx1])


overlap1 = LabelImages.keepLabels(l1, ids_c1)
overlap2 = LabelImages.keepLabels(l2, ids_c2)
overlap3 = LabelImages.keepLabels(l3, ids_c3)


overlap1.show() 
overlap2.show() 
overlap3.show()






