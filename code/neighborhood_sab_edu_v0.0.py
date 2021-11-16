# ============================================= README ======================================================
#
# Author: Gayathri Nadar, SCF MPI-CBG, nadar@mpi-cbg.de
#
# Version history:
# v0.0: 16.03.2021 counting EdU cells in an image which are within certain distance to each Sabgal cell in another image 
#
# Description:
#	This script takes two input label maps - sabgal and edu. For each label in sabgal, all cells within certain threshold distance are counted as neighbors. 
#
# Requirements:
#	Please install the following plugins on FIJI before running this script.
#	1. CSBDeep
#	2. IJPB- Plugins 
#	3. SCF-MPI-CBG 
#	4. Stardist 
#	On FIJI go to Help > Update and then "Manage update sites" and check the above mentioned plugins boxes. 
#	Click "close" and then "Apply changes". Then restart FIJI. 
#
# Steps:
# 	- Open images, get all labels 
#	- for the labels, get centroids 
#	- compute distance from each label in sabgal to all labels in edu 
#	- if a distance is below the threshold, this edu cell is a neighbor 
#	
# Usage:
#	- Open this script with FIJI and click run. 
#	- A window will pop-up, enter the necessary values and click ok to run the script.
#	- Hover over to uncover more details on the params.
#
# Parameters to be set:
#	- Distance threshold in pixels: max distance allowed between centroids of sabgal and edu cells 
#
# Output:
#	- label images with only neighbors 
#	- table with label no of sabcell and corresponding neighbor labels in edu
#	- original image with overlay of neighbors - yellow = sabgal, red = edu neighbor
#
# =============================================================================================================



# =============================================================================================================
# Don't modify anything below here without good reason
# =============================================================================================================

#@File (label= "(Cropped) SABGAL image", style= File) sab_
#@File (label= "SABGAL label map", style= File, description = "Label map of the same image as above") sabfile 
#@File (label= "(Cropped) EdU image", style= File) edu_
#@File (label= "EdU label map", style= File, description = "Label map of the same image as above") edufile 
#@Float (label= "Distance threshold between centroids (pixels)", value= 30) dist_threshold

import os, os.path 
import glob, math
from glob import glob 
from os.path import isfile, join, basename

from ij import ImagePlus, IJ
from ij.io import FileSaver
from ij.plugin import Duplicator
from ij.plugin.frame import RoiManager
from ij.gui import Overlay, Roi, TextRoi
from ij.gui import Line, ProfilePlot
from loci.plugins import BF
from loci.formats import ImageReader
from loci.formats import MetadataTools
from loci.plugins.prefs import OptionsList
from loci.plugins.in import ImporterOptions
from java.awt import Color, Font 

from de.mpicbg.scf.fijiplugins.ui.roi import LabelMapToRoiManagerPlugin;
from de.mpicbg.scf.fijiplugins.ui.labelmap import ThresholdLabelingPlugin
from net.imglib2.img.display.imagej import ImageJFunctions;

from inra.ijpb.label import LabelImages
from ij.measure import ResultsTable
from inra.ijpb.binary import BinaryImages
from inra.ijpb.measure.region2d import Centroid
from inra.ijpb.measure import IntensityMeasures
from inra.ijpb.measure.region2d import BoundingBox

rm = RoiManager.getInstance()
if not rm:
	rm = RoiManager()

def main():
	# open original and label maps
	sabimp = openImageWithBF(sab_.getPath(), virtual= False)
	eduimp = openImageWithBF(edu_.getPath(), virtual= False)

	sab_label = openImageWithBF(sabfile.getPath(), virtual= False)
	edu_label = openImageWithBF(edufile.getPath(), virtual= False)

	impname = os.path.splitext(sab_label.getTitle())[0]
	name = impname.split("labels_")[1]

	# change LUT 
	applyGlasbeyLUT(sab_label)
#	sab_label.show() 
	applyGlasbeyLUT(edu_label)
#	edu_label.show()

	# find all labels in sab and edu 
	sab = LabelImages.findAllLabels(sab_label)
	edu = LabelImages.findAllLabels(edu_label)

	# compute neigbors 
	neighbor_pairs, no_edu_within_distance, sab_neighbors, edu_neighbors = getNoNeighbors(sab_label, edu_label, distance_threshold= dist_threshold)
#	sab_neighbors.show() 
#	edu_neighbors.show()

	# get rois for overlay 
	overlay = Overlay()
	roiarray_sabN = getRois(sab_neighbors)
	for r in roiarray_sabN:
		r.setStrokeColor(Color.yellow)
		overlay.add(r)
	roiarray_eduN = getRois(edu_neighbors)
	for r in roiarray_eduN:
		r.setStrokeColor(Color.red)
		overlay.add(r)

	sabimp.setDisplayMode(IJ.COMPOSITE)
	sabimp.setOverlay(overlay)
	sabimp.show() 
	eduimp.setOverlay(overlay)
	eduimp.show()

	# list to store edu neighbor labels for each sab label (for which there are neighbors) 
	edu_labelN = []
	sab_labelN = LabelImages.findAllLabels(sab_neighbors)

	# show which label has how many neighbors 
	for l, n in zip(sab, no_edu_within_distance):
		print "SAB label {} has EdU neighbors= {}".format(l, n)
		if n != 0:
			# extract neighbors for each sab where no neighbor > 0
			edus = [x[1] for x in neighbor_pairs if x[0] == l]            # for eg for neighbor pairs (2, 5), (2, 10), output = [5, 10]
			edu_labelN.append(edus)  

	table = ResultsTable()
	for i, j in zip(sab_labelN, edu_labelN):
		table.incrementCounter();
		table.addValue("Sab label", i)
		table.addValue("EdU neighbor labels", str(j))
		table.addValue("No of EdU neighbors", len(j)) 
		table.addValue("Distance threshold (pixels)", dist_threshold)

	table.show(name + " - Sab-edu neighbors")

	# save dir for output 
	results_path = join(os.path.dirname(sabfile.getPath()), 'neighborhoodAnalysis_' + name)
	if not os.path.exists(results_path):
		os.makedirs(results_path)

	# save output 
	FileSaver(sab_neighbors).saveAsTiff(join(results_path, "sabgal_withneighbors_" + name + ".tif"))
	FileSaver(edu_neighbors).saveAsTiff(join(results_path, "edu_neighbors_" + name + ".tif"))
	FileSaver(eduimp).saveAsTiff(join(results_path, "edu_withoverlay_" + name + ".tif"))
	FileSaver(sabimp).saveAsTiff(join(results_path, "sab_withoverlay_" + name + ".tif"))
	table.saveAs(join(results_path, "Sab-edu neighbors -" + name + ".csv"));


def getNoNeighbors(src_label, target_label, distance_threshold= 5.0):
	"""
	Find labels in src which has a neighbor in target within the distance threshold. 

	params: src and target label maps, min distance between centroids 
	returns: filtered src label map with only labels which have a neighbor in target label map 
	"""
	# find all labels and their centroids
	labels1 = LabelImages.findAllLabels(src_label) 
	labels2 = LabelImages.findAllLabels(target_label)
	centroids_1 = Centroid().centroids(src_label.getProcessor(), labels1)
	centroids_2 = Centroid().centroids(target_label.getProcessor(), labels2)

	# array to store no of neighbors and neighbor pairs
	no_labels_withinDist = [0] * len(labels1)
	neighbor_pairs = []
	
	for i1, (c1, l1) in enumerate(zip(centroids_1, labels1)):
		for i2, (c2, l2) in enumerate(zip(centroids_2, labels2)):
			dist = getDistance(c1 , c2)
			if dist < distance_threshold:
				print "Neighbor pairs: label1 = {}, label2 = {}". format(l1, l2)
				neighbor_pairs.append((l1, l2))
				no_labels_withinDist[i1] = no_labels_withinDist[i1] + 1

	# src labels with only labels which have neighbors in trg 
	src_label_with_neighbors = list(set([i[0] for i in neighbor_pairs]))

	# trg labels with neighbors of src 
	trg_neighbors = [i[1] for i in neighbor_pairs]

	# create label map from labels 
	neighbor_src = LabelImages.keepLabels(src_label.getProcessor(), src_label_with_neighbors)
	neighbor_trg = LabelImages.keepLabels(target_label.getProcessor(), trg_neighbors) 
	
	src_n = ImagePlus("src neighbors", neighbor_src)
	trg_n = ImagePlus("target neighbors", neighbor_trg)

	applyGlasbeyLUT(src_n)
	applyGlasbeyLUT(trg_n)

	return neighbor_pairs, no_labels_withinDist, src_n, trg_n

		
def getDistance(c1 , c2):
	"""
	Compute euclidean distance between 2 coordinates.

	params: c1, c2 where c1 = [x1, y1] and c2 = [x2, y2]
	retruns: distance 
	"""
	dist = ( ((c2[0] - c1[0])**2) + ((c2[1] - c1[1])**2) )

	return math.sqrt(dist)


def applyGlasbeyLUT(imp):
	"""Applies Glasbey on Dark LUT on label maps"""
	IJ.run(imp, "glasbey_on_dark", "");


def getRois(labelimp):
	"""
	Converts a label map into ROIs and adds it to ROI manager.

	params: label map 
	returns: Array of rois from roi manager
	"""
	LabelMapToRoiManagerPlugin.apply(labelimp)
	roiarray = rm.getRoisAsArray();

	try:
		return roiarray
	finally:
		rm.reset();


def openImageWithBF(path, virtual= True):
	# set options to open image - use virtual for quick loading
	options = ImporterOptions()
	options.setColorMode(ImporterOptions.COLOR_MODE_DEFAULT)
	options.setAutoscale(True)
	options.setStackFormat("Hyperstack")
	options.setVirtual(virtual)
	options.setId(path)
	
	imp = BF.openImagePlus(options)[0]	

	return imp 


main() 
print("Done")
