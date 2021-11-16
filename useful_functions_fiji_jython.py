#@ CommandService command 

from de.csbdresden.stardist import StarDist2D
from ij import IJ
from ij import Prefs
from ij.io import FileSaver
from ij.io import Opener
from ij import ImagePlus
from ij import ImageStack
from ij.plugin import Duplicator
from ij.plugin import HyperStackConverter
from ij.plugin import ZProjector
from ij.plugin import Concatenator
from jarray import array
from loci.plugins import BF
from loci.formats import ImageReader, FilePattern
from loci.formats import MetadataTools
from loci.plugins.prefs import OptionsList
from loci.plugins.in import ImporterOptions
from loci.common import Region
from ij.gui import GenericDialog, NonBlockingGenericDialog, WaitForUserDialog, Overlay, Roi, PointRoi
from de.mpicbg.scf.spotcoloc import SpotProcessor, SpotVisualization
from trainableSegmentation import WekaSegmentation
from trainableSegmentation.utils import Utils; 
import os, math, re
import fnmatch as fnm
from os.path import isfile, join
from os import listdir 
import glob
import time
from de.mpicbg.scf.fijiplugins.ui.roi import LabelMapToRoiManagerPlugin;
from de.mpicbg.scf.fijiplugins.ui.labelmap import ThresholdLabelingPlugin
from net.imglib2.img.display.imagej import ImageJFunctions;
from inra.ijpb.label import LabelImages
from inra.ijpb.binary import BinaryImages
from inra.ijpb.measure.region2d import Centroid
from inra.ijpb.measure import IntensityMeasures
from inra.ijpb.measure import IntrinsicVolumes2D
from inra.ijpb.measure.region2d import BoundingBox
from ij.process import ImageStatistics
from ij.process import StackStatistics
from ij.gui import Overlay, Roi
from ij.process import FloodFiller
from ij.measure import ResultsTable
from java.awt import Color
import time, datetime, sys

def openImageWithBF(path, virtual= True, groupfiles = False, seriesdata = True, openseries = 1):
	"""
	set options to open image using bio-formats- use virtual for quick loading

	params:
	path to file 
	virtual: bool set True to load image as virtual stack (for big data)
	groupfiles: bool set True to load images from folder having similar name pattern (stored in the metadata of first image)
	seriesdata: bool set True if image contains series 
	openseries: int series ID no to be opened 

	returns: 
	imp: ImagePlus 
	"""
	options = ImporterOptions()
	options.setColorMode(ImporterOptions.COLOR_MODE_DEFAULT)
	options.setAutoscale(True)
	options.setStackFormat("Hyperstack")
	options.setVirtual(virtual)
	options.setGroupFiles(groupfiles) 
	
	if seriesdata: 
		reader = ImageReader()
		reader.setId(path)
		seriesCount = reader.getSeriesCount()
		reader.close()
		print "Found series data. Image series count =", seriesCount
		print "Reading series ID ", openseries, "\n"
		options.setOpenAllSeries(True)

	options.setId(path)
	allimps = BF.openImagePlus(options)

	if seriesdata:
		imp = allimps[openseries - 1]

	else:
		imp = allimps[0]
	
	return imp

def getStackThreshold(imp):
	"""
	Compute manually a threshold value for the image based on the image histogram. 
	Works for 2D and 3D stack. 
	For stack, the entire stack histogram is considered to compute threshold (NOT slice by slice!). 
	
	Threshold = IMaxFreq + threshMultiplyer * ISigma
	Where:
	IMaxFreq = the intensity where histogram value is maximum
	threshMultiplyer = to be set by user, higher values tend to provide smoother outlines 
	ISigma = std dev of low intensity half of the histogram peak (till IMaxFreq)

	params:
	imp 
	
	returns: 
	threshold, IMaxFreq, ISigma
	"""
	stats = StackStatistics(imp)
	hist = stats.getHistogram()

	hMin = stats.histMin
	hMax = stats.histMax
	nBins = len(hist);

	# find the maximum of the histogram and peak variance
	maxFreq = 0;
	binIdxMaxFreq=0
	binIdx=0

	for freq in hist:
		if freq > maxFreq:
			maxFreq = freq
			binIdxMaxFreq = binIdx
		binIdx = binIdx+1;

	IMaxFreq = hMin + binIdxMaxFreq*(hMax-hMin)/nBins;

	# measure the peak variance. use only the low intensity half of the histogram peak
	sumdI2 = 0.0001
	sumFreq = 0.0001 # avoid division by zero 

	for binIdx in range(binIdxMaxFreq): 
		I = hMin + binIdx * (hMax - hMin)/nBins;
		freq = hist[binIdx]
		sumdI2 = sumdI2 + freq * math.pow(I - IMaxFreq, 2)
		sumFreq = sumFreq + freq

	ISigma = math.sqrt(sumdI2/sumFreq)

	# define a threshold
	threshold = IMaxFreq + threshMultiplyer*ISigma

	return threshold, IMaxFreq, ISigma

def getNRegion(mask):
	"""
	Manual implementation of connected components labelling. 

	params: 
	mask: binary image, objects have values 255, bg is 0 

	returns:
	nRegions: no of connected components 
	labelimage: labelmap of the connected components 
	"""
	imp = mask.duplicate()
	width = imp.getWidth();
	
	ip = imp.getProcessor();
	ip = ip.duplicate().convertToShortProcessor()
	ff = FloodFiller(ip);

	# start count from 255 
	count = 255;
	nRegions = 0
	pix_data = ip.getPixels();

	for i in range(0, len(pix_data)):
		# find the first foreground object and its location on image 
		if (pix_data[i]==255):
			y = i/width
			x = i%width

			# increase count by 1, fill the entire region 
			count += 1;
			ip.setValue(count);
			ff.fill8(x,y);

	labelmap = ImagePlus("labelmap", ip)
	IJ.run(labelmap, "glasbey_on_dark", "");

	# as we start from 255
	nRegions = count - 255

	return nRegions, labelmap 


def getFilesFromDir(path, extension=None):
	"""
	Function to get list of files in a directory.

	params:
	path: str data dir path
	extension: str file extension for e.g. '.tif', '.lsm'
	
	Returns: 
	filelist: list of files with entire path 
	"""
	if extension is None:
		filelist = [join(path,f) for f in sorted(listdir(path)) if isfile(join(path, f))]
	else:
		filelist = [join(path,f) for f in sorted(listdir(path)) if isfile(join(path, f)) and f.endswith(extension)]

	return filelist

def getUserPoints(impt, roi_array):
	"""
	Function which allows user to add more points 
	to the spot detection. 
	params: z-stack at timeframe= t, spots roi array 
	returns: array of 2D points added by user 
	"""
	imp = impt.duplicate() 
	spot_ov = Overlay()

	for rr in roi_array:
		rr.setPosition(rr.getZPosition())
		spot_ov.add(rr)

	imp.setOverlay(spot_ov);
	imp.show(); 

	IJ.setTool("multipoint");
	gd = NonBlockingGenericDialog("Action required")
	gd.addMessage("Select additional points and click ok")
	gd.showDialog()

	if gd.wasCanceled():
		imp.close()
		return 

	else:
		points = imp.getRoi()
	
		try:
			return points
		finally:
			imp.close()
			IJ.setTool("rectangle");

def getRois(labelimp):
	"""
	Converts a label map into ROIs and adds it to ROI manager.

	params: label map 
	returns: Array of rois from roi manager
	"""
	LabelMapToRoiManagerPlugin.apply(labelimp)
	roiarray = rm.getRoisAsArray();
	
	return roiarray
	
	
def runStardistGetLabel(imp):
	"""
	Run Stardist on an imp using default params
	adapted from: https://gist.github.com/maweigert/8dd6ef139e1cd37b2307b35fb50dee4a

	params: imp 
	returns: label map 
	"""

	res = command.run(StarDist2D, False,
		"input", imp, "modelChoice", "Versatile (fluorescent nuclei)",
		"outputType", "Label Image", "nTiles", 1, "excludeBoundary", 0).get()

	label = res.getOutput("label")
	labelimp = getImpfromImg(label, "label_" + os.path.basename(imp.getTitle()))
	applyGlasbeyLUT(labelimp)
	IJ.run(labelimp, "Enhance Contrast", "saturated=0.35");
	
	return labelimp


def getImpfromImg(img, title="image"):
	"""
	Convert IJ2 img into IJ1 imp
	
	params: IJ2 img, image title 
	returns: imp
	"""
	imp = ImageJFunctions.wrap(img, title);
	
	return imp 

def applyGlasbeyLUT(imp):
	"""Applies Glasbey on Dark LUT on label maps"""
	IJ.run(imp, "glasbey_on_dark", "");

def doMultiStackReg(imp0):
	"""
	Multistackreg for 2 channel- timelapse. 
	Align ch1 (translation). Use ch1 t=1 as reference
	Align ch2 to ch1 
	 
	params: imp 2 channel- timelapse
	returns: drift-corrected composite imp 
	"""
	imp = imp0.duplicate()
	
	imps = ChannelSplitter().split(imp)
	impch1 = imps[0]
	impch2 = imps[1]
	
	impch1_name = impch1.getTitle()
	impch2_name = impch2.getTitle()

	impch1.setT(1);
	impch1.show()
	impch2.show()

	# set 1st-frame as reference frame to register to 
	impch1.setT(1);

	# run multistackreg: align ch1 with translation and use it as reference to register ch2
	IJ.log("Starting drift correction...this might take time! Do not close the images!")
	IJ.run("MultiStackReg", "stack_1="+ impch1_name +" action_1=Align file_1=[] stack_2="+ impch2_name +" action_2=[Align to First Stack] file_2=[] transformation=Translation");
	
	IJ.run("Merge Channels...", "c1="+ impch1_name +" c2="+ impch2_name +" create");
	composite = IJ.getImage();
	composite.setTitle("drift corrected - " + imp0.getTitle())
	IJ.log("Done with drift correction!")

	return composite

def getIntensityMeasures(mask, imp):
	"""
	Get intensity measurements using mask
	
	params: imp, mask
	returns: arrays - mean, total, min, max intensity of all labels (from imp)
	"""
	label = BinaryImages.componentsLabeling(mask, 4, 16)
	label_ids = LabelImages.findAllLabels(label)	

	label.show()

	# preprocess before intensity measurements
#	IJ.run(imp, "Subtract Background...", "rolling=50 stack");
	im = IntensityMeasures(imp, label)
	
	rt = im.getMean()
	mean = rt.getColumn(0)

	rt = im.getMedian()
	median = rt.getColumn(0)
			
	rt = im.getMin() 
	minI= rt.getColumn(0)
			
	rt = im.getMax() 
	maxI = rt.getColumn(0)

	# this is the rawIntDen from imageJ - tested manually
	rt = im.getSumOfVoxels() 
	total = rt.getColumn(0)

	rt.reset();

	return label_ids, mean, median, minI, maxI, total 


def getGeometricalMeasurements(mask, imp):
	"""
	Get area measurements using mask
	
	params: imp, mask
	returns: arrays - area
	"""
	label = BinaryImages.componentsLabeling(mask, 4, 16)
	label_ids = LabelImages.findAllLabels(label)

	cal = imp.getCalibration();
	areas_array = IntrinsicVolumes2D.areas(label.getProcessor(), label_ids, cal)

	return areas_array

def getWekaSegmentation(imp_process, modelPath):
	"""
	Apply weka model to image

	returns: probability map of object class 
	"""
	weka = WekaSegmentation(imp_process) 
	getProbs = True 
	weka.loadClassifier(modelPath.getPath())
	probmap = weka.applyClassifier(imp_process, 0, getProbs)

	# object class = 1 
	# background class = 2
	obj_class = Duplicator().run(probmap, 1, 1, 1, 1, 1, 1)

	return obj_class

def getMaskfromProbMap(probability_map, prob=0.9):
	"""
	Threshold probability image to get 8-bit mask 
	Threshold value = prob 

	params: prob map 32 bit, threshold value 
	returns: 8-bit mask 
	"""
	impp = probability_map.duplicate()
	impp.show()
	title = impp.getTitle() 
	IJ.selectWindow(title);
	
	IJ.setThreshold(prob, 1.0000);
	IJ.run("Convert to Mask"); 
	impp.hide()

	return impp

def getNeighborLabels(src_label, target_label, distance_threshold= 5.0):
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

	ids_c1 = []
	ids_c2 = []
	used = []
	
	# for each in src find neighbor in target
	for i1, (c1, l1) in enumerate(zip(centroids_1, labels1)):
		close = closest(c1, centroids_2, distance_threshold)		
		
		if close is not None:
			idx = centroids_2.index(close)
			if labels2[idx] not in used:
				ids_c2.append(labels2[idx])
				ids_c1.append(l1)
				used.append(labels2[idx])

	neighbor_src = LabelImages.keepLabels(src_label, ids_c1)
	neighbor_target = LabelImages.keepLabels(target_label, ids_c2)

	return neighbor_src, neighbor_target, ids_c1, ids_c2

def getNoNeighbors_new(src_label, target_label, distance_threshold= 5.0):
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

def scaleImp(imp, factor = 0.5):
	"""
	scale image by factor defined in factor 
	"""
	height = imp.getHeight(); 
	width = imp.getWidth(); 

	newHeight = height * factor 
	newWidth = width * factor 

	imp_resized = imp.resize(int(newWidth), int(newHeight), "Bilinear");
	imp_resized.setCalibration(imp.getCalibration());

	imp_resized.setTitle("downscaled_" + imp.getTitle())

	return imp_resized

def getOverlapLabelUsingJaccardIndex(labelImage1, labelImage2):
	"""
	Use Jaccard index to find overlapping labels in 2 label images. 
	For overlapping labels Jaccard index = 1 

	params: label image 1, label image 2 
	returns: label image with labels from label image 1 which have overlap in label image 2
	"""

	rt = ResultsTable() 
	rt.reset()

	# get jaccard index 
	rt = LabelImages.getJaccardIndexPerLabel(labelImage1, labelImage2)
	j_index = rt.getColumn(0)

	labels1 = LabelImages.findAllLabels(labelImage1)
	labels_keep = []

	# keep labels for which jaccard index = 1 i.e. which have overlap in label image 2
	for i, j in zip(labels1, j_index):
		if j != 0.0:
			labels_keep.append(i)

	overlap = LabelImages.keepLabels(labelImage1, labels_keep)

	return overlap, labels_keep

def doNucleiSpotsOverallMeasurement(impnuc_sumT, impspot_sumT, nuclei_roi, impName, cal):
	"""
	Main function to get overall measurements. 
	Spots are segmented in the image impspot_sumT inside the nuclei roi and measurements are obtained for these spots. 
	Results are populated into the overall measurement table. 
	The spots label image is saved. 

	Note: the number of spots could vary from the tracking output; this measure is just meant to get an approximate idea 
	of ratio of spots vs nucleoplasm
	
	params: imp nuclei and spots channel sum-projected over time (2D + t), nuclei roi, image title, calibration 
	returns: None
	"""
	impnuc_sumT.setRoi(nuclei_roi)
	IJ.run(impnuc_sumT, "Clear Outside", "stack"); 

	impspot_sumT.setRoi(nuclei_roi)
	IJ.run(impspot_sumT, "Clear Outside", "stack"); 

	# get nuclei stats
	nuclei_area, nuclei_mean, nuclei_total = getNucleiMeasurements(impnuc_sumT, nuclei_roi)  

	# get spots roi and label 
	spots_label, spots_roi = getSpotsRoiAndLabel(impspot_sumT, nuclei_roi)

	# spots stats 
	labels, meanI, minI, maxI, totalI, diameters, areas = getSpotsMeasurements(impspot_sumT, spots_label, cal)
	
	# table 1 for overall measurements 
	table_overall.incrementCounter()
	table_overall.addValue("Image", impName) 
	table_overall.addValue("Nuclei area (" + cal.getUnits() + "^2)", nuclei_area)
	table_overall.addValue("Nuclei mean", nuclei_mean)
	table_overall.addValue("Nuclei total", nuclei_total)
	table_overall.addValue("No spots (approx)", len(labels))
	table_overall.addValue("Spots area (" + cal.getUnits() + "^2)", sum(areas))
	table_overall.addValue("Spots mean", sum(meanI))
	table_overall.addValue("Spots total", sum(totalI))
	table_overall.addValue("Nucleoplasm (nuclei - spots) total", nuclei_total - sum(totalI))
	table_overall.addValue("Spots total/ nucleoplasm total", sum(totalI) / (nuclei_total - sum(totalI)))

	if displayimage:
		spots_label.show()

	FileSaver(spots_label).saveAsTiff(join(results_path.getPath(), impName + "-spots_label.tif"))

def segmentNuclei(imp0):
	"""
	Segment the nuclei from the nuclei channel imp 
	params: imp with nuclei 
	returns: single nucleus label map and selection
	"""
	imp_ = imp0.duplicate() 
	imp = ZProjector().run(imp_, "max all");
	
	# get mask of nuclei
	IJ.run(imp, "Gaussian Blur...", "sigma=3");
	IJ.setAutoThreshold(imp, "Huang" + " dark");
	IJ.run(imp, "Convert to Mask", "");
	IJ.run(imp, "Watershed", ""); 

	# label map
	label0 = BinaryImages.componentsLabeling(imp, 4, 16)
	label = LabelImages.sizeOpening(label0, 200)
	labels = LabelImages.findAllLabels(label)

	# ask user to keep only one label for processing
	while len(labels) > 1:
		label.show()
		IJ.setTool("multi-point");
		WaitForUserDialog("Action required", "KEEP A SINGLE nucleus to process. Please click on the other nuclei to DISCARD and press ok.").show()
	
		roi = label.getRoi()
		LabelImages.removeLabels(label, roi, True)
		IJ.run(label, "glasbey_on_dark", ""); 
		labels = LabelImages.findAllLabels(label)
		
	label.killRoi()

	# remap labels so that the remaining label has value = 1
	LabelImages.remapLabels(label)
	label.hide()

	# convert the single nucleus label map into mask again
	mask = label.duplicate()
	LabelImages.replaceLabels(mask, [1], 255)
	IJ.run(mask, "8-bit", "");
	IJ.run(mask, "Grays", "");
#	mask.show()
	
	IJ.run(mask, "Create Selection", "");
	nuclei_roi = mask.getRoi()

	return label, nuclei_roi

def getNucleiMeasurements(imp, nuclei_roi):
	"""
	Get area, mean and total intensity of nuclei from roi drawn by user 
	params: imp, roi 
	returns: area, mean and total intensity of the roi in imp 
	"""

	# preprocess before measurement - recommended 
	IJ.run(imp, "Subtract Background...", "rolling=50 stack");

	imp.setRoi(nuclei_roi);
	stats = imp.getStatistics(Measurements.ALL_STATS)
	
	area = stats.area 
	meanI = stats.mean 

	# this is the rawIntDen measurement from imageJ 
	# check: https://forum.image.sc/t/intden-vs-rawintden/5147 - nice comparison
	totalI = float(stats.pixelCount * stats.mean)

	return area, meanI, totalI 

def getSpotsRoiAndLabel(imp0, nuclei_roi):
	"""
	Get spots label, and spots roi 
	params: imp (2D)
	returns: spots label map, spots roi (selection from mask of all spots)
	"""
	spots_mask = imp0.duplicate()
	spots_mask.setRoi(nuclei_roi)

	# get spots mask 
	IJ.run(spots_mask, "Gaussian Blur...", "sigma=1");
	IJ.setAutoThreshold(spots_mask, "Li" + " dark");
	IJ.run(spots_mask, "Convert to Mask", "");
	IJ.run(spots_mask, "Watershed", "");
#	spots_mask.show()

	# get spots roi
	spots_mask.killRoi()
	IJ.run(spots_mask, "Create Selection", "");
	spots_roi = spots_mask.getRoi()
	spots_roi.setStrokeColor(Color.red)
	spots_mask.killRoi()

	# get spots label map
	spots_label = BinaryImages.componentsLabeling(spots_mask, 4, 16)
	applyGlasbeyLUT(spots_label)
	spots_label.setTitle(os.path.splitext(imp0.getTitle())[0] + "_spotsLabels")
#	spots_label.show()

	return spots_label, spots_roi


def getSpotsMeasurements(imp, label, calibration):
	"""
	Get spots measurements using labels 
	
	params: imp, label imp, calibration 
	returns: arrays - area, diameter (from label) and mean, total, min, max intensity of all labels (from imp)
	"""
	label_ids = LabelImages.findAllLabels(label)	

	# max-feret diameter
	pairs = MaxFeretDiameter().analyzeRegions(label.getProcessor(), label_ids, calibration)
	diameters = [p.diameter() for p in pairs]
	
	areas = IntrinsicVolumes2D.areas(label.getProcessor(), label_ids, calibration)

	# preprocess before intensity measurements
	IJ.run(imp, "Subtract Background...", "rolling=50");
	im = IntensityMeasures(imp, label)
	
	rt = im.getMean()
	mean = rt.getColumn(0)
			
	rt = im.getMin() 
	minI= rt.getColumn(0)
			
	rt = im.getMax() 
	maxI = rt.getColumn(0)

	# this is the rawIntDen from imageJ - tested manually
	rt = im.getSumOfVoxels() 
	total = rt.getColumn(0)

	rt.reset();

	return label_ids, mean, minI, maxI, total, diameters, areas 

def getRoisByParticleAnalysis(mask, measures, minsize, maxsize, mincirc, maxcirc):
	"""
	Function to get array of rois from input mask.
	Input: mask imp, measurement options, min and max size, min and max circ
	Returns: array of rois
	"""
	roi_manager = RoiManager.getInstance()
	if not roi_manager:
		roi_manager = RoiManager()
		
	options = PA.ADD_TO_MANAGER + PA.SHOW_NONE 
	rtt= ResultsTable()
	
	pa = PA(options, measures, rtt, minsize, maxsize, mincirc, maxcirc)
	pa.setHideOutputImage(True)
	pa.analyze(mask)
	rois_array = roi_manager.getRoisAsArray();
	roi_manager.reset()

	return rois_array

def removeAdjacentLabelsAndReturnMask(labelmap):
	"""
	Function to discard labels touching each other from the label map
	and return a mask of the filtered label map.
	"""
	pairs = RegionAdjacencyGraph.computeAdjacencies(labelmap)
	new = Duplicator().run(labelmap)

	for p in pairs:
		rm_labels = [p.label1, p.label2]
		LabelImages.replaceLabels(new.getProcessor(), rm_labels, 0.0)

	ipb = LabelImages.labelBoundaries(new.getProcessor())
	mask = ImagePlus("bin", ipb)
	IJ.run(mask, "Fill Holes", "");
	for i in range(1, 4):
		IJ.run(mask, "Dilate", "");
		
	return mask

def getminmaxnormalized(imp):
	'''
	Every intensity is subtracted with min intensity in image and divided by maxI - minI in image.
	In the end, we multiply image again with 255 to preserve range between 0-255.
	'''
	IJ.run(imp, "32-bit", "");
	
	stats = imp.getStatistics(Measurements.ALL_STATS)
	minI = stats.min
	maxI = stats.max
	Irange = maxI-minI

	img = ImageJFunctions.wrap(imp)
	img_norm = ops.math().subtract(img, minI)
	img_norm = ops.math().divide(img, Irange)
	img_norm = ops.math().multiply(img, 255)

	imp_norm = ImageJFunctions.wrap(img_norm, "minmax_norm - " + os.path.splitext(imp.getTitle())[0])
	imp_norm.setCalibration( imp.getCalibration() )
	imp_norm.setTitle("minmax_norm - " + os.path.splitext(imp.getTitle())[0])

	return imp_norm

def getminmaxnormalized_withpercentiles(imp, minI, maxI):
	'''
	Similar to min max but instead of using image min and max, we use the min and max percentiles of image histogram.
	The percentiles of values set by user is taken example 5% and 95% or 2% and 98%.
	'''
	IJ.run(imp, "32-bit", "");

	Irange = maxI-minI

	img = ImageJFunctions.wrap(imp)
	img_norm = ops.math().subtract(img, minI)
	img_norm = ops.math().divide(img, Irange)
	img_norm = ops.math().multiply(img, 255)

	imp_norm = ImageJFunctions.wrap(img_norm, "minmaxwithperc_norm - " + os.path.splitext(imp.getTitle())[0])
	imp_norm.setCalibration( imp.getCalibration() )
	imp_norm.setTitle("minmaxwithperc_norm - " + os.path.splitext(imp.getTitle())[0])

	IJ.run(imp_norm, "Enhance Contrast...", "saturated=0.1 normalize");
	IJ.run(imp_norm, "Multiply...", "value=255.000");

	return imp_norm

def normalizeWithPercentiles( imp0, scale, offset ):
	'''
	Scale and offset are computed by using the min and max percentile values.
	New image is obtained by multiplying with scale and adding the offset.
	'''
	IJ.run(imp0, "32-bit", "");
	imp = Duplicator().run(imp0)
	img = ImageJFunctions.wrap(imp)

	# scale the image by multiplication and add the offset
	img = ops.math().multiply(img, scale)
	img = ops.math().add(img, offset)

	imp1 = ImageJFunctions.wrap(img, "percentile_norm - " + os.path.splitext(imp0.getTitle())[0])
	imp1 = Duplicator().run(imp1)
	imp1.setCalibration( imp0.getCalibration() )
	imp1.setTitle("percentile_norm - " + os.path.splitext(imp0.getTitle())[0])

	IJ.run(imp1, "Enhance Contrast...", "saturated=0.1 normalize");
	IJ.run(imp1, "Multiply...", "value=255.000");

	return imp1

def percentile(imp, fractions):
	# fractions = [5.0, 95.0]

	# get stats for stack or single imp
	if imp.getStackSize()>1 :
		stats = StackStatistics(imp) 
		#stats = imp.getStack().getProcessor(1).getStatistics()
	else:
		stats = imp.getProcessor( ).getStatistics()

	# get histogram and its min, max, bins
	hist = stats.getHistogram()
	
	hMin = stats.histMin
	hMax = stats.histMax
	nBins = len(hist);

	# get sum of all values in histogram
	histSum = [hist[0]]
	count=1
	for val in hist[1:] :
		histSum.append( hist[count] + histSum[count-1] )
		count = count+1

	# get max
	xmax= float(histSum[-1])

	# get percentages instead of values 
	histSum = [x/xmax*100 for x in histSum]

	percentiles = []
	count=0

	# f is min of fractions i.e. 5%
	f = fractions[0]
	for i,f1 in enumerate(histSum) :

		# if % is greater than 5%
		if f1 >= f :
			# get % at i - 1
			f0 = histSum[i-1]
			ifrac = ( f - f0 ) / (f1-f0) + i - 1

			# get actual value and append to empty array
			percentiles.append( hMin + (hMax-hMin)*ifrac/nBins )

			count = count+1
			# restrict count to just len(fractions) so that we have just two values, the min and max wrt percentiles
			if count >= len(fractions) :
				break

			# change the value of f i.e. now it is 95%
			f = fractions[count] 
			
	return percentiles;