# ============ README ================= 
# 
# Author: Gayathri Nadar, SCF MPI CBG nadar@mpi-cbg.de 
# Date: 20.07.2021 
#
# Description: 
#	This script applies a trained WEKA model on a set of images 2D SINGLE CHANNEL in a folder and generates a probability map. 
#	The probability map is threshold to get a mask of the foreground objects. 
#	Lastly, we measure area and mean intensity from the objects detected. 
#
#
# Usage:
#	Open FIJI and open an image you want to process. 
#	Launch Trainable Weka Segmentation plugin and train a classifier. 
#	Save the classifier file using the Save Classifier option. This will save the .model file.
#	Save the other images you want to apply the classifier to in a folder.
#	Open this script in FIJI and click run.
#	Specify the image folder to process and the choose the .model file. 
#	Click OK. 
#
#	Output: Results table
#
# Note: 
#	The following plugins are needed to run this script:
#	1. Bio-formats 
#	2. MorpholibJ (plugin name: IJPB plugins)
#	3. LOCI
#
#	In case you do not have these plugins then please install it by going to Help > Update...
#	- Then in window that comes up, click on "Manage Update Sites". Click the box next to mentioned name. 
#	- Click "close" and then "Apply changes".
#
# 
# =============================================================================================================



# =============================================================================================================
# Don't modify anything below here without good reason
# =============================================================================================================


#@ File (label= "Select image folder to batch-process", style= "directory", required= True) processFolder 
#@ File(label="Weka model", description="Select the Weka model to apply", required= True) modelPath

from ij import IJ
from ij import ImagePlus
from ij.plugin import Duplicator
from ij.measure import Measurements, ResultsTable
import loci.plugins
from loci.plugins import BF
from loci.plugins.in import ImporterOptions 
from loci.formats import ImageReader
import os
from os.path import isfile, join
from trainableSegmentation import WekaSegmentation
from trainableSegmentation.utils import Utils;
from inra.ijpb.binary import BinaryImages
from inra.ijpb.label import LabelImages
from inra.ijpb.measure import IntensityMeasures 
from inra.ijpb.measure import IntrinsicVolumes2D

# specify file extension
ext = ".tif"

def main():
	# find all files with tif
	imgFolder = processFolder.getPath()
	files = sorted([ f for f in os.listdir(imgFolder) if isfile(join(imgFolder, f)) and f.endswith(ext)])

	table = ResultsTable()

	for f in files: 
		imp = openImageWithBF(join(imgFolder, f), False, False) 
		imp.show();
		impname = os.path.splitext(imp.getTitle())[0]
		
		# apply trained weka model 
		print "Applying weka model..this might take time!"
		imp_probabilitymap = getWekaSegmentation(imp, modelPath) 

		# threshold prob map to get mask 
		# probability map has values between 0 and 1 
		# using threshold we set all values above the threshold to 255 and those below to 0 
		# output is a binary mask of the objects of interest
		imp_mask = getMaskfromProbMap(imp_probabilitymap, prob=0.5)

		# process mask to get connected components and get measurements for example using particle analysis 
		# here we use MorpholibJ plugin 
		# intensity measurements for all objects from original image
		label_ids, mean, median, minI, maxI, total = getIntensityMeasures(imp_mask, imp)       # mean, min, etc. are all lists (array)

		# geometry measurements: area 
		areas = getGeometricalMeasurements(imp_mask, imp)

		# add to table 
		# populate tables 
		for idx, l in enumerate(label_ids):
			table.incrementCounter(); 
			table.addValue("Imagename", impname)
			table.addValue("Obj ID", l)
			table.addValue("Mean", mean[idx])
			table.addValue("Total", total[idx]) 
			table.addValue("Area", areas[idx])

	table.show("Results:All Images")

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

	try:
		return imp 
	finally:
		imp.close();
		del imp 


main() 
print "Done"
		