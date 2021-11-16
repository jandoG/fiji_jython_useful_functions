# ============================================= README ======================================================
#
# Author: Gayathri Nadar, SCF MPI-CBG, nadar@mpi-cbg.de
#
# Version history:
# 04.02.2021 v0.0: segment droplets and measure size  
# 19.04.2021 v0.1: remove overlay from drift corrected image before saving
#				   filter droplets by size before adding to results table
#
# Description:
#	This script measures size and intensity of droplets. Data expected: 2D time lapse with 2 channels ch1: TAMRA ch2: FAM. All measurements done from FAM channel.
#
# Steps:
#	For each file in the list:
#	- split channels, correct for drift using multistackreg.  
#	- max-project across time ch1 image and use stardist 2D to get labels. Remove labels touching boundaries. 
#	- get rois from label map
#	- set roi, statistics - area, diameter, mean intensity 
#
# Note: 
#	The following plugins are needed to run this script:
# 	1. SCF MPI CBG 
#	2. MorpholibJ (plugin name: IJPB plugins)
#	3. CSBDeep and StarDist 
#	4. Big-EPFL 
#	
#	In case you do not have these plugins then please install it by going to Help > Update...
#	- Then in window that comes up, click on "Manage Update Sites". Click the box next to mentioned names. 
#	- Click "close" and then "Apply changes".
#
#   5. MultiStackReg: download MultiStackReg (https://github.com/miura/MultiStackRegistration/releases)
#	Copy paste it into the plugins folder in FIJI ---> Right click FIJI > Show Package contents > plugins
#      
# Usage:
#	Arrange your data into a folder for processing. 
# 	Open this script with FIJI and click run.
# 	A window will pop-up, enter the necessary values and click ok to run the script.
#
# Results:
# 	Results are stored in a folder called "_results" within the image folder. 
# 	Results saved:
#	1. CSV file with measurements from all the files selected to process
# 	2. Drift corrected image: no overlay 
#	3. Label image with droplet IDs of filtered droplets
#
# Input parameters:
#	- Select 1 or multiple images to process: Select a single file or multiple files to process. 
#	  Note on "Add folder content":
#		- This will add all the .tif from the chosen directory including from sub-folders excluding the _results folder. So make sure to have your data well organized.  
#	
#   - The diameter of droplets that need to be filtered. We select all droplets within range diameter - + 10
#
# =============================================================================================================



# =============================================================================================================
# Don't modify anything below here without good reason
# =============================================================================================================


#@ String (value = "v0.1", visibility = "MESSAGE", required = false) label
#@ File[] (label= "Select 1 or more images", style="files") listfiles
#@ Integer (label= "Approximate diamater (pixels) of droplets to be filtered out?", value = 30) filter_diameter
#@ CommandService command 

from de.csbdresden.stardist import StarDist2D
import os, os.path 
import glob, math
from glob import glob 
from os.path import isfile, join, basename
from ij import ImagePlus, IJ
from ij.io import FileSaver
from ij.plugin import Duplicator, ChannelSplitter, ZProjector
from ij.plugin.frame import RoiManager
from ij.measure import Measurements, ResultsTable
from ij.gui import Overlay, Roi
from loci.plugins import BF
from loci.formats import ImageReader
from loci.formats import MetadataTools
from loci.plugins.prefs import OptionsList
from loci.plugins.in import ImporterOptions
from java.awt import Color
from ij.plugin.filter import Analyzer
from de.mpicbg.scf.fijiplugins.ui.roi import LabelMapToRoiManagerPlugin;
from de.mpicbg.scf.fijiplugins.ui.labelmap import ThresholdLabelingPlugin
from net.imglib2.img.display.imagej import ImageJFunctions;
from inra.ijpb.label import LabelImages
from inra.ijpb.binary import BinaryImages
from inra.ijpb.measure.region2d import Centroid
from inra.ijpb.measure import IntensityMeasures
from inra.ijpb.measure.region2d import BoundingBox

rm = RoiManager.getInstance()
if not rm:
	rm = RoiManager()
rm.reset();
table = ResultsTable() 
ext = '.tif'

def main():
	# get list of files ending with ext defined 
	sorted_list = sorted([f.getPath() for f in listfiles if f.isFile() and f.getPath().endswith(ext) and not os.path.dirname(f.getPath()).endswith("_results")])

	# create dir to save results 
	results_path = join(os.path.dirname(sorted_list[0]), os.path.basename(os.path.dirname(sorted_list[0])) + '_results')
	if not os.path.exists(results_path):
		os.makedirs(results_path)

	# process each file in the list
	for imgfile in sorted_list:
		imp = openImageWithBF(imgfile, virtual= False)
		impName = os.path.splitext(imp.getTitle())[0]
		print "Processing image =", impName 

		# get drift corrected imp, and stardist output from ch1 
		corrected_imp, ch1, ch2, ch1_label = processImp(imp)
		label_ids = LabelImages.findAllLabels(ch1_label)
#		ch1_label.show()

		# get rois from label map 
		LabelMapToRoiManagerPlugin.apply(ch1_label)
		roi_array = rm.getRoisAsArray()

		cal = imp.getCalibration(); 
		overlay = Overlay()

		# empty list to store filtered ids
		keeplabels = [] 
		
		for t in range(1, ch2.getNFrames() + 1):
			ch2.setT(t)
			ip = ch2.getProcessor()
			for labelid, roi in zip(label_ids, roi_array):
				# give a name to roi labels 
				idx = rm.getRoiIndex(roi)
				rm.rename(idx, "droplet-" + str(labelid));

				# get measurements 
				ch2.setRoi(roi)
				stats = ch2.getStatistics(Measurements.ALL_STATS)
				diameter = float((stats.minor + stats.major)/2) / cal.pixelWidth   # in pixels 

				# check if the measured diameter is within range  
				# continue measurements on only the filtered droplets
				if (filter_diameter - 10) <= diameter and diameter <= (filter_diameter + 10):
					overlay.add(roi);
					keeplabels.append(labelid)

					table.incrementCounter();
					table.addValue("image", impName);
					table.addValue("frame no", t);
					table.addValue("droplet id", rm.getName(idx));
					table.addValue("measurement ch", 2)
					table.addValue("area (" + cal.getUnits() + "^2)", stats.area);
					table.addValue("diamater (" + cal.getUnits() + ")", float((stats.minor + stats.major)/2));
					table.addValue("mean intensity", stats.mean);

#		corrected_imp.setOverlay(overlay) # overlay not required 16.04.2021
		corrected_imp.setDisplayMode(IJ.COLOR)
		corrected_imp.show()
		outputName = 'drift_corrected_' + impName + '.tif'
		FileSaver(corrected_imp).saveAsTiff(join(results_path, outputName))

		# save image with droplet IDs
		label = ch1_label.duplicate()
		label.killRoi()
		label_keep = LabelImages.keepLabels(label, list(set(keeplabels)))
		label_keep.setTitle("dropletIDs_filtered_" + impName)
		applyGlasbeyLUT(label_keep)
		label_keep.show()
		FileSaver(label_keep).saveAsTiff(join(results_path, "dropletIDs_filtered_by_size_" + impName + ".tif"))

		# clear roi manager for next file
		rm.reset();

	table.show("Results: All files")
	tableName = 'Droplet_size_measurements_FAM_channel_all_files.csv'
	table.saveAs(join(results_path, tableName))


def processImp(imp0):
	"""
	Main processing

	returns: ch1, ch2 (drift corrected)
	label imp of ch1
	"""
	imp = imp0.duplicate()

	# drift correction 
	drift_corrected = doMultiStackReg(imp)
#	drift_corrected.show()
	drift_corrected.hide()

	ch1 = Duplicator().run(drift_corrected, 1, 1, 1, 1, 1, imp.getNFrames())
	ch2 = Duplicator().run(drift_corrected, 2, 2, 1, 1, 1, imp.getNFrames())

	# duplicate ch1 and get its label map
	ch1_projected = ZProjector().run(ch1, "max");
	label = runStardistGetLabel(ch1_projected)
	LabelImages.removeBorderLabels(label)
	
	return drift_corrected, ch1, ch2, label

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

def openImageWithBF(path, virtual= False):
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
print "Done"
roiManager = RoiManager.getRoiManager();
roiManager.close();	