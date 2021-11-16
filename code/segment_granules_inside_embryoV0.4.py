################# README #######################

# Author Gayathri Nadar, SCF MPI-CBG, nadar@mpi-cbg.de
# Version 0.0: segment granules inside c-elegans embryo over time and measure different parameters
# Version 0.1: get orientation angle for every time point and rotate it individually slice-wise
# Version 0.2: rotate the embryo, then ask user if the anterior-posterior position is flipped
# Version 0.3: improved preprocessing pipeline for segmenting granules
#				option to save results table
#				option to specify minsize of the granule to be segmented
#				added granules count to measurements
# Version 0.4: instead of creating multiple granule roi, we create one selection of all granules together per slice
#				added background subtraction to the image from which intensity is measured, to avoid adding background intensity to the measurements.

'''
# Description:

- This script loops through time lapse movie of c-elegans embryo, segments the granules (bright spots) and quantifies them.
- Data is originally z-stack which is max-projected to get time lapse image.
- In first step, we segment the embryo and get its orientation by fitting an ellipse. 
- Then the embryo is rotated using the orientation angle such that it is parallel to x-axis.
- We ask the user if the AP position of the embryo is flipped, is yes, we flip the image and proceed to segment granules, if not we process the max-projected, rotated image.
- Next, in a loop, we do:
	-- segment the embryo to get its roi
	-- segment the granules inside this roi
	-- measure mean I, total I, area, count, etc.
- Output:
	-- result image with detected granules 
	-- results table

# Usage:
- Open the script with FIJI
- Click run
- Enter the required parameters in pop up window
		1. file to be processed
		2. option to save results to data folder.
- Click ok to run the script
- Please select yes or no in the pop-up window comes up to ask the user to flip the AP position. If yes, image will be flipped if not, script continues to measure in the original rotated image.

'''
################################################

#@File (label= "Image to process") processFile
#Float (label= "min size of granule in radius^2 (um)") minsize
#@boolean (label="Save results as csv (checked = yes, unchecked = no)") saveResults

from ij import IJ
from ij import ImagePlus
from ij import ImageStack
from ij import WindowManager
from ij.text import TextWindow
from ij.plugin import Duplicator, FolderOpener
from ij.plugin import ZProjector
from ij.plugin.filter import ParticleAnalyzer as PA
from ij.plugin.filter import Analyzer
from ij.plugin.filter import MaximumFinder
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable, Measurements, Calibration
from ij.gui import ShapeRoi, Roi, Overlay, PolygonRoi, WaitForUserDialog, YesNoCancelDialog
from ij.process import StackStatistics
from ij.plugin.filter import ThresholdToSelection;
from ij.io import FileSaver
from java.awt import Color, Frame
import loci.plugins
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import os, math
from os.path import isfile, join

global ov
global thresholding_method
thresholding_method = 'Li'
ov = Overlay()

roi_manager = RoiManager.getInstance()
if not roi_manager:
	roi_manager = RoiManager()

def main():
	imgFile= processFile.getPath()
	imp= BF.openImagePlus(imgFile)[0]

	global info
	global imp_name
	info = imp.getOriginalFileInfo()
	imp_name = (os.path.splitext(imp.getTitle())[0])
	
	imp1 = Duplicator().run(imp)
	imp_max_proj = ZProjector.run(imp1,"max all");
#	imp_max_proj.show()
	
	rotatedimp = rotatestack_at_every_timepoint(imp_max_proj)
	imp_embryo1 = Duplicator().run(rotatedimp, 1, 1, 1, 1, 5, 5)
	IJ.run(imp_embryo1, "Enhance Contrast", "saturated=0.3");
	imp_embryo1.show()

	dialogbox = Frame()
	isflipped = YesNoCancelDialog(dialogbox, "AP axis", "Is anterior - posterior flipped? Click yes or no..") 
	imp_embryo1.close()
	if not isflipped.cancelPressed() and isflipped.yesPressed():
		IJ.run(rotatedimp, "Flip Horizontally", "stack");
		segmentgranules(rotatedimp)	
	elif isflipped.cancelPressed():
		IJ.log("Please select yes or no!")
		isflipped.show()
		if isflipped.yesPressed():
			IJ.run(rotatedimp, "Flip Horizontally", "stack");
#			rotatedimp.show()
			segmentgranules(rotatedimp)	
		else:
			print "No, there is no flipping."
			segmentgranules(rotatedimp)
	else:
		print "No flip!"
		segmentgranules(rotatedimp)


def getEmbryoOrientation(imp):
	imp_embryo = Duplicator().run(imp, 1, 1, 1, 1, 1, 1)
	embryo_roi = getembryoroi(imp_embryo)
	imp_embryo.setRoi(embryo_roi)
	stats = imp_embryo.getStatistics(Measurements.ALL_STATS)
	angle = stats.angle

	return angle

def getembryoroi(imp_embryo):
	IJ.run(imp_embryo, "Gaussian Blur...", "sigma=4");
	IJ.setAutoThreshold(imp_embryo, thresholding_method + " dark");
	IJ.run(imp_embryo, "Convert to Mask", "");
	IJ.run(imp_embryo, "Median...", "radius=15");
	IJ.run(imp_embryo, "Fill Holes", "");
	IJ.run(imp_embryo, "Create Selection", "");
	embryo_roi = imp_embryo.getRoi()

	return embryo_roi

def rotatestack_at_every_timepoint(imp_timelapse):
	# get angle and rotate at every time!
	imp = Duplicator().run(imp_timelapse)
	frames = imp_timelapse.getNFrames()
	for t in range(1,frames +1):
		imp.setT(t)
		angle = getEmbryoOrientation(Duplicator().run(imp, 1, 1, 1, 1, t, t))
		IJ.run(imp, "Rotate... ", "angle="+str(angle)+" grid=1 interpolation=Bilinear slice");

	return imp

def segmentgranules(embryo_stack):
#	embryo_stack.show()
	frames = embryo_stack.getNFrames()
	print "No of frames= ", frames
	imp_result = Duplicator().run(embryo_stack)
	IJ.run(imp_result, "Subtract Background...", "rolling=100 stack");

	resultsTableName = "Granules"
	rt = getResultsTable(resultsTableName)

	for t in range(1, frames + 1):
		embryo_stack.setT(t)

		# get embryo roi
		imp_proc = Duplicator().run(embryo_stack, 1, 1, 1, 1, t, t)
		roi = getembryoroi(Duplicator().run(imp_proc))
		roi.setStrokeColor(Color.yellow)
		roi.setImage(imp_result)
		roi.setPosition(t)
		ov.add(roi)

		# add background subtraction before measurements
		IJ.run(imp_proc, "Subtract Background...", "rolling=100 slice");
		
		imp_proc.setRoi(roi)
		embryo_stats = imp_proc.getStatistics(Measurements.ALL_STATS)
		embryo_total_intensity = float(embryo_stats.pixelCount * embryo_stats.mean)
		embryo_integrated_intensity = float(embryo_stats.area * embryo_stats.mean)
		embryo_area = embryo_stats.area
		imp_proc.killRoi()

		granules_mask = Duplicator().run(imp_proc)
		granules_mask.setRoi(roi)

		# test 1 = >90%
		IJ.setAutoThreshold(granules_mask, "Triangle" + " dark");
		IJ.run(granules_mask, "Convert to Mask", "");
		IJ.run(granules_mask, "Median...", "radius=1");
		IJ.run(granules_mask, "Watershed", "");
#		granules_mask.show()

		# maxima detection to get granules count
		ip= granules_mask.getProcessor()
		maxPolygon = MaximumFinder().getMaxima(ip, 0, False)
		spots_count = maxPolygon.npoints
		print "No of granules =", spots_count

		# create selection, get roi, set selection, measure, add to table
		IJ.run(granules_mask, "Create Selection", "");
		granules_roi = granules_mask.getRoi()
		granules_roi.setStrokeColor(Color.red)
		granules_roi.setImage(imp_result)
		granules_roi.setPosition(t)
		ov.add(granules_roi)
		granules_mask.killRoi()

		imp_proc.setRoi(granules_roi)
		stats= imp_proc.getStatistics(Measurements.ALL_STATS)
		spot_total_raw_intensity = float(stats.pixelCount * stats.mean)
		spot_integrated_intensity = float(stats.area * stats.mean)
		spots_area = stats.area

		volume_fraction_granules_to_embryo = spots_area / embryo_area
		volume_fraction_embryo_to_granules = embryo_area / spots_area
		intensity_fraction_granules_to_embryo = spot_total_raw_intensity / embryo_total_intensity
		intensity_fraction_embryo_to_granules = embryo_total_intensity / spot_total_raw_intensity

		rt.incrementCounter()
		rt.addValue('Frame', t)
		rt.addValue('No of granules', spots_count)
		rt.addValue('Area of granules', spots_area)
		rt.addValue('Area of embryo', embryo_area)
		rt.addValue('Intensity fraction granules/embryo', intensity_fraction_granules_to_embryo)
		rt.addValue('Intensity fraction embryo/granules', intensity_fraction_embryo_to_granules)
		rt.addValue('Area fraction granules/embryo', volume_fraction_granules_to_embryo)
		rt.addValue('Area fraction embryo/granules', volume_fraction_embryo_to_granules)
		rt.addValue('Embryo sum of all pixels', embryo_total_intensity)
		rt.addValue('Embryo integrated intensity: area x mean', embryo_integrated_intensity)
		rt.addValue('Granules sum of all pixels', spot_total_raw_intensity)
		rt.addValue('Granules integrated intensity: area x mean', spot_integrated_intensity)

	roi_manager.reset()

	imp_result.setOverlay(ov)
	IJ.run(imp_result, "Enhance Contrast", "saturated=0.3");
	imp_result.setTitle("Final result")
	imp_result.show()
	rt.show(resultsTableName)

	if saveResults:
		savePath = join(info.directory, 'results')
		if not os.path.exists(savePath):
			os.makedirs(savePath)
		rtName = imp_name + '_results.csv'
		rt.save(join(savePath, rtName))
		
def getResultsTable(name, reset=False):
	win = WindowManager.getWindow(name)
	if ( (win!=None) & isinstance(win,TextWindow) ) :	
		rt = win.getTextPanel().getResultsTable();
		if reset:
			rt.reset()
	else:
		rt = ResultsTable();
	return rt

main()
print "DONE"
IJ.selectWindow("ROI Manager");  
IJ.run("Close"); 

























		