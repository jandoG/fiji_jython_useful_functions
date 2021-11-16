################# README #######################

# Author Gayathri Nadar, SCF MPI-CBG, nadar@mpi-cbg.de
# V0.0: Script for measuring the degree of polarization of the granules, i.e. at which pole are they accumulated more.
# V0.1: Added additinal measures -
#		1. total I of granules in a region / total I of granules in all regions 
#		2. total area of granules in a region / total area of granules in all regions
#		3. Count in each region

'''
# Description:
- This script loops through time lapse movie of c-elegans embryo, segments the granules (bright spots) and quantifies them.
- Data is originally z-stack which is max-projected to get time lapse image.
- In first step, we segment the embryo and get its orientation by fitting an ellipse. 
- Then the embryo is rotated using the orientation angle such that it is parallel to x-axis.
- We ask the user if the AP position of the embryo is flipped, is yes, we flip the image and proceed to segment granules, if not we process the max-projected, rotated image.
- Next, in a loop, we do:
	-- segment the embryo to get its roi
	-- divide its rectangular bounding box into no of regions defined by the user along the length (aka width in IMAGEJ)
	-- segment the granules inside each of these regions in a loop
	-- measure total I and mean I in each of the region and get the ratio
- Output:
	-- result image with detected granules 
	-- results table with volume fraction in each region at each time point 
	
# Usage:
- Open the script with FIJI
- Click run
- Enter the required parameters in pop up window
		1. file to be processed
		2. no of regions the embryo is to be divided into
		3. option to save results 
- Click ok.
- Please select yes or no in the pop-up window comes up to ask the user to flip the AP position. If yes, image will be flipped if not, script continues to measure in the original rotated image.

'''
################################################

#@File (label= "Image to process") processFile
#@Integer (label= "No of divisions of the embryo to compute polarization") no_regions
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
from java.awt import Font
from ij.gui import TextRoi
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
	IJ.run(imp_embryo1, "Enhance Contrast", "saturated=0.35");
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

def rotatestack_at_every_timepoint(imp_timelapse):
	# get angle and rotate at every time!
	imp = Duplicator().run(imp_timelapse)
	frames = imp_timelapse.getNFrames()
	for t in range(1,frames +1):
		imp.setT(t)
		angle = getEmbryoOrientation(Duplicator().run(imp, 1, 1, 1, 1, t, t))
		IJ.run(imp, "Rotate... ", "angle="+str(angle)+" grid=1 interpolation=Bilinear slice");

	return imp

def getembryoroi(imp_embryo):
	IJ.run(imp_embryo, "Gaussian Blur...", "sigma=4");
	IJ.setAutoThreshold(imp_embryo, thresholding_method + " dark");
	IJ.run(imp_embryo, "Convert to Mask", "");
	IJ.run(imp_embryo, "Median...", "radius=15");
	IJ.run(imp_embryo, "Fill Holes", "");
	IJ.run(imp_embryo, "Create Selection", "");
	embryo_roi = imp_embryo.getRoi()

	return embryo_roi

def getembryoroi_in_parts(imp_embryo, divisions):
	# get embryos entire roi
	IJ.run(imp_embryo, "Gaussian Blur...", "sigma=2");
	IJ.setAutoThreshold(imp_embryo, thresholding_method + " dark");
	IJ.run(imp_embryo, "Convert to Mask", "");
	IJ.run(imp_embryo, "Median...", "radius=15");
	IJ.run(imp_embryo, "Fill Holes", "");
	IJ.run(imp_embryo, "Create Selection", "");
	embryo_roi = imp_embryo.getRoi()

	imp_embryo.killRoi()

	# get its rectangular bounds
	rectangle = embryo_roi.getFloatBounds()
	width = rectangle.width
	height = rectangle.height
	x = rectangle.x
	y = rectangle.y

	# array to store multiple roi regions 
	division_rois = []

	# first rectangular region from left
	first_roi_partitioned = Roi(x, y, width/divisions, height)
	division_rois.append(first_roi_partitioned)

	# subsequent regions by changing x-coordinate and width of total width/no of divisions
	for j in range(1, divisions):
		roi_partitioned = Roi(x + ((j * width)/divisions), y, width/divisions, height)
		division_rois.append(roi_partitioned)

	return division_rois

def segmentgranules(embryo_stack):
#	embryo_stack.show()
	frames = embryo_stack.getNFrames()
	print "No of frames= ", frames
	imp_result = Duplicator().run(embryo_stack)

	resultsTableName = "Granules: degree of polarization"
	rt = getResultsTable(resultsTableName)

	for t in range(1, frames + 1):
		embryo_stack.setT(t)

		# get embryo roi regions
		imp_proc = Duplicator().run(embryo_stack, 1, 1, 1, 1, t, t)
		roi_parts = getembryoroi_in_parts(Duplicator().run(imp_proc), no_regions)

		# array to store to I and area of granules from all regions 
		all_granule_integrated_intensity_in_fractions = []
		all_granule_area_in_fractions = []
		for idx, rr in enumerate(roi_parts):
			roi_embryo_part, granules_roi, embryo_total_intensity_in_part, embryo_integrated_intensity_in_part, embryo_area_in_part, spots_count = getGranulesRoiAll(imp_proc, rr)

			# get stats from the granules 
			imp_proc1 = Duplicator().run(imp_proc)
			IJ.run(imp_proc1, "Subtract Background...", "rolling=100 slice");
		  	imp_proc1.setRoi(granules_roi)
			stats= imp_proc1.getStatistics(Measurements.ALL_STATS)
			spot_total_raw_intensity_in_part = float(stats.pixelCount * stats.mean)
			spot_integrated_intensity_in_part = float(stats.area * stats.mean)
			spots_area_in_part = stats.area 
			imp_proc1.killRoi()

			# fill the array defined before loop
			all_granule_integrated_intensity_in_fractions.append(spot_integrated_intensity_in_part)
			all_granule_area_in_fractions.append(spots_area_in_part)

			# measurements 
			volume_fraction_granules_to_embryo = spots_area_in_part / embryo_area_in_part
			volume_fraction_embryo_to_granules = embryo_area_in_part / spots_area_in_part
			intensity_fraction_granules_to_embryo = spot_total_raw_intensity_in_part / embryo_total_intensity_in_part
			intensity_fraction_embryo_to_granules = embryo_total_intensity_in_part / spot_total_raw_intensity_in_part

			# results table 
			rt.incrementCounter()
			rt.addValue('Frame', t)
			rt.addValue('Region no from left to right', idx +1)
			rt.addValue('Count', spots_count)
			rt.addValue('Area of granules in part', spots_area_in_part)
			rt.addValue('Area of embryo in part', embryo_area_in_part)
			rt.addValue('Intensity fraction granules/embryo', intensity_fraction_granules_to_embryo)
			rt.addValue('Intensity fraction embryo/granules', intensity_fraction_embryo_to_granules)
			rt.addValue('Area fraction granules/embryo', volume_fraction_granules_to_embryo)
			rt.addValue('Area fraction embryo/granules', volume_fraction_embryo_to_granules)

			# adding to overlay
			rr.setStrokeColor(Color.yellow)
			rr.setImage(imp_result)
			rr.setPosition(t)
			ov.add(rr)
			roi_embryo_part.setStrokeColor(Color.cyan)
			roi_embryo_part.setImage(imp_result)
			roi_embryo_part.setPosition(t)
			ov.add(roi_embryo_part)
			granules_roi.setStrokeColor(Color.red)
			granules_roi.setImage(imp_result)
			granules_roi.setPosition(t)
			ov.add(granules_roi)
			text = "Regions: left ---> right"; 
		  	font = Font("SansSerif", Font.PLAIN, 18); 
		  	roitext = TextRoi(10, 10, text, font); 
		  	roitext.setStrokeColor(Color.black); 
		  	roitext.setFillColor(Color.white); 
		  	roitext.setImage(imp_result)
		  	roitext.setPosition(t)
		  	ov.add(roitext)

#		print all_granule_integrated_intensity_in_fractions
#		print sum(all_granule_integrated_intensity_in_fractions)
#		print all_granule_area_in_fractions
#		print sum(all_granule_area_in_fractions)

		# computing 1 region: total wrt intensity and area of granules
		for j in range(len(all_granule_area_in_fractions)):
			rt.addValue("Granules Integrated intensity ratio:region " +str(j + 1) + " to all", float(all_granule_integrated_intensity_in_fractions[j] / sum(all_granule_integrated_intensity_in_fractions)))
			rt.addValue("Granules Area ratio:region " +str(j + 1) + " to all", float(all_granule_area_in_fractions[j] / sum(all_granule_integrated_intensity_in_fractions)))

		roi_manager.reset()

	# add overlay, display results
	imp_result.setOverlay(ov)
	IJ.run(imp_result, "Enhance Contrast", "saturated=0.2");
	imp_result.setTitle("Final result")
	imp_result.show()
	rt.show(resultsTableName)

	if saveResults:
		savePath = join(info.directory, 'results')
		if not os.path.exists(savePath):
			os.makedirs(savePath)
		rtName = imp_name + '_results.csv'
		rt.save(join(savePath, rtName))


def getGranulesRoiAll(imp_proc, roi_rect):
	imp_proc.killRoi()
	
	# get rectangular roi and embryo roi, an AND operation gives embryo part inside the rect region
	roi_embryo = getembryoroi(Duplicator().run(imp_proc))
	roi_embryo1 = ShapeRoi(roi_embryo)
	roi_rect1 = ShapeRoi(roi_rect)
	roi_new = roi_embryo1.and(roi_rect1)  # embryo roi inside rect roi 

	imp = Duplicator().run(imp_proc)
	IJ.run(imp, "Subtract Background...", "rolling=100 slice");

	# embryo stats in the fraction
	imp.setRoi(roi_new)
	embryo_stats = imp.getStatistics(Measurements.ALL_STATS)
	embryo_total_intensity =  float(embryo_stats.pixelCount * embryo_stats.mean)
	embryo_integrated_intensity = float(embryo_stats.area * embryo_stats.mean)
	embryo_area = embryo_stats.area
	imp.killRoi()

	# get granules inside the embryo in the fraction
	granules_mask = Duplicator().run(imp)
	granules_mask.setRoi(roi_new)

	# test 1 = >90%
	IJ.setAutoThreshold(granules_mask, "Triangle" + " dark");
	IJ.run(granules_mask, "Convert to Mask", "");
	IJ.run(granules_mask, "Median...", "radius=1");
	IJ.run(granules_mask, "Watershed", "");
#	granules_mask.show()

	# particle analysis to get granule rois inside the fraction (create selection not suited as it selects all granules in entire granules mask)
	options = PA.ADD_TO_MANAGER \
			+ PA.SHOW_NONE \
			+ PA.SHOW_PROGRESS 
			
	measurements = 0
						
	rt_final = ResultsTable()
	minSize = 1
	maxSize = float('inf')
	min_circularity = 0.0
	max_circularity = 1.0

	# do particle analysis and get roi
	pa = PA(options, measurements, rt_final, minSize, maxSize, min_circularity, max_circularity)
	pa.setHideOutputImage(True)
	Analyzer.setRedirectImage(imp);
	granules_mask.setRoi(roi_new)
	pa.analyze(granules_mask);

	# get no of granules
	spots_count = roi_manager.getCount()

	# use combine function to all granule rois as a single selection
	roi_manager.runCommand(granules_mask, "Select All")
	roi_manager.runCommand(granules_mask,"Combine");
	roi_manager.addRoi(granules_mask.getRoi())
	noRois = roi_manager.getCount()
	combined_roi_index = noRois - 1
	granules_roi = roi_manager.getRoi(combined_roi_index)
	roi_manager.reset()
	
	return roi_new, granules_roi, embryo_total_intensity, embryo_integrated_intensity, embryo_area, spots_count

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


		