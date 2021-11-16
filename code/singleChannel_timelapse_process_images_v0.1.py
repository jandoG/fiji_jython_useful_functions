# ============================================= README ======================================================
#
# Author Gayathri Nadar, SCF MPI-CBG, nadar@mpi-cbg.de
#
# Version history:
# 	16.02.2021 Version 0.0: single channel 2D/ 3D z-stack segment spots inside nuclei and measure different parameters
#   10.03.2021 v0.1: extend script to process time-lapse images. Measurements at each time point 
#
# Description:
# 	EXPECTED DATA: SINGLE CHANNEL 2D or 3D Z-STACK DATA with time-lapse 
#   This script asks user to outline a roi and tracks the proteins (bright spots) inside this roi and quantifies them
#
# Steps:
#	For each file in the list:
#		- check if it's single channel - if not skip to next  
#		- check if 2D or z-stack - if z-stack, we do a sum projection and proceed for analysis in 2D
#		- ask user to draw a nuclei roi in 2D image 
#		- get spots inside them - track spots and get their measurements 
#		- sum project over time: get spots inside labels and roi 
#		- get measurements for the nuclei roi = area, total intensity (sum of pixel values), mean 
#		- using spots label get measurements for the spots = area, diameter, mean, total intensity (sum of pixel values) 
#		
# Output:
#	For each file:
#		- result image with detected spots and nuclei roi drawn by user as overlay + tracks 
#		- label image of the detected spots from sum projection of time-lapse 
#	For all files:
#		- table with overall measurements - nuclei, all spots, nucleoplasm (nuclei - spots), spots vs nucleoplasm 
#		- table with spots measurements - individual spots diameter, area, etc 
# 
# Note: 
#	The following plugins are needed to run this script:
# 	1. SCF MPI CBG 
#	2. MorpholibJ (plugin name: IJPB plugins)
#
#	In case you do not have these plugins then please install it by going to Help > Update...
#	- Then in window that comes up, click on "Manage Update Sites". Click the box next to mentioned names. 
#	- Click "close" and then "Apply changes".
#
# Usage:
#	- Open the script with FIJI
#	- Click run
#	- Enter the required parameters in pop up window
#		1. file/files to be processed: 2D/3D single channel time-lapse files 
#		2. option to display result images 
#		3. spot radius in um 
#       4. quality threshold for spot detection 
#	- Click ok to run the script
#
# Results:
# 	Results are stored in a folder selected by user 
# 	Results saved:
#	1. Two CSV files with overall and spots measurements from all the files selected to process
# 	2. For each file, images with overlay of the spots and nuclei segmentation: yellow = nuclei, magenta = spots, yellow = tracks 
#	3. For each file, label image with spots IDs from sum projection of time-lapse 
#
# =============================================================================================================


# =============================================================================================================
# Don't modify anything below here without good reason
# =============================================================================================================


#@ File[] (label= "Select 1 or more images (2D or stack + time) with 1 CHANNEL: choose images with similar spot size", style="files") listfiles
#@ File (label= "Select folder to save results", style="directory", required= True) results_path
#@ String (label= "File extension", choices= {".tif", ".tiff"}) ext 
#@ Float (label= "Approximate RADIUS of the spot in um (spots smaller will be discarded)", value= 0.4) radius 
#@ Float (label= "Quality threshold: higher the value, less spots detected", value= 50) quality 
#@ boolean (label="Show image results (checked = yes, unchecked = no)") displayimage

from ij import IJ
from ij import ImagePlus
from ij import ImageStack
from ij.plugin import Duplicator, ZProjector
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable, Measurements, Calibration
from ij.gui import Roi, Overlay, PolygonRoi, WaitForUserDialog
from ij.io import FileSaver
from java.awt import Color
import loci.plugins
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from inra.ijpb.label import LabelImages
from inra.ijpb.binary import BinaryImages
from inra.ijpb.binary.conncomp import FloodFillComponentsLabeling
from inra.ijpb.measure import IntensityMeasures, IntrinsicVolumes2D
from inra.ijpb.measure.region2d import MaxFeretDiameter
import os, math, sys
from os.path import isfile, join

from fiji.plugin.trackmate import Spot as Spot
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking.sparselap import SparseLAPTrackerFactory
from fiji.plugin.trackmate.tracking.kalman import KalmanTrackerFactory
from fiji.plugin.trackmate.tracking import LAPUtils
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer as HyperStackDisplayer
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from fiji.plugin.trackmate.providers import SpotAnalyzerProvider
from fiji.plugin.trackmate.action import ExportStatsToIJAction
from fiji.plugin.trackmate.gui import TrackMateGUIController
from fiji.plugin.trackmate.action import CaptureOverlayAction
import java.util.ArrayList as ArrayList 


# ======= 

table_spots = ResultsTable() 
table_overall = ResultsTable() 

def main():
	# get list of files ending with ext defined 
	sorted_list = sorted([f.getPath() for f in listfiles if f.isFile() and f.getPath().endswith(ext) and not os.path.dirname(f.getPath()).endswith("_results")])

	# process each file in the list
	for imgfile in sorted_list:
		imp = openImageWithBF(imgfile, virtual= False)
		impName = os.path.splitext(imp.getTitle())[0]
		print "Processing image =", impName 

		# overlay for each file
		overlay = Overlay()
		cal = imp.getCalibration() 
		dims = imp.getDimensions()

		if dims[2] > 1:
			print "Image is multichannelled, skipping current image."
			continue

		# if z-stack, sum-project and process 
		if dims[3] > 1:
			imp_max_proj = ZProjector().run(imp, "sum all");
			imp_max_proj.setCalibration(cal)
			imp_max_proj.setTitle(impName + "_sumprojection")
			
			## Individual spots measurements: track spots and get measurements 
			nuclei_imp, nuclei_roi = trackSpotsMeasure(imp_max_proj, impName, cal, overlay)
			
			## For overall measurements: sum project the time-lapse and measure in 2D  
			print "Sum projecting time-lapse for overall measurements..."
			impsum = nuclei_imp.duplicate()
			imp_sumT = ZProjector().run(impsum, "sum all");

			doNucleiSpotsOverallMeasurement(imp_sumT, nuclei_roi, impName, cal)

		# if 2D, process directly 
		else:
			print "Processing 2D!"
			## Individual spots measurements: track spots and get measurements 
			nuclei_imp, nuclei_roi = trackSpotsMeasure(imp, impName, cal, overlay)

			## For overall measurements: sum project the time-lapse and measure in 2D  
			print "Sum projecting time-lapse for overall measurements..."
			impsum = nuclei_imp.duplicate()
			imp_sumT = ZProjector().run(impsum, "sum all");

			doNucleiSpotsOverallMeasurement(imp_sumT, nuclei_roi, impName, cal)
				 

	table_overall.show("Table: All files overall measurements (from sum-projection over time)")
	table_spots.show("Table: All files spots measurements")

	# save results tables 
	tableName_overall = 'SingleCh_time-lapse_all_files_overall_measurements_sumproj_over_time.csv'
	table_overall.saveAs(join(results_path.getPath(), tableName_overall))
	tableName_spots = 'SingleCh_time-lapse_all_files_spots_measurements.csv'
	table_spots.saveAs(join(results_path.getPath(), tableName_spots))


#===== helper 

def trackSpotsMeasure(imp_max_proj, impName, cal, overlay):
	"""
	Main function to get track the spots over time and get its measurements. 
	The spots are detected in the image and tracked using trackmate. The measurements from TrackMate are populated into the spots table. 
	The image with tracks and detected spots is saved. 
	
	params: images (2d + t), image title, calibration, overlay object 
	returns: nuclei image (image cleared outside the nuclei roi), nuclei roi 
	"""
	# get nuclei roi and imp with outside of nuclei roi = cleared 
	nuclei_roi, nuclei_imp = getNucleiRoi(imp_max_proj)
	nuc = nuclei_imp.duplicate()

	# trackmate 
	nuclei_imp.killRoi()
	trackmate = create_trackmate(nuclei_imp, quality_threshold=quality, initial_spot_filtering_threshold = 1000)
	
	ok = process(trackmate)
	if not ok:
		sys.exit(str(trackmate.getErrorMessage()))

	# display overlay of tracks
	model = trackmate.getModel()
	selectionModel = SelectionModel(model)
	displayer =  HyperStackDisplayer(model, selectionModel, nuclei_imp)
	displayer.render()
	displayer.refresh() 
	capture = CaptureOverlayAction.capture(nuclei_imp, 1, nuclei_imp.getNFrames(), Logger.IJ_LOGGER)
	overlay.add(nuclei_roi);
	capture.setOverlay(overlay)
	captureName = "tracks_overlay_for_" + impName
	capture.setTitle(captureName) 

	# Get track model
	tm = model.getTrackModel()
	trackIDs = tm.trackIDs(True) 

	# populate results table for each spot in each track 
	for id in trackIDs:
		spots = tm.trackSpots(id) 
		# Let's sort them by frame.
		ls = ArrayList(spots);
		
		# ls.sort(Spot.frameComparator) # gives error https://forum.image.sc/t/trackmate-error-after-updating-fiji/11137/13 
		ls_sorted = sorted(ls, key=lambda x: x.getFeature('FRAME'), reverse=False) 

		for spot in ls_sorted:
			sid = spot.ID()
			t = spot.getFeature('FRAME')
			x = spot.getFeature('POSITION_X')
			y = spot.getFeature('POSITION_Y')
			total = spot.getFeature('TOTAL_INTENSITY')
			mean = spot.getFeature('MEAN_INTENSITY')
			max_intensity = spot.getFeature('MAX_INTENSITY')
			diam = spot.getFeature('ESTIMATED_DIAMETER') 
			
			# table 2 for individual spots measurements 
			table_spots.incrementCounter(); 
			table_spots.addValue("Image", impName)
			table_spots.addValue("Track ID", id)
			table_spots.addValue("Frame no", t + 1)
			table_spots.addValue("Diameter (" + cal.getUnits() + ")", diam)
			table_spots.addValue("Area (" + cal.getUnits() + "^2)", math.pi * (diam/2 ** 2))
			table_spots.addValue("Mean", mean)

	if displayimage:
		capture.show()
		
	FileSaver(capture).saveAsTiff(join(results_path.getPath(), impName + "-tracks.tif"))
	
	return nuc, nuclei_roi


def doNucleiSpotsOverallMeasurement(imp_sumT, nuclei_roi, impName, cal):
	"""
	Main function to get overall measurements. 
	Spots are segmented in the image imp_sumT inside the nuclei roi and measurements are obtained for these spots. 
	Results are populated into the overall measurement table. 
	The spots label image is saved. 

	Note: the number of spots could vary from the tracking output; this measure is just meant to get an approximate idea 
	of ratio of spots vs nucleoplasm
	
	params: imp sum-projected over time (2D + t), nuclei roi, image title, calibration 
	returns: None
	"""
	# get nuclei stats
	nuclei_area, nuclei_mean, nuclei_total = getNucleiMeasurements(imp_sumT, nuclei_roi)  

	# get spots roi and label 
	spots_label, spots_roi = getSpotsRoiAndLabel(imp_sumT, nuclei_roi)

	# spots stats 
	labels, meanI, minI, maxI, totalI, diameters, areas = getSpotsMeasurements(imp_sumT, spots_label, cal)
	
	# table 1 for overall measurements 
	# nuclei = roi drawn by user (with spots inside)
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


def create_trackmate(imp0, quality_threshold=30.0, initial_spot_filtering_threshold = 2000):
	"""
	Initiate trackmate with params:

	imp: image to track 
	quality: higher the value, less spots detected 
	initial spot filtering threshold: for intensity measurements 

	returns: trackmate object 
	"""
	cal = imp0.getCalibration()
	imp = imp0.duplicate()

	# Model
	model = Model()
	model.setLogger(Logger.IJ_LOGGER)
	model.setPhysicalUnits(cal.getUnit(), cal.getTimeUnit())
	
	# Settings
	settings = Settings()
	settings.setFrom(imp)

	# Configure detector
	settings.detectorFactory = LogDetectorFactory()
	settings.detectorSettings = { 
	    'DO_SUBPIXEL_LOCALIZATION' : True,
	    'RADIUS' : float(radius),
	    'TARGET_CHANNEL' : 1,
	    'THRESHOLD' : float(quality_threshold),
	    'DO_MEDIAN_FILTERING' : True,
	}  

	# Configure spot filters
	filter1 = FeatureFilter('MEAN_INTENSITY', initial_spot_filtering_threshold, True)
	settings.addSpotFilter(filter1)

	# Configure tracker 
	settings.trackerFactory = SparseLAPTrackerFactory()
	settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap()
	settings.trackerSettings[ 'LINKING_MAX_DISTANCE' ] 		= 3.0
	settings.trackerSettings[ 'GAP_CLOSING_MAX_DISTANCE' ]	= 3.0
	settings.trackerSettings[ 'MAX_FRAME_GAP' ]				= 1

	# Add spots analyzer 
	spotAnalyzerProvider = SpotAnalyzerProvider()
	for key in spotAnalyzerProvider.getKeys():
		settings.addSpotAnalyzerFactory(spotAnalyzerProvider.getFactory(key))

	# Create the TrackMate instance
	trackmate = TrackMate(model, settings)
	trackmate.getModel().getLogger().log(settings.toStringFeatureAnalyzersInfo())
	trackmate.computeSpotFeatures(True)

	settings.initialSpotFilterValue = 0
	
	return trackmate


def process(trackmate):
	"""
	Execute the full process
	"""
	# Check settings.
	ok = trackmate.checkInput()
	print ok

  	# Process 
	ok = trackmate.process()

	# Initial filtering
	print( 'Spot initial filtering' )
	ok = ok and trackmate.execInitialSpotFiltering()

	# Filter spots 
	print( 'Filtering spots' )
	ok = ok and trackmate.execSpotFiltering( True )
	
	# Track spots
	print( 'Tracking' )
	ok = ok and trackmate.execTracking()

	# Compute spot features.
	print( 'Computing spot features' )
	ok = ok and trackmate.computeSpotFeatures( True ) 

	return ok
	

def getNucleiRoi(imp0):
	"""
	Get nuclei roi and image with no spots about the drawn roi  
	params: imp (2D)
	returns: nuclei roi (input by user), image with everything cleared outside of drawn roi 
	"""
	imp = imp0.duplicate()

	# get nuclei roi from user
	imp.show();
	IJ.setTool("freehand");
	WaitForUserDialog("Action required", "Please draw a ROI around the nuclei of interest and press ok.").show() 
	nuclei_roi = imp.getRoi();
	imp.hide();
	imp.killRoi(); 

	# clear outside
	nuclei_mask = Duplicator().run(imp)
	nuclei_mask.setRoi(nuclei_roi)
	IJ.run(nuclei_mask, "Clear Outside", "stack");
	IJ.run(nuclei_mask, "Gaussian Blur...", "sigma=1 stack"); 

	return nuclei_roi, nuclei_mask 

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

def applyGlasbeyLUT(imp):
	"""Applies Glasbey on Dark LUT on label maps"""
	IJ.run(imp, "glasbey_on_dark", "");


def openImageWithBF(path, virtual= True):
	"""
	set options to open image using bio-formats- use virtual for quick loading
	"""
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




























			