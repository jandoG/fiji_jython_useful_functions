# ============================================= README ======================================================
#
# Author: Gayathri Nadar, SCF MPI-CBG, nadar@mpi-cbg.de
#
# Read the word document (README_TrackMate-StarDist_usage_worklow.docx) for more details on running these scripts and to get more info on the parameters to be set!
#
# Version history:
# v0.0: 05.11.2021 Segmenting and tracking droplets using Trackmate v7's StarDist detector and measure their mean intensity at each time point 
#
# Input data:
#	2D single or multichannel time-lapse movie 
#
# Description/Workflow:
#	This script using the StarDist detector in Trackmate to detect and segment the droplets. Next the spots are tracked along time and intensity is measured at each time point.
#	
#	Workflow steps:
#	- Read image
#	- Create trackmate model using StarDist detector 
#	- Detect droplets, filter them based on area, diameter and circularity
#	- Tracking using simple LAP tracker 
#	- Filter tracks: keep only those whose track duration is > 50% of the entire time-lapse duration (to remove smaller tracks which are not of interest)
#	- For the spots inside each of the tracks, measure mean intensity and other properties.
#
# Usage:
# 	Open this script with FIJI (drag and drop) and click run.
# 	A window will pop-up, enter the necessary values and click ok to run the script. Default values have already been added.
#	
#	Params to be set:
#		- Mode of processing - Current image or folder 
#		- Directory containing single time-point images (2D) if "Folder mode"
#		- Circularity threshold (ones with > threshold will be kept), default value = 0.5 
#		- Maximum diameter of droplet in um, default = 50 
#		- (optional) filter by diameter threshold (ones with > max diameter will be removed)
#	
# 	"Current image mode":
#		- Open a time-lapse movie and select "current image mode" 
#		- Set the parameters as before, and dir to save results 
#		- Click ok. 
#
#	Other params to be set:
#		- Frame interval in seconds: frame interval is read from the metadata but in case it is 0.0, the value set by the user is used. IMPORTANT TO BE SET! Default = 1.0 sec 
#		- Channel to process 
#       - Track duration threshold as percentage of entire duration: keep tracks which have duration > x% of total duration. 
#			duration threshold is computed as (%x/100.0) * (no of frames * frameInterval). Tracks having duration higher than this value will be kept. 
#			This is to discard smaller tracks. Smaller the percentage more tracks will be kept
#
# Results:
#	All results will be found in the folder "results_imagename" outside of image folder (in folder mode of processing). In case of current image, it will be saved to 
#	folder specified by the user.
# 	Saved:
#		CSV file with measurements of spots intensity over time imagename_droplets_intensity.csv 
#		RGB image with tracks overlayed imagename_tracks_overlay.tif
#		Raw image saved as movie imagename_timelapse.tif (in folder processing mode)
#
# Note:
# 	Required update sites: "IJPB-plugins" (morpholibj), "Bio-Formats", "SCF-MPI-CBG", "CSBDeep", "StarDist", "TrackMate-StarDist"
#	Launch FIJI. Go to Help > Update... and then go to Manage Update Sites and check the boxes next to the names mentioned above. Click Close. 
#	Click Apply Changes. Relaunch FIJI. 
#
# =============================================================================================================



# =============================================================================================================
# Don't modify anything below here without good reason
# =============================================================================================================


#@ String (value=" Track droplets using Trackmate-Stardist version 0.0", visibility=MESSAGE) _m1_
#@ String (label= "Processing mode", choices=("Current time-lapse image", "Folder")) processMode
#@ File (label= "Save dir for results (if 'Current image' mode)", style= "directory", required= False) save_dir
#@ String (value=" -------- Folder mode --------", visibility=MESSAGE) _m2_
#@ File (label= "Input directory containing time images (if 'Folder' mode)", style= "directory") processFolder 
#@ String (value=" -------- Additional params --------", visibility=MESSAGE) _m3_
#@ Float (label= "Maximum droplet diameter in um (default=50)", description= "Droplets having area higher than (pi*diameter^2)/4 will be discarded", value= 50.0) diameter_threshold
#@ Boolean (label= "Filter droplets based on diameter?", description= "Additional filtering option. Remove droplets which have diameter > max diameter set", required = False) filter_by_diameter
#@ Float (label= "Minimum droplet circularity (default=0.5)", description= "Set value between 0 and 1, 1 being perfect circle. Droplets with lower values will be discarded", value= 0.5) circularity_threshold
#@ Integer (label= "Minimum track duration as % of entire duration (default=50)", description= "Tracks with duration > x% of total will be kept", value= 50) duration_threshold
#@ String (value=" -------- Image params --------", visibility=MESSAGE) _m4_
#@ Float (label= "Frame interval in seconds (default=1)", value= 1.0) frame_interval
#@ Integer (label= "Channel to process (default=1)", value=1) channel_to_process

from ij import IJ, Prefs
from ij import ImagePlus
from ij.plugin import Duplicator
from ij.measure import Measurements, ResultsTable
from ij.io import FileSaver 
import loci.plugins
from loci.plugins import BF
from loci.plugins.in import ImporterOptions 
import os, math, sys 
from os.path import isfile, join 
from datetime import datetime as dt
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.util import LogRecorder;
from fiji.plugin.trackmate.tracking.sparselap import SparseLAPTrackerFactory
from fiji.plugin.trackmate.tracking import LAPUtils
from fiji.plugin.trackmate.util import TMUtils
from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate.stardist import StarDistDetectorFactory
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettings
from fiji.plugin.trackmate.gui.displaysettings import DisplaySettingsIO
from fiji.plugin.trackmate.action import CaptureOverlayAction
import java.util.ArrayList as ArrayList

reload(sys)
sys.setdefaultencoding('utf-8')

show_output = True
ext = ".TIF"
sigmaBG = 30
maxArea = ((math.pi * diameter_threshold *  diameter_threshold) / 4) 
maxArea = int(math.ceil(maxArea / 1000.0)) * 1000
print "Computed max possible area for droplets =", maxArea

def main():
	# clean start 
	Prefs.setBlackBackground=True
	Prefs.padEdges=True
	IJ.setBackgroundColor(0,0,0) # for 'clear outside'
	IJ.run("Collect Garbage");
	IJ.log("\\Clear");
	
	if processMode == "Current time-lapse image":
		print "\n \nProcess mode = Current time-lapse image"
		imp00 = IJ.getImage();
		imagename = os.path.splitext(imp00.getTitle())[0]

		print "Processing image =", imagename, "\n\n"
		
		savedir = save_dir.getPath()

		print "-------Preprocessing image...\n"
		imp0 = preprocessImage(imp00, sigmaBG)
		imp0.setTitle(imagename)

	else:
		print "\n \nProcess mode = Folder"
		imgFolder= processFolder.getPath() 
		files = sorted([ f for f in os.listdir(imgFolder) if isfile(join(imgFolder, f)) and f.endswith(ext)])
		imp00 = openImageWithBF(join(imgFolder, files[0]), True, True)
		name = os.path.splitext(imp00.getTitle())[0]
		imagename = name.split("_t<")[0]

		print "Processing image =", imagename, "\n\n"

		savedir = join(os.path.dirname(imgFolder), "results_" +  imagename)
		if not os.path.exists(savedir):
			os.makedirs(savedir)

		FileSaver(imp00).saveAsTiff(os.path.join(savedir, imagename + "_timelapse.tif"))

		print "-------Preprocessing image...\n"
		imp0 = preprocessImage(imp00, sigmaBG)
		imp0.setTitle(imagename)

	print "-------Detecting droplets and tracking...\n"
	print "This might take time...."
	trackmate = create_trackmate(imp0) 
	
	ok = trackmate.checkInput()
	if not ok:
	    print( str( trackmate.getErrorMessage() ) )
	    return
	    
	ok = trackmate.process()
	if not ok:
	    print( str( trackmate.getErrorMessage() ) )
	    return

	# Get track model
	model = trackmate.getModel()
	tm = model.getTrackModel()
	fm = model.getFeatureModel()
	trackIDs = tm.trackIDs(True)
	table = ResultsTable() 

	print "-------Computing spot paramaters and creating results...\n"

	# populate results table for each spot in each track 
	for id in trackIDs:
		# Fetch the track feature from the feature model.
		t_duration = fm.getTrackFeature(id, 'TRACK_DURATION')
		v = fm.getTrackFeature(id, 'TRACK_DISPLACEMENT')
		print "Tracks details:"
		print "Track ID: {}, has duration {} secs and displacement {} um".format(id, t_duration, v)
	    
		# get spots and sort them by frame.
		spots = tm.trackSpots(id) 
		ls = ArrayList(spots);
		ls_sorted = sorted(ls, key=lambda x: x.getFeature('FRAME'), reverse=False)

		for spot in ls_sorted:
			sid = spot.ID()
#			print(spot.getFeatures())
			t = spot.getFeature('FRAME')
			x = spot.getFeature('POSITION_X')
			y = spot.getFeature('POSITION_Y')
			total = spot.getFeature('TOTAL_INTENSITY_CH1')
			mean = spot.getFeature('MEAN_INTENSITY_CH1')
			max_intensity = spot.getFeature('MAX_INTENSITY')
			radius = spot.getFeature('RADIUS')
			diam = radius * 2.0
			area = spot.getFeature('AREA')
		
			table.incrementCounter();
			table.addValue("Track ID", id);
			table.addValue("Spot ID", sid);
			table.addValue("Frame no", t + 1);
			table.addValue("Channel no", channel_to_process)
			table.addValue("Mean intensity", mean);
			table.addValue("Diameter", diam);
			table.addValue("Area", area);

	tablename = "Droplets: Intensities over time"
	table.show(tablename);

	if show_output:
		model = trackmate.getModel()
		selectionModel = SelectionModel(model)
		ds = DisplaySettings()
		ds = DisplaySettingsIO.readUserDefault()
		ds.spotDisplayedAsRoi = True
		displayer =  HyperStackDisplayer( model, selectionModel, imp0, ds)
		displayer.render()
		displayer.refresh()

	model.getLogger().log(str(model))
	logger = model.getLogger()
	settings = trackmate.getSettings()

	# capture overlay 
	image = trackmate.getSettings().imp 
	capture = CaptureOverlayAction.capture(image, -1, imp0.getNFrames(), logger)
#	capture.show()

	print "\n-------Writing to XML file...\n"
	saveFile = TMUtils.proposeTrackMateSaveFile(settings, logger)
	writer = TmXmlWriter(saveFile,logger)
	writer.appendLog(logger.toString())
	writer.appendModel(trackmate.getModel())
	writer.appendSettings(trackmate.getSettings())
	writer.writeToFile();
	print "XML file saved to: " + saveFile.toString();
	print "NOTE: if XML is not saved in the same folder as image, move and save it next to the image being processed!!"
	print "Use 'Load a Trackmate file' in FIJI to load XML file and explore options like Trackscheme and plots.\n"

	print "-------Saving results to folder =", savedir + "\n"
	table.saveAs(os.path.join(savedir, imagename + "_droplets_intensity.csv"));
	FileSaver(capture).saveAsTiff(os.path.join(savedir, imagename + "_tracks_overlay.tif"))
	


# ===================== HELPER ===============================

def create_trackmate(imp):
	"""
	Creates a TrackMate instance configured to operated on the specified
	ImagePlus imp. 
	
	params: imp
	returns: trackmate object with all the settings 
	"""
	
	# if frame interval is 0 then set it manually
	cal = imp.getCalibration()
	if cal.frameInterval == 0.0:
		cal.frameInterval = frame_interval

	model = Model()
	model.setLogger(Logger.IJ_LOGGER)

	settings = Settings(imp)
	setup = settings.toStringImageInfo()
	       
	# Configure StarDist default detector.
	settings.detectorFactory = StarDistDetectorFactory()
	settings.detectorSettings = {
	    'TARGET_CHANNEL' : channel_to_process
	}

	# Configure tracker
	settings.trackerFactory = SparseLAPTrackerFactory()
	settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap()
	settings.trackerSettings[ 'LINKING_MAX_DISTANCE' ] 		= 10.0
	settings.trackerSettings[ 'GAP_CLOSING_MAX_DISTANCE' ]	= 10.0
	settings.trackerSettings[ 'MAX_FRAME_GAP' ]				= 3
	settings.initialSpotFilterValue = -1.

	# analyzers 
	settings.addAllAnalyzers()

	# Add some filters for tracks/spots 
	# filter on duration = keep tracks > threshold % of total duration 
	maxduration = (duration_threshold/100.0) * (imp.getNFrames() * cal.frameInterval)                
	filter1_track = FeatureFilter('TRACK_DURATION', maxduration, True)
	settings.addTrackFilter(filter1_track)

	# filter on spots 
	filter1_spot = FeatureFilter('AREA', maxArea, False)						# discard spots above max area 
	filter2_spot = FeatureFilter('RADIUS', (diameter_threshold/2.0), False)     # discard spots above max diameter 
	filter3_spot = FeatureFilter('CIRCULARITY', circularity_threshold, True)	# keep droplets above min circ 
	settings.addSpotFilter(filter1_spot)
	if filter_by_diameter:
		settings.addSpotFilter(filter2_spot)
	settings.addSpotFilter(filter3_spot)

	print "Spot filters added = ", settings.getSpotFilters()
	print "Track filters added = ", settings.getTrackFilters(), "\n"
	
	# trackmate 	    
	trackmate = TrackMate(model, settings)
	trackmate.computeSpotFeatures(True)
	trackmate.computeTrackFeatures(True)

	return trackmate 


def preprocessImage(timelapseImp, sigmaBG = 30):
	"""
	Subtracts the backgound (one value per slice)
	Value per slice is needed because some slices are different: completely black

	params: timelapse imp, sigma for blurring (default=30)
	returns: background subtracted timelapse imp
	"""
	timelapse_imp = timelapseImp.duplicate()
	background = getBackgroundPerFrame(timelapse_imp, sigmaBG)
	frames = timelapse_imp.getNFrames()
	
	for f in range(1, frames + 1):
		timelapse_imp.setT(f)
#		print "Value subtracted at t = {}, is {}".format(f, background[f-1])
		IJ.run(timelapse_imp, "Subtract...", "value="+str(background[f-1])+" slice"); 

	return timelapse_imp
		

def getBackgroundPerFrame(imp0, sigma):
	"""
	Computes the background values per frame (as minimum of a blurred image) for a timelapse image. Input must be single channel
	note: some frames may have value=0 everywhere -> don't average over full stack for background but do per frame.
	
	params: imp0 timelapse image, sigma for blurring background before taking minvalue. default:30
	returns: returns an array with background intensity per frame
	"""
	imp = imp0.duplicate()
	IJ.run(imp, "Gaussian Blur...", "sigma="+str(sigma)+" stack");
	
	background = [] 
	frames = imp.getNFrames() 
	
	for frame in range(1, frames + 1):
		imp.setT(frame)
		stat = imp.getStatistics()
		background.append(stat.min)

	return background 


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



main() 
print "\n\nDone"