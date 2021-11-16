# ============================================= README ======================================================
#
# Author: Gayathri Nadar, SCF MPI-CBG, nadar@mpi-cbg.de
# 
# THIS SCRIPT HAS BEEN ADAPTED FROM THE SCRIPT DROPLET_ANALYSIS_0.4 BY BENOIT LOMBARDOT (former SCF MPI CBG). FOR ANALYSIS FROM 2021 USE THIS UPDATED SCRIPT! 
# WHAT IS NEW:
#	- Open images with Bio-formats - can read many proprietary file formats 
#	- Option to process a series from multiple series in an image
#	- Option to set calibration, frame, channel number to process  
#	- Option to compute background from raw image if not available. 
#	- Save results automatically in folder processing mode. 
#	- Documentation 
#
#
# Version history from 2021:
# v0.0: droplet_analysis_2021_v0.0.py 29.09.2021 Generic droplet analysis script for 2D droplet images or slice by slice images of 3D stack - can process an open image of a folder of images 
#				   Count the number of droplets and measure the ratio of signal sum in the droplet to signal sum in the media (area outside of droplets)
# v0.1: 26.10.2021 Additional options added 
#					1. Project images to 2D - if yes, images will be projected using either max or sum projection as set by user. Note: this image will be saved as well.
#					2. Display output - if kept unchecked, no image pop-ups occur
#
# 
# Version history till 2020: 
# v0  : 2016-02-11
# v0.1: 2016-03-06 (add the ratio of average intensity to the results table with name ratioAvgIDropletToAvgIMedia, request from Jie)
# v0.2: 2016-12-16 add possibility to process a whole folder (by Robert Haase, SCF, MPI CBG, rhaase@extern.mpi-cbg.de)
# v0.3
# v0.4: Dec 2019 Gayathri Nadar SCF MPI-CBG nadar@mpi-cbg.de
# 		0.3 was selecting the background instead of the droplets. Added line to invert the selection.
#		Added additional column to table dropletFraction= ratioIDropletsToIMedia/(1+ratioIDropletsToIMedia)
#		Feb 2020 changed dropletFraction= ratioIDropletsToIMedia/1+ratioIDropletsToIMedia to just dropletFraction. The measure still is the same.
#
# 
# Description:
#	This script segments droplets in 2D images of 2D slices of z-stack and counts the number of droplets and measures the ratio of signal sum in the droplet to signal sum in the media (area outside of droplets)
#   Measurements are corrected for background value (mean of background image). The background mean is subtracted from all measurements. 
#	If background image is not available, it is computed from the raw image using 'Subtract background' function in ImageJ 
#	Can process:
#		2D image or 
#		Series of independant images in a stack. 
#		Only single channel or time point! 
#	Background if available:
#		Should be 2D single channel, frame
# 
#
# Steps for 2D image:
#		- Get background image if available or compute background image 
#		- Compute threshold for the image using custom function
#		- Threshold image and get mask of droplets
#		- Do measurements by computing stats for droplets and stats for background or media 
#		- Create results: image with overlay and results table. 
#
#
# Measures in the table: 
#		- Image name 
#		- Time stamp
#		- Slice no processed 
#		- Droplets count (nDroplets): total number of droplets found 
#	    - MinI background: min intensity of background image 
#	    - SumINorm droplets (sumIDroplets): sum intensity of droplets normalized wrt area (sum droplets/ area of image)
#	    - SumINorm media (sumIMedia): sum intensity of the media normalized wrt area (sum media/ area of image)
#		- AvgI droplets (IAvgDroplets): avg intensity droplets (minus bg)
#		- AvgI media (IAvgMedia): avg intensity media (minus bg)
#		- Ratio: SumINorm droplets/ SumINorm media (ratioIDropletsToIMedia): sum intensity droplets / sum intensity media (use unnormalized values)		
#		- Ratio: AvgI droplets/ AvgI media (ratioAvgIDropletsToAvgIMedia): avg intensity droplets / avg intensity media 
#		- Droplet fraction: ratioIDropletsToIMedia / (1 + ratioIDropletsToIMedia) 
#
#
# Usage:
#	- Open this script with FIJI and click run. 
#	- A window will pop-up, enter the necessary values and click ok to run the script.
#
#
# Params to be set: 
#	- Process mode: Current image or Folder 
#		- current image mode: open an image and run the script. 
#		- Save directory: choose directory to save results from processing current image. 
#	- Options for folder processing 
#		- Input directory containing images to process 
#		- Series data - set "Yes" or "No"
#		- Series ID to open - set series no to open and process 
#		- File extenstion: for example .vsi, .tiff  
#	- Background 
#		- Background image available: set "Yes" or "No" (if no, it will be computed from image)
#		- Background image file: If set "Yes", choose the background image (2D)
#	- Advanced 
#		- Process as 2D (image stacks will be projected): yes, no check box 
#		- Type of projection (if 'yes' to project as 2D): sum or max 
#		- Display output: if yes, output images will be displayed with overlay of detected droplets. If unchecked, this will not displayed, but saved to output folder.
#		- Channel, frame number to process (in case of hyperstack)
#		- Calibration details: pixel width, unit, etc. 
#		- Threshold multiplier parameter: Threshold = I0+A*StdDev, Set param A i.e. threshold multiplier. I0 is the intensity where histogram value is maximum. Higher value, smoother segmentation
#
#
# Output:
#	- Images: 
#		Yellow overlay = droplets segmented
#		Max_imagename of Sum_imagename = when process as 2D is selected
#	- Table
#	
#	For current image mode, results saved in folder set by user. 
#	For folder mode, results saved in folder results_folderProcessing which lies outside of data folder set by user to process.
#
#
# =============================================================================================================


# =============================================================================================================
# Don't modify anything below here without good reason
# =============================================================================================================

#@ String (value=" Droplet Analysis 2021 v0.0", visibility=MESSAGE) _m1_
#@ String (label= "Processing mode", choices=("Current image", "Folder")) processMode
#@ File   (label= "Save directory (if 'Current image' mode)", style= "directory") saveDir
#@ String (value=" -------- Options to choose for 'Folder' mode --------", visibility=MESSAGE) _m2_
#@ File   (label= "Input directory containing images (if 'Folder' mode)", style= "directory") processDir
#@ String (label= "Series data", choices=("Yes", "No"), required=True) hasseries
#@ Integer (label= "Series ID to open (if 'Series') (default = 1)", value=1) openseriesId 
#@ String (label="File extenstion (default = .tif)", default=".tif") ext
#@ String (value=" -------- Background (will be computed if not available) --------", visibility=MESSAGE) _m3_
#@ String (label= "Background image available", choices=("Yes", "No"), required=True) bgAvailable
#@ File   (label= "Background image file (if 'yes' above)", style= "file", required=False) bgFile
#@ String (value=" -------- Advanced options --------", visibility=MESSAGE) _m4_
#@ Boolean (label= "Process as 2D (image stacks will be projected)", required=False) projectImp 
#@ String  (label= "Type of projection (if 'yes' above)", choices= ("Max", "Sum")) projectionType 
#@ Boolean (label= "Display output", description= "Lots of popups if checked. If unchecked, they will be only saved to output folder", required=False) displayRes
#@ Integer (label= "Channel number (default = 1)", value=1) channel 
#@ Integer (label= "Frame number (default = 1)", value=1) frame  
#@ Float   (label="Pixel width (default = 0.1625)", value=0.1625) pixelWidth
#@ String  (label="Pixel unit (default = um)", value="um") pixelUnit
#@ Float   (label="Threshold multiplier parameter (default = 20)", value=20, description= "Threshold = I0+A*StdDev, Set param A i.e. threshold multiplier. I0 is the intensity where histogram value is maximum. Higher value, smoother segmentation.") threshMultiplyer

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
import os, math, re
import fnmatch as fnm
from os.path import isfile, join
from os import listdir 
import glob
import time


from ij.process import ImageStatistics
from ij.process import StackStatistics
from ij.gui import Overlay, Roi
from ij.process import FloodFiller
from ij.measure import ResultsTable
from java.awt import Color
import time, datetime, sys

# clean start 
Prefs.setBlackBackground=True
Prefs.padEdges=True
IJ.setBackgroundColor(0,0,0) 

# will be used if background image is not available
rolling_ball_radius_for_bg_creation = 100

# =========== main =============== 

def main():	
	IJ.log("\\Clear")
	
	if processMode == "Folder":
		IJ.log("Process mode: Folder")
		imgfolder = processDir.getPath()
		
		IJ.log("Processing folder: " + imgfolder + "\n")
		processFolder(imgfolder)
	
	else:
		IJ.log("Process mode: Current image")
		imp = IJ.getImage(); 
		
		if bgAvailable == "Yes":
			impbg = openImageWithBF(bgFile.getPath(), virtual= True, groupfiles = False, seriesdata = False)
		else:
			impbg = computeBackgroundImage(imp, rolling_ball_radius=rolling_ball_radius_for_bg_creation)

		IJ.log("Processing image: " + imp.getTitle())
		processImage(imp, impbg, saveDir.getPath()) 



# =========== helper functions =============== 


def processFolder(fileDir):
	"""
	Main function to process directory of image files. 
	Each file is processed using function 'processImage'

	params: 
	fileDir: str path to directory containing files. 

	Output: 
	A folder called results_folderProcessing is created in the same location as the selected data dir. 
	This folder contains N CSV files and N image files where N is no of files processed. 
	"""
	# get files with set ext from dir 
	fileList = getFilesFromDir(fileDir, extension=ext)
	nFile = len(fileList)

	# create save dir 
	savedir = join(os.path.dirname(fileDir), "results_folderProcessing")
	if not os.path.exists(savedir):
		os.makedirs(savedir)

	# find out if series data and get series ID to open 
	if hasseries == "Yes":
		seriesdata = True
		seriesno = openseriesId
	else:
		seriesdata = False
		seriesno = None

	# process each file 
	for idx, imageFile in enumerate(fileList):
		IJ.log("Processing " + imageFile + " (" + str(idx + 1) + "/" + str(nFile) + ")")
		imp0 = openImageWithBF(imageFile, virtual= True, groupfiles = False, seriesdata = seriesdata, openseries = seriesno) 

		if bgAvailable == "Yes":
			impbg = openImageWithBF(bgFile.getPath(), virtual= True, groupfiles = False, seriesdata = False)
		else:
			impbg = computeBackgroundImage(imp0, rolling_ball_radius=rolling_ball_radius_for_bg_creation)
		
		processImage(imp0, impbg, savedir)

	return


def processImage(imp0, impBg, savedir):
	"""
	Main function to process image. Works for both 2D and 3D stack. In case of 3D stack, the processing happens slice by slice. 
	BUT ONLY SINGLE CHANNEL OR FRAME! 
	The processing occurs in the channel and frame set by user in case of processing hyperstack. 
	Raw image is duplicated to have channel, frame set by user and then processed. 
	
	For each image/slice, compute threshold, using threshold create mask of droplets. Get roi of droplets.
	Do measurements for droplets and bg (media) and populate results table. 

	If "Process as 2D" is checked, the 3D stacks will be projected into 2D either by sum or max projection (as set). 
	For the rest, the processing happens as before in specified frame and channel. 

	params: 
	imp0: image to processs (could be 2D or 3D stack)
	imgBg: background image 
	savedir: str path to directory where results will be saved. 

	returns: 
	None 

	Output:
	One results table and one result image with overlay of droplets is saved to the savedir. 
	One table for each image. Table has measurements per slice if z-stack 
	
	"""

	imagename = os.path.splitext(imp0.getTitle())[0]
	
	# get time stamp
	timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d-%H_%M_%S')

	# set measurements 
	IJ.run("Set Measurements...", "area mean redirect=None decimal=5"); 

	# get bg stats 
	Ibg = StackStatistics(impBg).min

	# get image metadata 
	imp00 = imp0; 
	cal = imp00.getCalibration()
	cal.pixelHeight = pixelWidth
	cal.pixelWidth = pixelWidth
	cal.setUnit(pixelUnit)
	dims = imp00.getDimensions()

	# process duplicate 
	if projectImp:
		if projectionType == "Sum":
			impz = ZProjector.run(imp00,"sum");
			IJ.saveAsTiff(impz, join(savedir, "SUM__" + imagename + ".tif"))
			imp0 = Duplicator().run(impz, channel, channel, 1, 1, frame, frame)
		else:
			impz = ZProjector.run(imp00,"max");
			IJ.saveAsTiff(impz, join(savedir, "MAX__" + imagename + ".tif"))
			imp0 = Duplicator().run(impz, channel, channel, 1, 1, frame, frame)
	else:
		imp0 = Duplicator().run(imp00, channel, channel, 1, dims[3], frame, frame)
	
	imp0.setCalibration(cal)

	# check if hyperstack 
	isHyperstack = imp0.isHyperStack()
	if isHyperstack :
		nSlices = imp0.getNSlices() # number of z slices only
	else :
		nSlices = imp0.getStackSize()

	# step 1: compute threshold for image 
	threshold, IMaxFreq, ISigma = getStackThreshold(imp0)

	# create empty stack for results 
	resultStack = ImageStack(imp00.getWidth(), imp00.getHeight() );
	ovResults = Overlay();

	# ensure one table for all images in case of folder processing 
	rt = ResultsTable();
	rt.reset()

	for z in range(1, nSlices+1) :
		imp = Duplicator().run(imp0, channel, channel, z, z, frame, frame)
		imp.setCalibration(cal)

		# step 2: analyze image
		measures, roiDroplets, roiMedia= analyzeImage(imp, threshold, Ibg)

		# fill a result table with img name, slice, ratio, sumIdrop, sumIMedia,
		rt.incrementCounter();
		rt.addValue("Image name", 								imagename );
		rt.addValue("Time", 									timestamp);
		rt.addValue("Slice", 									z );
		if isHyperstack:
			rt.addValue("Channel", 								channel)
			rt.addValue("Frame", 								frame)
		if projectImp:
			rt.addValue("Projection", 							projectionType)
		rt.addValue("Droplets count", 							measures["nDroplets"] );
		rt.addValue("MinI background", 							Ibg)
		rt.addValue("SumINorm droplets", 							measures["sumIDroplets"] );
		rt.addValue("SumINorm media", 								measures["sumIMedia"] );
		rt.addValue("AvgI droplets", 							measures["IAvgDroplets"] );
		rt.addValue("AvgI media", 								measures["IAvgMedia"] );
		rt.addValue("Ratio: SumI droplets/ SumI media", 		measures["ratioIDropletsToIMedia"] );
		rt.addValue("Ratio: AvgI droplets/ AvgI media", 		measures["ratioAvgIDropletsToAvgIMedia"] );
		rt.addValue("Droplet fraction",          				float(measures["ratioIDropletsToIMedia"]/ (1 + measures["ratioIDropletsToIMedia"])))

		# fill a stack with analyzed slice
		resultStack.addSlice(imp.getProcessor())

		# fill an overlay with rois of segmented regions
		if roiDroplets != None:
			# droplets 
			roiDroplets.setStrokeColor(Color.yellow)
			roiDroplets.setStrokeWidth(1);
			roiDroplets.setPosition(z)
			ovResults.add(roiDroplets)

			# bg or media   
			roiMedia.setStrokeColor(Color.magenta)
			roiMedia.setStrokeWidth(0.2);
			roiMedia.setPosition(z)
#			ovResults.add(roiMedia)    # skip

	# step 3 show table, output image and save results 
	IJ.log("Saving results to folder =" + savedir)
	rtname = "results__" + imagename
	rt.saveAs(join(savedir, rtname + ".csv"));
	
	impResultName = "results__" + imagename + ".tif"
	impResults = ImagePlus(impResultName, resultStack)
	impResults.setOverlay(ovResults)
	impResults.setCalibration(cal)
	IJ.saveAsTiff(impResults, join(savedir, impResultName));

	if displayRes:
		rt.show(rtname)
		impResults.show()


def analyzeImage(imp, threshold, Ibg):
	"""
	Main function for 2D image analysis
	Computes various params for droplet analysis.
	 
	Sum and mean values for 
	droplet (entire region covered by droplets) and
	media (entire region outside of droplets) are corrected by subtracting mean of the background image. 

	Measurements:
	- Image name 
	- Time stamp
	- Slice no processed 
	- Droplets count (nDroplets): total number of droplets found 
    - MinI background: min intensity of background image 
    - SumINorm droplets (sumIDroplets): sum intensity of droplets normalized wrt area (sum droplets/ area of image)
    - SumINorm media (sumIMedia): sum intensity of the media normalized wrt area (sum media/ area of image)
	- AvgI droplets (IAvgDroplets): avg intensity droplets (minus bg)
	- AvgI media (IAvgMedia): avg intensity media (minus bg)
	- Ratio: SumINorm droplets/ SumINorm media (ratioIDropletsToIMedia): sum intensity droplets / sum intensity media (use unnormalized values)		
	- Ratio: AvgI droplets/ AvgI media (ratioAvgIDropletsToAvgIMedia): avg intensity droplets / avg intensity media 
	- Droplet fraction: ratioIDropletsToIMedia / (1 + ratioIDropletsToIMedia) 

	All set to 0 if no droplets found.

	params: 
	imp: image to analyze (2D)
	threshold: float threshold value 
	Ibg: Int/Float mean intensity of background image 

	returns:
	measures: dict with all measurements mentioned above  
	roiDropletsafterWatershed: roi of droplets  
	roiMedia: roi of region outside of droplets (inverse of roiDropletsafterWatershed)
	
	"""
	# droplet mask creation
	impMask = createDropletMask(imp, threshold)
#	impMask.show()

	# get droplet selection (without watershed, here selecting the entire region covered by droplets)
	IJ.run(impMask, "Create Selection", "")
	roiDroplets = impMask.getRoi()
	impMask.killRoi()

	statsImage = imp.getStatistics()

	if roiDroplets != None:
		# get SumIdroplet, SumImedia (norm by surface), make the ratio
		imp.setRoi(roiDroplets)
		statsDroplets = imp.getStatistics()
		imp.killRoi();

		# droplet measures 
		IAvgDroplets = statsDroplets.mean - Ibg
		sumIDroplets = IAvgDroplets * statsDroplets.area

		# media measures 
		sumIMedia = (statsImage.mean - Ibg) * statsImage.area - sumIDroplets
		IAvgMedia = sumIMedia/(statsImage.area - statsDroplets.area)

		# normalized wrt area 
		normFactor = statsImage.area
		sumIDropletsNormalized = sumIDroplets / normFactor
		sumIMediaNormalized = sumIMedia / normFactor

		# ratio 
		ratioDropletsMedia = sumIDroplets / sumIMedia
		ratioAvgIDropletsToAvgIMedia = IAvgDroplets / IAvgMedia

		# do watershed to get accurate number of droplets 
		impWatershed = Duplicator().run(impMask)
		IJ.run(impWatershed, "Watershed", "")
		nDroplets, labelimage = getNRegion(impWatershed)

		# create droplet selection after watershed 
		IJ.run(impWatershed, "Create Selection", "")
		roiDropletsafterWatershed = impWatershed.getRoi()
		impWatershed.killRoi()

		# create background selection outside of droplets 
		IJ.run(impWatershed, "Create Selection", "")
		IJ.run(impWatershed, "Make Inverse", "");
		roiMedia = impWatershed.getRoi()
		impWatershed.killRoi()

	else:
		# no droplets found, set everything to zero
		nDroplets = 0
		
		sumIMediaNormalized = statsImage.mean-Ibg
		sumIDropletsNormalized = 0

		ratioDropletsMedia = 0
		ratioAvgIDropletsToAvgIMedia = 0

		IAvgMedia = 0
		IAvgDroplets = 0
		
		ratioAvgIDropletsToAvgIMedia = 0

		roiDropletsafterWatershed = None 
		roiMedia = None 
		

	measures = {"nDroplets":nDroplets, "sumIMedia":sumIMediaNormalized, "sumIDroplets":sumIDropletsNormalized, "ratioIDropletsToIMedia":ratioDropletsMedia, \
		"IAvgMedia":IAvgMedia,"IAvgDroplets":IAvgDroplets,"ratioAvgIDropletsToAvgIMedia":ratioAvgIDropletsToAvgIMedia}

	return measures, roiDropletsafterWatershed, roiMedia


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


def computeBackgroundImage(imp, rolling_ball_radius=20):
	"""
	Create background image using ImageJ's Subtract Background function.
	Works for both 2D image or 3D stack image. 
	For 3D, bg is computed slice by slice. 

	params: 
	imp, rolling ball radius 

	returns:
	impbg: background image 
	"""
	impbg = imp.duplicate() 

	IJ.run(impbg, "Subtract Background...", "rolling="+str(rolling_ball_radius)+" create disable stack");
	
	return impbg 


def createDropletMask(imp, threshold):
	"""
	Manual thresholding of an image using the "threshold" computed previously. 

	params:
	imp, threshold value 

	returns: 
	impMask: mask image of foreground objects (LUT inverted)
	
	"""
	impMask = Duplicator().run(imp)
	IJ.setThreshold(impMask, threshold, 65535, None)
	IJ.run(impMask, "Convert to Mask", "")
	IJ.run(impMask, "Median...", "radius=2")

#	impMask.show()

	return impMask 


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



# ===========================

main()
IJ.log("\nDONE")




