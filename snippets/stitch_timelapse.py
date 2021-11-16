# ==============================================================================================
# VERSION 1.1
# ==============================================================================================
# Description
# ==============================================================================================
# date: 2019-07-07
# author: Noreen Walker, Scientific Computing Facility, MPI-CBG
# 
# versions
# v1.1 (2019-07-07): save log files, added slight blur of initial tiles to avoid blinking in zmax projections
# v1.0 (2019-05-06): initial version
#
# === Summary ===
# This script stitches a time-lapse movie where each time points consists of several tiles (each tile a zstack). 
# It was written for .lif files (multi-photon microscopes) as input.
# The preprocessing step (tif conversion) also blurs the tiles slightly (as they are very noisy).
# The core functionality used for stitching is BigStitcher.
# The resulting stitched images are saved automatically, together with movies of the max-projected stitched images.
#
# ==== Usage: ====
# Drag the script into ImageJ & Click Run.
# * Input file: Select the .lif file. File should be a series with time-points, tiles, and a zstack per tile
# * Results folder: all results will be stored in this folder. Choose different folders for different datasets.
# * TilesX/Y: number of tiles in each direction. assumes Snake right&down pattern
# * downsampling: decreases resolution of final image. significant speedup though.
# * advanced options to skip steps: use them to restart an interrupted session at an intermediate point.
# 
# Required packages: 
# * Update Sites: BigStitcher, BioFormats (see here for instructions on how to add  an update site: https://imagej.net/Following_an_update_site)
#
# === Algorithm description ===
# 1) Conversion to tif and slight blurring. Each tile is converted to a tif file, blurred with sigma=1, and stored in the folder /1_TiffConverted/.
# 2) Stitching with BigStitcher. 
# 		* stitching all time points at the same time was buggy, therefore the script often loops over time points and handles them individually.
#		* the stitched zstack images of single time points are stored in /3_StitchedSingleTimePoints/
# 3) Movie creation: Does max projection of stitched images. Additionally creates a 2nd movie with gaussian blur (sigma=2) and then applies zmax projection. Saves both movies as tif and avi (8 bit rescaled!!) 
#		to the results directory.
#		It is recommended to continue with the .tif movie files for further analysis.
#
# === Important to know ===
# * The assumed tile pattern is: Snake right & down 
# * The assumed rough tile overlap is: 10%
# * The mysterious additional 2D tile is omitted during conversion to tif.
# * The transformations applied to the image tiles are translations (i.e. no ICP).
#
# ================================================================================================



# ================================================================================================
# Don't modify anything below here
# ================================================================================================


#@ String (visibility=MESSAGE, value="Stitching a time-lapse movie.   version 1.1") msg
#@File (label= "Input file to be stitched (typically .lif)", description="single file containing a time series with tiles, each tile being a zstack" , style= "file") inputFileName
#@File (label= "Results folder (use a different one per dataset)" , style= "directory") saveDirRoot
#@Integer (label="Number of tiles in x direction", description="assumes pattern: snake right and down") tilesX
#@Integer (label="Number of tiles in y direction", description="assumes pattern: snake right and down") tilesY
#@Integer (label="Downsampling factor for stitched image",description="decreases the x,z,y resolution of stitched images. big influence on speed. use =1 as default. Workflow suggestion: =4 for first run. Check if all good. Then =1 with all Skip boxes below checked",value=1) downsampleFused
#@ String (visibility=MESSAGE, value="   ") msg2
#@ String (visibility=MESSAGE, value="Advanced options - use them only to restart an interrupted stitching session. Hover over the fields for extra info. Default: all unchecked") msg3
#@ Boolean (label="Skip step 1: Conversion to .tif and slight blurring", description="Prerequisite: Files were converted to tif and all .tif tiles can be found in folder 1_TiffConverted",value=False) skipTiffConversion
#@ Boolean (label="Skip step 2: Conversion to BigStitcher format", description="Prerequisite: dataset.h5 and dataset.xml exist in folder 2_StitchingIntermediates",value=False) skipHdf5Conversion
#@ Boolean (label="Skip step 3: Calculation of tile-shifts", description="Prerequisite: Not easy to check. One way: all tiles have a field Stitching Transform in the xml",value=False) skipStitchingCalculation
#@ Boolean (label="Skip step 4: Calculation of stitched image", description="Prerequisite: All time points are fused and saved as .tif in 3_StitchedSingleTimePoints. If checked, will only do step 5: create movie",value=False) skipStitchingFusion


# === Further parameters ====
tileOverlap=10 # overlap of tiles (in x & y direction) in percent, initial guess for stitching. taken from microscope settings
downsampleXY=4 # or 4 # for pairwise shift calculation
mincorr=0.5 # filter pariwise shifts: min correlation required
maxdz=10 # filter pairwise shifts: max shift in z in um that is allowed. good tiles gave shifts of 2-4.5 um
sigmaPreprocess=1 # gaussian blur of the raw tiles during tif conversion
sigmaMovie=2 # additional Gaussian blur before creating the blurred max-projected movie
# ============================


from ij import IJ
from loci.plugins import BF
from ij.plugin import FolderOpener
from ij.plugin import ZProjector
from loci.formats import ImageReader
from loci.plugins.in import ImporterOptions
from ij import WindowManager
import os


def convertSeriesToTiff(filename, saveDir, prefix="", exclude2DImage=False, sigma=0):
	"""Converts all images of a series (e.g. .lif, or .czi with masterfile) into tif tiles and saves them.
	filename: (str) full path+filename of the image file (if series stored in single file), otherwise the filename of the masterfile
	saveDir: (str). Will be created if not existent
	prefix: optional prefix to the file save names. If left empty: default save names are Tile_01.tif etc.
	exclude2DImage: workaround to deal with .lif files (from multiphoton). Skips the spurious 2D image (single-plane, no time-lapse) and only converts the real tiles.
		Note: if the real tiles are also 2D this would have to be changed to a more elaborate dimension checking)
	sigma: if >0: image is blurred before being resaved
	"""

	if not os.path.exists(saveDir):
		os.makedirs(saveDir)

	# check how many files belong to the series
	bfreader=ImageReader()
	bfreader.setId(filename)
	seriesCount=bfreader.getSeriesCount()
	printLog("Number of files in series: seriesCount="+str(seriesCount))
	del bfreader

	# load images and convert them
	counter=1
	for idx in range(1,seriesCount+1): # one-based
		IJ.log("Loading series element "+str(idx));
		imp=loadOneImageOfSeries(filename,idx)
		
		if (exclude2DImage and imp.getNSlices()==1):
			IJ.log("... skipping resaving. Image is 2D and not a real tile.")
		else:
			if sigma>0:
				IJ.run(imp, "Gaussian Blur...", "sigma="+str(sigma)+" stack")

			IJ.saveAsTiff(imp, os.path.join(saveDir,prefix+"Tile_"+"{:02d}".format(counter)));
			IJ.log("... resaved image as tif.");
			counter+=1


def loadOneImageOfSeries(filename, seriesIdx, verbose=False):
	"""Loads one image of a series format (e.g. .lif, or .czi with masterfile), using the bioformats importer.
	params:
	filename: full path+filename of the image file (if series stored in single file), otherwise the filename of the masterfile
	seriesIdx: (int). which image to load from series. idx is one-based to match the names displayed in the GUI importer ("Series 1", "Series 2", etc)
	returns: ImagePlus
	"""
	# check how many files belong to the series
	bfreader=ImageReader()
	bfreader.setId(filename)
	seriesCount=bfreader.getSeriesCount()
	if verbose:
		IJ.log("Number of files in series: seriesCount:"+str(seriesCount))
	
	# intialize import options
	bfoptions=ImporterOptions()
	bfoptions.setId(filename)
	
	# specify the image to be loaded
	if seriesIdx in range(1,seriesCount+1):
		# shift from one to zero based (internal series indices are zero based, but displayed as one-based in the GUI)
		bfoptions.setSeriesOn(seriesIdx-1,True)
	else: 
		IJ.log ("Not a valid image index: "+str(seriesIdx)+ ". Skipping ...")
	
	# load image(s)
	imp=BF.openImagePlus(bfoptions)[0]

	return imp


def getFileList(datadir, extension=None, verbose=False):
	"""returns sorted list with all files (full paths) in datadir. Filtered by extension if it is not None."""
	from os import listdir
	from os.path import isfile, join
	
	if extension is None:
		filelist = [os.path.join(datadir,f) for f in sorted(listdir(datadir)) if isfile(join(datadir, f))]
	else:
		filelist = [os.path.join(datadir,f) for f in sorted(listdir(datadir)) if isfile(join(datadir, f)) and f.endswith(extension)]
	
	if verbose:
		IJ.log("Found nr of files: "+str(len(filelist)))
	
	return filelist


def getIntegerSortedFileList(datadir, extension=None, verbose=False):
	"""Similar to getFileList but it does a special sorting (needed for stitcher output filenames) by the increasing integer values within the filename.
	It takes thereby care of missing 0 padding of the filenames.
	["f1_ch0","f2_ch0","f10_ch0"] is returned as ['f1_ch0', 'f2_ch0', 'f10_ch0']
	"""

	from os import listdir
	from os.path import isfile, join

	if extension is None:
		basenames = [str(f) for f in (listdir(datadir)) if isfile(join(datadir, f))]
	else:
		basenames = [str(f) for f in (listdir(datadir)) if isfile(join(datadir, f)) and f.endswith(extension)]

	#print basenames
	basenames.sort(key=lambda f: int(filter(str.isdigit, f)))

	filelist=[os.path.join(datadir,f) for f in basenames]

	if verbose:
		IJ.log("Found nr of files: "+str(len(filelist)))
	
	return filelist


def addSlice(impStack, impTimePoint):
	"""Adds a 2D slice to an image stack (width,height,slices). Note: the third dim is later converted to time axis
	"""
	if impStack is None: # first iteration
		impStack=impTimePoint
	else:
		stack = impStack.getStack() 
		stack.addSlice(impTimePoint.getProcessor())		
		impStack.setStack(stack);
	return impStack


def printLog(string):
	"""prints logging info to two outputs"""
	print string
	IJ.log(string);



def createProjectionMovies(inputDir, saveDir, imageName="",sigma=2, specialIntegerSorting=False):
	"""Creates a movie from  a folder of (stitched) zstack images with single time point per image file. The movie is a zmax projection, and a 
	second flavor is first blurred with isgma, and then max projected.
	The movies are saved in saveDir as tif & avi, and are also returned.

	Params:
	inputDir: path (string) of the single time point zstacks. typically they should be fused tiles. images must be tif.
	saveDir: movie save location (string)
	imageName: suffix for save name of movies
	sigma: gaussian blur
	specialIntegerSorting: if true: sorts the filenames according to increasing integer in name. naeeded for bigstitcher. if false: standard sorting
	
	Returns:
	movieZMax
	movieBlurred	
	"""

	if specialIntegerSorting:
		imageFiles= getIntegerSortedFileList(inputDir,"tif")
	else:
		imageFiles=getFileList(inputDir,"tif")
		
	movieZMax=None
	movieBlurred=None
	
	for idx in range(len(imageFiles)):
		imp=BF.openImagePlus(imageFiles[idx])[0] 
		
		imp2=imp.duplicate();

		impZMax = ZProjector.run(imp,"max")
		impZMax.killRoi() 

		IJ.run(imp2, "Gaussian Blur...", "sigma="+str(sigma)+" stack")
		impZMaxBlurred=ZProjector.run(imp2,"max")
		
		# collect time-frames as z-stack
		movieZMax=addSlice(movieZMax, impZMax)
		movieBlurred=addSlice(movieBlurred, impZMaxBlurred)
		
	# convert slices to frames
	movieZMax.setDimensions(1, 1, movieZMax.getNSlices()) # nChannels, nSlices, nFrames
	movieBlurred.setDimensions(1, 1, movieBlurred.getNSlices())
	
	# save movies
	IJ.run(movieZMax, "AVI... ", "compression=None frame=5 save=["+os.path.join(saveDir,"MaxProjected_Stitched_8bit_"+imageName+".avi")+"]");
	IJ.saveAsTiff(movieZMax,os.path.join(saveDir,"MaxProjected_Stitched_"+imageName))
	IJ.run(movieBlurred, "AVI... ", "compression=None frame=5 save=["+os.path.join(saveDir,"MaxProjected_Stitched_Blurred_8bit_"+imageName+".avi")+"]");
	IJ.saveAsTiff(movieBlurred,os.path.join(saveDir,"MaxProjected_Stitched_Blurred_"+imageName))

	return movieZMax, movieBlurred



# ====================================================================================================================================
def main():
	# == define & create save directories ==
	baseName=os.path.splitext(os.path.basename(inputFileName.getPath()))[0] # base name without extension
	
	saveDirTiff=os.path.join(saveDirRoot.getPath(),"1_TiffConverted")
	saveDirStitchingIntermediates=os.path.join(saveDirRoot.getPath(),"2_StitchingIntermediates")
	saveDirStitchedTimpoints=os.path.join(saveDirRoot.getPath(),"3_StitchedSingleTimePoints")
	saveDirFinal=saveDirRoot.getPath()
	
	if not os.path.exists(saveDirTiff):
			os.makedirs(saveDirTiff)
	if not os.path.exists(saveDirStitchingIntermediates):
			os.makedirs(saveDirStitchingIntermediates)
	if not os.path.exists(saveDirStitchedTimpoints):
			os.makedirs(saveDirStitchedTimpoints)
	if not os.path.exists(saveDirFinal):
			os.makedirs(saveDirFinal)
	
	
	
	################################
	# ==preprocessing: convert tiles to tif ==
	if not skipTiffConversion:
		printLog("\nConverting series "+inputFileName.getPath()+" to tiff (slow).")
		convertSeriesToTiff(inputFileName.getPath(), saveDirTiff,  exclude2DImage=True, sigma=sigmaPreprocess)
		printLog("\nFinished tif conversion.")
	
	################################
	# extract number of timepoints from a tiff file
	if not (skipHdf5Conversion and skipStitchingCalculation and skipStitchingFusion):
		fn=getFileList(saveDirTiff,"tif")[0]
		imp=BF.openImagePlus(fn)[0] 
		numTimePoints=imp.getNFrames();
		del imp, fn
		
	# == stitching ==
	if not skipHdf5Conversion:
		printLog("\nStarting stitching. \nRunning dataset conversion for BigStitcher (slow)")
		
		# --a) convert dataset==
		IJ.run("Define dataset ...", "define_dataset=[Automatic Loader (Bioformats based)] project_filename=dataset.xml path=["+saveDirTiff+"] exclude=10 "+\
		"pattern_0=Tiles move_tiles_to_grid_(per_angle)?=[Move Tile to Grid (Macro-scriptable)] grid_type=[Snake: Right & Down      ] tiles_x="+str(tilesX)+\
		" tiles_y="+str(tilesY)+" tiles_z=1 overlap_x_(%)="+str(tileOverlap)+" overlap_y_(%)="+str(tileOverlap)+" overlap_z_(%)=10 keep_metadata_rotation how_to_load_images=[Re-save as multiresolution HDF5] "+\
		"dataset_save_path=["+saveDirStitchingIntermediates+"] check_stack_sizes subsampling_factors=[{ {1,1,1}, {2,2,1}, {4,4,1} }] hdf5_chunk_sizes=[{ {32,32,4}, {32,32,4}, {16,16,16} }] "+\
		"timepoints_per_partition=1 setups_per_partition=0 use_deflate_compression export_path=["+os.path.join(saveDirStitchingIntermediates,"dataset")+"]");

		printLog("Finished dataset conversion.")

		
	if not skipStitchingCalculation:
		# --b) pairwise shifts: process each time point individually--
		printLog("\nCalculating pairwise shifts")
		
		# loop manually over time points
		for timePoint in range(numTimePoints): # zero based
			IJ.log("\ntime point: "+str(timePoint))
			IJ.run("Calculate pairwise shifts ...", "select=["+os.path.join(saveDirStitchingIntermediates,"dataset.xml")+"] process_angle=[All angles] "+\
			"process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[Single Timepoint (Select from List)] "+\
			"processing_timepoint=[Timepoint "+str(timePoint)+"] method=[Phase Correlation] downsample_in_x="+str(downsampleXY)+" downsample_in_y="+str(downsampleXY)+\
			" downsample_in_z=1");
		
		# NOTE This should work but does not (see my submitted github issue): compute for all time points and treat them individually:
		# does not work: IJ.run("Calculate pairwise shifts ...", "select=["+os.path.join(saveDirStitchingIntermediates,"dataset.xml")+"] process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[All Timepoints] method=[Phase Correlation] downsample_in_x="+str(downsampleXY)+" downsample_in_y="+str(downsampleXY)+" downsample_in_z=1");
		
		
		# --c) filter pairwise shifts--
		printLog("Filtering pairwise shifts")
		dummyMaxDxy=100000 # random high value
		IJ.run("Filter pairwise shifts ...", "select=["+os.path.join(saveDirStitchingIntermediates,"dataset.xml")+"] filter_by_link_quality min_r="+str(mincorr)+" max_r=1 "+\
		"filter_by_shift_in_each_dimension max_shift_in_x="+str(dummyMaxDxy)+" max_shift_in_y="+str(dummyMaxDxy)+" max_shift_in_z="+str(maxdz)+" max_displacement=0");
		
		# --d) global optimization--
		printLog("Doing global optimization")
		for timePoint in range(numTimePoints): # zero based
			IJ.log("\ntime point"+str(timePoint))
			IJ.run("Optimize globally and apply shifts ...", "select=["+os.path.join(saveDirStitchingIntermediates,"dataset.xml")+"] process_angle=[All angles] "+\
			"process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[Single Timepoint (Select from List)] "+\
			"processing_timepoint=[Timepoint "+str(timePoint)+"] relative=2.500 absolute=3.500 global_optimization_strategy=[Two-Round using Metadata to align unconnected Tiles] "+\
			"fix_group_"+str(timePoint)+"-0");
		
		# --e) icp--
		# ICP: For fine-alignment. Deactivated because ICP is a full-affine Transformation and introduced rotations of a few degrees to some tiles. 
		# Rotations could be an artifact of high noise/low signal features which are accidentally matched. Unlikely to have really rotated tiles.
		# Note: ICP is interest-point based
		# Update: ICP Expert Mode: option to choose the transformation model (translation, rigid, affine)
		#doICP=False
		#if doICP:
		#	print "doing ICP";
		#	for timePoint in range(numTimePoints): 
		#		IJ.run("ICP Refinement ...", "select=["+os.path.join(saveDirStitchingIntermediates,"dataset.xml")+"] process_angle=[All angles] process_channel=[All channels] "+\
		#		"process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[Single Timepoint (Select from List)] "+\
		#		"processing_timepoint=[Timepoint "+str(timePoint)+"] icp_refinement_type=[Simple (tile registration)] downsampling=[Downsampling 8/8/4] interest=[Average Threshold] "+\
		#		"icp_max_error=[Normal Adjustment (<5px)]");

		printLog("Finished shift calculations (pairwise shifts and global optimization).\n")

		# save logfile
		logdata=IJ.getLog();
		if logdata is not None:
			IJ.saveString(logdata, os.path.join(saveDirStitchingIntermediates,"logfile_1_postshiftcalc.txt"));

	if not skipStitchingFusion:
		# --f) fuse and export--
		printLog("\nExporting stitched image (slow)")
		for timePoint in range(numTimePoints):
			IJ.run("Fuse dataset ...", "select=["+os.path.join(saveDirStitchingIntermediates,"dataset.xml")+"] process_angle=[All angles] process_channel=[All channels] "+\
			"process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[Single Timepoint (Select from List)] "+\
			"processing_timepoint=[Timepoint "+str(timePoint)+"] bounding_box=[All Views] downsampling="+str(downsampleFused)+" pixel_type=[16-bit unsigned integer] "+\
			"interpolation=[Linear Interpolation] image=[Precompute Image] interest_points_for_non_rigid=[-= Disable Non-Rigid =-] blend preserve_original "+\
			"produce=[Each timepoint & channel] fused_image=[Save as (compressed) TIFF stacks] output_file_directory=["+saveDirStitchedTimpoints+"]");

		# save logfile
		logdata=IJ.getLog();
		if logdata is not None:
			IJ.saveString(logdata, os.path.join(saveDirStitchingIntermediates,"logfile_2_postfusion.txt"));
	
		printLog("Finished exporting of stitched image.\n")
		
	################################
	# == postprocessing: create projection movies ==
	printLog("\nCreating projection movies")
	suffix=baseName
	movieZMax,movieZMaxBlurred=createProjectionMovies(saveDirStitchedTimpoints,saveDirFinal,imageName=suffix,sigma=sigmaMovie, specialIntegerSorting=True)
	
	movieZMax.show()
	movieZMaxBlurred.show()
			
	printLog("\nDone.")



# ======================================================================================
main()
