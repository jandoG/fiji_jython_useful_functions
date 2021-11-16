# ==============================================================================================
# VERSION 1.0
# ==============================================================================================
# Description
# ==============================================================================================
# date: 2019-09-03
# author: Noreen Walker, Scientific Computing Facility, MPI-CBG
# 
# versions:
# v1.0 2019-09-03: initial version
#
# === Summary ===
# This scripts concatenates all stitched single timepoints into a timelapse movie (3D).
# It should be run after stitch_timelapse.py
# Note: the stitch_timelapse.py also creates a movie, but that movie is a max-projection (optionally blurred) while here the 3D+t movie is created.
#
# === Usage ===
# Drag the script into ImageJ & Click Run.
# * Hover over the text fields to get additional tips.
#
# ================================================================================================


# ================================================================================================
# Don't modify anything below
# ================================================================================================
#@ String (visibility=MESSAGE, value="Concatenate stitched time points to movie (3d+t).   version 1.0") msg
#@File (label= "Raw input file that was stitched (typically .lif)", description="only needed for proper output file name", style= "file") rawFile
#@File (label= "Folder of stitched single timepoints", description="typically /somepath/3_StitchedSingleTimePoints/" , style= "directory") inputDirTimePoints
#@File (label= "Folder to save results", description="typically parent folder of the single-time-points folder" , style= "directory") saveDir

extension="tif"


from ij import IJ
from loci.plugins import BF
from ij.plugin import FolderOpener
from ij.plugin import HyperStackConverter
import os



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





# ====================================================================================================================================
def main():

	# == define & create save directories ==
	baseName=os.path.splitext(os.path.basename(rawFile.getPath()))[0] # base name without extension

	assert os.path.exists(inputDirTimePoints.getPath())
		
	if not os.path.exists(saveDir.getPath()):
			os.makedirs(saveDir.getPath())

	# extract number of slices
	fn=getIntegerSortedFileList(inputDirTimePoints.getPath(),extension)[0]
	imp=BF.openImagePlus(fn)[0] 
	nSlices=imp.getNSlices()
	del imp

	# load movie as virtual stack (3dims)	
	imp = FolderOpener.open(inputDirTimePoints.getPath()+os.sep, "virtual"); # needs explicit filesep

	# re-order to 3d+t
	nSlicesTmp=imp.getNSlices();
	nFrames=int(nSlicesTmp/nSlices)
	imp = HyperStackConverter.toHyperStack(imp, 1, nSlices, nFrames, "Color");

	print "Starting saving"
	IJ.saveAsTiff(imp, os.path.join(saveDir.getPath(),"Movie3D_Stitched_"+baseName+".tif"));

	imp.show();
	print "Done"
			


# execute
main()