############## README #############

# Author Gayathri Nadar, SCF MPI-CBG, nadar@mpi-cbg.de
# V0.0: Scripting "register virtual stack slices" plugin to register images inside a folder
# V0.1: Switching the reference to the unstained image in order to know what transformation is needed to move from stained->unstained
# V0.2: Use bioformats to open images and save them to disk again in proper calibrated format. Then use this folder as source folder to register images.
#		Option to also provide feature selection model. More options for feature extraction and registration.

"""
# Description:
 	- The script uses a reference image selected by the user inside a folder and registers all the other images in the folder 
 	to this reference image. 
 	- For our case, reference image = unstained image

# Parameters:
	- Feature extraction model: The expected transformation model finding inliers (i.e. correspondences or landmarks between images) 
		in the feature extraction: translation, rigid, similarity or affine.
	- Registration model. The image transformation model: translation, rigid, similarity, affine, elastic or moving least squares.

# Usage:
	- Create a folder with all the images to be registered along with the reference image.
	- Open this script with FIJI and click run.
	- In the pop-up window, 
		1. select the folder which contains the images to be registered and 
		2. reference image 
		and click ok.

# Comments:
	- Make sure you have "Register virtual stack slices" plugin installed in FIJI.
	- It comes with FIJI and can be found under Plugins > Registration > Register virtual stack slices
	- For some reason if it is not found, update FIJI by using Help > Update...

"""
##################################

#@File (label= "Source folder containing stained and unstained image to process" , style= "directory", required= "True") source
#@File (label= "Reference file to register to (use unstained one)", style= "extension:tif", required= "True") processFile
#@String (label="Feature extraction model? 0=TRANSLATION, 1=RIGID, 2=SIMILARITY, 3=AFFINE", choices={"0", "1", "2", "3"}, value = "2") fr_choice
#@String (label="Transformation model? 0=TRANSLATION, 1=RIGID, 2=SIMILARITY, 3=AFFINE, 4=ELASTIC, 5=MOVING_LEAST_SQUARES", choices={"0", "1", "2", "3", "4", "5"}, value = "2") tr_choice

from register_virtual_stack import Register_Virtual_Stack_MT
from ij import IJ
from ij.io import FileSaver
from ij.plugin import ChannelSplitter

from loci.plugins import BF
from loci.plugins.prefs import OptionsList
from loci.plugins.in import ImporterOptions

import os
from os.path import isfile, join

output_subdir= "output/" # aligned images. dir will be created
transforms_subdir = "transforms/" # transformations in xml format. dir will be created

extension="tif"

IJ.log("\\Clear")

def main():
	# close windows
	dummyimp=None 
	IJ.run(dummyimp, "Close All", "")

	options = ImporterOptions()
	options.setColorMode(ImporterOptions.COLOR_MODE_COLORIZED)
	options.setAutoscale(True)
	options.setStackFormat("Hyperstack")

	# read images in folder using bio-formats, get ch1 and save them inside another folder
	source_dir = source.getPath()

	IJ.log("Source folder =" + source_dir)
	
	source_process_dir = join(source_dir, 'to_register', '') # '' to add trailing slash
	if not os.path.exists(source_process_dir):
		os.makedirs(source_process_dir)

	files = [f for f in os.listdir(source_dir) if isfile(join(source_dir, f)) and f.endswith(extension)]
	for f in files:
		full_path = join(source_dir, f)
		imp = BF.openImagePlus(full_path)[0]
		impName = os.path.splitext(imp.getTitle())[0]
		
		nchannels = imp.getNChannels();
		if nchannels > 1:
			imps = ChannelSplitter().split(imp)
			imp_to_save = imps[0]
			IJ.run(imp_to_save, "Grays", "");

			outputName = "C1-" + impName + '.tif'
			FileSaver(imp_to_save).saveAsTiff(join(source_process_dir, outputName))

		else:
			IJ.run(imp, "Grays", "");
			outputName = "C1-" + impName + '.tif'
			FileSaver(imp).saveAsTiff(join(source_process_dir, outputName))

	output_dir=os.path.join(source_process_dir,output_subdir)
	transforms_dir=os.path.join(source_process_dir,transforms_subdir)
		
	# initialize savedirs
	if not os.path.exists(output_dir): 
		os.makedirs(output_dir)
	if not os.path.exists(transforms_dir): 
		os.makedirs(transforms_dir)

	refFile= processFile.getPath()
	IJ.log("Reference image to register to (unstained image) =" + refFile);
	
	reference_name= refFile.split(os.sep)[-1] # basename
	ref_name = "C1-" + reference_name
#	print ref_name
	print "Processing folder: ", source_process_dir," Aligning to reference: ", ref_name

	use_shrinking_constraint = 0
	p = Register_Virtual_Stack_MT.Param()

	p.featuresModelIndex=int(fr_choice) # featureModelIndex: 0=TRANSLATION, 1=RIGID, 2=SIMILARITY, 3=AFFINE
	p.registrationModelIndex=int(tr_choice) # registrationModelIndex: (0=TRANSLATION, 1=RIGID, 2=SIMILARITY, 3=AFFINE, 4=ELASTIC, 5=MOVING_LEAST_SQUARES)
	
	for i, fr in zip(range(4), ['TRANSLATION', 'RIGID', 'SIMILARITY', 'AFFINE']):
		if int(fr_choice) == i:
			IJ.log("Feature model used =" + fr)
	
	for j, tr in zip(range(6), ['TRANSLATION', 'RIGID', 'SIMILARITY', 'AFFINE', 'ELASTIC', 'MOVING_LEAST_SQUARES']):
		if int(tr_choice) == j:
			IJ.log("Transform model used =" + tr)

	# do registration
	Register_Virtual_Stack_MT.exec(source_process_dir, output_dir, transforms_dir, ref_name, p, use_shrinking_constraint)

	IJ.log("Aligned images found in =" + output_dir);
	IJ.log("Image transform xml files found in =" + transforms_dir);

	logdata=IJ.getLog();
	if logdata is not None:
		IJ.saveString(logdata, os.path.join(source_dir,"logfile_registration.txt"));
	
	print "Done"


main()