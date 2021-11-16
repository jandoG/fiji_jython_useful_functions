# ============================================= README ======================================================
#
# Author: Gayathri Nadar, SCF MPI-CBG, nadar@mpi-cbg.de
#
# Version history:
# v0.0: 02.03.2021 FRAP analysis for 2D time-lapse movie 
#				   modified script from the version available on imagej wiki: https://imagej.net/Analyze_FRAP_movies_with_a_Jython_script
# v0.1: 10.03.2021 add diffusion coeffient measure = D (length^2/time-unit) = 0.224 * w^2 / t-half  w = radius of frap roi 
#				   table with intensities before normalized 
#
# Description:
#	This script does FRAP analysis for user-drawn FRAP and normalizing ROI.
#	Assumptions: Image is 2D + time 
#
# Steps:
#	- The image is opened with Bio-formats and displayed. 
#	- The user is asked for draw a Roi around FRAP object and normalizing object. 
#	- The intensities are obtained for these two Rois 
#	- Curve is fit, and FRAP params are computed.
#
# Usage:
#	- Open this script with FIJI and click run. 
#	- A window will pop-up, enter the necessary values and click ok to run the script.
#
# Parameters to be set:
#	- File to process 
#	- Last time point to process: The movie will be processed only till this frame.
#
# Output:
#	Results are stored in folder imagename_results. 
#	Results saved:
#		- image with overlay of rois
#		- log window as txt file 
#		- results tables as csv file - frap measures and intensities before normalization
#		- plot as png file
#
# Frap Measures notes: http://frapbot.kohze.com/features.html
#
# =============================================================================================================



# =============================================================================================================
# Don't modify anything below here without good reason
# =============================================================================================================


#@ File (label= "Select an image to process", style="file") impFile
#@ Integer (label= "Process until time frame") n_slices0

import java.awt.Color as Color
from ij import WindowManager as WindowManager
from ij.gui import Roi, PointRoi, WaitForUserDialog, Overlay
from ij.plugin.frame import RoiManager as RoiManager
from ij.process import ImageStatistics as ImageStatistics
from ij.measure import Measurements as Measurements
from ij.io import FileSaver
from ij import IJ as IJ
from ij.measure import CurveFitter as CurveFitter
from ij.gui import Plot as Plot
from ij.gui import PlotWindow as PlotWindow
from ij.measure import ResultsTable
from ij.text import TextWindow;
import loci.plugins
from loci.plugins import BF
from loci.plugins.in import ImporterOptions 
from loci.formats import ImageReader
from loci.formats import MetadataTools
import math, os, sys 
from os.path import isfile, join


resultsTableName = "frap measures"
rt = ResultsTable()
rt1 = ResultsTable()   # for intensity measures 
overlay = Overlay()
IJ.log("\\Clear")

rm = RoiManager.getInstance()
if not rm:
	rm = RoiManager()
rm.reset();

def main(n_slices):
	current_imp = openImageWithBF(impFile.getPath(), virtual= False, groupfiles = True)
	impname = os.path.splitext(current_imp.getTitle())[0]

	if current_imp.getNSlices() > 1 or current_imp.getNFrames() == 1:
		IJ.log("Image is not 2D or has only 1 time point. Required image = 2D time-lapse. Try with a different image!")
		sys.exit()		

	results_path = join(os.path.dirname(impFile.getPath()), impname + '_frap_results')
	if not os.path.exists(results_path):
		os.makedirs(results_path)	

	# Get current image plus and image processor
	stack        = current_imp.getImageStack()
	calibration  = current_imp.getCalibration()

	# Specify up to what frame to fit and plot
	if n_slices > current_imp.getStackSize():
		n_slices = current_imp.getStackSize()

	# Get ROIs
	current_imp.killRoi()
	current_imp.show()
	
	WaitForUserDialog("Action required", "Please draw a ROI around FRAP object and press ok.").show() 
	rm.addRoi(current_imp.getRoi());
	current_imp.killRoi()
	
	WaitForUserDialog("Action required", "Please draw a ROI around normalizing object and press ok.").show()
	rm.addRoi(current_imp.getRoi()); 
	current_imp.killRoi()

	current_imp.hide()

	# rename rois for convenience 
	roi_list = rm.getRoisAsArray()
	rm.rename(0, "frap");
	rm.rename(1, "normalizing");
	
	# We assume first one is FRAP roi, the 2nd one is normalizing roi.
	roi_FRAP    = roi_list[0];
	roi_norm    = roi_list[1];

	# add to overlay 
	overlay.add(roi_FRAP);
	overlay.add(roi_norm);

	# save the rois 
	rm.runCommand("Save", join(results_path, impname + "_RoiSet.zip"));
 
	# Collect intensity values (add these to table)
	If = []  # Frap values 
	In = []  # Norm values
	 
	# Loop over each slice of the stack
	for i in range(0, n_slices):
	  
	    # Get the current slice
	    ip = stack.getProcessor(i+1)
	  
	    # Put the ROI on it
	    ip.setRoi(roi_FRAP)
	  
	    # Make a measurement in it
	    stats = ImageStatistics.getStatistics(ip, Measurements.MEAN, calibration);
	    mean  = stats.mean
	  
	    # Store the measurement in the list
	    If.append( mean  )
	 
	    # Do the same for non-FRAPed area
	    ip.setRoi(roi_norm)
	    stats = ImageStatistics.getStatistics(ip, Measurements.MEAN, calibration);
	    mean = stats.mean
	    In.append( mean  )

	# get time params 
	frame_interval = calibration.frameInterval
	timeValues = [i * frame_interval for i in range( n_slices ) ]
	time_units = calibration.getTimeUnit()
	IJ.log('Image name = ' + current_imp.getTitle())
	IJ.log('Time interval = ' + str(frame_interval) + ' ' + time_units)

	# get image params 
	image_units = calibration.getUnit()
	  
	# Find minimal intensity value in FRAP and bleach frame
	min_intensity = min( If )
	bleach_frame = If.index( min_intensity )
	IJ.log('FRAP frame = ' + str(bleach_frame+1) + ' at t = ' + str(timeValues[bleach_frame]) + ' ' + time_units )
	  
	# Compute mean pre-bleach intensity
	mean_If = 0.0
	mean_In = 0.0
	for i in range(bleach_frame):         # will loop until the bleach time
	    mean_If = mean_If + If[i]
	    mean_In = mean_In + In[i]
	mean_If = mean_If / bleach_frame
	mean_In = mean_In / bleach_frame
	  
	# Calculate normalized curve
	normalized_curve = []
	for i in range(n_slices):
	    normalized_curve.append( (If[i] - min_intensity) / (mean_If - min_intensity)   *   (mean_In - min_intensity) / (In[i]- min_intensity) )
	     
	# plot 
	x = timeValues
	y = normalized_curve
	 
	xtofit = x[ bleach_frame : n_slices ]
	x0 = xtofit[0]
	xtofit = [ xt-x0 for xt in xtofit]
	ytofit = normalized_curve[ bleach_frame : n_slices ]
	  
	# Fitter
	fitter = CurveFitter(xtofit, ytofit)
	fitter.doFit(CurveFitter.EXP_RECOVERY_NOOFFSET)
	IJ.log("Fit FRAP curve by " + fitter.getFormula() )
	param_values = fitter.getParams()
	IJ.log( fitter.getResultString() )
	  
	# Overlay fit curve, with oversampling (for plot)
	#xfit = [ (t / 10.0  + bleach_frame) * frame_interval for t in range(10 * len(xtofit) ) ]
	xfit = x[ bleach_frame : n_slices]
	yfit = []
	for xt in xfit:
	    yfit.append( fitter.f( fitter.getParams(), xt - xfit[0]) )
	 
	plot = Plot("Normalized FRAP curve for " + impname, "Time ("+time_units+')', "Normalized Intensity (a.u.)")
	plot.setLimits(0, max(x), 0, 1.2 );
	plot.setLineWidth(2)
	 
	plot.setColor(Color.BLACK)
	plot.addPoints(x, y, Plot.LINE)
	plot.addPoints(x,y,PlotWindow.X);
	  
	plot.setColor(Color.RED)
	plot.addPoints(xfit, yfit, Plot.LINE)
	 
	plot.setColor(Color.black);
	plot_window =  plot.show()
	 
	# Output FRAP parameters
	thalf = math.log(2) / param_values[1]
	mobile_fraction = param_values[0]

	# diffusion coefficient D length^2/time-unit = 0.224 * w^2 / t-half  w = radius of frap roi 
	current_imp.setRoi(roi_FRAP)
	stats = current_imp.getStatistics(Measurements.ALL_STATS);
	diameter = float((stats.minor + stats.major)/2)
	radius = float(diameter/ 2)
	diffusion_coeff = 0.224 * (radius**2) / thalf 
	current_imp.killRoi()
	 
	str1 = ('Half-recovery time = %.2f ' + time_units) % thalf
	IJ.log( str1 )
	str2 = "Mobile fraction = %.1f %%" % (100 * mobile_fraction)
	IJ.log( str2 )

	# results table
	rt.incrementCounter()
	rt.addValue("Image name", current_imp.getTitle() )
	rt.addValue("half time", thalf )
	rt.addValue("mobile fraction (%)", 100 * mobile_fraction )
	rt.addValue("diffusion coefficient " + image_units + "^2/" + time_units, diffusion_coeff)
#	rt.addValue("R^2", fitter.getFitGoodness() )
	rt.addValue("time unit", time_units)
	rt.addValue("image unit", image_units)
	rt.show( resultsTableName )
	rt.saveAs((join(results_path, impname + "_" + resultsTableName + ".csv")));

	# intensity results 
	for idx, (i, j) in enumerate(zip(If, In)):
		rt1.incrementCounter() 
		rt1.addValue("Image name", current_imp.getTitle())
		rt1.addValue("Frame no", idx + 1)
		rt1.addValue("Frap roi mean intensity", i)
		rt1.addValue("Norm roi mean intensity", j)
	rt1.show( "intensities before normalization" )
	rt1.saveAs((join(results_path, impname + "_" + "intensitiesBeforeNorm" + ".csv")));
	
	# save log 
	IJ.selectWindow("Log");
	IJ.saveAs("Text", join(results_path, impname + "_log.txt"));

	# results images 
	imp_result = current_imp.duplicate()
	imp_result.setOverlay(overlay)
	imp_result.setTitle(impname + "-rois")
	imp_result.show()
	FileSaver(imp_result).saveAsTiff(join(results_path, imp_result.getTitle() + ".tif"))
	FileSaver(plot.getImagePlus()).saveAsTiff(join(results_path, impname + "_normalized_FRAP_curve.png"))


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

			
main(n_slices0)

print("Done")



