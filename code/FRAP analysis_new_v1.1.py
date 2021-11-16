#@ Boolean (label="extract image timestamp") extractTimeStamp
#@ Integer (label="proces until frame")n_slices0

import java.awt.Color as Color
from ij import WindowManager as WindowManager
from ij.plugin.frame import RoiManager as RoiManager
from ij.process import ImageStatistics as ImageStatistics
from ij.measure import Measurements as Measurements
from ij import IJ as IJ
from ij.measure import CurveFitter as CurveFitter
from ij.gui import Plot as Plot
from ij.gui import PlotWindow as PlotWindow
import math, os


# modified script from the version available on imagej wiki: https://imagej.net/Analyze_FRAP_movies_with_a_Jython_script
# version 1.1: 
#   + modify the normalisation to take into account the minimum signa measured at the frap point
#   - add the possibility to measure the time from slice timestamp
#   - create a result table to export results

resultsTableName = "frap measures"


def main( n_slices ) :
	# Get ROIs
	roi_manager = RoiManager.getInstance()
	roi_list    = roi_manager.getRoisAsArray()
	 
	# We assume first one is FRAP roi, the 2nd one is normalizing roi.
	roi_FRAP    = roi_list[0];
	roi_norm    = roi_list[1];
	 
	 
	# Get current image plus and image processor
	current_imp  = WindowManager.getCurrentImage()
	stack        = current_imp.getImageStack()
	calibration  = current_imp.getCalibration()

	# Specify up to what frame to fit and plot.
	if n_slices > current_imp.getStackSize() :
		n_slices = current_imp.getStackSize()
	print n_slices
	
	#############################################
	 
	# Collect intensity values
	 
	# Create empty lists of number
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
	  
	# Gather image parameters
	if extractTimeStamp :
		timeValues, time_units = getTimeStamps(current_imp)
		timeValues = timeValues[0:n_slices]
		print timeValues
		print n_slices
	else :
		frame_interval = calibration.frameInterval
		timeValues = [i * frame_interval for i in range( n_slices ) ]
		time_units = calibration.getTimeUnit()
		
		
	#IJ.log('For image ' + current_imp.getTitle() )
	#IJ.log('Time interval is ' + str(frame_interval) + ' ' + time_units)
	  
	# Find minimal intensity value in FRAP and bleach frame
	min_intensity = min( If )
	bleach_frame = If.index( min_intensity )
	IJ.log('FRAP frame is ' + str(bleach_frame+1) + ' at t = ' + str(timeValues[bleach_frame]) + ' ' + time_units )
	  
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
	 
	  
	plot = Plot("Normalized FRAP curve for " + current_imp.getTitle(), "Time ("+time_units+')', "NU", [], [])
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
	 
	str1 = ('Half-recovery time = %.2f ' + time_units) % thalf
	IJ.log( str1 )
	str2 = "Mobile fraction = %.1f %%" % (100 * mobile_fraction)
	IJ.log( str2 )


	rt = getResultsTable( resultsTableName )
	rt.incrementCounter()
	rt.addValue("Image name", current_imp.getTitle() )
	rt.addValue("half time", thalf )
	rt.addValue("mobile fraction (%)", 100 * mobile_fraction )
	rt.addValue("R^2", fitter.getFitGoodness() )
	rt.addValue("time unit", time_units)

	rt.show( resultsTableName )
			


from ij.measure import ResultsTable
from ij import WindowManager;
from ij.text import TextWindow;

def getResultsTable(name, reset=False):
	win = WindowManager.getWindow(name)
	if ( (win!=None) & isinstance(win,TextWindow) ) :	
		rt = win.getTextPanel().getResultsTable();
		if reset:
			rt.reset()
	else:
		rt = ResultsTable();
	return rt





from loci.formats import ImageReader
from loci.formats import MetadataTools

def getTimeStamps(imp0):
	fileInfo = imp0.getOriginalFileInfo()
	fileName = None
	timeUnit = None
	if not fileInfo == None :
		if (not fileInfo.directory==None) & (not fileInfo.fileName==None):
			fileName= os.path.join(fileInfo.directory, fileInfo.fileName)
		print fileName
	
	else :
		print 'the time stamp could not be extracted'
		print 'default time interval will be used'
		nSlice = imp0.getStackSize()
		dt = imp0.getCalibration().frameInterval
		timeUnit = imp0.getCalibration().getTimeUnit()
		if dt == 0 :
			dt =1
			timeUnit = 'frame'
		timevalues = range(0,nSlices*dt)
		
		return timevalues, timeUnit
	
		
	reader = ImageReader()
	omeMeta = MetadataTools.createOMEXMLMetadata()
	reader.setMetadataStore(omeMeta)
	reader.setId(fileName)
	reader.close()

	
	
	nPlane = omeMeta.getPlaneCount(0) # we assume there is only one time series

	print nPlane
	
	timeValues = []
	for plane in range(0,nPlane):
		fr = omeMeta.getPlaneTheT(0,plane)
		ch = omeMeta.getPlaneTheC(0,plane)
		if ch.getNumberValue() == 0 :
			timeValues.append( omeMeta.getPlaneDeltaT(0,plane).value() )
			timeUnit = omeMeta.getPlaneDeltaT(0,plane).unit().getSymbol()

	if timeUnit is None:
		timeUnit = imp0.getCalibration().getTimeUnit()
	
	#print timeValues[]

	return timeValues, timeUnit



main( n_slices0 )



