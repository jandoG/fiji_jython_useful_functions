#@File (label="Frap data with FRAP, Non-FRAP and Normalized intensities") table
#@Boolean (label="Compute normalized values (checked), get it from file (unchecked)") compute_normValues

# Author: Gayathri Nadar, SCF MPI_CBG nadar@mpi-cbg.de
# Version: V0.0 19.07.18

# Description:
# This scripts fits a curve to the FRAP data usig the CurveFitter in FIJI. 


# Usage:
# Please make sure your data has . for decimal instead of comma. 
# In MacBook this can be done by going to System preferences > Language & Region > advanced..> General and selecting , in Grouping and . in Decimal
# Save "only the required" data as txt file or csv file, since this script uses tab/comma seperated files only.
# The first lines of the code define the required column names for plotting the data. Enter the names accordingly and click run. 

time_column = "Time"
frap_column = "RIMZONE_FRAP"
nonfrap_column = "UNFRAP/REFERENCE ROI"
normalize_column = "RIM ZONE_FRAP _NORMALISED"
time_units = "secs"

from ij import IJ
from ij.measure import ResultsTable
from ij.io import Opener
import java.awt.Color as Color
from ij import WindowManager 
from ij.measure import CurveFitter 
from ij.gui import Plot 
from ij.gui import PlotWindow 
import math, os
from os.path import join

rt = ResultsTable()
table = rt.open(table.getPath())
table.show("Frap Data")

headings = table.getColumnHeadings()

If = table.getColumn(table.getColumnIndex(frap_column)) # Frap values
In = table.getColumn(table.getColumnIndex(nonfrap_column)) # non Frap values
timeValues = table.getColumnAsDoubles(table.getColumnIndex(time_column)) # time values

frame_interval = timeValues[1] - timeValues[0]

min_intensity = min( If )
bleach_frame = If.index( min_intensity )
IJ.log('FRAP frame is ' + str(bleach_frame+1) + ' at t = ' + str(timeValues[bleach_frame]) + ' ' + time_units )

n_slices = len(If)
 
# Calculate normalized curve
if compute_normValues:
	# Compute mean pre-bleach intensity
	mean_If = 0.0
	mean_In = 0.0
	for i in range(bleach_frame):         # will loop until the bleach time
	    mean_If = mean_If + If[i]
	    mean_In = mean_In + In[i]
	mean_If = mean_If / bleach_frame
	mean_In = mean_In / bleach_frame
	normalized_curve = []
	for i in range(n_slices):
	    normalized_curve.append( (If[i] - min_intensity) / (mean_If - min_intensity)   *   (mean_In - min_intensity) / (In[i]- min_intensity) )
else:
	normalized_curve = table.getColumnAsDoubles(table.getColumnIndex(normalize_column))  # Norm values

x = timeValues
y = normalized_curve 
xtofit = x[ bleach_frame : n_slices ]
x0 = xtofit[0]
xtofit = [ xt-x0 for xt in xtofit]
ytofit = normalized_curve[ bleach_frame : n_slices ]

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

plot = Plot("Normalized FRAP curve ", "Time ("+time_units+")", "NU", [], [])
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

rt.incrementCounter()
rt.addValue("half time", thalf )
rt.addValue("mobile fraction (%)", 100 * mobile_fraction )
rt.addValue("R^2", fitter.getFitGoodness() )
rt.addValue("time unit", time_units)
rt.show("Results from analysis")

