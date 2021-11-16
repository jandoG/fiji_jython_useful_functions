#@ImagePlus imp

from ij import IJ
from ij.plugin.frame import RoiManager
from ij.gui import Overlay, Roi, TextRoi 
from ij.measure import ResultsTable, Measurements
from ij.process import StackStatistics, ImageStatistics

rm = RoiManager.getInstance()
if not rm:
	rm = RoiManager() 


IJ.run("Set Measurements...", "area mean centroid stack display redirect=None decimal=3"); 

ip = imp.getProcessor()
 
# Instead of multimeasure command in RoiManager, get an array of ROIs.
ra = rm.getRoisAsArray()
 
# loop through ROI array and do measurements. 
# here is only listing mean intensity of ROIs
# if you want more, see 
# http://rsbweb.nih.gov/ij/developer/api/ij/process/ImageStatistics.html

for t in range(1, imp.getNFrames() + 1):
	imp.setT(t)
	print t
	for r in ra:
		idx = rm.getRoiIndex(r)
		print rm.getName(idx)
		imp.setRoi(r)
		istats = imp.getStatistics(Measurements.ALL_STATS)
	#	print istats
		print istats.xCentroid
#		print istats.yCentroid

