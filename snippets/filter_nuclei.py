
from ij import IJ, Prefs, ImagePlus, WindowManager
from ij.process import ImageProcessor, StackStatistics
from ij.plugin import ZProjector, Duplicator, Scaler, ImageCalculator
from ij.plugin.filter import ThresholdToSelection;
from ij.gui import Overlay, WaitForUserDialog
from ij.plugin.frame import RoiManager
from ij.measure import ResultsTable, Measurements
from de.csbdresden.stardist import StarDist2D 
from net.imglib2.img.display.imagej import ImageJFunctions
from inra.ijpb.morphology import LabelImages
from inra.ijpb.binary import BinaryImages
from inra.ijpb.plugins import AnalyzeRegions3D
from java.awt import Color
import os 
import glob
from java.io import File 

frame=5

IJ.selectWindow("max");
maxImp=IJ.getImage();
IJ.selectWindow("stardist");
nucleiLabels=IJ.getImage()
print nucleiLabels

def removeDarkNuclei(maxImp2Dt,nucleiLabels, fraction=0.5):
	"""Removes all segmented nuclei which are too dark.
	Selection condition: Nuclei are kept if >50% (fraction=0.5) of their interior is detected with a 
			default thresholding operation,
	maxImp2Dt: 2D+t
	nucleiLabels: 2D 
	returns: updated nucleiLabels
	"""

	maskT=maxImp2Dt.duplicate()
	IJ.run(maskT, "Gaussian Blur...", "sigma=1 stack");
	IJ.setAutoThreshold(maskT, "Default dark stack");
	IJ.run(maskT, "Convert to Mask", "method=Default background=Dark black");

	mask=Duplicator().run(maskT,1,1,1,1, frame,frame)

	# keep all bright nuclei: must be at least 50% detected in a thresholded mask
	nucleiLabels.show()
	mask.show()
	
	# measure foreground fraction for each nuclei
	rt=ResultsTable.getResultsTable()
	if rt is not None:
		rt.reset()
	IJ.run(nucleiLabels, "3D Intensity Measure", "objects=["+nucleiLabels.getTitle()+"] signal=["+mask.getTitle()+"]");
	rt=ResultsTable.getResultsTable()
	intensities=rt.getColumn(rt.getColumnIndex("Average"))
	foregroundFracs=[i/float(255) for i in intensities]# foreground fraction

	# remove dark nuclei 
	for labelID in range(1,len(foregroundFracs)+1):
		if foregroundFracs[labelID-1]<fraction:
			print nucleiLabels
			LabelImages.replaceLabels(nucleiLabels,[labelID],0)
	
		
	# clean up
	IJ.selectWindow("Results")
	IJ.run("Close")
	mask.hide()
	# nucleiLabels.hide()

	LabelImages.remapLabels(nucleiLabels)

	return nucleiLabels 


nucleiLabels=removeDarkNuclei(maxImp,nucleiLabels, fraction=0.5)
nucleiLabels.show()
