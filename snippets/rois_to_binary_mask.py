from ij import IJ, Prefs, ImagePlus, WindowManager
from ij.plugin.frame import RoiManager
import os 
import glob
from ij.process import ImageProcessor

#Prefs.setForegroundColor(255,255,255)
#Prefs.setBackgroundColor(0,0,0)

rm=RoiManager.getRoiManager()
nrois=rm.getCount()

imp=IJ.getImage();
width, height, nChannels, nSlices, nFrames=imp.getDimensions();

mask=IJ.createImage("mask", "8-bit black", width, height, 1)
mask.show()
# add calibration after roi processing! (reason: use make band in px)
for idx in range(nrois):
	roi=rm.getRoi(idx)
	mask.setRoi(roi)
	IJ.run("Make Band...", "band=1.0");	
	bandroi=mask.getRoi()
	
	ip=mask.getProcessor()
	ip.setValue(255) # fill roi with white
	ip.fill(roi) 
	ip.setValue(0) # draw a black outline (avoid merging regions)
	ip.fill(bandroi)

# TODO: copy the calibration from 2D
calib=imp.getCalibration()
mask.setCalibration(imp.getCalibration())
print calib




	
	