from ij import IJ 
from de.mpicbg.scf.spotcoloc import SpotProcessor, SpotVisualization 
from ij import IJ, ImagePlus, ImageStack 
from ij.gui import Overlay, Roi, PointRoi 
from java.awt import Color 
from ij.plugin import Duplicator, Concatenator
from ij.measure import ResultsTable 
from ij.gui import GenericDialog, NonBlockingGenericDialog, WaitForUserDialog, Overlay, Roi, PointRoi
from ij.plugin.frame import RoiManager 

overlay = Overlay()

imp = IJ.getImage(); 
cal = imp.getCalibration();
rm = RoiManager.getInstance()
if not rm:
	rm = RoiManager();
rm.reset();
rm.runCommand("Associate", "true");

doSubpixel=True
doMedian=False
radius = 0.5
quality_threshold = 8000

spotProcessor = SpotProcessor(imp)
spots = spotProcessor.detectSpots(1, radius, quality_threshold, doSubpixel, doMedian)

for peak in spots:
	print peak.getDoublePosition(0), peak.getDoublePosition(1), peak.getDoublePosition(2)
	roi = PointRoi(peak.getDoublePosition(0) / cal.pixelWidth, peak.getDoublePosition(1) / cal.pixelHeight)
	roi.setPosition(int(round(peak.getDoublePosition(2) / cal.pixelDepth))+1)
	roi.setPointType(0)
	overlay.add(roi);

imp1 = imp.duplicate();
imp1.setOverlay(overlay);
imp1.show()

gd = NonBlockingGenericDialog("Action required")
gd.addMessage("Select additional points and click ok")
gd.showDialog()

roinew = imp1.getRoi()

for pp in roinew:
	print pp
	p1 = PointRoi(pp.x, pp.y)
	print p1.getPosition()
	rm.addRoi(p1)
	
rm.runCommand(imp1, "measure");

	