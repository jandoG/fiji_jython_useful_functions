from de.mpicbg.scf.spotcoloc import SpotProcessor, SpotVisualization 
from ij import IJ, ImagePlus, ImageStack 
from ij.gui import Overlay, Roi, PointRoi 
from java.awt import Color 
from ij.plugin import Duplicator, Concatenator
from ij.measure import ResultsTable 

imp = IJ.getImage(); 

doSubpixel=True
doMedian=False
radius = 0.5
quality_threshold = 8000

overlay = Overlay()

for t in range(1, imp.getNFrames() + 1):
	impt = Duplicator().run(imp, 1, 1, 1, imp.getNSlices(), t, t)

	spotProcessor = SpotProcessor(impt)
	spots = spotProcessor.detectSpots(1, radius, quality_threshold, doSubpixel, doMedian)
	IJ.log("No of spots at t" + str(t) + " =" + str(len(spots)))
	
	ov = SpotVisualization.createOverlayOfSpots(impt, spots, Color.cyan)
	roi_array = ov.toArray()
	
	for roi in roi_array:
		centroids = roi.getContourCentroid()
		pr = PointRoi(centroids[0], centroids[1])
		pr.setPosition(0, 0, t);
		pr.setPointType(2)
		overlay.add(pr.clone())
#		roi.setPosition(0, 0, t);
#		overlay.add(roi.clone())

	del roi_array

imp1 = imp.duplicate();
imp1.setOverlay(overlay)
imp1.show()




