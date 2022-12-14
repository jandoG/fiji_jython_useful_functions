from ij import ImagePlus, IJ
from ij.gui import Roi, PolygonRoi, Overlay

imp = IJ.getImage()
#roi = imp.getRoi()
#print roi 
#
#roicvx = PolygonRoi(roi.getFloatConvexHull(), Roi.TRACED_ROI)
#print roicvx
#
#width_um = 5
#
#IJ.run("Make Band...", "band="+str(width_um));

stats = imp.getStatistics();
count = int(stats.max)

rois = []
overlay = Overlay()

for t in range(1, count + 1):
	mask = imp.duplicate()
	IJ.setAutoThreshold(imp, "Default");
	IJ.setThreshold(mask, t, t);
	IJ.run(mask, "Convert to Mask", "");
	IJ.run(mask, "Create Selection", "")
	roi = mask.getRoi()
	rois.append(roi)

	overlay.add(roi);
	
imp.setOverlay(overlay);

	

	 