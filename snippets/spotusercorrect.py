from de.mpicbg.scf.spotcoloc import SpotProcessor, SpotVisualization 
from ij import IJ, ImagePlus, ImageStack 
from ij.gui import Overlay, Roi, PointRoi 
from java.awt import Color 
from ij.plugin import Duplicator, Concatenator
from ij.measure import ResultsTable 
from ij.gui import GenericDialog, NonBlockingGenericDialog, WaitForUserDialog, Overlay, Roi, PointRoi
from ij.plugin.frame import RoiManager


#def discard(roi_pts, roi_array):
#
#	new_pt = []
#	new_rois = []
#
#	for roi in roi_array:
##		print roi 
#		for p in roi_pts:
##			print p
#			if not roi.contains(p.x, p.y):
#				print "does not contain"
#				print p, roi 
#				new_pt.append(p)
#				new_rois.append(roi)
#
#	print (new_pt), (new_rois)
#	print set(new_pt), set(new_rois)
#
#	return set(new_pt), set(new_rois) 


imp = IJ.getImage(); 
IJ.run(imp, "Select None", "");

IJ.setTool("rectangle");
overlay = Overlay()
rm = RoiManager.getInstance()
if not rm:
	rm = RoiManager();
rm.reset()

doSubpixel=True
doMedian=False
radius = 0.5
quality_threshold = 8000

for t in xrange(1, imp.getNFrames() + 1, 2):
	impt = Duplicator().run(imp, 1, 1, 1, imp.getNSlices(), t, t)

	if t == 1:
		spotProcessor = SpotProcessor(impt)
		spots = spotProcessor.detectSpots(1, radius, quality_threshold, doSubpixel, doMedian)
		ov = SpotVisualization.createOverlayOfSpots(impt, spots, Color.cyan)

		pt_roi = ov.toArray();
		for r in pt_roi: rm.addRoi(r)
		
		imptt = impt.duplicate() 
		imptt.show()
		rm.runCommand(imptt, "Show all");

		WaitForUserDialog("Please delete false-positive spots from roi manager and click ok").show()
		positives = rm.getRoisAsArray();
		for rp in positives:
			rp.setPosition(0, 0, t)
			overlay.add(rp.clone());
		
		IJ.setTool("multipoint");
		gd = NonBlockingGenericDialog("Action required at T=" + str(t))
		gd.addMessage("Select additional points and click ok")
		gd.showDialog()

		points = imptt.getRoi()
		for idx, p in enumerate(points.iterator()):
			proi = PointRoi(p.x, p.y)
			print points.getPointPosition(idx)
			proi.setPosition(points.getPointPosition(idx))
#			proi.setPosition(0, points.getPointPosition(idx), t)
			overlay.add(proi.clone())
			rm.addRoi(proi);

		print "No spots at ", str(t), "=", len(rm.getRoisAsArray())

		imptt.close()
		rm.reset(); 

	else:
		spotProcessor = SpotProcessor(impt)
		spots = spotProcessor.detectSpots(1, radius, quality_threshold, doSubpixel, doMedian)
		ov = SpotVisualization.createOverlayOfSpots(impt, spots, Color.cyan)

		pt_roi = ov.toArray();
		for r in pt_roi: rm.addRoi(r)
		
		imptt = impt.duplicate() 
		imptt.show()
		rm.runCommand(imptt, "Show all");

		WaitForUserDialog("Please delete false-positive spots from roi manager and click ok").show()
		positives = rm.getRoisAsArray();
		for rp in positives:
			rp.setPosition(0, 0, t)
			overlay.add(rp.clone());
		
		IJ.setTool("multipoint");
		gd = NonBlockingGenericDialog("Action required at T=" + str(t))
		gd.addMessage("Select additional points and click ok")
		gd.showDialog()

		points = imptt.getRoi()
		for idx, p in enumerate(points.iterator()):
			proi = PointRoi(p.x, p.y)
			print points.getPointPosition(idx)
			proi.setPosition(points.getPointPosition(idx))
#			proi.setPosition(0, points.getPointPosition(idx), t)
			overlay.add(proi.clone())
			rm.addRoi(proi);

		print "No spots at ", str(t), "=", len(rm.getRoisAsArray())

		imptt.close()
		rm.reset(); 

#	del roi_array
#	break

imp1 = imp.duplicate();
imp1.setOverlay(overlay)
imp1.show()

print "Done"