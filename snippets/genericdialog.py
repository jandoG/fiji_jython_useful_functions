from ij import IJ 
from ij.gui import GenericDialog, NonBlockingGenericDialog, WaitForUserDialog, Overlay
from ij.measure import ResultsTable

imp = IJ.openImage("http://wsr.imagej.net/images/Spindly-GFP.zip");
#imp = IJ.getImage();

imp.show()
imp.killRoi(); 
imp.setT(1);

WaitForUserDialog("Action required", "Draw a rectangular ROI around organoid of interest and press ok").show() 
roi = imp.getRoi();

rt = ResultsTable()
ov = Overlay()

print imp.getStackSize();

for t in xrange(1, imp.getNFrames(), 10):
	imp.setT(t);
	rt.incrementCounter()

	if t ==1:
		imp.setRoi(roi);
		stats = imp.getStatistics();
		rt.addValue("T", t)
		rt.addValue("Mean", stats.mean)
#		imp.setPosition(imp.getC(), imp.getSlice(), t);
#		imp.setRoi(roi);
#		roi.setPosition(0, 0, t);
		ov.add(roi)

	else:
		gd = NonBlockingGenericDialog("Action required at T=" + str(t))
		gd.addMessage("Move the ROI to your organoid being processed by clicking inside the rectangle")
		gd.addMessage("Addtionally, check the box if you want to stop here")
		gd.addMessage("Click ok")
		gd.addCheckbox("Stop here?", False)
		gd.showDialog()
		
		roi_new = imp.getRoi()
		imp.setRoi(roi_new)
		stats = imp.getStatistics();
		rt.addValue("T", t)
		rt.addValue("Mean", stats.mean)

#		roi_new.setPosition(0, 0, t);
		ov.add(roi_new)
	
		if gd.wasCanceled():
			pass
	
		if gd.getNextBoolean() == True:
			print("Time to end")
			break


print("Out of loop")
rt.show("Final")

imp_final = imp.duplicate();
imp_final.setOverlay(ov)
imp_final.show()