from ij import IJ, ImagePlus 
from ij.measure import Calibration
from ij.gui import GenericDialog, NonBlockingGenericDialog, WaitForUserDialog, Overlay, Roi, PointRoi
from ij.plugin import Duplicator, ZProjector
from ij.measure import ResultsTable 

imp0 = IJ.getImage();

IJ.setTool("rectangle");
ov = Overlay() 
rt = ResultsTable()
global roi 

for t in xrange(1, imp0.getNFrames() + 1, 10):

	rt.incrementCounter()
	
	imp_z = Duplicator().run(imp0, 1, 1, 1, imp0.getNSlices(), t, t)

	imp = ZProjector.run(imp_z,"max");
	imp.killRoi(); 
	imp.show();

	# if t= 1, ask user to draw roi 
	if t == 1:
		WaitForUserDialog("Action required", "Draw a rectangular ROI around organoid of interest and press ok").show() 

		roi = imp.getRoi() 

		imp_process = imp.duplicate();
		imp_process.setRoi(roi);
		
		IJ.run(imp_process, "Gaussian Blur...", "sigma=3");
		IJ.setAutoThreshold(imp_process, "Li dark");
		IJ.run(imp_process, "Convert to Mask", "");
		IJ.run(imp_process, "Create Selection", "");
		roi_organoid = imp_process.getRoi() 

		imp.setRoi(roi_organoid);
		IJ.setTool("brush"); 

		WaitForUserDialog("Action required", "Adjust the selection and press ok").show() 
		roi_organoid1 = imp.getRoi();
		roi_organoid1.setPosition(t)
		ov.add(roi_organoid1.clone())

		stats = imp.getStatistics(); 
		print "Area= ", stats.area
		rt.addValue("Frame", t)
		rt.addValue("Area", stats.area)
		
	else: 
		imp.setRoi(roi);
		IJ.setTool("rectangle");
		gd = NonBlockingGenericDialog("Action required at T=" + str(t))
		gd.addMessage("Move the ROI to your organoid being processed by clicking inside the rectangle")
		gd.addMessage("Addtionally, check the box if you want to stop here")
		gd.addMessage("Click ok")
		gd.addCheckbox("Stop here?", False)
		gd.showDialog()

		roi_new = imp.getRoi();

		imp_process = imp.duplicate();
		imp_process.setRoi(roi_new)
		IJ.run(imp_process, "Gaussian Blur...", "sigma=3");
		IJ.setAutoThreshold(imp_process, "Li dark");
		IJ.run(imp_process, "Convert to Mask", "");
		IJ.run(imp_process, "Create Selection", "");
		roi_organoid = imp_process.getRoi() 

		imp.setRoi(roi_organoid);
		
		IJ.setTool("brush"); 
		WaitForUserDialog("Action required", "Adjust the selection and press ok").show() 
		roi_organoid1 = imp.getRoi();
		roi_organoid1.setPosition(t)
		ov.add(roi_organoid1.clone())

		stats = imp.getStatistics(); 
		print "Area= ", stats.area

		rt.addValue("Frame", t)
		rt.addValue("Area", stats.area)
		
		# if cancel was pressed, exit loop
		if gd.wasCanceled():
			imp.close()
			break

		# if box is checked, exit the loop
		if gd.getNextBoolean() == True:
			print "Exiting analysis, last time point processed= ", t 
			imp.close()
			break

	imp.close()

imp_result = ZProjector.run(imp0.duplicate(), "max all")
imp_result.setOverlay(ov)
imp_result.show()

rt.show("Results")
		
		
		

		

		

