from ij import IJ 
from ij.gui import GenericDialog, NonBlockingGenericDialog, WaitForUserDialog, Overlay

# open image and display
#imp = IJ.openImage("http://wsr.imagej.net/images/Spindly-GFP.zip");
imp = IJ.getImage();
imp.show()
imp.killRoi(); 
imp.setT(1);

WaitForUserDialog("Action required", "Draw a rectangular ROI and press ok").show() 
roi = imp.getRoi();

ov = Overlay()

for t in xrange(1, imp.getNFrames() + 1, 10):
	imp.setT(t);
        
    # at t=1, set roi and add roi to overlay
	if t ==1:
		imp.setRoi(roi)

        # which function to use here?
        # goal: set roi on all Z, C but only on current T

#		roi.setPosition(0, 0, t);                          # overlay seen on last t stopped, all z, c
#		roi.setPosition(imp.getC(), imp.getSlice(), t);    # no overlay seen
#		roi.setPosition(imp.getSlice())                    # overlay seen on z=1, t=1 and all c
#		ov.add(roi)										   # last roi drawn on all c, z, t
 
		# this works
		roi.setPosition(0, 0, t);
		ov.add(roi.clone())                                        
      
    # at t=11, 21, .. ask user to move roi 
    # get the new roi and add to overlay
	else:
		gd = NonBlockingGenericDialog("Action required at T=" + str(t))
		gd.addMessage("Move the ROI to your area of interest by clicking inside the rectangle")
		gd.addMessage("Addtionally, check the box if you want to stop here")
		gd.addMessage("Click ok")
		gd.addCheckbox("Stop here?", False)
		gd.showDialog()
		
		roi_new = imp.getRoi()
		imp.setRoi(roi_new)
       
        # same goal as above, set roi on all Z, C but only on current T
        
#		roi_new.setPosition(0, 0, t);                       # overlay seen on last t stopped, all z, c
#		roi.setPosition(imp.getC(), imp.getSlice(), t);     # no overlay seen
#		roi.setPosition(imp.getSlice())						# overlay seen on z=1, t=1 and all c
#		ov.add(roi_new)                                     # last roi drawn on all c, z, t

		# this works
		roi_new.setPosition(0, 0, t);  
		ov.add(roi_new.clone())
		
		if gd.wasCanceled():
			pass

		# if box is checked, exit the loop
		if gd.getNextBoolean() == True:
			print("Exiting loop")
			break


print("Out of loop")

imp_final = imp.duplicate();
imp_final.setOverlay(ov)
imp_final.show()