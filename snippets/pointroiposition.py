from ij import IJ, ImagePlus, ImageStack
from ij.gui import Roi, PointRoi 
from ij.measure import ResultsTable 
from ij.gui import WaitForUserDialog
from ij.plugin.frame import RoiManager

#imp = IJ.openImage("http://imagej.nih.gov/ij/images/mri-stack.zip"); 
#imp.show()
imp = IJ.getImage();
IJ.run(imp, "Select None", "");

rm = RoiManager.getInstance()
if not rm:
	rm = RoiManager();
rm.reset()

IJ.setTool("multipoint");
WaitForUserDialog("Action required", "Mark some points and click ok").show() 

# get points and add them to roi manager
points = imp.getRoi(); 

for idx, p in enumerate(points.iterator()):
	proi = PointRoi(p.x, p.y)
	print points.getPointPosition(idx)
	proi.setPosition(points.getPointPosition(idx))
	rm.addRoi(proi);


# measure --> slice number is always the last slice marked
#rm.runCommand(imp, "measure");

## here z is always 0
roi_array = rm.getRoisAsArray();
for roi in roi_array:
	print "Roi Z position =", roi.getZPosition()