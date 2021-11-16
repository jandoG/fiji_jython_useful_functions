from ij import IJ, ImagePlus, ImageStack
from ij.gui import Roi, PointRoi 
from ij.measure import ResultsTable 
from ij.gui import WaitForUserDialog
from ij.plugin.frame import RoiManager

IJ.createImage("Untitled","8-bit black",200,200,1).show()
rm = RoiManager(False)
roi = PointRoi(10, 10)
roi.setPosition(2, 2, 1)
rm.addRoi(roi)
roi2 = rm.getRoi(0)
print(roi2.getCPosition(), roi2.getZPosition(), roi2.getTPosition());