#@ ImagePlus imp0 

from ij import IJ
from ij import ImagePlus
from ij.gui import Roi, Overlay, PolygonRoi, PointRoi, WaitForUserDialog
from inra.ijpb.label import LabelImages
from inra.ijpb.binary import BinaryImages
from inra.ijpb.binary.conncomp import FloodFillComponentsLabeling


imp = imp0.duplicate()
IJ.run(imp, "Gaussian Blur...", "sigma=3");
IJ.setAutoThreshold(imp, "Huang" + " dark");
IJ.run(imp, "Convert to Mask", "");
IJ.run(imp, "Watershed", ""); 

label0 = BinaryImages.componentsLabeling(imp, 4, 16)
label = LabelImages.sizeOpening(label0, 200)
#label.show() 
#
#IJ.setTool("multi-point");
#WaitForUserDialog("Action required", "Please click on nuclei to discard and press ok.").show()
#
#roi = label.getRoi() 
#print roi
#
#LabelImages.removeLabels(label, roi, True)
#label.killRoi()
#IJ.run(label, "glasbey_on_dark", ""); 

labels = LabelImages.findAllLabels(label)

while len(labels) > 1:
	label.show()
	IJ.setTool("multi-point");
	WaitForUserDialog("Action required", "Keep only one nuclei to process. Please click on nuclei to DISCARD and press ok.").show()

	roi = label.getRoi()
	LabelImages.removeLabels(label, roi, True)
	IJ.run(label, "glasbey_on_dark", ""); 

	labels = LabelImages.findAllLabels(label)
	
label.killRoi()
LabelImages.remapLabels(label)

label1 = label.duplicate()
LabelImages.replaceLabels(label1, [1], 255)
IJ.run(label1, "8-bit", "");
IJ.run(label1, "Grays", "");

label1.show()

IJ.run(label1, "Create Selection", "");
nuclei_roi = label1.getRoi()
print nuclei_roi



