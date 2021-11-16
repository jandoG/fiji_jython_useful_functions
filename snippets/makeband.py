from ij import IJ
from ij import ImagePlus 
from ij.plugin import Duplicator, ImageCalculator, ChannelSplitter
from ij.gui import Roi, Overlay
from ij.plugin.frame import RoiManager
from inra.ijpb.binary import BinaryImages;
import math
from java.awt import Color

def getMask(imp0):
	"""
	Creates mask of the image for ch1 and ch3.
	Input: imp
	Output: mask of imp
	"""
	imp = Duplicator().run(imp0)
	IJ.run(imp, "Subtract Background...", "rolling=50 disable");
	IJ.run(imp, "Gaussian Blur...", "sigma=2");
	IJ.setAutoThreshold(imp, "Li dark no-reset");
	IJ.run(imp, "Convert to Mask", "")
	imp_mask_opened = BinaryImages.sizeOpening(imp, 300)
	IJ.run(imp_mask_opened, "Fill Holes", "");

	return imp_mask_opened


def getRoi(mask):
	"""
	Creates a selection in the mask and returns it as Roi.
	"""
	IJ.run(mask, "Create Selection", "");
	IJ.run(mask, "Make Inverse", "");
	roi = mask.getRoi()

	return roi


imp0 = IJ.getImage();
imp0.killRoi()

imp = Duplicator().run(imp0)
imps = ChannelSplitter().split(imp)

c1 = imps[0]
c2 = imps[1]
c3 = imps[2]

c1mask = getMask(c1)
c3mask = getMask(c3)

#c1mask.show()
#c3mask.show()

ic = ImageCalculator()
hepatocytes_only = ic.run("AND create", c1mask, c3mask)
IJ.run(hepatocytes_only, "Fill Holes", "");
IJ.run(hepatocytes_only, "Watershed", "");
hepatocytes_only = BinaryImages.sizeOpening(hepatocytes_only, 300)
#hepatocytes_only.show()
heparoi = getRoi(hepatocytes_only)
imp0.setRoi(heparoi)

rm = RoiManager.getInstance()
if not rm:
	rm = RoiManager()

rm.reset()
rm.addRoi(heparoi)
rm.runCommand(imp0,"Split");
rm.select(0);
rm.runCommand(imp,"Delete");

roi_array = rm.getRoisAsArray()
bands = []
rm.reset()

for r in roi_array:
	imp0.setRoi(r)
	IJ.run("Make Band...", "band=1.2");
	bandroi = imp0.getRoi()
	bands.append(bandroi)
#	rm.addRoi(bandroi)

imp0.killRoi()

c2.show()

dark = []
inter = []
bright = []
ov = Overlay()


for roi, bandroi in zip(roi_array, bands):
	c2.setRoi(roi)
#	IJ.run(c2,"Measure", "");
	statsroi = c2.getStatistics()
	mean_roi = statsroi.mean


	c2.setRoi(bandroi)
#	IJ.run(c2,"Measure", "");
	statsband = c2.getStatistics()
	mean_band = statsband.mean

	if mean_roi < mean_band:
		dark.append(roi)
		print "dark", abs(mean_roi - mean_band)
#		print "mean band", mean_band
		roi.setStrokeColor(Color.white);
		ov.add(roi)
		
	elif abs(mean_roi - mean_band) <= 800:
		inter.append(roi)
		print "inter", abs(mean_roi - mean_band)
		roi.setStrokeColor(Color.cyan);
		ov.add(roi)
		
	else:
		print "bright", mean_roi
		bright.append(roi)
		roi.setStrokeColor(Color.yellow);
		ov.add(roi)

c2.killRoi()
c2.setOverlay(ov)
	




















