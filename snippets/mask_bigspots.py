#@ImagePlus imp0

from ij import IJ, WindowManager
from ij import ImagePlus
from ij import ImageStack
from ij.plugin import Duplicator, ChannelSplitter, ImageCalculator

def maskBigSpots(imp_t):
	"""
	Masks the big prominent spots in the and 
	fills them with background value.

	Input: single time frame image 
	Output: image with bigger spots masked
	"""

	# get the mak of the big spots on the image 
	imp = imp_t.duplicate()
	IJ.run(imp, "Subtract Background...", "rolling=30");
	IJ.run(imp, "Grays", "");
	IJ.setAutoThreshold(imp, "Huang dark");
	IJ.run(imp, "Convert to Mask", "");

	# run analyze particles to filter out the smaller spots 
	IJ.run(imp, "Analyze Particles...", "size=2-10000 show=Masks exclude");
	imp_ = IJ.getImage();
	imp_mask = imp_.duplicate()
	imp_.close()
	IJ.run(imp_mask, "Invert LUT", "");
	IJ.run(imp_mask, "Dilate", "");
	IJ.run(imp_mask, "Dilate", "");
	IJ.run(imp_mask, "Dilate", "");

	# creates a roi of the big spots left on the mask
	IJ.run(imp_mask, "Create Selection", "");
	roi = imp_mask.getRoi()

	# creates a roi of the background and gets its minimum value
	imp_t.setRoi(roi)
	IJ.run(imp_t, "Make Inverse", "");
	stats = imp_t.getStatistics()
	minI = stats.min 

	# fill the big spots with the min value + buffer value 
	ip = imp_t.getProcessor()
	ip.setRoi(roi);
	ip.setValue(minI + 150.0);
	ip.fill(roi);

	# remove the roi and return the image 
	IJ.run(imp_t, "Select None", "");
	
	return imp_t





#============== M A I N ====================

imp1 = imp0.duplicate()
imps = ChannelSplitter().split(imp1)
impch1 = imps[0]
impch2 = imps[1]

imp_process = impch1.duplicate()

new_ch1 = []

for t in range(1, imp_process.getNFrames() + 1):
	imp_process.setT(t)

	imp = Duplicator().run(imp_process, 1, 1, 1, 1, t, t)
#	imp.show()

	imp_smallSpots = maskBigSpots(imp)
#	imp_smallSpots.show()

	imp_smallSpots.setTitle("t_" + str(t))
	new_ch1.append(imp_smallSpots)

#	if t==3:
#		break

stack = ImageStack()
for tt in new_ch1:
	print tt.getTitle()

	stack.addSlice(tt.getProcessor())

imp_ch1_new = ImagePlus("ch1 new", stack)
imp_ch1_new.show()
	
	




	
	
	
	