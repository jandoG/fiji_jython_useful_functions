from ij import IJ, ImagePlus
from inra.ijpb.label import LabelImages as li 

pred = IJ.getImage();

label1 = li.cropLabel(pred, 1, 5)

label11 = ImagePlus("label1", ipp)
label11.show()
label1.show()
