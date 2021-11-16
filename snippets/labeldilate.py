from ij import IJ
from ij import ImagePlus
from ij import ImageStack
from ij.plugin import Duplicator, ZProjector
from inra.ijpb.binary import BinaryImages
from inra.ijpb.plugins import DistanceTransformWatershed3D
from inra.ijpb.label import LabelImages 


imp_label = IJ.getImage(); 

labeldilate = LabelImages.dilateLabels(imp_label, 5);

labeldilate.show() 

print("Done")

