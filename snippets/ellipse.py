from inra.ijpb.binary import BinaryImages
from inra.ijpb.measure import IntensityMeasures
from inra.ijpb.label import LabelImages
from inra.ijpb.measure.region2d import Centroid, InertiaEllipse
from inra.ijpb.measure import IntensityMeasures 

from ij import IJ 


imp = IJ.getImage();
labels = LabelImages.findAllLabels(imp) 

ellipses = InertiaEllipse.inertiaEllipses(imp.getProcessor(), labels, imp.getCalibration())

radius1 = [e.radius1() for e in ellipses]
radius2 = [e.radius2() for e in ellipses]

print radius1, radius2

avgdiameters = [r1 + r2 for r1, r2 in zip(radius1, radius2)] 
print avgdiameters 

keep = []
for l,d in zip(labels, avgdiameters):
	print l, d 
	if d > 8.0:
		keep.append(l)


print keep 

imp_keep = imp.duplicate() 
imp_keep = LabelImages.keepLabels(imp, keep)

imp_keep.show()