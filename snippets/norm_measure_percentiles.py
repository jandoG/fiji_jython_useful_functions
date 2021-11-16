#@File[] (label= "Select 1 or multiple images to process", style="files") listfiles
#@Integer (label="DAPI channel no", value= 4) dapi_ch
#@Integer (label="HNF4 channel no", value= 2) hepatocyte_ch
#@Integer (label="Measure intensities from channel no", value= 1) measure_intensity_ch
#@String (label= "Select thresholding method", choices={"Triangle", "Otsu", "Li"}, style="listBox", value= "Triangle") thresholding_method
#@Float (label= "Factor to threshold computation", value=0.7) threshMultiplyer
#@boolean (label="Show segmentation results? (checked = yes, unchecked = no)") show
#@OpService ops


from ij import IJ
from ij.io import FileSaver
from ij import ImagePlus
from ij.process import ImageStatistics
from ij.process import StackStatistics
from ij.measure import Calibration, ResultsTable 
from ij.gui import Roi, Overlay
from ij.plugin.filter import ParticleAnalyzer as PA
from ij.plugin import Duplicator, ImageCalculator, ChannelSplitter
from ij.plugin.frame import RoiManager
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
from inra.ijpb.binary import BinaryImages;
from net.imglib2.img.display.imagej import ImageJFunctions
import math
import os, re, glob
from os.path import isfile, join

ext = '.tif'
fractions = [1, 100]

def main():

	# get all tif files from the list of files (if "Add Folder content" is selected, the files from _results folder is excluded)
	sorted_list = sorted([f.getPath() for f in listfiles if f.getPath().endswith(ext) and not os.path.dirname(f.getPath()).endswith("_results")])

	# create dir to save results 
	results_path = join(os.path.dirname(sorted_list[0]), '_results')
	if not os.path.exists(results_path):
		os.makedirs(results_path)

	# boilerplate
	table = ResultsTable()
	table.reset();
	rm = RoiManager.getInstance()
	if not rm:
	  rm = RoiManager()
	rm.reset()

	for imgfile in sorted_list:
		imp = openImageWithBF(imgfile, virtual= False)
		impName = os.path.splitext(imp.getTitle())[0]
		print "Processing image =", impName
#		imp.show();

		measure_ch = Duplicator().run(imp, measure_intensity_ch, measure_intensity_ch, 1, imp.getNSlices(), 1, imp.getNFrames())
#		measure_ch.show()

		""" 
		Better than idea 2 but negative values need to be taken care of 
		Idea 3: create mask of cytoplasm only in channel 3, 
		Set it on channel 1, get min, max 
		Then get hepatocyte nuclei mask, measure min, max in channel 1. 
		Scale the hepatocyte values to match range of cytoplasm 
		"""
		cyto_m = Duplicator().run(imp, 3, 3, 1, imp.getNSlices(), 1, imp.getNFrames())

		# get cytoplasm roi from channel 3
		threshold = getThreshold(cyto_m)
		cytomask = getCytoMask(cyto_m, threshold)
#		cytomask.show()
		IJ.run(cytomask, "Create Selection", "");
		cyto_roi = cytomask.getRoi()

		# get 2 and 98 percentiles of cytoplasm in ch1 
		cyto_min, cyto_max = percentile(measure_ch, cyto_roi, fractions)
		measure_ch.killRoi()
		
		# get mask of the hepatocytes: objects common to both DAPI and HNF4 channel
		hepa_mask = processImp(imp)
#		hepa_mask.show()

		# get roi of all hepatocytes 
		IJ.run(hepa_mask, "Create Selection", "");
		hepa_roi = hepa_mask.getRoi() 

		# get 2 and 98 percentiles of hepatocytes 
		hepa_min, hepa_max = percentile(measure_ch, hepa_roi, fractions)

		scale = (cyto_max - cyto_min)/ float(hepa_max - hepa_min)
		offset = cyto_min - hepa_min * scale 

		ch1_hepaNorm = normalizeToRef(measure_ch, scale, offset, hepa_roi)
		ch1_hepaNorm.show()
		ch1_hepaNorm.setRoi(hepa_roi)
		measure_ch.setRoi(hepa_roi)
		print "ch1_hepa mean =", measure_ch.getStatistics().mean 
		print "ch1_hepaNorm mean =", ch1_hepaNorm.getStatistics().mean 
		measure_ch.killRoi()
				

		table.incrementCounter();
		table.addValue("Name", impName);
		table.addValue("Nuclei min", hepa_min);
		table.addValue("Nuclei max", hepa_max);
		table.addValue("Cytoplasm min", cyto_min);
		table.addValue("Cytoplasm max", cyto_max);
#		table.addValue("Nuclei mean old", nuclei_oldMean);
#		table.addValue("Nuclei mean new scaled wrt cyto", nuclei_newMean);
#		table.addValue("Nuclei mean new scaled wrt cyto new formula", nuclei_newMean1);

		table.show("Results: All files");

































#==================================== Helper ================================================

def percentile(imp0, roi, fractions):
	# fractions = [5.0, 95.0]

	imp = imp0.duplicate()

	# get stats for stack or single imp
	imp.setRoi(roi);
	
	if imp.getStackSize()>1 :
		stats = StackStatistics(imp) 
	else:
		stats = imp.getProcessor( ).getStatistics()

	# get histogram and its min, max, bins
	hist = stats.getHistogram()
	hMin = stats.histMin
	hMax = stats.histMax
	nBins = len(hist);

	# get sum of all values in histogram
	histSum = [hist[0]]
	count=1
	for val in hist[1:] :
		histSum.append( hist[count] + histSum[count-1] )
		count = count+1

	# get max
	xmax= float(histSum[-1])

	# get percentages instead of values 
	histSum = [x/xmax*100 for x in histSum]

	percentiles = []
	count=0

	# f is min of fractions i.e. 5%
	f = fractions[0]
	for i,f1 in enumerate(histSum) :

		# if % is greater than 5%
		if f1 >= f :
			# get % at i - 1
			f0 = histSum[i-1]
			ifrac = ( f - f0 ) / (f1-f0) + i - 1

			# get actual value and append to empty array
			percentiles.append( hMin + (hMax-hMin)*ifrac/nBins )

			count = count+1
			# restrict count to just len(fractions) so that we have just two values, the min and max wrt percentiles
			if count >= len(fractions) :
				break

			# change the value of f i.e. now it is 95%
			f = fractions[count] 
			
	return percentiles;

def normalizeToRef( imp0, scale, offset, hepatocyte_roi ):
	imp = Duplicator().run(imp0)
#	IJ.run(imp, "32-bit", "");
	print hepatocyte_roi
		
	imp.setRoi(hepatocyte_roi);
	print "Before norm =", imp.getStatistics().mean 
	
	IJ.run(imp, "Clear Outside", "");
	img = ImageJFunctions.wrap(imp)

	# scale the image by multiplication and add the offset
	img = ops.math().multiply(img, scale)
	img = ops.math().add(img, offset )

	imp1 = ImageJFunctions.wrap(img, "norm - " + imp0.getTitle())
	imp1 = Duplicator().run(imp1)
	imp1.show()
	imp1.setCalibration( imp0.getCalibration() )
	imp1.setTitle("norm - " + imp0.getTitle())

#	IJ.run(imp1, "Multiply...", "value=65535.000");

	return imp1


def openImageWithBF(path, virtual= False):
	"""
	opens an image using bio-formats options. 
	params: path to image, virtual = set to true to open an virtual stack 
	"""
	# set options to open image - use virtual for quick loading
	options = ImporterOptions()
	options.setColorMode(ImporterOptions.COLOR_MODE_DEFAULT)
	options.setAutoscale(True)
	options.setStackFormat("Hyperstack")
	options.setVirtual(virtual)
	options.setId(path)
	
	imp = BF.openImagePlus(options)[0]	

	return imp 


def processImp(imp0):
	"""
	Main function to process imp.
	"""
	dapi = Duplicator().run(imp0, dapi_ch, dapi_ch, 1, imp0.getNSlices(), 1, imp0.getNFrames())
	hepatocyte = Duplicator().run(imp0, hepatocyte_ch, hepatocyte_ch, 1, imp0.getNSlices(), 1, imp0.getNFrames())
	
	dapi_mask = getMask(dapi)
	hepa_allmask = getMask(hepatocyte)

	hepatocytes_mask = getHepatocytesOnlyMask(dapi_mask, hepa_allmask)

	return hepatocytes_mask


def getMask(imp0):
	"""
	Returns mask of an image.
	"""
	imp = Duplicator().run(imp0)
	IJ.run(imp, "Subtract Background...", "rolling=50");
	IJ.run(imp, "Gaussian Blur...", "sigma=2");
	IJ.setAutoThreshold(imp, thresholding_method+" dark");
	IJ.run(imp, "Convert to Mask", "")
	IJ.run(imp, "Fill Holes", "")

	return imp


def getHepatocytesOnlyMask(dapi_mask, hepa_allmask):
	"""
	Returns intersection of two masks.
	hepa_allmask contains hepatocytes, dapi_mask shows all cell nuclei
	Hence, multiplying the two masks results in a mask having only hepatocytes.
	"""
#	dapi_mask.show()
#	hepa_allmask.show()
	ic = ImageCalculator()
	hepatocytes_only = ic.run("AND create", dapi_mask, hepa_allmask)
	IJ.run(hepatocytes_only, "Fill Holes", "");
	IJ.run(hepatocytes_only, "Watershed", "");
	hepatocytes_only = BinaryImages.sizeOpening(hepatocytes_only, 160)

	return hepatocytes_only

def getThreshold(imp):
	stats = StackStatistics(imp)
	hist = stats.getHistogram()
	hMin = stats.histMin
	hMax = stats.histMax
	nBins = len(hist);
	#print nBins
	#print hist
	#print hMin
	#print hMax

	# find the maximum of the histogram and peak variance
	maxFreq = 0;
	binIdxMaxFreq=0
	binIdx=0
	for freq in hist:
		if freq>maxFreq:
			maxFreq = freq
			binIdxMaxFreq = binIdx
		binIdx = binIdx+1;
	IMaxFreq = hMin + binIdxMaxFreq * (hMax-hMin) / nBins;

	# measure the peak variance. use only the low intensity half of the histogram peak
	sumdI2 = 0
	sumFreq = 0
	#Ilist = [hMin + i*(hMin-hMax)/nBins for i in range(binIdxMaxFreq+1)]
	for binIdx in range(binIdxMaxFreq): # I,freq in zip(Ilist,hist) :
		I = hMin + binIdx*(hMax-hMin)/nBins;
		freq = hist[binIdx]
		sumdI2 = sumdI2 + freq*math.pow(I-IMaxFreq,2)
		sumFreq = sumFreq + freq
	
	if sumFreq!= 0:
		ISigma = math.sqrt(sumdI2/sumFreq)

	else:
		ISigma = 1

	# define a threshold
	threshold = IMaxFreq - threshMultiplyer * ISigma  
	# original script: IMaxFreq + threshMultiplyer * ISigma, here we subtract, 
	# since we are interested to get mask of cytoplasm 
	
#	print threshold, IMaxFreq, ISigma

	return threshold

def getCytoMask(imp, threshold):

	# droplet mask creation
	
	impMask = Duplicator().run(imp)
	IJ.run(impMask, "Gaussian Blur...", "sigma=3");
	IJ.setRawThreshold(impMask, threshold, 65535, None)
	IJ.run(impMask, "Convert to Mask", " ")

	return impMask


main()
print("Done")

roiManager = RoiManager.getRoiManager();
roiManager.close();
