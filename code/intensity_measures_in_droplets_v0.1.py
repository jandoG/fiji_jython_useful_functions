#################README#####################
# Author Gayathri Nadar, SCF MPI-CBG, nadar@mpi-cbg.de
# V0.0: segmentation of GUV membrane in 2- channel images, ch1 = membrane, ch2 = GUV. Segment ch1, measure intensity and background in ch2
# V0.1: option to set min diameter of the GUV to be detected to filter out small false positives.

"""
Description:
	- This script segments the GUV membrane in one channel and uses this roi to measure intensity from GUV in another channel
	- Functions of the script:
		1. Get GUV and membrane channel numbers from user.
		2. Segment the membrane channel image to get GUV rois.
		3. Filter out rois which are below the circularity specified by user.
		4. For remaining rois, measure intensity from GUV channel and output table.

Usage:
	- Open script with FIJI, click run.
	- Set the required parameters.
		1. processFile: file to be processed
		2. min_circ: minimum circularity of the GUV, used to filter objects which are not GUV. 
		All objects with circ below this will be rejected.
		3. min_diam: minimum diameter of the GUV in microns, used to filter objects which are not GUV. 
		All objects with area below diam^2 below will be rejected.
		4. membrane_ch: channel no for membrane
		5. guv_ch: channel no for GUV
		6. saveResults: option to save results.
	- Click ok.

Output:
	- Results table
	- Image with overlay of detected GUV

Remark:
	- Please install plugin MorpholibJ required to run this script. (https://imagej.net/MorphoLibJ). Do this in FIJI by going to Help > Update... 
	- Then select "Manage Update Sites" and select "IJPB plugins" from the list. Click close to update and then restart FIJI.
"""

############################################

#@File(label= "File to process") processFile
#@Float (label="Min circularity required to be identified as GUV (0.0 - 1.0)", value = 0.8) min_circ
#@Float (label="Min diameter required to be identified as GUV in image units", value= 3.0) min_diam
#@Integer (label= "Channel which shows membrane") membrane_ch
#@Integer (label= "Channel which shows GUVs") guv_ch
#@boolean (label="Save image and table results (checked = yes, unchecked = display)") saveResults


from ij import IJ
from ij import ImagePlus
from ij import ImageStack
from ij import WindowManager
from ij.text import TextWindow
from ij.io import FileSaver
from ij.gui import PolygonRoi, Roi, ShapeRoi, Overlay
from ij.plugin import Duplicator, RGBStackMerge
from ij.plugin import HyperStackConverter
from ij.plugin import ImageCalculator
from ij.plugin.filter import ParticleAnalyzer
from ij.plugin.filter import MaximumFinder
from ij.plugin.frame import RoiManager
from ij.process import ImageStatistics
from ij.measure import ResultsTable, Measurements, Calibration

from inra.ijpb.morphology import MinimaAndMaxima;
from inra.ijpb.binary import BinaryImages;
from inra.ijpb.morphology import Morphology;
from inra.ijpb.watershed import Watershed;
from inra.ijpb.plugins import MorphologicalSegmentation
from inra.ijpb.plugins.MorphologicalSegmentation import ResultMode
from inra.ijpb.label import LabelImages

from de.mpicbg.scf.fijiplugins.ui.roi import LabelMapToRoiManagerPlugin;

from loci.plugins import BF
from loci.formats import ImageReader
from loci.formats import MetadataTools
from loci.plugins.prefs import OptionsList
from loci.plugins.in import ImporterOptions

from java.awt import Color;
import math
import os
from os.path import isfile, join

global options

def main():
	options = ImporterOptions()
	options.setColorMode(ImporterOptions.COLOR_MODE_DEFAULT)
	options.setAutoscale(True)
	options.setStackFormat("Hyperstack")
	options.setOpenAllSeries(True)

	imgFile= processFile.getPath()
	filename = os.path.splitext(os.path.basename(imgFile))[0]

	if saveResults:
		results_path = join(os.path.dirname(imgFile), 'results', filename)
		subfolders = ("images_overlay", "csv")
		if not os.path.exists(results_path):
			os.makedirs(results_path)
		for sb in subfolders:
			if not os.path.exists(join(results_path, sb)):
				os.makedirs(join(results_path, sb))

	# get no of positions/ series
	reader = ImageReader()
	reader.setId(imgFile)
	seriesCount = reader.getSeriesCount()
	reader.close()

	if (seriesCount==1):
		imp= BF.openImagePlus(imgFile)[0]
		impName = os.path.splitext(imp.getTitle())[0]

		print "Processing file: ", impName
		
		imp_final, rt_final = processImp(imp)
		
		if saveResults:
			rtName = impName + '_guv_measures.csv'
			rt_final.saveAs(join(results_path, subfolders[1], rtName))
			imp_final_name = impName + '_results_guv_overlay.tif'
			FileSaver(imp_final).saveAsTiff(join(results_path, subfolders[0], imp_final_name))

		else:
			imp_final.setTitle(impName + "_final_result")
			imp_final.show()
			rt_final.show(impName + ": guv measures")


	else:
		options.setId(imgFile)
		imp_series= BF.openImagePlus(options)

		all_imps = []
		for im in imp_series:
			all_imps.append(im)
			
#		test = 0
#		if test:
#			for i in range(3):
#				imp_final, rt_final = processImp(all_imps[i])
#				imp_final.show()
#				rt_final.show("guv measures for series no: " + str(i + 1))

		else:
			for i in range(len(all_imps)):
				impName = os.path.splitext(all_imps[i].getTitle())[0]
				print "Processing position: ", str(i+1), "/", str(seriesCount), " ", impName
				
				imp_final, rt_final = processImp(all_imps[i])
				
				if saveResults:
					rtName = impName + "_guv measures for series no: " + str(i + 1) + ".csv"
					rt_final.saveAs(join(results_path, subfolders[1], rtName))
					imp_final_name = impName + '_results_guv_overlay_series_no_' + str(i + 1)+ '.tif'
					FileSaver(imp_final).saveAsTiff(join(results_path, subfolders[0], imp_final_name))
	
				else:
					imp_final.setTitle(impName + "_final_result")
					imp_final.show()
					rt_final.show(impName + " series no " + str(i + 1) + ": guv measures")

		
def processImp(imp0):
	"""
	Main function to process image. 
	Input: image to process
	Returns: results table with intensity measures, image with GUV detection overlay
	"""
	try:
		# get background intensity
		imp0.setC(guv_ch)
		bg_estimation = imp0.getStatistics(ImagePlus.MEDIAN).median

		imp_guv = Duplicator().run(imp0, guv_ch, guv_ch, 1, 1, 1, 1)
		imp_guv.setTitle("imp guv")
#		imp_guv.show()
		imp_membrane = Duplicator().run(imp0, membrane_ch, membrane_ch, 1, 1, 1, 1)
#		imp_membrane.show()

		# preprocess membrane channel to remove noise, droplets inside vesicles 
		imp = Duplicator().run(imp0, membrane_ch, membrane_ch, 1, 1, 1, 1)
		IJ.run(imp, "Subtract Background...", "rolling=10 disable");
		IJ.run(imp, "Median...", "radius=5");

		# use MorpholibJ morpholical segmentation to get labelled image of GUV membranes 
		tolerance = 2
		conn = 4
		dams = True
		ip = imp.getProcessor()
		regionalMinima = MinimaAndMaxima.extendedMinima(ip, tolerance, conn);
		imposedMinima = MinimaAndMaxima.imposeMinima(ip, regionalMinima, conn);
		labeledMinima = BinaryImages.componentsLabeling(regionalMinima, conn, 32);
		result_image= Watershed.computeWatershed(imposedMinima, labeledMinima, conn, dams);
		result_imp = ImagePlus("watershed", result_image);
		result_imp.setCalibration(imp.getCalibration());
		IJ.run(result_imp, "glasbey", "")
#		result_imp.show()

		# convert label map to rois to be added to roi manager
		LabelMapToRoiManagerPlugin.apply(result_imp);
	
		roi_manager = RoiManager.getInstance()
		if not roi_manager:
			roi_manager = RoiManager()
		ov = Overlay()
		rt_final = getResultsTable("GUV measures", True)

		# get rois from roi manager and clear it
		rois = roi_manager.getRoisAsArray()
		roi_manager.reset()

		# get area in calibrated values
		min_area = min_diam * min_diam
		
		# set roi and filter out wrt circularity and area, for remaining rois get intensity from GUV channel
		for i, roi in enumerate(rois):
			imp_membrane.setRoi(roi)
			IJ.run("Set Measurements...", "area shape redirect=None decimal=3");		# need to go this way since API doesnt have way to get circ directly
			IJ.run(imp_membrane, "Measure", "");
	
			rtt = ResultsTable()
			rt = rtt.getResultsTable()
			circ = rt.getValue("Circ.", i)
			area = rt.getValue("Area", i)
			
			if circ > min_circ:
				if area > min_area:
					ov.add(roi)
					imp_guv.setRoi(roi)
					guv_stats = imp_guv.getStatistics()
					rt_final.incrementCounter()
					rt_final.addValue("circularity", circ)
					rt_final.addValue("mean intensity from GUV channel", guv_stats.mean)
					rt_final.addValue("background intensity from GUV channel", bg_estimation)

		imp_membrane.killRoi()
		imp_guv.killRoi()
		grab_results = rt.getResultsWindow()  #close results table from IJ.run(measure)
		grab_results.close(False)

		# prepare final output
		imp_result = RGBStackMerge().mergeChannels([Duplicator().run(imp_membrane), Duplicator().run(imp_guv)], True)
		imp_result.setOverlay(ov)
		
		for c in range(imp_result.getNChannels() + 1):
			imp_result.setC(c + 1);
			IJ.run(imp_result, "Enhance Contrast", "saturated=0.50");

		return imp_result, rt_final

	finally:
		logg = WindowManager.getWindow("Log")
		if logg:
			logg.close()
		rmm = WindowManager.getWindow("ROI Manager")
		if rmm:
			rmm.close()
			
	
def getResultsTable(name, reset=False):
	win = WindowManager.getWindow(name)
	if ( (win!=None) & isinstance(win,TextWindow) ) :	
		rt = win.getTextPanel().getResultsTable();
		if reset:
			rt.reset()
	else:
		rt = ResultsTable();
	return rt
	

	
main()
print "DONE"