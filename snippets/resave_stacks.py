# File (label = "Enter the directory with the files", style = directory, required = True) processdir 

processdir = "/Users/nadar/Documents/IA_Projects/Hyman_lab/hari_hyman_cell_area_quantifications/data/Montage/All channel Snap11_20210228_84105 PM"

import os 
from os import path 
from os.path import isfile, join
from ij import IJ 
from ij import ImagePlus
from ij.io import FileSaver
from ij.measure import Calibration
from ij.plugin import Duplicator
from ij.plugin import ChannelSplitter, RGBStackMerge, Concatenator, ZProjector
from loci.plugins import BF
from loci.formats import ImageReader
from loci.formats import MetadataTools
from loci.plugins.prefs import OptionsList
from loci.plugins.in import ImporterOptions

def main():
	dirpath = processdir
	print "Processing folder =", dirpath

	results_path = join(os.path.dirname(dirpath), 'results_' + os.path.basename(dirpath))
	if not os.path.exists(results_path):
		os.makedirs(results_path)

	# get list of tif files in the folder 
	ext = ".tif"
	filelist = sorted([f for f in os.listdir(dirpath) if f.endswith(ext)])

	tiles = []
	slices = []
	impname = None

	for f in filelist:
		fname = os.path.splitext(f)[0]
		strs = fname.split("_")
		tiles.append(strs[1])
		slices.append(strs[2])
		impname = strs[0]


	no_tiles = len(set(tiles))
	no_slices = len(set(slices))
	tilenames = sorted(list(set(tiles)))

	print(no_tiles, no_slices)

	newfilelist = (list(chunks(filelist, no_slices)))
	print len(newfilelist)
	print tilenames

	z_projections = []

	for l, tilename in zip(newfilelist, tilenames):
		sequence = []
		n_ch = None 
		n_frames = None
		cal = None
		
		
		for im in l: 
			imp = openImageWithBF(join(dirpath, im), virtual= True) 
			n_ch = imp.getNChannels();
			n_frames = imp.getNFrames();
			cal = imp.getCalibration();
			sequence.append(imp) 

		cc = Concatenator()
		zstack = cc.concatenate(sequence, False)
		zstack.setCalibration(cal)
		zstack.setDimensions(n_ch, no_slices, n_frames)
		imp_max_proj = ZProjector().run(zstack, "max");
		imp_max_proj.setTitle(impname + "_" + tilename)
		z_projections.append(imp_max_proj)


	for imp in z_projections[:2]:
		imp.show();

	


	

		





def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in xrange(0, len(lst), n):
        yield lst[i:i + n]


def openImageWithBF(path, virtual= True):
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























main()
print("Done")