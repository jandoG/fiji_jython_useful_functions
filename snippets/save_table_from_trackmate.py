################# README #######################

# Author Gayathri Nadar, SCF MPI-CBG, nadar@mpi-cbg.de
# Save results after running TrackMate 
#
# Description:
#	This script will save the results tables acquired from TrackMate GUI. 
#	The script takes the name of the active image, renames the tables accordingly and saves to a folder specified by the user.
#	
#
# Usage:
#	- Open this script with FIJI and run it after running TrackMate. 
#	- Ideally, one should have only the current image and tables from TrackMate displayed.
#	- If not, at least make sure to have the correct image activated (click on the image to activate) before running this script.
#	- Also, do not forget to click the Analysis button from TrackMate to get the tables.
#
#
# ===============================================

#@String (visibility=MESSAGE, value="Select active image before running this script!", required=false) msg
#@File (label= "Folder to save results to", style="directory", required= "True") save_folder

from ij import IJ
from ij import WindowManager
from ij.io import FileSaver 
from ij.measure import ResultsTable
from ij.text import TextWindow
from ij.plugin import Duplicator
import os
from os.path import isfile, join

results_path = save_folder.getPath()
global wm
wm = WindowManager 

# get active image 
imp = IJ.getImage(); 
impname = os.path.splitext(imp.getTitle())[0]

# remove C1 if present in image 
if "C1-" in impname:
	finalname = impname.split("C1-")[1]
else:
	finalname = impname 
	
print "Image name to be used for saving: ", finalname

# get all non-image titles 
titles = wm.getNonImageTitles()
#print(titles)

# loop through titles, get windows from trackmate and save 
for t in titles:
	if "Spots" in t or "Links" in t or "Track" in t:
		print "Found window: ", t
		win = wm.getWindow(t)
	
		if win is not None:
			txtpanel = win.getTextPanel()
	
			if "Spots" in t:
				tablename = finalname + '__spots_stats.csv'
				txtpanel.saveAs(join(results_path, tablename))
	
			elif "Links" in t:
				tablename = finalname + '__links_stats.csv'
				txtpanel.saveAs(join(results_path, tablename))
	
			elif "Track" in t:
				tablename = finalname + '__tracks_stats.csv'
				txtpanel.saveAs(join(results_path, tablename))
	
	else:
		IJ.log("No results window found from TrackMate; run analysis and try again!")

print("Done")
