import os 
from os import path 

processdir = "/Users/nadar/Documents/IA_Projects/salma_eaton_cilia_intensity_measurements/data_202105/210503-AB-optmize-set6-ARL13B-Gli2/B02" 


dirpath = processdir
print "Processing folder =", dirpath

# get list of tif files in the folder 
ext = ".tif"
filelist = sorted([f for f in os.listdir(dirpath) if f.endswith(ext)])

i = 0
for f in filelist:
#	if f.endswith("T0001F001L01A01Z01C01.tif"):
#		print(f)
#
#	if f.endswith('C01.tif'):
#		print(f)
#		i+=1

	if 'F002' in f and 'C04' in f:
		print(f)
		i+=1

#	if 'A02Z01C01' in f:
#		print(f)

#	else:
#		print "File not found"

print "out", i

regex_pattern = "(?P<barcode>.*)_(?P<well>[A-P][0-9]{2})_.*(?P<field>F[0-9]{3}).*(?P<channelNumber>C[0-9]{2}).tif$"