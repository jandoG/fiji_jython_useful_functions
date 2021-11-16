import os
from os import path
from os import listdir
from os.path import isfile, join
from java.lang import String
from ij.gui import GenericDialog


#filename = "/Users/nadar/Documents/IA_Projects/isabel_honigmann_protein_concentration_analysis/data/res.csv"
filename = "/Users/nadar/Documents/IA_Projects/isabel_honigmann_protein_concentration_analysis/data/results_folderProcessing/results__multi_position_FUSdrops_488_CSU_F8_05.csv"

separator = ","

gd = GenericDialog("Parameter")

types = ["8-bit", "16-bit", "32-bit", "RGB"]

#gd.addChoice​("Choose one param:", types, "RGB");
#gd.showDialog();


#for t in types:
#	print type(t)


#gd.addChoice("Type:", types, "RGB");
#gd.showDialog();
#choice = gd.getNextChoice();
#print choice

if (isfile(filename)):
	content = open(filename)

	lines = [line.rstrip('\n') for line in content]

	if lines:
		firstline = lines[0].split(separator)
		print firstline
		column_values = [str(firstline[x]) for x in range(0, len(firstline))] 
		for c in column_values:
			print type(c)
			print c
			
		gd.addChoice("Type:", column_values, column_values[0]);
		gd.showDialog();
		choice = gd.getNextChoice();
		print choice