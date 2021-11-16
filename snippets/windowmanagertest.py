from ij import ImagePlus
from ij import WindowManager
from ij.measure import ResultsTable

def getTableName(titles_array):
	for t in titles_array:
		if "all" in t:
			return t

			
wm = WindowManager
titles = wm.getNonImageTitles()
print titles

#particle_m = wm.getWindow(titles[1])
#particle_m.setVisible(False)


name = getTableName(titles)
print name

all = wm.getWindow(name)

print all
all.rename("Results")

rt = ResultsTable()
grab_results = rt.getResultsWindow()
pt = grab_results.getResultsTable()

totaloverlap = pt.getColumn(1)
jaccard = pt.getColumn(2)
dice_coeff = pt.getColumn(3)
volume_similarity = pt.getColumn(4)

grab_results.close(False)



