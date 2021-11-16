#@ImagePlus imp0
#@String 	(label="Analyze from" , choices={"semi-manual data", "image metadata"}) analysisType
#@boolean 	(label="plot") plotFit
#@boolean 	(label="reset results tables") resetResultsTable
#@boolean 	(label="save results") saveResults

# Author: Benoit Lombardot, Scientific Computing, Facility, MPI-CBG
# 2018-06-15: version 1.0 MultiFrap script created from frap_analysis_multipoint script. 
#							perform the exponential fit plus a linear fit, 
#							normalize by the whole droplet thus bleaching
#							normalisation by external droplet is not necessary
# 2018-06-15: version 1.1 update normalisation signal to be 'NotFrapped' region, use signal right after frap as 
#							reference signal and change the exponential model to  1-a.exp(-bt)+c to be a.(1-exp(-bt))
#
# 2018-08-22: version 1.2 (Gayathri Nadar) take time unit and frame interval from image itself instead of script parameters.
# 2018-11-08: version 1.3 (Gayathri Nadar) 
#							- use CurveFitter.EXP_RECOVERY_NOOFFSET instead of custom curve model
#							- instead of using frames use time values, i.e in image time units from the calibration 
#							- for fitting the curve use actual time values and not frameInterval x Intensity 
#							- add background subtraction to pre-processing for better segmentation

# previous script history (written for Shamba Hyman lab)
#2016-10-11: version 0.9
#2016-10-7 : version 1.1  	Add results table output, correct issue with roi to have floating point position,
#							allow semi-manual selection of frapped droplets
#2016-10-17: version 1.2 	Cancel the subtraction of background value from the intensity measures (i.e. 'signalFrap',
#							'signalNotFrap')
#2017-08-25: version 1.3	correct to account for a change of the  jython api ( replace 'is' by '==')
#							added support for pixel size (important for tracking)
#							added support for frap region inside droplet (not only on the border of droplet anymore)
#							(analyze  particle does not detect holes so I took the intersection between analyze particle and create selection)
#							change thresholding method from Li to Default (frap region not detected otherwise)



# usage:
# 	to analyze in semi manual mode: use multipoint roi to point the frapped droplet at the frap time.

from ij.gui import Roi
from ij.gui import PolygonRoi
from ij.gui import ShapeRoi
from ij.gui import PointRoi
from ij.gui import Overlay
from ij import ImagePlus
from ij.process import ImageStatistics
from ij.plugin import Duplicator
from ij.plugin import HyperStackConverter
from ij import IJ
from ij.plugin.frame import RoiManager
from ij.gui import Wand
from ij.measure import ResultsTable
from ij import WindowManager
from ij.text import TextWindow
from ij.io import FileSaver
from ij import Prefs
from ij.measure import CurveFitter

from java.awt import Color

import os, re
import math
import time, datetime

from loci.formats import ImageReader
from loci.formats import MetadataTools

def getTimeStamps(imp0):
	fileInfo = imp0.getOriginalFileInfo()
	fileName = None
	if not fileInfo == None :
		if (not fileInfo.directory==None) & (not fileInfo.fileName==None):
			fileName= os.path.join(fileInfo.directory, fileInfo.fileName)
		print "Processing... ", fileName
	
	else :
		print 'the time stamp could not be extracted'
		print 'default time interval will be used'
		nSlice = imp0.getStackSize()
		dt = imp0.getCalibration().frameInterval
		timeUnit = imp0.getCalibration().getTimeUnit()
		if dt == 0 :
			dt =1
			timeUnit = 'frame'
		timevalues = range(0,nSlices*dt)
		
		return timevalues, timeUnit
	
		
	reader = ImageReader()
	omeMeta = MetadataTools.createOMEXMLMetadata()
	reader.setMetadataStore(omeMeta)
	reader.setId(fileName)
	reader.close()
	nPlane = omeMeta.getPlaneCount(0) # we assume there is only one time series
#	print "Plane count= ", nPlane
	timeValues = []
	for plane in range(0,nPlane):
		fr = omeMeta.getPlaneTheT(0,plane)
		ch = omeMeta.getPlaneTheC(0,plane)
		if ch.getNumberValue() == 0 :
			timeValues.append( omeMeta.getPlaneDeltaT(0,plane).value() )
			timeUnit = omeMeta.getPlaneDeltaT(0,plane).unit().getSymbol()
	
	print "Time values= ",timeValues

	return timeValues, timeUnit

global t_val
t_val, t_unit = getTimeStamps(imp0)
dt = imp0.getCalibration().frameInterval
if dt == 0:
	print "Frame interval is zero, will be set to 1"
	dt = 1
	t_unit = 'frame'

print "frame interval=", dt
print "timeunit= ", t_unit

dims = imp0.getDimensions()
nFrame = max(dims[3],dims[4])
imp0 = HyperStackConverter.toHyperStack(imp0, dims[2], min(dims[3],dims[4]), nFrame);
pixelWidth = imp0.getCalibration().pixelWidth

###############################################################
# parse properties to get frap position and frap region shape #
infoRaw = imp0.getProperty("Info")

info = infoRaw.split('\n')
info = [row.strip() for row in info]

#############################################
# get path name
fileInfo = imp0.getOriginalFileInfo();
filePath = ""
if not fileInfo==None :
	if (not fileInfo.directory==None and not fileInfo.fileName==None ):
		filePath = os.path.join(fileInfo.directory, fileInfo.fileName )
		
if analysisType=="image metadata" :

	# get frap frame 
	infoTime = [row for row in info if row.startswith('Repeat T - ')]
#	print "Time info=", infoTime
	
	mTime = re.match('^Repeat T - (\d*) times.*$', infoTime[0])
	frapTime = int(mTime.group(1))+1 # frame numbering starts at 1
	
	# get frap region of interest 
	infoFrap = [row for row in info if row.startswith('_FRAPPA') and 'NumberOfPoints' in row]
	
	frapRegions = []
	count=0
	for row in infoFrap :
		count = count+1
		roiType = row.split('\t')[1]
		posInfo = row.split('\t')[2]
		nPoints = int( posInfo.split(':')[0].replace( 'NumberOfPoints( ','' ).replace( ') ', '' )  )
		posInfo = posInfo.split(':')[1]
		pos = []
		m = re.findall('\(\s?(\d+)\s?,\s?(\d+)\s?\)', posInfo)
		p1 = [int(m[0][0]) , int(m[0][1]) ]
		p2 = [int(m[1][0]) , int(m[1][1]) ]
		roi = Roi(p1[0], p1[1], p2[0]-p1[0], p2[1]-p1[1])
		frapRegionData = {	'roi': 		roi,								\
							'roiType': 	roiType,							\
							'center': 	[(p1[0]+p2[0])/2,(p1[1]+p2[1])/2],	\
							'width':	p2[0]-p1[0], 						\
							'height':	p2[1]-p1[1], 						\
							'origin': 	p1, 								\
							'name':		'frap '+str(count)					\
						 }
		frapRegions.append(frapRegionData)

else : # if analysisType=="semi-manual data" :
	
	roiManager = RoiManager.getInstance()
	
	frapRois = roiManager.getRoisAsArray()[0]
	xPoints = frapRois.getFloatPolygon().xpoints
	yPoints = frapRois.getFloatPolygon().ypoints
	frapRegions = []
	count=0
	for x,y in zip(xPoints,yPoints) :
		count += 1
		x = int(x)
		y = int(y)
		frapRegionData = {	'roi': 		PointRoi(x,y),			\
							'roiType': 	'Point',				\
							'center': 	[x,y],					\
							'width':	0, 						\
							'height':	0, 						\
							'origin': 	[x,y], 					\
							'name':		'frap '+str(count)		\
						 }	
		frapRegions.append(frapRegionData)
	
	frapTime = frapRois.getPosition()
	

		
#############################################
# get background value 
ip = imp0.getStack().getProcessor(frapTime)
bgValue = ImageStatistics().getStatistics(ip,ImageStatistics.MEDIAN, imp0.getCalibration()).median


#############################################
# get segmentation of particles 

def segment(imp):

	IJ.run(imp, "Subtract Background...", "rolling=50 disable");
	IJ.run(imp, "Gaussian Blur...", "sigma=0.70")

	oldPref = Prefs.blackBackground
	Prefs.blackBackground = False # important to do a proper region selection

	IJ.setAutoThreshold(imp, "Default dark")
	IJ.run(imp, "Convert to Mask", "")
	
	# get a roi of the background
	IJ.run(imp, "Create Selection", "")
	rois_refined = imp.getRoi()
	rois_refined = ShapeRoi( rois_refined )
	imp.killRoi()
	
	Prefs.blackBackground = oldPref # revert parameter to older state

	roiManager = RoiManager.getInstance()
	if not roiManager==None :
		roiManager.reset()
	else:
		roiManager = RoiManager()
		
	IJ.run(imp, "Analyze Particles...", "size=40-Infinity pixel add");

	rois = roiManager.getRoisAsArray()
	roiManager.reset()
	particles = []

	
	for roi in rois:
		rectangle = roi.getBounds()
		origin = [rectangle.x, rectangle.y]
		center = [rectangle.x + rectangle.width/2, rectangle.y+rectangle.height/2]

		roi2 = ShapeRoi(roi)
		roi2 = roi2.and( rois_refined )
		
		particleData = {'roi': roi2, 'center': center , 'origin':origin }
		particles.append(particleData)

	return particles

				
impBeforeFrap = Duplicator().run(imp0,1,1,1,1,frapTime-1,frapTime-1)
particlesBeforeFrap = segment(impBeforeFrap)
impAfterFrap = Duplicator().run(imp0,1,1,1,1,frapTime,frapTime)
particlesAfterFrap = segment(impAfterFrap)


######################################################################
# associate region after frap to region before frap if any 
def getDistance(p1, p2):
	d = math.sqrt( math.pow(p1[0]-p2[0],2) + math.pow(p1[1]-p2[1],2) )
	return d

def matchRegions2(regionList1, regionList2, name ):
	n1 = len(regionList1)
	n2 = len(regionList2)
	matchList = []
	for reg1, id1 in zip(regionList1, range(n1)) :
		reg1['roi'+name] = None
		for reg2, id2 in zip(regionList2, range(n2)) :
			d = getDistance(reg1['center'], reg2['center'])
			matchList.append([d, id1, id2])
	
	matchList2 = sorted(matchList, key=lambda match: match[0])
	matchedAlready1 = [False]*n1
	matchedAlready2 = [False]*n2
	for match in matchList2 :
		d = match[0]
		if d>100:
			break
		id1 = match[1]
		id2 = match[2]
		#print d,id1,id2
		if (not matchedAlready1[id1]) and (not matchedAlready2[id2]) :
			regionList1[id1]['roi'+name] = regionList2[id2]['roi']
			#print id2
			matchedAlready1[id1]=True
			matchedAlready2[id2]=True
	#print '---'		
	return regionList1

def getLargestRoi(rois):
	areaMax = 0
	roiMax=None
	for roi in rois :
		area = roi.getStatistics().area
		if area>areaMax :
			areaMax = area
			roiMax = roi
	return roiMax	


particlesBeforeFrap = matchRegions2(particlesBeforeFrap, frapRegions, 'Frap')
particlesBeforeFrap = matchRegions2(particlesBeforeFrap, particlesAfterFrap, 'AfterFrap')

particles = particlesBeforeFrap

for particle in particles :
	if not particle['roiFrap']==None :
		
		roiFrap2    = ShapeRoi(particle['roi'].clone())
		roiNotFrap2 = ShapeRoi(particle['roi'].clone())
		roiAF = particle['roiAfterFrap']
		if roiAF==None :
			roiFrap2 = particle['roi'].clone()
			roiNotFrap2 = None
		else:	
			roiAF = ShapeRoi(roiAF)
			roiFrap2.not(roiAF)
			roiNotFrap2 = roiNotFrap2.not(roiFrap2) #.shapeToRoi() does not work if multiple rois
			
			#roiFrap2 = roiFrap2.shapeToRoi()
			#roiNotFrap2 = roiNotFrap2.shapeToRoi() #ok if only on roi in the region
			roiNotFrap2 = getLargestRoi( roiNotFrap2.getRois() )
			roiFrap2 = getLargestRoi( roiFrap2.getRois() )
				
		particle['roiFrap2'] = roiFrap2
		particle['roiNotFrap2'] = roiNotFrap2
		
		# roiNotFrapped = roiBF-roiFrap2
		
	else :
		particle['roiFrap2'] = None
		particle['roiNotFrap2'] = particle['roi'].clone()
		

#####################################################################
# track a particle as reference particle displacement 

def trackRegion( p0, sequence ):
	nFrame = sequence.getDimensions()[4]
	p=p0
	trajectory = []
	for frame in range(nFrame) :
		ip = sequence.getStack().getProcessor(frame+1)
		wand = Wand(ip)
		wand.autoOutline(int(p[0]),int(p[1]),255,255, Wand.FOUR_CONNECTED)
		if not wand.npoints == 0 :
			roi = PolygonRoi(wand.xpoints, wand.ypoints, wand.npoints, PolygonRoi.POLYGON)
			#print roi.getBounds()
			ip.setRoi(roi)
			ipStats = ImageStatistics().getStatistics(ip,ImageStatistics.SHAPE_DESCRIPTORS, sequence.getCalibration())
			p = [ipStats.xCentroid/pixelWidth , ipStats.yCentroid/pixelWidth ] 
			trajectory.append(p)
		else:
			#print "None"
			trajectory.append(None)
	#print trajectory
#	print "-----------------"
	return trajectory

# create mask of mparticles in the sequence	
sequenceMask = Duplicator().run(imp0)
IJ.run(sequenceMask, "Subtract Background...", "rolling=50 disable stack");
IJ.run(sequenceMask, "Gaussian Blur...", "sigma=0.70 stack")

IJ.run(sequenceMask, "Options...", "iterations=1 count=1 black do=Nothing") # import to set background black before to use the binary tools
IJ.run(sequenceMask, "Convert to Mask", "method=Default background=Dark calculate black")
IJ.run(sequenceMask, "Watershed", "stack")
#IJ.run(sequenceMask, "StackReg", "transformation=[Rigid Body]");

#sequenceMask.show()

#for particle in particles :
#	print particle

# track the particles that are not
for particle in particles :
	
	if particle['roiFrap']==None:
		trajectory = trackRegion(particle['center'], sequenceMask)
		if all(pos is not None for pos in trajectory) :
			particle['trajectory'] = trajectory
			#print "valid trajectory"

count=0
deltaTrajectory = [[0,0]]*nFrame
for particle in particles :
#	print particle.keys()
#	print str([key == 'trajectory' for key in particle.keys()] )
	if any(key == 'trajectory' for key in particle.keys() ):
		count+=1
#		print(count)
		trajectory = particle['trajectory']
		c = trajectory[frapTime-1]

		#for p in trajectory :
		#	print "dx:" + str(p[0]-c[0]) + " ; x:" + str(p[0]) + " ; x0:" + str(c[0])
			
		auxTraj = [  [ p[0]-c[0] , p[1]-c[1] ]  for p in trajectory ]
		#print auxTraj
		deltaTrajectory = [ [dpAcc[0]+dp[0], dpAcc[1]+dp[1]]  for dpAcc, dp in zip(deltaTrajectory, auxTraj )]
		
		#print deltaTrajectory

#print count
deltaTrajectory = [ [dp[0]/count, dp[1]/count] for dp in deltaTrajectory ]

#print deltaTrajectory

#####################################################################
# do the measures 

# area of the measure region and initialisation of signal sequence
for particle in particles :
	particle['signalFrap'] = [0]*nFrame
	particle['signalNotFrap'] = [0]*nFrame
	particle['signalAll'] = [0]*nFrame
	particle['signalNorm'] = [0]*nFrame
	roiF = particle['roiFrap2']
	roiNF =  particle['roiNotFrap2']
	roiF_area = 0
	roiNF_area = 0
	if not roiF==None :
		roiF_area = roiF.getStatistics().area
	if not roiNF==None :
		roiNF_area = roiNF.getStatistics().area
	particle['roiFrap2_Area'] = roiF_area
	particle['roiNotFrap2_Area'] = roiNF_area

# measures signal sequence for each region
for frame, dp in zip(range( nFrame), deltaTrajectory[:]) :
	
	ip = imp0.getStack().getProcessor(frame+1)
	
	for particle in particles :
		roiF = particle['roiFrap2']
		roiNF =  particle['roiNotFrap2']
		
		if not roiF==None:
			bounds = roiF.getBounds()
			roiF.setLocation(bounds.x+dp[0], bounds.y+dp[1])
			ip.setRoi(roiF)
			stats = ip.getStatistics()
			particle['signalFrap'][frame] = (stats.mean)
			roiF.setLocation(bounds.x, bounds.y)
			
		if not roiNF==None:
			bounds = roiNF.getBounds()
			roiNF.setLocation(bounds.x+dp[0], bounds.y+dp[1])
			ip.setRoi(roiNF)
			stats = ip.getStatistics()
			particle['signalNotFrap'][frame] = (stats.mean)
			roiNF.setLocation(bounds.x, bounds.y)
			if not roiF==None:
				particle['signalAll'][frame] = ( particle['signalFrap'][frame]*particle['roiFrap2_Area'] \
										+ particle['signalNotFrap'][frame]*particle['roiNotFrap2_Area'] )\
										/ ( particle['roiFrap2_Area'] + particle['roiNotFrap2_Area'] )
				#particle['signalNorm'][frame] =  (   particle['signalFrap'][frame] - bgValue ) \
				#								 / ( particle['signalAll'][frame]  - bgValue )

for particle in particles :
	if not particle['roiFrap2'] == None :
		sNF0   = sum(particle['signalNotFrap'][:frapTime])/frapTime
		sAll0   = sum(particle['signalAll'][:frapTime])/frapTime
		sFrap0  = sum(particle['signalFrap'][:frapTime])/frapTime
		sFrapMin =  particle['signalFrap'][frapTime]
		particle['signalNorm'] = [ (sFrap-sFrapMin)/(sNF-sFrapMin)*sNF0/sFrap0 for sFrap,sNF in zip(particle['signalFrap'],particle['signalNotFrap']) ]
		
		
# fit each particle with an exponential recovery
# in the range [frap time, frap time + 0.33*tau] do a linear fit

def curve_fit(x, y, frame_interval, curve_model, params0):
	
#	x = [i * frame_interval for i in range( len(y) ) ]
	curvefitter = CurveFitter(x, y)
#	curvefitter.doCustomFit(curve_model, params0, False)
	curvefitter.doFit(CurveFitter.EXP_RECOVERY_NOOFFSET)
	params = curvefitter.getParams()
	params.append(curvefitter.getFitGoodness())
	
	y_fit = [ curvefitter.f(params, xval - x[0]) for xval in x]
	return params, y_fit;


from ij.gui import Plot
from ij.gui import PlotWindow
def plot_data_and_fit(x, y, yfit_exp, yfit_lin , frame_interval, title="test", Yrange =[0, 1.45] ):
	
#	x = [i * frame_interval for i in range( len(y) ) ]

#	print "y ", y
	plot = Plot( title, "Time ("+t_unit+')', "NU", [], [])
	plot.setLimits(0, max(x), Yrange[0], Yrange[1] );
#	plot.setLimitsToFit(True)
	plot.setLineWidth(2)
	
	plot.setColor(Color.BLUE)
	plot.addPoints(x, y, Plot.LINE)
	plot.addPoints(x,y,PlotWindow.X);
	
	plot.setColor(Color.RED)
	plot.addPoints(x, yfit_exp, Plot.LINE)
	plot.setColor(Color.GREEN)
	plot.addPoints(x, yfit_lin, Plot.LINE)

	plot.updateImage()
	plot.show()


stopTime = nFrame
frac = 0.25
count = 0
for particle in particles :
	if not(particle['roiNotFrap2']==None) and not(particle['roiFrap2']==None):
	
		I_norm = particle['signalNorm'][(frapTime):stopTime]
		I_norm_to_plot = particle['signalNorm'][(0):stopTime]
		curve_model = "y = a*(1-exp(-b*x))" #"y = 1-a*exp(-b*x)+a-1" # CurveFitter.EXP_RECOVERY_NOOFFSET;
		#thalf = len(I_norm)*dt/4
		params0 = [ 1 , I_norm[10]/(10*dt)]

		print "fraptime", frapTime
		x = t_val
#		print len(x)
#		print len(I_norm)
		xtofit = x[ frapTime : len(x) ]
		print len(xtofit)
		x0 = xtofit[0]
		xtofit = [ xt-x0 for xt in xtofit]
		
		param, I_fit = curve_fit(xtofit, I_norm, dt, curve_model, params0)

		slope_exp = param[1]*param[0]
		range_exp = param[0]
		print "slope_exp", slope_exp

#		print param
		nFrameLin = max( 10 , int(frac/slope_exp/float(dt)) )
#		print '--'
#		print 'nframe lin: ', int(nFrameLin)
		I_norm2 = particle['signalNorm'][(frapTime-1):(frapTime-1 + nFrameLin)]
		curve_model2 = "y = a*x+b "
		params0 = [param[0] , 1-param[0]-param[1]]

		xtofit1 = x[ (frapTime-1):(frapTime-1 + nFrameLin) ]
		x0 = xtofit1[0]
		xtofit1 = [ xt-x0 for xt in xtofit1]
		
		print len(xtofit1)
		print len(I_norm2)
		param2, I_fit2 = curve_fit(xtofit1, I_norm2, dt, curve_model2, params0)
		slope_lin = float(param2[0])
		
		particle['1/slope (exp)'] = 1/slope_exp
		particle['1/slope (lin)'] = 1/slope_lin
		particle['1/slope norm (exp)'] = range_exp/slope_exp
		particle['1/slope norm (lin)'] = range_exp/slope_lin
		particle['range (exp)'] = range_exp
		particle['GOF (exp)'] = param[-1]
		particle['GOF (lin)'] = param2[-1]

#		print param
#		print param2
		count += 1
		if plotFit :
			plot_data_and_fit( xtofit, I_norm_to_plot, I_fit, I_fit2 , dt, title="vesicle "+str(count) )
			

###############################################################
# visualisation  

def addRoiToOverlay(roi0, ov, color, frame, dp):
	if not roi0==None :
		roi = roi0.clone()
		roi.setStrokeColor(color)
		roi.setPosition(frame)
		bounds = roi.getBounds()
		roi.setLocation(bounds.x+dp[0], bounds.y+dp[1])
		ov.add(roi)
	return ov

	
ov = Overlay()
for particle in particles:
	roiBF = particle['roi']
	roiAF = particle['roiAfterFrap']
	roiF0 = particle['roiFrap']
	roiF  = particle['roiFrap2']
	roiNF = particle['roiNotFrap2']
	
	for frame, dp in zip(range(nFrame) , deltaTrajectory ) :
		ov = addRoiToOverlay(roiF , ov, Color.red , frame+1, dp)
		ov = addRoiToOverlay(roiNF , ov, Color.green , frame+1, dp)
		
			
imp0.setOverlay(ov)

# for each analyzed droplet
# add each analyzed droplet to the to the roi manager
roiManager = RoiManager.getInstance()
count=0
for particle in particles:
	if (not particle['roiFrap']==None) : # and (not particle['roiNotFrap2']==None) :
		count+=1
		roi = particle['roi']
		name = 'vesicle_'+str(count)
		roi.setName(name)
		roi.setPosition(frapTime)
		roiManager.addRoi(roi)
		particle['name'] = name
	
# 2 result table for excel analysis
# one with cell summary info one with time sequence (on column per sequence)
# image name, cell name, area not fraped, area frapped , bg Value, cell type
# sequence data: col name imgName_dropletName_F and imgName_dropletName_NF


def getResultsTable(name, reset=False):
	win = WindowManager.getWindow(name)
	if ( (win!=None) & isinstance(win,TextWindow) ) :	
		rt = win.getTextPanel().getResultsTable();
		if reset:
			rt.reset()
	else:
		rt = ResultsTable();
	return rt


timeStamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y_%m_%d-%H_%M_%S') # if the file is analyzed multiple times

dropletInfoName = 'Vesicles informations'
rt = getResultsTable(dropletInfoName, reset=resetResultsTable)
for particle in particles:
	if not particle['roiFrap2']==None :
		rt.incrementCounter();
		rt.addValue("Sequence name", 				imp0.getTitle() 			);
		rt.addValue("Vesicle name", 				particle['name']			);
		rt.addValue("frapped area", 				particle['roiFrap2_Area']	);
		rt.addValue("not frapped area", 			particle['roiNotFrap2_Area']);
		rt.addValue("background", 					bgValue						);
		rt.addValue("frap frame", 					frapTime-1					);
		rt.addValue("time stamp", 					timeStamp					);
		rt.addValue("1/slope norm (exp)", 			particle.get('1/slope norm (exp)'));
		rt.addValue("1/slope norm (lin)", 			particle.get('1/slope norm (lin)') );
		rt.addValue("1/slope (exp)", 				particle.get('1/slope (exp)')      );
		rt.addValue("1/slope (lin)", 				particle.get('1/slope (lin)')	);
		rt.addValue("range (exp)", 					particle.get('range (exp)')		);
		rt.addValue("GOF (exp)", 					particle.get('GOF (exp)')		);
		rt.addValue("GOF (lin)", 					particle.get('GOF (lin)')		);
		
#rt.show(dropletInfoName)
dropletSequenceName = 'Vesicles signal sequence'
rt2 = getResultsTable(dropletSequenceName, reset=True)	

for particle in particles[:] :
	if not particle['roiFrap2']==None :
		seqF  = particle['signalFrap']
		seqA = particle['signalAll']
		seqN = particle['signalNorm']
		colF = particle['name']+"_Frap_"+imp0.getTitle()
		colA = particle['name']+"_All_"+imp0.getTitle()
		colN = particle['name']+"_Norm_"+imp0.getTitle()
		
		for valF, valA, valN, frame in zip(seqF, seqA, seqN, range(nFrame)):
			rt2.setValue( colF, frame, valF )
			rt2.setValue( colA, frame, valA )
			rt2.setValue( colN, frame, valN )
			
			
#rt2.show(dropletSequenceName)		

# saveResults
if saveResults :
	IJ.log("Saving results ...")
	basename = os.path.split(filePath)[1]
	basename = os.path.splitext(basename)[0]
	savePath = os.path.split(filePath)[0]
	
	saveNameRoiSet = os.path.join( savePath, basename+"_Roiset.zip" )
	roiManager = RoiManager.getInstance()
	roiManager.runCommand(imp0,"Deselect")
	roiManager.runCommand("Save", saveNameRoiSet )

	saveName = os.path.join( savePath, basename+"_vesiclesInfo.xls" )
	rt.saveAs(saveName)

	saveName = os.path.join( savePath, basename+"_vesiclesSequences.xls" )
	rt2.saveAs(saveName)

	savename = os.path.join( savePath, basename+"_results.tif" )
	FileSaver(imp0).saveAsTiff( savename )
	

# 1 result image with overlay of detectected region at all time

# 1 result tables (for R analysis):
# image name, droplet name, area frapped, area not frapped, I mean frapped, I mean not frapped, frame, normalizing factor , Ibg, timestamp


print "DONE"

