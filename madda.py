 
import os
import matplotlib.pyplot as plt
from astropy.table import Column
import matplotlib.patches as mpatches
import pywcsgrid2
import img_scale
from mpl_toolkits.axes_grid1.axes_rgb import imshow_rgb
import numpy as np
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from matplotlib import rc
import numpy.ma as ma
from pyds9 import DS9
from astropy.convolution import convolve, Gaussian2DKernel, Tophat2DKernel
import matplotlib.colors as colors
import pyregion
import glob
from astropy.modeling.models import Gaussian2D
from astropy.convolution import CustomKernel, Model2DKernel 
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset,inset_axes
from numpy import mean, sqrt, square
import pywcsgrid2  
from astropy.wcs import WCS 
import corner
from astropy.table import Table
from matplotlib.patches import Rectangle,FancyArrowPatch

from astropy.io import fits
from matplotlib import rcParams

from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
from mpl_toolkits.axes_grid1 import Grid
from matplotlib import gridspec

from myPYTHON import *

from matplotlib.colors import LogNorm
from progressbar import * 
from imports.gaiaTB import  GAIATB

from imports.gaiaDIS import GAIADIS

from  testPowerLaw.powerLaw import myPowerLaw

from scipy  import special
from mpl_toolkits.axes_grid1 import make_axes_locatable

from astropy.coordinates import SkyCoord

from astropy import units as u
from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
import mpl_toolkits.axes_grid1.axes_grid as axes_grid
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from starlinkTest import myCore
from  distanceTB import disTB



from astropy.modeling import models, fitting

import astropy.units as u
from astropy.utils import data
from spectral_cube import SpectralCube

from HIITB import myHII

from cmfTB import myCMF
#GlobalVariables

core_PIDENT="PIDENT"
core_Peak1="Peak1"
core_Peak2="Peak2"
core_Peak3="Peak3"
core_Cen1="Cen1"
core_Cen2="Cen2"
core_Cen3="Cen3"
core_Size1="Size1"
core_Size2="Size2"
core_Size3="Size3"
core_Sum="Sum"
core_Peak="Peak"
core_Volume="Volume"




class Velo(object):
    def __init__(self, header):
        wcs =  WCS(header)
        self.wcs_vel = wcs.sub([3])

    def to_vel(self, p):
        v = self.wcs_vel.wcs_pix2world([[p]], 0)
        return int(round(v[0][0]) )

    def to_pixel(self, p):
    	 
			
			v = self.wcs_vel.wcs_world2pix([[p]], 0)
			return int(round(v[0][0]) )


def box( centerL,centerB,lSize,bSize,dummy=0):

	"""
	return lRange and B Range
	"""

	lSize=lSize/3600.
	bSize=bSize/3600.
	
	
	return [centerL-lSize/2.,  centerL+lSize/2.   ], [centerB-bSize/2., centerB+bSize/2. ],[centerL,centerB],[lSize,bSize]
	
def circle( centerL,centerB,radius ):

	"""
	return lRange and B Range
	"""
 
	
	return  centerL,centerB,radius/3600.
	




# crop G216

#lRange_Maddalena=[213.5,219.74] #deg 
#bRange_Maddalena=[-4.95 ,-1.25 ]  #deg
#vRange_Maddalena= [16.604, 38.6874] #km/s

#myFITS.cropFITS("mosaic_U_-30_70.fits",outFITS="G216_12.fits",Vrange=[14,38],Lrange=[219.75,213.5],Brange=[-5,-1.25] )



class myG210:
	
	
	#dataPath="./data/"

	dataPath="/home/qzyan/WORK/projects/maddalena/data/"

	
	CO12="mosaic_U.fits"
	CO13="mosaic_L.fits"
	CO18="mosaic_L2.fits"
	
	FUGINCO12=dataPath+"FUGING210M0.fits" 
	
	
	DAMECO12=dataPath+"Wco_DHT2001_C.fits"

	
	
	singleChannelNoise12=0.5 #k
	RMS13=0.25
	RMS12=0.5
 
	
	SMCOFile12=dataPath+"SMCO12.fits"
	
	COFile12= dataPath+ CO12
	COFile13= dataPath+ CO13
	COFile18= dataPath+ CO18

	vRange=[0,70] #kms
	
	
	
	#absolutPath=
	
	YSOTBFile="/home/qzyan/WORK/projects/maddalena/data/YSOMaddalena.fit"
	
	
	figurePath="/home/qzyan/WORK/projects/maddalena/figures/"
	
	
	lCol="_Glon"
	bCol="_Glat"
	classYSO="Class"
	doFITS=myFITS()
	
	WISEHIICatName="/home/qzyan/WORK/projects/maddalena/data/wise_hii_V2.0.csv"
	
	IRASHIIName="/home/qzyan/WORK/projects/maddalena/data/IRAS_HII_candidates.fit"
	
	tempFilsPath="/home/qzyan/WORK/projects/maddalena/tempFiles/"

	figurePath="./home/qzyan/WORK/projects/maddalena/figures/"
	WISE22resolution=0.00038189999999999996 #deg
	WISE22="WISE 22"
	WISE12="WISE 12"	
	
	drawColors={}
	drawColors["K"]='red'
	drawColors["C"]='cyan'
	drawColors["Q"]='yellow'
	drawColors["G"]='green'

	WISEFITS="WISE22.fits"
	
	coint="coint" #actuall ,this is also IRAS
	
	outflowFigurePath="/home/qzyan/WORK/projects/maddalena/outflowFigures/"
	
	
	WISEcropPath="/home/qzyan/WORK/projects/maddalena/wisecrop/"
	
	
	distance= 2450. # pc
	#W345distance= 1950. # pc
	
	G214Distance=2250.
	

	
	W5distance=2300. #pc
	W3distance= 1950. # pc
	W345SplitDL= 135.6  # distance split  

	parsecToMeter= 3.0857e16 #m 
	
	lRange_Maddalena=[213.5,219.74] #deg 
	bRange_Maddalena=[-4.95 ,-1.25 ]  #deg
	vRange_Maddalena= [16.604, 38.6874] #km/s
 


	lRange_S287,bRange_S287,aaaa,aaaa = box(217.7942744,-0.0231650,6669.560,7421.875,0)

	lRange_G211,bRange_G211,aaaa,aaaa =  box(211.7835238,2.4512112,6610.974,6921.311,0)

	lRange_G214,bRange_G214,aaaa,aaaa = box(214.6240124,-1.8017578,4638.337,3909.934,0)
	
	G214L= max(lRange_G214) #less than this 
	G214B= min(bRange_G214) #larger than this
	
	
	
	#lRange_G214,bRange_G214,a1,b1 = box(214.6240124,-1.8017578,4638.337,3909.934,0)
	#lRange_Maddalena=[210.8 ,213]
	#bRange_Maddalena=[-1.94 , 0.3 ]

	lRange_S284=[211.5,212.7]
	bRange_S284=[-1.8,-0.5]
	vRangeb_S284=[38,52]


	


	V_OutArm=38. #obsolaed
	V_PerseusArm=15.

	
	G216FITS="G216.fits"
	S284FITS12= "S284_CO12.fits"
	S284FITS13= "S284_CO13.fits"

	vResolution13= 0.166
	singleChannelRMS13= 0.3

	vResolution12=  0.159
	singleChannelRMS12= 0.5


	W345FITS12="/home/qzyan/WORK/projects/maddalena/data/cropW345CO12.fits"
	W345FITS13="/home/qzyan/WORK/projects/maddalena/data/cropW345CO13.fits" 


	#MWFITS="./data/MWGray.fits"
	
	MWFITS="/home/qzyan/WORK/projects/maddalena/data/agMap.fits"

	

	G210LRange=[209.75,219.75]
	G210BRange=[-5,5]

	allCoreTBFile='/home/qzyan/WORK/projects/maddalena/data/G200220DuchampCO13.fit'
	
	
	#WISEYSOTB=Table.read("/home/qzyan/WORK/projects/maddalena/WISEYSOG210.fit")
	
 	

	outMask12="outMask12.fits"
	perMask12="perMask12.fits"

	outMask13="outMask13.fits"
	perMask13="perMask13.fits"


	extraNameList=[ "G214.8+04.0", "G214.4-04.3", "G211.1-02.1"   ]
	extraLBVList= { "G214.8+04.0":[214.8387749,4.0087410,12.3],  "G214.4-04.3":[214.4089628,-4.3045862,21.6], "G211.1-02.1":[211.1411458,-2.1028260,10.2] } #



	def __init__(self ): 
		
		#read three fits
		
		
		self.CO12HDU= fits.open(self.COFile12)[0]
		
		self.CO13HDU= fits.open(self.COFile13)[0]

		
		
		
		self.CO12DATA,self.COHEADER=self.doFITS.readFITS(self.COFile12)
		self.CO12DATASM,dummy=self.doFITS.readFITS(self.SMCOFile12)
 
 
 
		self.W345CO12DATA,self.W345COHEADER=self.doFITS.readFITS(self.W345FITS12)

		self.WCSW345CO12=WCS(self.W345COHEADER)
 
 
		#self.WCSCO12=WCS(self.convert4DTo3D(self.COHEADER)  )
		self.WCSCO12=WCS( self.COHEADER   )

		self.CO13DATA,dummy=self.doFITS.readFITS(self.COFile13)
		#del self.COHEADER["CTYPE4"]
		#del self.COHEADER["CDELT4"]
		#del self.COHEADER["CRPIX4"]
		#del self.COHEADER["CROTA4"]
		#del self.COHEADER["CRVAL4"] 
		
		
		#aa,self.CO12DATA=myFITS.readFITS(self.CO12)
		self.YSOTB=Table.read(self.YSOTBFile)
	 
		#read wise catalog
		
		self.WISEHIICat=Table.read(self.WISEHIICatName)
		self.IRASHIICat=Table.read(self.IRASHIIName)

		self.gaiaDisDo=GAIADIS()
		self.doPowerLaw=myPowerLaw(2.3,0.5,1000)


	def getLBrange(self,centerLB,size):
		"""
		
		The unit is degree
		
		"""
		
		lRange=[centerLB[0] - size, centerLB[0] + size  ]
		bRange=[centerLB[1] - size, centerLB[1] + size  ]

		return lRange, bRange

	def getWISE(self,centerLB,size):
		"""
		this function is used to to crop wise fits and returen

		"""
		
		outFITSName="cropWISE_{}deg.fits".format(size)
		lRange,bRange= self.getLBrange(centerLB,size) # [centerLB[0] -    ]
		myFITS.cropFITS2D(self.WISEFITS, outFITS=outFITSName,Lrange=lRange,Brange=bRange,overWrite=True)
		return myFITS.readFITS(outFITSName)


	def getAverageSpectraByBox(self,FITSname,centerLB,size):
		"""
		get average spectra by this box

		FITSname indicate it is 12CO or 13CO
		"""
		lRange,bRange= self.getLBrange(centerLB,size) # [centerLB[0] -    ]
		
		outFITSName='cropFITS_avgSpec.fits'
		
		myFITS.cropFITS(FITSname, outFITS=outFITSName,Lrange=lRange,Brange=bRange,overWrite=True)
		cropedFITSData,cropedFITSHead=myFITS.readFITS(outFITSName)

		#np.mean(cropedFITSData,a)
		Nz,Ny,Nx=cropedFITSData.shape

		sumY=np.sum(cropedFITSData,axis=1)
		
		sumXY=np.sum(sumY,axis=1)

		averageSpec=sumXY/1./Nx/Ny
		
		vSIndex=np.arange(Nz)*1.

		WCScrop=WCS(cropedFITSHead)

 
		a,b,velocities=WCScrop.wcs_pix2world(0,0, vSIndex,0 ) 


		return   averageSpec, velocities/1000. 

	def getNoise12(self,vRange):
		return np.sqrt(abs(vRange[0]-vRange[1])/self.vResolution12)*self.singleChannelNoise12*self.vResolution12

	def cutArray(self,inputArray):
		pass

	def drawThreeColor(self,processFile,outFileName,doMoment=True,trimMix=0,trimMax=10):
		
		"""
		draw three colors of 12CO fits
		"""
		drawCO13=False
		
		if "13" in outFileName:
			drawCO13=True
		
		rRange=[0,self.V_PerseusArm] #local arm 
		gRange=[self.V_PerseusArm,self.V_OutArm] #perseus
		bRange=[self.V_OutArm, 70] # out arm
		
		CO12RFile=outFileName+"_R.fits"
		CO12GFile=outFileName+"_G.fits"
		CO12BFile=outFileName+"_B.fits"
		
		
		if doMoment: #do moment
			
			
			
			Rdata,Rhead=doMadd.doFITS.momentFITS(processFile, rRange,0,outFITS=CO12RFile)
			Gdata,Ghead=doMadd.doFITS.momentFITS(processFile, gRange,0,outFITS=CO12GFile)
			Bdata,Bhead=doMadd.doFITS.momentFITS(processFile, bRange,0,outFITS=CO12BFile)

 		else:
			Rdata,Rhead=doMadd.doFITS.readFITS( CO12RFile)
			Gdata,Ghead=doMadd.doFITS.readFITS( CO12GFile)
			Bdata,Bhead=doMadd.doFITS.readFITS( CO12BFile)

		rData=Rdata[0]
		gData=Gdata[0]
		bData=Bdata[0]
		
		
		#draw RGB$^{12}CO ()
			
		img = np.zeros((bData.shape[0], bData.shape[1], 3), dtype=float)

		img[:,:,2]=img_scale.sqrt(rData, scale_min=trimMix, scale_max=trimMax)
		img[:,:,1]= img_scale.sqrt(gData, scale_min=trimMix, scale_max=trimMax)
		img[:,:,0]=img_scale.sqrt(bData, scale_min=trimMix, scale_max=trimMax)
			
		plt.clf()
		#print dir(pywcsgrid2)
		#plt.rcParams['font.sans-serif'] = ['Helvetica']
		fig = plt.figure()
		ax=pywcsgrid2.subplot(111,header=WCS(Bhead))
		#ax.imshow(bData, origin="lower")
		#from astropy.visualization import make_lupton_rgb
		#image = make_lupton_rgb(rData, gData, bData, stretch=0.5)
		#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		ax.set_ticklabel_type("absdeg", "absdeg")
		
		#ax.set_xlabel(r'$l$')
		
		#rc('text.latex', preamble='\usepackage{xcolor}')
		#plt.rcParams['font.sans-serif'] = ['Helvetica']
		#ax.imshow(bData,origin="lower",cmap="jet")zz
		#trimMax=10
		#trimMix=0
		#i = imshow_rgb(ax, self.doFITS.mytrim(rData,trimMix,trimMax) ,self.doFITS.mytrim(gData,trimMix,trimMax) , self.doFITS.mytrim(bData,trimMix,trimMax) , origin="lower", interpolation="none",vmin=trimMix,vmax=trimMax )
		i = imshow_rgb(ax,  img[:,:,0]  , img[:,:,1]  ,  img[:,:,2]  , origin="lower", interpolation="none" )

		#ax.imshow(img[:,:,0],interpolation="none",vmin=0.3,vmax=3 )
 
		#i = imshow_rgb(ax, img[:,:,0] ,img[:,:,1] ,img[:,:,2] , origin="lower", interpolation="none",vmin=2,vmax=15 )
		#axins = zoomed_inset_axes(ax,   loc=3)
		#axins=inset_axes(ax, width="32%", height="25%", loc=1)

		if "12" in outFileName:
			patchLine=mpatches.Patch(color='None',label=r"$^{12}$CO ($J=1\rightarrow 0$)")
		
		else:
			patchLine=mpatches.Patch(color='None',label=r"$^{13}$CO ($J=1\rightarrow 0$)")

		patchRed=mpatches.Patch(color='blue',label=r"the Local Arm (0-15 km s$^{-1}$)")
		#patchGreen=mpatches.Patch(color='green',label=r"14-38 km s$^{-1}$ (the Perseus Arm)")
		#patchBlue=mpatches.Patch(color='blue',label=r"38-70 km s$^{-1}$ (the Outer Arm)")
		patchGreen=mpatches.Patch(color='green',label=r"the Perseus Arm")
		patchBlue=mpatches.Patch(color='red',label=r"the Outer Arm")
		
		patches=[patchLine,patchRed,patchGreen,patchBlue]
		
		
		l=ax.legend(handles=patches,loc=2,handlelength=0.5,handleheight=0.3,fontsize=7,framealpha=0.6)
		
		cs=['black','blue','green','darkred'  ]
		for countN,text in enumerate(l.get_texts()):
			text.set_color(cs[countN])
		
 
		#plt.subplots_adjust(left=0.02, right=1.05, top=0.99, bottom=0.05)
		
		ax.axis[:].major_ticks.set_color("w")
		plt.tight_layout(pad=0)
 




		#lRange_S287,bRange_S287,aaaa,aaaa = box(217.7942744,-0.0231650,6669.560,7421.875,0)
	
		#lRange_G211,bRange_G211,aaaa,aaaa =  box(211.7835238,2.4512112,6610.974,6921.311,0)
	
 
	
		#lRange_G214=[211.5,212.7]
		#bRange_G214=[-1.8,-0.5]
		#vRangeb_G214=[38,52]
	
	





		#draw S287
		#lRange_S287=[216.9,218.5]
		#bRange_S287=[-1 ,1]
		
		#lRange_S287,bRange_S287,a1,b1 = box(217.7942744,-0.0231650,6669.560,7421.875,0)
		if 0:
			self.drawBox(ax,self.lRange_S287,self.bRange_S287 ,color="g")
			#ax["gal"].text(max(self.lRange_S287)-0.1,max(self.bRange_S287)-0.3,"Sh 2-287 (S287)",color="g",fontsize=6)
			
			ax["gal"].text(max(self.lRange_S287)-0.1,max(self.bRange_S287)-0.3,r" G217.7-00.2 (S287)",color="g",fontsize=6)
	
			
			#draw G211
			
			#lRange_G211,bRange_G211,a1,b1 =  box(211.7835238,2.4512112,6610.974,6921.311,0)
	 
			
			self.drawBox(ax,self.lRange_G211,self.bRange_G211 ,color="b")
			ax["gal"].text(max(self.lRange_G211)-0.1,max(self.bRange_G211)-0.25,r"G211.6+02.3",color="b",fontsize=6)
			
			
			
			#draw G214
			
			#lRange_G214,bRange_G214,a1,b1 = box(214.6240124,-1.8017578,4638.337,3909.934,0)
			self.drawBox(ax,self.lRange_G214,self.bRange_G214 ,color="g")
			ax["gal"].text(max(self.lRange_G214)-0.15,max(self.bRange_G214)-0.2,r"G214.6-01.8",color="g",fontsize=6)
			
	 
			
			#draw Maddlena
	 
	 
			self.drawBox(ax,self.lRange_Maddalena,self.bRange_Maddalena,color="g")
			ax["gal"].text(max(self.lRange_Maddalena)-0.1,min(self.bRange_Maddalena)+0.1,r"G216.2-02.5 (G216, Maddalena)",color="g",fontsize=6)
	 
			# SNR G211.7-01.1
			#draw Maddlena
			#lRange_S284=[210.8 ,213]
			#bRange_S284=[-1.94 , 0.3 ]
			self.drawBox(ax,self.lRange_S284,self.bRange_S284,color="r")
			#ax["gal"].text(max(lRange_Maddalena)-0.1,max(bRange_Maddalena)-0.2,r"G212-1.3 (\mbox{H\,\textsc{ii}}  region)",color="green",fontsize=6)
			#ax["gal"].text(max(self.lRange_S284)-0.1,min(self.bRange_S284)+0.1,r"Sh 2-284 (S284)",color="b",fontsize=6)
			ax["gal"].text(max(self.lRange_S284)-0.1,min(self.bRange_S284)+0.1,r"G212.2-00.9",color="r",fontsize=6) #(S284)
	
			#13CO cut myFITS.cropFITS("mosaic_L.fits",Vrange=[17.6,38],Lrange=[210.8 ,213],Brange=[-1.94 , 0.3 ],outFITS="maddCO13.fits")



		orderDict=doMadd.getNDict()
		fontsize=7
		if 1: # draw contours of clouds 
			
			#cIDs=[0, 2, 4, 5, 6, 8, 10, 14, 18, 22, 27, 28, 32, 39, 40, 44, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56]
			cIDs=[0, 2, 4, 5, 6, 33, 21 , 10, 14, 18, 22, 27, 28, 32, 39, 40, 44, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56]

			for ID in cIDs:
				maskPath="/home/qzyan/WORK/projects/maddalena/distanceData/allFigureComplete/"
				
				doneTxt=False
				
				cloudName= doMadd.getCloudNameByID(ID) 
				
				N=orderDict[cloudName]
				
				l,b,v=self.getLBVByID(ID)
				
				cloudID='Cloud{}'.format(ID)
				
				maskFITS=		glob.glob(maskPath+"*"+cloudID  +"_*mask.fits")[0]
				intFITS=		glob.glob(maskPath+"*"+cloudID  +"_*int.fits")[0]
				
				#contour
				
				
				maskData,maskHead=myFITS.readFITS(maskFITS)
				if not drawCO13:
					ax.contour(maskData, colors= "white", linewidths=0.25,origin="lower",levels=[1],alpha=0.8)

				offsetL=0
				offsetB=0



				if N==5:
					offsetL=0.1
					offsetB= -0.2
					ax["gal"].text(l+offsetL,b+offsetB,str(N),color='white',fontsize=fontsize)
					##ax["gal"].arrow(l+offsetL,b+offsetB,-offsetL , -offsetB  , head_width=0.05, head_length=0.05, fc='w', ec='w',lw=0.2,capstyle="butt")
					
					#ax["gal"].annotate('local max', xy=(2, 1), xytext=(100, 100), arrowprops=dict(facecolor='white', shrink=0.05),color="white" )
					#ax["gal"].annotate(str(N) , xy=(l  , b), xytext=(l+offsetL, b+offsetB), arrowprops=dict(facecolor='white', shrink=0.01,headwidth=0.3, width =0.1),color="white",horizontalalignment='center'   )



					doneTxt=True

				else:
					offsetL=0.1
					offsetB= -0.05

				if N==7:

					offsetB= -0.1


				if not doneTxt:

					ax["gal"].text(l+offsetL,b+offsetB,str(N),color='white',fontsize=fontsize)


		if 1: #draw extra files
			
			for eachExtra in self.extraNameList:
				
				serachName= eachExtra.replace("+0","+")
				
				serachName= serachName.replace("-0","-")

				intFITS="/home/qzyan/WORK/dataDisk/G210/distance/"+serachName+".fits"

				dodis=disTB()
				
				
				
				extraRow=dodis.getRowByName( eachExtra )
				
				l=extraRow["l"]
				b=extraRow["b"]

				#print extraRow["l"]
				
				maskData,data,head=dodis.getMask( eachExtra)
				
				conterData=maskData*data
				
				drawContour= 4. # K km/s
				if not drawCO13:
					ax.contour(conterData, colors= "white", linewidths=0.25,origin="lower",levels=[drawContour],alpha=0.8)

				N=orderDict[eachExtra]







				ax["gal"].text(l,b,str(N),color='white',fontsize=fontsize , weight='ultralight')

				
			
			#extraNameList  extraLBVList


		fig.savefig(outFileName+'.pdf', bbox_inches="tight")
		fig.savefig(outFileName+'.png', bbox_inches="tight",dpi=300)





	def drawYSOSpectraOne(self):
		"""
		Draw 12CO YSO Spectra one 
		
		"""
		
		
		testRow=self.YSOTB[0]
 
		spectrum,velo=myFITS.getSpectraByLB(self.CO12DATA,self.COHEADER,testRow[self.lCol],testRow[self.bCol])
 
		f=plt.figure(figsize=(8,6))
		ax=f.add_subplot(111)
 
		ax.step(velo,spectrum)
		
		plt.savefig("spectraTestOne.eps",bbox_inches='tight')

	def drawYSOSpectraAll(self,drawClass='Proto',smooth=True):
		"""
		Draw 12CO YSO Spectra one 
		
		"""
  
		useData=self.CO12DATASM
 		
 		#for eachRow in  self.YSOTB:
 
		N=33
		namePrefix="Protostar"
		if drawClass=="Disk":
			N=41
			namePrefix="Young Star"
		
		nameSurfix="SM"
		
		if not smooth:
			nameSurfix="NOSM"
			useData=self.CO12DATA
		
		fig=plt.figure(figsize=(8,130))
		
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		grid=ImageGrid(fig,111,nrows_ncols=(N,1), axes_pad=0.1,share_all=True)
		
		
		axIndex=0
		for eachRow in self.YSOTB:
		#testRow=self.YSOTB[0]
			#print eachRow[self.classYSO],drawClass
 			if str.strip(eachRow[self.classYSO])!=drawClass:
 				continue
  
 			 
			spectrum,velo=myFITS.getSpectraByLB(useData,self.COHEADER,eachRow[self.lCol],eachRow[self.bCol])
			spectrum13,velo13=myFITS.getSpectraByLB(self.CO13DATA,self.COHEADER,eachRow[self.lCol],eachRow[self.bCol])

 
		#ax=f.add_subplot(3111)
			#print axIndex,eachRow[self.lCol],eachRow[self.classYSO]
			grid[axIndex].step(velo,spectrum,lw=0.5, color='blue',label=r"12CO" )
			grid[axIndex].step(velo13,spectrum13,lw=0.5, color='red',label=r"13CO" )
			
			
			grid[axIndex].plot( velo,velo*0,lw=1,color='black',label="{} {}".format(namePrefix,axIndex+1))
			
			grid[axIndex].legend(loc=1)
			grid[axIndex].set_ylabel(r'$\rm T_{\rm mb}$ (K)')
			axIndex=axIndex+1
			
			
		grid[-1].set_xlabel(r"$\rm V_{\rm LSR}$ (km s$^{-1}$)")	
			
		#ax=f.add_subplot(3112)
		#ax.step(velo,spectrum,lw=0.5)
		
		 
		
		plt.savefig("YSO_{}_{}.pdf".format(drawClass,nameSurfix),bbox_inches='tight')



	def filterTBByRange(self,TB,colname,valueRange):
		
		#No none value are concerned
		
		processTB=TB.copy()
		processTB.add_index(colname)
		returnTB=processTB.loc[colname,min(valueRange):max(valueRange)]
		try:
			returnTB.remove_indices(colname)
		except:
			pass
		return returnTB





 	def drawHIIRegions(self,backGroundFITS):
		"""
		"""
		Bdata,Bhead=doMadd.doFITS.readFITS( backGroundFITS)
		
		if len(Bdata.shape)==3:
			Bdata=Bdata[0]
		
		
		WISEcat=self.WISEHIICat.copy()
		
		colL="GLong<br>(deg.)"
		colB="GLat<br>(deg.)"
		
		lMin,lMax=min( self.G210LRange), max( self.G210LRange), 
		
		
		maddaHII=self.filterTBByRange(WISEcat,colL, [lMin,lMax ])

		maddaHII=self.filterTBByRange(maddaHII,colB,[-5.,5.])

		
		MaddaIRASHII=self.filterTBByRange(self.IRASHIICat,"glon",[lMin-360,lMax-360])

		MaddaIRASHII=self.filterTBByRange(MaddaIRASHII,"glat",[-5.,5.])


		print "Total number of HII reions?",len( MaddaIRASHII)+len(MaddaIRASHII)
		
		
		#darw my IRAS HII region
		
		
		#for eachHII in maddaHII:
			
			#print eachHII[colL]
		
		#filter with l, b
		
		
		
		fig = plt.figure()
		
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		ax=pywcsgrid2.subplot(111,header=WCS(Bhead))

		rmsnoise=2.9

		ax.imshow(Bdata,origin='lower',cmap="bone",vmin=1,vmax=20,interpolation='none')


		print maddaHII.colnames

		
		from astropy.coordinates import SkyCoord
		from astropy import units as u
		
		wiseCatalog=SkyCoord( maddaHII[colL]*u.deg,  maddaHII[colB]*u.deg,frame="galactic")
		
		
		Ls=[];Bs=[]
		
		for eachIRASHII in MaddaIRASHII:

			Ls.append(eachIRASHII["glon"])
			Bs.append(eachIRASHII["glat"])
			#print eachIRASHII["glon"],eachIRASHII["glat"]
			
		c=SkyCoord(np.array(Ls)*u.deg,np.array(Bs)*u.deg,frame="galactic")
		
		#c.separaton(wiseCatalog)<0.1*u.deg
		
		#scalarc = SkyCoord(1*u.deg, 2*u.deg)
		idxc, idxcatalog, d2d, d3d= wiseCatalog.search_around_sky(c, 3./60*u.deg)  
			
		print max(idxc),max(idxcatalog)
		
		
		for i in range(len(MaddaIRASHII)):
			if i not in list(set(idxc)):
				nameFeature="{:.2f}_{:.2f}_{:.1f}_{:.1f}".format(Ls[i],Bs[i],0.5,0.5)
				print "WISE22_{}.fits".format(nameFeature), Ls[i]+360,Bs[i]
		
		
		#print len(c),len(wiseCatalog)
			
			#doMadd.drawPVslice([eachIRASHII["glon"],eachIRASHII["glat"]],[0.5,0.5],45,0.1)
		
		
		#draw WISE
		if 1:
			
			for eachWISEHII in  maddaHII:
				ax["gal"].scatter(eachWISEHII[colL],eachWISEHII[colB],color= "blue",s=10,facecolors='None',lw=0.4)

				
				#ax["gal"].scatter(eachWISEHII[colL],eachWISEHII[colB],color=self.drawColors[eachWISEHII["Catalog"]],s=eachWISEHII["Radius<br>(arcsec.)"]/5.,facecolors='None',lw=0.5)
		 
				
		#if 1: #draw IRAS
			
 
			#ax["gal"].scatter(  MaddaIRASHII["glon"],MaddaIRASHII["glat"],color='blue',facecolors='None',lw=0.5,s=10 )

		if 1: #draw my own IRAS
			for i in range(len(MaddaIRASHII)):
				if i not in list(set(idxc)):
					#nameFeature="{:.2f}_{:.2f}_{:.1f}_{:.1f}".format(Ls[i],Bs[i],0.5,0.5)
					#print "WISE22_{}.fits".format(nameFeature), Ls[i]+360,Bs[i]
					ax["gal"].scatter(  Ls[i]+360 ,Bs[i],color='red',facecolors='None',lw=0.4,s=10 )
					#ax["gal"].text(  Ls[i]+360 ,Bs[i], "({},{})".format(Ls[i]+360,Bs[i]),color='red',fontsize=2)

		
		#Examine the WISE image of all the HII regions, to
		if 0:
			for eachWISEHII in maddaHII:
			
				doMadd.drawPVslice([eachWISEHII[colL],eachWISEHII[colB]],[0.5,0.5],45,0.1)
			for eachIRASHII in MaddaIRASHII:
			
				doMadd.drawPVslice([eachIRASHII["glon"],eachIRASHII["glat"]],[0.5,0.5],45,0.1)
			 
		#ax["gal"].set_xlim([ 219.5,214])
		#ax["gal"].set_ylim([ -5,2])
		#ax.set_xlim([ 0,750])
		#ax.set_ylim([ 0,750]) 
		
		self.drawBox(ax,self.lRange_S287,self.bRange_S287 ,color="g")
		ax["gal"].text(max(self.lRange_S287)-0.1,max(self.bRange_S287)-0.3,"Sh 2-287 (S287)",color="g",fontsize=6)
		
 
		
		#draw G211
		
		#lRange_G211,bRange_G211,a1,b1 =  box(211.7835238,2.4512112,6610.974,6921.311,0)
		
		
		self.drawBox(ax,self.lRange_G211,self.bRange_G211 ,color="r")
		ax["gal"].text(max(self.lRange_G211)-0.1,max(self.bRange_G211)-0.3,"G211",color="r",fontsize=6)
		
		
		
		#draw G214
		
		#lRange_G214,bRange_G214,a1,b1 = box(214.6240124,-1.8017578,4638.337,3909.934,0)
		self.drawBox(ax,self.lRange_G214,self.bRange_G214 ,color="g")
		ax["gal"].text(min(self.lRange_G214)+0.53,max(self.bRange_G214)-0.2,"G214",color="g",fontsize=6)
		
 
		
		#draw Maddlena
 
 
		self.drawBox(ax,self.lRange_Maddalena,self.bRange_Maddalena,color="g")
		ax["gal"].text(max(self.lRange_Maddalena)-0.1,min(self.bRange_Maddalena)+0.1,"G216-2.5 (G216, Maddalena)",color="g",fontsize=6)
 
		# SNR G211.7-01.1
		#draw Maddlena
		#lRange_S284=[210.8 ,213]
		#bRange_S284=[-1.94 , 0.3 ]
		self.drawBox(ax,self.lRange_S284,self.bRange_S284,color="b")
		#ax["gal"].text(max(lRange_Maddalena)-0.1,max(bRange_Maddalena)-0.2,r"G212-1.3 (\mbox{H\,\textsc{ii}}  region)",color="green",fontsize=6)
		ax["gal"].text(max(self.lRange_S284)-0.1,min(self.bRange_S284)+0.1,r"Sh 2-284 (S284)",color="b",fontsize=6)
		
		#13CO cut myFITS.cropFITS("mosaic_L.fits",Vrange=[17.6,38],Lrange=[210.8 ,213],Brange=[-1.94 , 0.3 ],outFITS="maddCO13.fits")
		





		ax.axis[:].major_ticks.set_color("w")
		plt.tight_layout(pad=0)
		fig.savefig( 'HIIMadda.pdf', bbox_inches="tight")
		#fig.savefig( 'YSOMadda.eps', bbox_inches="tight")


	def drawBox(self,axFITS,Lrange,Brange,color='w'):
		
		
			colorStr='{}--'.format(color)
			axFITS["gal"].plot([Lrange[0],Lrange[0]],Brange,colorStr,lw=0.4)
			axFITS["gal"].plot([Lrange[1],Lrange[1]],Brange,colorStr,lw=0.4)
			
			axFITS["gal"].plot(Lrange,[Brange[0],Brange[0]],colorStr,lw=0.4)
			axFITS["gal"].plot(Lrange,[Brange[1],Brange[1]],colorStr,lw=0.4)
		


	def extractProtoStarSpectra(self,backGroundFITS):
		"""
		default spectra is 12CO, because C13O, and C18O is weak
		"""
		
		
		#readBackgoundFITS and draw.
		
		
		Bdata,Bhead=doMadd.doFITS.readFITS( backGroundFITS)
 
		
		
		
		fig = plt.figure()
		
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		ax=pywcsgrid2.subplot(111,header=WCS(Bhead))

		ax.imshow(Bdata,origin='lower',cmap="bone",vmin=0,vmax=20,interpolation='none')

		
		#plot 
		
		for eachYSO in self.YSOTB:
			#Proto Disk
			#if eachYSO["Class"]!="Proto":
				#continue
				
			#else:
			ax["gal"].scatter(eachYSO[self.lCol],eachYSO["_Glat"],color='red',s=1 )
			 
		#ax["gal"].set_xlim([ 219.5,214])
		#ax["gal"].set_ylim([ -5,2])
		ax.set_xlim([ 0,750])
		ax.set_ylim([ 0,750]) 




		ax.axis[:].major_ticks.set_color("w")
		plt.tight_layout(pad=0)
		fig.savefig( 'YSOMadda.eps', bbox_inches="tight")
		#fig.savefig( 'YSOMadda.eps', bbox_inches="tight")


	def drawCMFByVrange(self,vRange,nameSuffix="Maddalena"):
		
		"""
		"""
		
		startV,endV=vRange
		
		coreTB=self.getClumpTB("/home/qzyan/WORK/projects/maddalena/data/clumps/smClumps/SMclumpFindCors.txt")



		#calculate the physical range of velocities
		
		wcs = WCS(self.COHEADER)
		startIndex=wcs.all_world2pix(215,0,startV*1000,0)[2]
		endIndex=wcs.all_world2pix(215,0,endV*1000,0)[2]
 
		coreTB.add_index(core_Cen3)
		

		drawTB=coreTB.loc[core_Cen3,startIndex:endIndex]
		
 		print "Drawing {} cores for {}".format(len(drawTB), nameSuffix)
		


		#filterTable According to vrange
		#print coreTB[40][core_Cen1],coreTB[40][core_Cen2],coreTB[40][core_Cen3]
		


		fakeMass=np.log(drawTB[core_Sum])
		
		hist,edges=np.histogram(fakeMass,bins=20)
		histCenter=(edges[1:]+edges[0:-1])/2.
		

		fig=plt.figure(figsize=(8,6))
		ax=fig.add_subplot(111)
		
		ax.step( histCenter , hist ,lw=1.5,c="blue",label= nameSuffix)
		ax.set_yscale('log')
		
		ax.legend(fontsize=13)
		
		ax.set_xlabel(r"Core mass (arbitrary unit)")	
		ax.set_ylabel(r"Number")	

		
		plt.savefig("CMF{}.eps".format(nameSuffix),bbox_inches='tight')

 

	def drawDistance(self):
		
		"""
		calculate the distance of Gaia 
		
		"""

		print "Draw the distance of maddalena"

		
		
		Lrange=[214,219]
		Brange=[-4.8,-1.5] #degree

		#calWithFITS(self,fitsName,Lrange,Brange,noise,signalLevel,noiseLevel) 

		noise=self.getNoise12(self.vRange)
		noise=2.9
		
		
		#print 2.9,self.getNoise12(self.vRange)
		signalLevel=3
		noiseLevel=1
		#use 12CO 
		fitsFileName="mosaic_U_M0.fits"

		self.gaiaDisDo.calWithFITS(fitsFileName,Lrange,Brange,noise,signalLevel,noiseLevel) 
		
		
		
		
		
		


	
	def getClumpTB(self,coreFile):
		"""
		read core table producd with Starlink
		"""
		
		 
		

 
		#startLineMark="BEGINTABLE"
		
		colNames=[core_PIDENT,core_Peak1,core_Peak2,core_Peak3,core_Cen1,core_Cen2, \
		 core_Cen3,core_Size1,core_Size2,core_Size3,core_Sum,core_Peak,core_Volume]
		
		f = open(coreFile, "r")
		
		lines = f.readlines()
		
		f.close()
		
		
		dataRows=[]
		
		for eachLine in lines:
			
			
			eachLine=eachLine.replace("D","e")
			lineSplit=eachLine.split()
			if lineSplit==[]:
				continue
				
			lineHead=str.strip(lineSplit[0])
			if not str.isdigit(lineHead):
				continue
			
			dataRows.append(map(float,lineSplit) )

		
		t=Table(rows=dataRows,names=colNames)

		return t
		
	def downLoadSkyView(self):
		"""
		
		a function used to download skyview data 
		"""
		
		
		myFITS.downLoadSkyview("WISE 22",[0,0],savePath="./")
		pass

	def drawPVslice(self, LB,sizeLB,pathAngle,pathWidth):
		"""
		A test of draing PPV diagram, based on a position,size,angle,with
		
		We need to dowload wise 24 um and crop 12CO fits,
		
		Do a specific example and generize it.
		
		position and size is in degree, while pathAngle is in the unit of degree...
		"""
		
		centerL,centerB=LB
		sizeL,sizeB=sizeLB
		
		lRange=[centerL-sizeL/2.,centerL+sizeL/2.]
		bRange=[centerB-sizeB/2.,centerB+sizeB/2.]
		vRange=[29,45] #km/sovided for myFITS

		#crop 12CO fits
		#tempFilsPath
		
		nameFeature="{:.2f}_{:.2f}_{:.1f}_{:.1f}".format(centerL,centerB,sizeL,sizeB)
		
		crop12COName="Crop_{}_12CO.fits".format(nameFeature)
		
		crop12COFile=self.tempFilsPath+crop12COName
		
		if os.path.isfile(crop12COFile):
			print "....12CO croped file exists...."
			pass
		else:
			print "....Croping 12CO fits....."
			
			myFITS.cropFITS( self.COFile12, outFITS=crop12COFile,Vrange=vRange,Lrange=lRange,Brange=bRange)				

 
		wise22Name="WISE22_{:.2f}_{:.2f}_{:.1f}_{:.1f}.fits".format(centerL,centerB,sizeL,sizeB)
		wise22File=self.tempFilsPath+wise22Name
		if os.path.isfile(wise22File):
			print "....WISE 22 micron fits file exists...."
			pass
		else:
			print "....Downloading WISE 22 fits....."
			myFITS.downLoadSkyview(self.WISE22,LB,sizeLB=sizeLB, \
				resolution=self.WISE22resolution,savePath=self.tempFilsPath,overWrite=False,saveName=wise22Name)

		#draw 12CO integration, and 
		
		#moment 12CO
		crop12COMoment=self.tempFilsPath+"Crop_{:.2f}_{:.2f}_{:.1f}_{:.1f}_12CO_M0.fits".format(centerL,centerB,sizeL,sizeB)

		myFITSDo=myFITS()

 
		myFITSDo.momentFITS(crop12COFile,[30,42],0,outFITS=crop12COMoment ,cutEdge=False) 
		
		COM0data,COM0head=doMadd.doFITS.readFITS( crop12COMoment)
		
		# draw 12CO, and wise
		fig = plt.figure()
		ax=pywcsgrid2.subplot(121,header=WCS(COM0head))
		#ax.imshow(bData, origin="lower")
		#from astropy.visualization import make_lupton_rgb
		#image = make_lupton_rgb(rData, gData, bData, stretch=0.5)
		#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		
		#draw CO =========================
		ax.set_ticklabel_type("absdeg", "absdeg")
		ax.imshow(COM0data[0],origin="lower",vmin=2,vmax=20,cmap="bone")
		
 
		
		ax["gal"].scatter(centerL,centerB,edgecolors='r',lw=0.5,s=30 ,facecolors="none")


		at = AnchoredText(r"$^{12}$CO($J\rightarrow 1-0$)", loc=1, frameon=True)
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		ax.add_artist(at)




		#draw WISE

		WISEData,WISEHead=doMadd.doFITS.readFITS( wise22File)
		axWISE=pywcsgrid2.subplot(122,header=WCS(WISEHead))
		axWISE.set_ticklabel_type("absdeg", "absdeg")
		
		#axWISE.imshow(WISEData,origin="lower",norm=LogNorm(vmin=160.353,vmax=166),cmap="bone" ,interpolation="none" )
		axWISE.imshow(WISEData,origin="lower", cmap="bone" ,interpolation="none" )

		
		
		axWISE["gal"].scatter(centerL,centerB,edgecolors='r',s=30,lw=0.5 ,facecolors="none")
		at = AnchoredText(r"WISE 22 $\mu$m", loc=1, frameon=True)
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axWISE.add_artist(at)
		if 0:#draw counter of 12CO,
			
			colevel=np.arange(4,30,3)
			
			axWISE[WCS(COM0head)].contour(COM0data[0],levels=colevel,linewidths=0.1,cmap="Reds")
		


		
		
		
		
		fig.savefig(self.figurePath+nameFeature+'.png', bbox_inches="tight",dpi=600)
		fig.savefig(self.figurePath+nameFeature+'.eps', bbox_inches="tight" )



	def examineXYDis(self,backGroundFITS):
		
		
		"""
		
		This function is used to examine the distribution of intensity with respect to x and y
		
		are they indenpended?
		
		"""
		
		#first, read the moment fits  
		Bdata,Bhead=doMadd.doFITS.readFITS( backGroundFITS)

		processData=Bdata[0]

		print processData.shape
		
 
		#get the noise level of the fits
		
		radius=1.1380966
		
		noiseLrange=[218.2568544-radius,218.2568544+radius]
		noiseBrange=[3.84812993-radius,3.84812993+radius]
		
		cutNoiseFITS="./noise.fits"
		
		myFITS.cropFITS2D( backGroundFITS, outFITS=cutNoiseFITS, Lrange=noiseLrange,Brange=noiseBrange,overWrite=True)				
 
		noiseData,noisehead=doMadd.doFITS.readFITS( cutNoiseFITS)

		sigma= myFITS.getRMS(noiseData)
		#noise got

		processData[processData<3*sigma]=0


 
		ny,nx=processData.shape
		
		XLists=[]
		YLists=[]
		
		for i in range(nx):
			
			XLists.append( np.sum(processData[:,i]) )
		for i in range(ny):
			
			YLists.append( np.sum(processData[i,:]) )
		
		
		#plotXLists
		n, bins, patches = plt.hist(XLists, 50, density=True, facecolor='g', alpha=0.75)
		n, bins, patches = plt.hist(YLists, 50, density=True, facecolor='b', alpha=0.75)

		plt.show()




	def getCointByLB(self,bgData,bgWCS,l,b):
		"""
		pass
		"""
		#COdata,COWCS,
		Ny,Nx=bgData.shape  # by default 2d
		
		x,y=bgWCS.all_world2pix(l,b,0)
		x=int(round(x))
		y=int(round(y))
		if x<0 or x>Nx-1 or y<0 or y>Ny-1:
			return -100
		return  bgData[y,x]



	def getOnsourceGaia(self,GaiaTB,bgData,bgHead ):
 
		GaiaTB= GaiaTB.copy()
		#print "processing {} gaia sourxes".format(len(GaiaTB))
		#cloudMark="cloudMark"
		
		if len(bgData.shape)==3:
			
			bgData=bgData[0]
 
		addCol=Column(name=self.coint,data=np.zeros(len(GaiaTB)))
		GaiaTB.add_column(addCol)
 
		bgWCS=WCS(bgHead)
		
		Ls=GaiaTB["l"]
		Bs=GaiaTB["b"]
 
		Xs,Ys=bgWCS.all_world2pix(Ls,Bs,0)
 		
 		Xs=map(round,Xs)
 		Ys=map(round,Ys)
 		Xs=map(int,Xs)
 		Ys=map(int,Ys)
 
		#[y,x]

		for i in range(len(Xs)):
			GaiaTB[i][self.coint]= bgData[Ys[i]][Xs[i]]  #self.getCointByLB(bgData, WCS(bgHead), eachStar["l"],eachStar["b"])
  
		return  GaiaTB 

	def convert3DTo2D(self,data,header):
		if len(data.shape)==3 and data.shape[0]==1:
		
			print "Reducing to 2D...."
			data=data[0]
			header["NAXIS"]=2
			
			del header["CRPIX3"]		
			del header["CDELT3"]		
			del header["CRVAL3"]		
			del header["CTYPE3"]	

		return data,header

	def convert4DTo3D(self, header):


		header=header.copy()
		print "Reducing to 3D...."
		
		header["NAXIS"]=3
		
		del header["CRPIX4"]		
		del header["CDELT4"]		
		del header["CRVAL4"]		
		del header["CTYPE4"]	

		return  header

	def smoothAg(self, list1, list2, smoothScale=5):

		#but found, the list1Error should be sorted together
		#every smoothScale
		
		#bin the data according to list1, which 1
		
		radius=smoothScale/2.
		
		list1Copy=list1.copy()
		
		NewList2=[]
		
		minDis=list1.min()
		maxDis=list1.max()
		
		
		
		#print "What happened?",minDis,maxDis+smoothScale,smoothScale
		
		newDisRange=np.arange(minDis,maxDis+smoothScale,smoothScale)

		newDisList=[]
		newExtictionList=[]
		
		error1=[]
		error2=[]
		
		for countN,eachD in enumerate(newDisRange[0:-1]):
			
			beginDis=eachD
			endDis=newDisRange[countN+1]
			
			#print beginDis,endDis
			cut1=list1Copy[list1Copy>=beginDis]
			cut2=list2[list1Copy>=beginDis]
			
			cut22=cut2[cut1<=endDis]
			cut12=cut1[cut1<=endDis]

			if len(cut12)==0:
				continue
			newDisList.append(np.mean(cut12))
			newExtictionList.append(np.mean(cut22))
			error1.append(np.var(cut12,ddof=1))
			error2.append(np.var(cut22,ddof=1))
		#print newDisList
		#print newExtictionList
			
		return np.array(newDisList),np.array(newExtictionList),error1,error2



	#

	def getGaiaDR2Stars(self,backgroundFITS, noiseLevel=3, smoothInterval=5):
		
		
		"""
		backgroundFITS should be a 2D fits, noiseLevel is the noise level of the backgroud
		
		default smooth is 5 pc
		
		
		
		"""
		gaiaDo=GAIATB()

		colParallax="parallax"
		
		colAg="a_g_val"

		bgData,bgHead= doMadd.doFITS.readFITS(backgroundFITS)
		
		
		#foundGaiaStar=gaiaDo.getByLBRange([210,219.5],[-5,5]) #better get the lb range according to bgFITS

		foundGaiaStar=gaiaDo.getByLBRange([214,218 ],[-4,1])


		##assign each star with a CO intensity
		bgData,bgHead=self.convert3DTo2D(bgData,bgHead)
		
		
		gaiaTBwithCO=self.getOnsourceGaia(foundGaiaStar,bgData, bgHead)
		
		
		gaiaTBwithCO.add_index(self.coint)
		
		
		onSourceGaia=gaiaTBwithCO.loc[self.coint,noiseLevel*2:]
		offSourceGaia=gaiaTBwithCO.loc[self.coint,:noiseLevel*0.5]
 

		sortDistanceON,sortAgON,a,b=self.smoothAg(1/onSourceGaia[colParallax]*1000, onSourceGaia[colAg],smoothScale=smoothInterval)
		sortDistanceOFF,sortAgOFF,a,b=self.smoothAg(1/offSourceGaia[colParallax]*1000, offSourceGaia[colAg],smoothScale=smoothInterval)

		
		
		
		if 1: #draw  Gaia sources
		
			f = plt.figure(figsize=(16.5,7))
			ax = f.add_subplot(121)
			ax.scatter(sortDistanceON,sortAgON,lw=0.3,facecolors='b',s=12, edgecolors='b',label="On-cloud stars" )
			ax.scatter(sortDistanceOFF,sortAgOFF,lw=0.3,facecolors='g',s=12, edgecolors='b',label="Off-cloud stars" )




			plt.savefig(  'TEST_extinctionGaiaAg.png', bbox_inches="tight",dpi=300)
			
			
			#ax.plot([sortDistance,sortAg],[0,max(sortCOExtention)],lw=0.8,color="black")



		
		print "Is this done?"
		
		
	def reduceTo3D(self,data,head):
		
		if len(data)==4:
			data=data[0] 

		del head["CTYPE4"]
		del head["CDELT4"]
		del head["CRPIX4"]
		del head["CROTA4"]
		
		return data,head
		
		
		
		
	def getPPVFits(self,fitsFile):
		"""
		"""
		
		ppvFITS13 ="FITSOFPVExtractor13.fits"
		data,head=self.doFITS.readFITS(fitsFile)
		
		# 
		wcs=WCS(head)
		
		Nz,Ny,Nx=data.shape
		beginP=[0, (Ny-1.)/2.]
		
		endP= [(Nx-1.) , (Ny-1.)/2.]
		#get pv diagrame 
		widthPix=Ny 
		from pvextractor import extract_pv_slice,Path
	
		endpoints = [beginP,endP]
		xy = Path(endpoints,width= widthPix )
	
		pv = extract_pv_slice(  self.CO13HDU, xy)
		
	
		
		pv.writeto(ppvFITS13)



	def drawPV(self,fitsFILE,beginP,endP,widthPix=5,saveName="G214PV.fits"):
		
		pvHDU=fits.open(fitsFILE)[0] 
		
		from pvextractor import extract_pv_slice,Path
		
		endpoints = [beginP,endP]
		xy = Path(endpoints,width= widthPix )
		
		pv = extract_pv_slice(  pvHDU, xy)
 
		os.system("rm "+saveName)
		pv.writeto(saveName)
		
	def drawPVL(self,fitsFile ):
		
		ppvFITS="myOwnPVDuchamp.fits"  
		if 1: #draw pv diagram
			fig = plt.figure(figsize=(8,3.5))
		
		
			rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
			#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
		
			rc('text', usetex=True)
			#readfits
			pvData,pvHead=myFITS.readFITS(ppvFITS) 
		
			grid = ImageGrid(fig, (1, 1, 1), nrows_ncols = (1, 1),
				cbar_mode="each", cbar_pad="0.5%",cbar_size="2%",
				cbar_location="right",
				axes_class=(pywcsgrid2.Axes, dict(header=WCS(pvHead))))
			
			
			axPV = grid[0]
		
			
			#fig,ax= plt.subplots()
			
			#axPV=pywcsgrid2.subplot(211, header=WCS(pvHead))
		
			imCO12=axPV.imshow(pvData, origin="lower",cmap="jet" , norm=LogNorm(vmin=2, vmax=500),    interpolation='none')  #draw color bar
			
			#axPV.set_facecolor('black')
			axPV.set_facecolor('silver')
	
			cb_axes = grid.cbar_axes[0] 
			cb_axes.colorbar(imCO12)
			cb_axes.set_ylabel("K Deg")
			#cb_axes.axis["right"].toggle(ticklabels=False)
			#fig.colorbar(imCO12, cax=cax, orientation='vertical')
			#print dir(axPV)
			#axPV.colorbar()
			
			#cb_axes.set_ticklabel1_type("manual",locs=[1,2,4,8,16,32,64,128] )
	
		
			
		
			if 1: #velocity
				vRange=[-30,70]
				vInterval=10 #km/s
				
				yLocs=np.arange(int(vRange[0]) , int(vRange[1])+vInterval,vInterval)
				
				yLabels= map(int, yLocs)
				yLabels=map(str,yLabels)
		 
				axPV.set_ticklabel2_type("manual",locs=yLocs*1000.,labels=yLabels)
		
				axPV.set_ylabel(r"Velocity ($\rm km \ s^{-1}$)")
		
			if 1: #longtitude
				lRange=[219,201]
				lInterval=-1 #degree
				
				xLocs=  np.arange( (lRange[0]) ,  (lRange[1])+lInterval,lInterval)    #np.arange( (lRange[0]) ,  (lRange[1])+lInterval,lInterval)
				
				#xLabels= map(int, xLocs)
				xLabels=map(str,xLocs)
		
				locsX=xLocs-360 #abs(xLocs-219.5)
		
				#print locsX
		
				axPV.set_ticklabel1_type("manual",locs=locsX,labels=xLabels)
		
				axPV.set_xlabel(r"$l$($^{\circ}$)")
		
		
			#replace the ticks
			#axPV.axis[:].major_ticks.set_color("w")
		
		
			fig.tight_layout()
			wcsPV=WCS(pvHead)
			
			if 0:
	
				aa,v0 =wcsPV.wcs_world2pix(0,self.V_OutArm*1000 ,0)
			
				Ny,Nx=pvData.shape
				axPV.plot([0,Nx-1],[v0,v0],'b--',lw=0.8,alpha=0.8)
			
			
				textOffsetX=15
				textOffsetY=20
			
				axPV.text(textOffsetX,v0+textOffsetY,"the Outer Arm",color='blue',alpha=0.8,)
			
				aa,v1 =wcsPV.wcs_world2pix(0,self.V_PerseusArm*1000 ,0)
			
			
				axPV.plot([0,Nx-1],[v1,v1],'g--',lw=0.8,alpha=0.8)
			
		
				axPV.text(textOffsetX,v1+textOffsetY,"the Perseus Arm",color='green',alpha=0.8,)
			
				axPV.text(textOffsetX,textOffsetY,"the Local Arm",color='red',alpha=0.8,)
			
			if 1:
				
				aa,v0 =wcsPV.wcs_world2pix(0,self.V_OutArm*1000 ,0)
				Ny,Nx=pvData.shape

				
				aa,v1 =wcsPV.wcs_world2pix(0,self.V_PerseusArm*1000 ,0)
			
			
				axPV.plot([0,Nx-1],[v1,v1],'g--',lw=0.8,alpha=0.8)
				
				
				
				
				
				#outArm Line

				x0,v0 =wcsPV.wcs_world2pix(219.75-360,55*1000 ,0)
				x1,v1 =wcsPV.wcs_world2pix(200-360,20*1000 ,0)

				#Ny,Nx=pvData.shape
				axPV.plot([x0,x1],[v0,v1],'--',lw=1,  color='black')
 
				axPV.set_xlim(x0,x1)
	
		
		
		
		
			at = AnchoredText(r"$^{12}\mathrm{CO}~(J=1\rightarrow0)$", loc=1, frameon=False,prop={"color":"black"} )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			axPV.add_artist(at)
		
		
		
		
			
		
		
		
		
			fig.savefig('PPVG200220.pdf', bbox_inches="tight",dpi=300)
	
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			return
		#draw pv diagram along the Galactic latitude,
		if 0: #
			data,head=self.doFITS.readFITS(fitsFile)
			
			# 
			wcs=WCS(head)
			
			Nz,Ny,Nx=data.shape
			
			#average the spectrum 
			
			
			#d#ata=data[:,Ny/4:Ny*3/4,:]
			
			PVL= np.sum(data,axis=1)/1./Ny
			
			print PVL.shape
			
			fig,ax= plt.subplots()
	 		
			
			ax.imshow(PVL ,origin="lower", interpolation="none" )
			
			
			#fig.savefig('pvtestMy.png', bbox_inches="tight",dpi=300)
			
			fig.savefig('pvtestMy.pdf', bbox_inches="tight",dpi=300)
	

		ppvFITS="FITSOFPVExtractor.fits"
		ppvFITS13 ="FITSOFPVExtractor13.fits"

		 
		if 0: #use ppv 12CO
			data,head=self.doFITS.readFITS(fitsFile)
			
			# 
			wcs=WCS(head)
			
			Nz,Ny,Nx=data.shape
			beginP=[0, (Ny-1.)/2.]
			
			endP= [(Nx-1.) , (Ny-1.)/2.]
			#get pv diagrame 
			widthPix=Ny 
			from pvextractor import extract_pv_slice,Path
	
			endpoints = [beginP,endP]
			xy = Path(endpoints,width= widthPix )
	
			pv = extract_pv_slice(  self.CO12HDU, xy)
			

			
			pv.writeto(ppvFITS)
			
			
		 
		if 0: #use ppv 13CO

			data,head=self.doFITS.readFITS(fitsFile)
			
			# 
			wcs=WCS(head)
			
			Nz,Ny,Nx=data.shape
			beginP=[0, (Ny-1.)/2.]
			
			endP= [(Nx-1.) , (Ny-1.)/2.]
			#get pv diagrame 
			widthPix=Ny 
			from pvextractor import extract_pv_slice,Path
	
			endpoints = [beginP,endP]
			xy = Path(endpoints,width= widthPix )
	
			pv = extract_pv_slice(  self.CO12HDU, xy)
			

			
			pv.writeto(ppvFITS13)
			
			
			
			
			

		if 1: #draw pv diagram
			fig = plt.figure(figsize=(8,3.5))
 

			rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
			#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	 
			rc('text', usetex=True)
			#readfits
			pvData,pvHead=myFITS.readFITS(ppvFITS) 

			grid = ImageGrid(fig, (1, 1, 1), nrows_ncols = (1, 1),
				cbar_mode="each", cbar_pad="0.5%",cbar_size="2%",
				cbar_location="right",
				axes_class=(pywcsgrid2.Axes, dict(header=WCS(pvHead))))
			
			
			axPV = grid[0]
 
 
			#fig,ax= plt.subplots()
	 		
			#axPV=pywcsgrid2.subplot(211, header=WCS(pvHead))

			imCO12=axPV.imshow(pvData, origin="lower",cmap="bone" ,vmin=0.02 ,vmax=0.25,   interpolation='none')  #draw color bar
			cb_axes = grid.cbar_axes[0] 
			cb_axes.colorbar(imCO12)
			cb_axes.set_ylabel("K Deg")
			#cb_axes.axis["right"].toggle(ticklabels=False)
			#fig.colorbar(imCO12, cax=cax, orientation='vertical')
			#print dir(axPV)
			#axPV.colorbar()
			
 
			

			if 1: #velocity
				vRange=[0,70]
				vInterval=10 #km/s
				
				yLocs=np.arange(int(vRange[0]) , int(vRange[1])+vInterval,vInterval)
				
				yLabels= map(int, yLocs)
				yLabels=map(str,yLabels)
		 
				axPV.set_ticklabel2_type("manual",locs=yLocs*1000.,labels=yLabels)
		
				axPV.set_ylabel(r"Velocity ($\rm km \ s^{-1}$)")
	
			if 1: #velocity
				lRange=[219,209]
				lInterval=-1 #degree
				
				xLocs=np.arange( (lRange[0]) ,  (lRange[1])+lInterval,lInterval)
				
				#xLabels= map(int, xLocs)
				xLabels=map(str,xLocs)

				locsX=abs(xLocs-219.5)

				#print locsX

				axPV.set_ticklabel1_type("manual",locs=locsX,labels=xLabels)
		
				axPV.set_xlabel(r"$l$($^{\circ}$)")


			#replace the ticks
			axPV.axis[:].major_ticks.set_color("w")
 
 
			fig.tight_layout()
 
			wcsPV=WCS(pvHead)
			
	 
			aa,v0 =wcsPV.wcs_world2pix(0,self.V_OutArm*1000 ,0)

			Ny,Nx=pvData.shape
			axPV.plot([0,Nx-1],[v0,v0],'b--',lw=0.8,alpha=0.8)


			textOffsetX=15
			textOffsetY=20

			axPV.text(textOffsetX,v0+textOffsetY,"the Perseus Arm",color='blue',alpha=0.8,)
 
			aa,v1 =wcsPV.wcs_world2pix(0,self.V_PerseusArm*1000 ,0)

 
			axPV.plot([0,Nx-1],[v1,v1],'g--',lw=0.8,alpha=0.8)


			axPV.text(textOffsetX,v1+textOffsetY,"the interarm region",color='green',alpha=0.8,)

			axPV.text(textOffsetX,textOffsetY,"the Local Arm",color='red',alpha=0.8,)
 

			at = AnchoredText(r"$^{12}\mathrm{CO}~(J=1\rightarrow0)$", loc=1, frameon=False,prop={"color":"w"} )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			axPV.add_artist(at)



 
			




			fig.savefig('PPVG210.pdf', bbox_inches="tight",dpi=300)




	def drawRegions(self):
		"""
		This figure is used to review the region of maddalena, and draw the name of sub regionss
		"""
		#use 12CO
		backGroundFITS="mosaic_U_M0.fits"
		Bdata,Bhead=doMadd.doFITS.readFITS( backGroundFITS)
		
		if len(Bdata.shape)==3:
			Bdata=Bdata[0]
		
		
		
		
		fig = plt.figure()
		
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		ax=pywcsgrid2.subplot(111,header=WCS(Bhead))



		ax.imshow(Bdata,origin='lower',cmap="bone",vmin=0,vmax=20,interpolation='none')

 


		ax.axis[:].major_ticks.set_color("w")
		plt.tight_layout(pad=0)
		fig.savefig( 'subregions.eps', bbox_inches="tight")
		#fig.savefig( 'YSOMadda.eps', bbox_inches="tight")


	def drawCompareCOFUGINDame(self):
		#compare our survey with FUGIN and Dame
		
		"""
		"""
		#

		ourFITS="mosaic_U_M0.fits"
		Bdata,Bhead=doMadd.doFITS.readFITS( ourFITS)
		
		if len(Bdata.shape)==3:
			Bdata=Bdata[0]

 
		FUGINdata,FUGINhead=doMadd.doFITS.readFITS( self.FUGINCO12)
		
		if len(FUGINdata.shape)==3:
			FUGINdata=FUGINdata[0]

		DAMEdata,DAMEhead=doMadd.doFITS.readFITS( self.DAMECO12)
		
		if len(DAMEdata.shape)==3:
			DAMEdata=DAMEdata[0]



		x = np.arange(0, 10, 0.2)
		y = np.sin(x)
		
		# plot it
		fig = plt.figure(figsize=(8, 6)) 
		
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1]) 
		#ax0 = plt.subplot(gs[0])
		
		
		axMy=pywcsgrid2.subplot(gs[0,0],header=WCS(Bhead))  
		axMy.imshow(Bdata,origin="lower",cmap="bone",vmin=3,vmax=30,interpolation='none')
		axMy.set_ticklabel_type("absdeg", "absdeg")
		axMy.axis[:].major_ticks.set_color("w")

		at = AnchoredText(r"MWISP", loc=2, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axMy.add_artist(at)

		
		
		axDAME=pywcsgrid2.subplot(gs[0,1],header=WCS(DAMEhead))  
		axDAME.imshow(DAMEdata,origin="lower",cmap="bone",vmin=1,vmax=10,interpolation='none')
		axDAME.set_ticklabel_type("absdeg", "absdeg")
		axDAME.axis[:].major_ticks.set_color("w")

		at = AnchoredText(r"Dame et al. (2001)", loc=2, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axDAME.add_artist(at)
		
		
		
		
		axFUGIN = pywcsgrid2.subplot(gs[1,:],header=WCS(FUGINhead))  
		axFUGIN.imshow(FUGINdata,origin="lower",cmap="bone",vmin=3,vmax=40,interpolation='none')
		axFUGIN.set_ticklabel_type("absdeg", "absdeg")
		#ax0.plot(x, y)
		axFUGIN.axis[:].major_ticks.set_color("w")
		axFUGIN.set_xlim( 215,4248)
		
		
		at = AnchoredText(r"Umemoto et al. (2017)", loc=1, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axFUGIN.add_artist(at)

		
 
		
		plt.tight_layout()
		plt.savefig('compareCOSurvey.pdf')
		plt.savefig('compareCOSurvey.png',dpi=600)





	#fitting function
	def calPowerLawProbLogcomponent1(self,theta, dataArraym ,minV,maxV  ):
		
		
		alpha1  =theta[0] #Mt is the turnover mass of molecular cores
		
 
 		
		beta1= alpha1-1 #1-alpha1
 
		if maxV==None:
 		#normalFactor1=beta1/( 0 -minV**beta1 )
			normalFactor1= beta1*minV**beta1   #beta1/(  -minV**beta1 )

		else:
			#aaaaaaaaaaaaaaaaaaaaa
			#normalFactor1=-1./beta1/( minV**(-beta1)-maxV**(-beta1)  )   #( alpha1 -1 )*( minV**(  alpha1 -1 ) ) #beta1/( maxV**beta1 -minV**beta1 )
			normalFactor1=  (1-alpha1)/(  maxV**(1-alpha1) - minV**(1-alpha1)   )

			
		return len(dataArraym)*np.log(normalFactor1)-alpha1*np.sum(np.log( dataArraym ))

		#print np.log(normalFactor1*dataArraym**(-alpha1))


		#return len(dataArraym)*( beta1*np.log(minV)- np.log(beta1))-alpha1*np.sum(np.log( dataArraym ))
 

	def fitPowerLawWithMCMCcomponent1(self,dataArray,sampleN=1000,burn_in=100,minV=None,maxV=None,thin=15):
		"""
		fit a power law with MCMC
		"""
		#if minV==None or maxV==None:
			#minV= min(dataArray)
			#maxV= max(dataArray)
		
		mass,massN=self.madeCount(dataArray)

		#print np.min(dataArray ),   max(dataArray)
		
		print "The minimum and maximum area are (in solar deg), in the MWISP cloud folder",minV,maxV
		
		
		np.random.seed()

		alpha1= np.random.uniform(1,4) #np.random.exponential(1)
 
		theta=[ alpha1 ]
		

		
		p0= self.calPowerLawProbLogcomponent1(theta,mass,minV,maxV)

		sampleK=[]
		sampleAlpha=[ alpha1 ]
		widgets = ['MCMCSmapleSlope: ', Percentage(), ' ', Bar(marker='>',left='|',right='|'),
		           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
		
		pbar = ProgressBar(widgets=widgets, maxval=sampleN+burn_in+1)
		pbar.start() 


 		recordSample=[]
		
		for i in range(100000):

			newAlpha1=  sampleAlpha[-1] +np.random.normal(0, 0.5 ) #np.random.exponential(1)

			if newAlpha1<1:
				continue


			theta=[ newAlpha1  ]
			
			p1= self.calPowerLawProbLogcomponent1(theta,mass,minV,maxV)

			if  np.isinf(p1) :
				continue
				
				
			randomR=np.random.uniform(0,1)



			if p1>=p0 or p1-p0>np.log(randomR):
				p0=p1; 

				sampleAlpha.append( newAlpha1 )
 

			else:
				sampleAlpha.append(  sampleAlpha[-1] )



			if i%thin==0:
				recordSample.append(sampleAlpha[-1]  )


			pbar.update(len(recordSample)) #this adds a little symbol at each iteration

			if len(recordSample)>sampleN+burn_in:
				break
		pbar.finish()
		#print mean( sampleAlpha[burn_in:] ), np.std(sampleAlpha[burn_in:]  )
		return np.array(  recordSample[burn_in:]  )

		#sample theta




	#fitting function
	def calPowerLawProbLog3p(self,theta, dataArraym ,minV,maxV  ):
		
		#maxV=maxV*1000
		
		alpha1,alpha2,Mt =theta #Mt is the turnover mass of molecular cores
		
		
		#should normallize them together
		
		part1=dataArraym[dataArraym<Mt]
		part2=dataArraym[dataArraym>=Mt]
		
		beta1=1-alpha1
		beta2=1-alpha2

 		normalFactor1=beta1/(  Mt**beta1-minV**beta1 )
 		#normalFactor2=beta2/(  maxV**beta2-Mt**beta2 )

 		normalFactor2=beta2/(  0-Mt**beta2 )


		return len(part1)*np.log(normalFactor1)-alpha1*np.sum(np.log( part1 ))  +len(part2)*np.log(normalFactor2)-alpha2*np.sum(np.log( part2 ))



	#fitting function
	def calPowerLawProbLog_lognormal(self,theta, dataArraym ,minV,maxV  ):
		
		
		isigma, mu, Mt, alpha =theta #Mt is the turnover mass of molecular cores
		
		
		
		

		
		#should normallize them together
		
		part1=dataArraym[dataArraym<Mt] #gaussian part 
		

		
		part2=dataArraym[dataArraym>=Mt] # 
		
		
		
 
		beta2=1-alpha
 		normalFactor2= 1./beta2*(  np.exp(beta2*maxV) -   np.exp(beta2*Mt)  ) #beta2/(  maxV**beta2-Mt**beta2 )

 		#
 		
 		#better use trucated gaussian
 		nominator=np.sqrt(2)
 		
		a= 1 #0.5*special.erf((Mt-mu)/nominator*isigma)+ special.erf((mu-minV)/nominator*isigma)  

 		prob1=     len(part1)*(-np.log(a)+ np.log(isigma)-0.5*np.log(2*np.pi)  ) -np.sum( 0.5*isigma**2*(part1-mu)**2)

 		
		if 0:#only gauss
			vCut=5.3
			dataArraym=dataArraym[dataArraym<vCut] #gaussian part 

			a= 1 #0.5*special.erf((vCut-mu)/nominator*isigma)+ special.erf((mu-minV)/nominator*isigma)  
	
			return  len(dataArraym)*(-np.log(a)+ np.log(isigma)-0.5*np.log(2*np.pi)  ) -np.sum( 0.5*isigma**2*(dataArraym-mu)**2)

		return    prob1 -len(part2)*np.log(normalFactor2) + np.sum( (1-alpha)*part2   )  #alpha*np.sum(np.log( part2 ))



 
	def getNextValues(self,runSamples,RMSs,paraIndex):
		
		dim=len(RMSs)

		currentValues=[]
		for i in range(dim):
			
			currentValues.append( runSamples[i][-1]  )
  
		currentValues[paraIndex]=currentValues[paraIndex]+np.random.normal(0,RMSs[paraIndex])
		
		return currentValues



	def fitPowerLawWithMCMC_lognormal(self,dataArray,sampleN=1000,burn_in=1000,thinning=10):
		"""
		fit a power law and a log normal with MCMC 
		
		#four parameters
		
		mu,sigma, alpha,mt
		
		
		#fitting a log normal 
		
		"""
		dataArray=np.log(dataArray) #convert to log normal 

		np.random.seed()
		
		minV= min(dataArray)
		maxV= max(dataArray)
		mass,massN=self.madeCount(dataArray) #sor the data,massN is useless

		
 
		isigma= np.random.exponential(2.5)
		
		Mt=  np.random.uniform(4,maxV)
		
		mu=np.random.uniform( minV ,Mt)



		alpha=np.random.exponential(3 )

		RMSs=[ 2., 2., 2., 2. ]

		
		theta0=[ isigma, mu, Mt, alpha] #four parameters

		p0= self.calPowerLawProbLog_lognormal(theta0,mass,minV,maxV)


		runSamples=[ [isigma], [mu], [Mt], [alpha] ]


		isigmaList=[]
		muList=[]
		MtList=[]
		alphaList=[]
		
		
		
		acceptN=0

		widgets = ['MCMC_Smaple_LogNormal: ', Percentage(), ' ', Bar(marker='>',left='|',right='|'),
		           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
		
		pbar = ProgressBar(widgets=widgets, maxval=sampleN+burn_in+1)
		pbar.start() 
		
		sampleIndex=0
		
		for i in range(20000):

 

			paraJ=0
			while paraJ<4:
 
				valueCand= self.getNextValues(runSamples,RMSs, paraJ) 
				
				
 
				isigma, mu, Mt, alpha =valueCand 
				
				
				
				

 
				#print  isigma, mu, Mt, alpha 

 
				if   Mt< mu  or Mt>maxV or mu< minV or isigma<0 :
					continue

				p1= self.calPowerLawProbLog_lognormal(valueCand,mass,minV,maxV)


				
				if  np.isinf(p1) or  np.isnan(p1) :
					continue 
				
				
 
				

				randomR=np.random.uniform(0,1)
				
				
				if p1>=p0 or p1-p0>np.log(randomR):
					#print disCloud,  imu1sigma,mu2,imu2sigma,"--->",processID
					
					
					
					acceptN=acceptN+1
					p0=p1;
					runSamples[paraJ].append( valueCand[paraJ] )
				
				else:
					runSamples[paraJ].append( runSamples[paraJ][-1] )
				paraJ=paraJ+1
			 
				runSamples=[ runSamples[0][-1:],runSamples[1][-1:], runSamples[2][-1:]  , runSamples[3][-1:] ]
				
	 
		
			if i%thinning==0:


				isigmaList.append( runSamples[0][-1] )
				muList.append(runSamples[1][-1] )
				MtList.append(runSamples[2][-1] )


				alphaList.append(runSamples[3][-1] )

				
				
				
				#alpha2List.append(runSamples[1][-1]  )
				#MtList.append(runSamples[2][-1] )
				#print runSamples
				
				print runSamples

			if len(alphaList)>burn_in+sampleN:
			
				break
			pbar.update(len(alphaList)) #this adds a little symbol at each iteration

					 
		pbar.finish() 

		return  np.array( [ isigmaList[burn_in+1:], muList[burn_in+1:], MtList[burn_in+1:] , alphaList[burn_in+1:] ] ) #runSamples # np.array(  sampleAlpha[burn_in:]  )


 







	def fitPowerLawWithMCMC3p(self,dataArray,sampleN=1000,burn_in=100,thinning=15):
		"""
		fit a power law with MCMC 
		
		#three parameters
		
		alpha1, alpha2, mass turnover
		
		"""
		minV= min(dataArray)
		maxV= max(dataArray)
		
		mass,massN=self.madeCount(dataArray)
 
		
		
		print "Moddling 3 parameters."
		
		#print "The minimum and maximum masses are (in solar mass)",minV,maxV
		
		
		np.random.seed()

		alpha1= np.random.normal(-1,1) #np.random.exponential(0 .1) 
		
		alpha2=  np.random.normal(3,1)  #np.random.exponential(3)

		Mt=np.random.uniform(minV,100)
		#Mt=np.random.exponential(500)

		if maxV>1000.:
			Mt=np.random.uniform(minV,300)

		print Mt,minV,maxV
		
		#Mt=  np.random.uni  #np.random.exponential(50)

		
		#print Mt
		theta0=[ alpha1,alpha2,Mt]
		
 
		
		p0= self.calPowerLawProbLog3p(theta0,mass,minV,maxV)

		sampleK=[]
		sampleAlpha=[]
		
		
		runSamples=[ [alpha1], [alpha2], [Mt] ]
		
 		
		alpha1List=[]
		alpha2List=[]
		MtList=[]
		
		RMSs=np.array( [0.5 ,0.5 , 10.] )
		
		
		acceptN=0

		widgets = ['MCMC_Smaple_3P: ', Percentage(), ' ', Bar(marker='>',left='|',right='|'),
		           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
		
		pbar = ProgressBar(widgets=widgets, maxval=sampleN+burn_in+1)
		pbar.start() 
		
		sampleIndex=0
		
		for i in range(200000):
			
 
			paraJ=0
			while paraJ<3:
 
				valueCand= self.getNextValues(runSamples,RMSs, paraJ) 
			
				alpha1, alpha2, Mt =valueCand 

				theta=[ alpha1,alpha2,Mt  ]


				if   Mt<minV or Mt>maxV:
					
					#print "What the hell?"
 
					continue

				p1= self.calPowerLawProbLog3p(theta,mass,minV,maxV)
 
				if   np.isinf(p1)  :
					
					#print "What the hell?"
 
					continue
				
			
				
				randomR=np.random.uniform(0,1)
				
				
				if p1>=p0 or p1-p0>np.log(randomR):
					#print disCloud,  imu1sigma,mu2,imu2sigma,"--->",processID
					acceptN=acceptN+1
					p0=p1;
					runSamples[paraJ].append( valueCand[paraJ] )
				
				else:
					runSamples[paraJ].append( runSamples[paraJ][-1] )
				paraJ=paraJ+1 #next parameter
				
			 
				runSamples=[ runSamples[0][-1:],runSamples[1][-1:], runSamples[2][-1:]]
				
	 
		
			if i%thinning==0:

				alpha1List.append(runSamples[0][-1] )
				alpha2List.append(runSamples[1][-1]  )
				MtList.append(runSamples[2][-1] )
				#print runSamples
				#print alpha1List[-1], alpha2List[-1]
			if len(alpha1List)>burn_in+sampleN:
			
				break
			pbar.update(len(alpha1List)) #this adds a little symbol at each iteration

					 
		pbar.finish() 


		#print len(MtList[burn_in+1:])

		return  np.array( [ alpha1List[burn_in+1:], alpha2List[burn_in+1:], MtList[burn_in+1:] ] ) #runSamples # np.array(  sampleAlpha[burn_in:]  )

 

 

	def madeCount(self,dataArray):
		
		mass=np.sort(dataArray)

		massN= len(mass)/1.- np.arange(len(mass))

		return mass,massN
 
	def clearCores(self,coreTB):
		
		#some cores are fake
		#['PIDENT', 'Peak1', 'Peak2', 'Peak3', 'Cen1', 'Cen2', 'Cen3', 'Size1', 'Size2', 'Size3', 'Sum', 'Peak', 'Volume']

 		emptyCoreTB=coreTB[0:1]

		emptyCoreTB.remove_row(0)
		
		for eachRow in coreTB:
			
			if eachRow["Peak3"]<=21:
				continue
				
			#if eachRow["Peak3"]<37 and (eachRow["Peak1"]>562.10765 and eachRow["Peak2"]<145.84456):
				#print eachRow["Peak1"], eachRow["Peak2"], eachRow["Peak3"]
			if eachRow["Peak3"]<37 and (eachRow["Peak1"]<562.10765 or eachRow["Peak2"]>145.84456):
				continue
				#print eachRow["Peak1"], eachRow["Peak2"], eachRow["Peak3"]
			
			if 1:# the peak vale to 6 sigma
				if eachRow["Peak"]<3:
					continue
			
			
			emptyCoreTB.add_row(eachRow)
		return emptyCoreTB


	def getTexFromT12(self,T12):
		
		return 5.53/np.log( 1.+ 5.53/(T12+0.819)  )



	def getMassForCore(self,coreRow):



		
		W345=False
		
		
		Cen1="Cen1"
		Cen2="Cen2"
		Cen3="Cen3"
		
		
		Size1="Size1"
		Size2="Size2"
		Size3="Size3"
		
		Sum="Sum"
		
		colnames=['PIDENT',
							'Peak1',
							'Peak2',
							'Peak3',
							'Cen1',
							'Cen2',
							'Cen3',
							'Size1',
							'Size2',
							'Size3',
							'Sum',
							'Peak',
							'Volume']
		
		
		
		if coreRow["Peak1"]< 160:
			W345=True
		
		
		degSize1=coreRow[Size1]/60./60.
		
		degSize2=coreRow[Size2]/60./60.
		
		pixSize=30./3600.
		

		#looks like useless
		core_lRange=[ coreRow[ Cen1]-degSize1 ,    coreRow[Cen1] + degSize1   ]
		
		core_bRange=[ coreRow[ Cen2]-  degSize2,    coreRow[Cen2] +  degSize2  ]

		core_vRange=[ coreRow[ Cen3]- coreRow[Size3] ,    coreRow[Cen3] + coreRow[Size3]    ]

		
		core_vRange[0]=core_vRange[0]/1000.
		core_vRange[1]=core_vRange[1]/1000.

		#crop 12COFITS, then average  
		
		coreL= coreRow[ Cen1]
		coreB= coreRow[ Cen2]

		LVB=[ coreRow[ Cen1], coreRow[ Cen2],coreRow[ Cen3]/1000. ]
		#print W345,"????????????"
		# better use peak values
		if not W345:
			T12=myFITS.getVoxValue(self.CO12DATA,self.WCSCO12,LVB ) 
		else:
			
			
			T12=myFITS.getVoxValue(self.W345CO12DATA,self.WCSW345CO12,LVB ) 

 
		
		Tex=self.getTexFromT12(T12)
		#this is wrong
		#ex= np.mean(coreData)
		#TexSTD=np.std(coreData,ddof=1)


 
		S=np.pi*degSize1*degSize2
		
		PartTex=1./(  1-np.exp(-5.29/Tex)  )
		
		
		
		NH2=1.49e20*PartTex*coreRow[Sum]*self.vResolution13 #0.158737644553 # cm-2

		#distance=2200 # pc
		
		#parsecToMeter= 3.0857e16 #m 
		#length1=np.radians(degSize1)*self.distance*self.parsecToMeter*100. #cm 
		#length2=np.radians(degSize2)*self.distance*self.parsecToMeter*100.


		if not 			W345: #G216
			
			
			G216D=self.getDisG216(coreB,coreL )
			#print G216D, coreB,coreL
 			length1=np.radians(30./60./60.)*G216D*self.parsecToMeter*100. #cm 
 			
 			
		else: #W345
			
			W345D=self.getDisW345(coreL )
			
 			length1=np.radians(30./60./60.)*W345D*self.parsecToMeter*100. #cm 

 
 
		mu=1.36
		
		Mh2=3.35e-27 #kg
		solarMass=1.9891e30 #kg
		#s=np.pi*length1*length2
		s= length1*length1

		coreSolar= s*NH2*Mh2*mu/solarMass

		#calculate virial mass
		Mvirial=0
		G=6.673e-11 #N m2/kg2
		
		k2=210. #use uniform 
		k1=5./3.
		if not 			W345: #G216
			lengthA=2.35*np.tan(np.radians(degSize1))*self.distance #*self.parsecToMeter*100. #cm 
			lengthB=2.35*np.tan(np.radians(degSize2))*self.distance #*self.parsecToMeter*100.

			sigmaV=coreRow[Size3]/1000.
			
 			R=np.sqrt( lengthA*lengthB)
 			
 			Mvirial=k2*R*(2.35*sigmaV)**2
 			
 			Mv1=k1*3*coreRow[Size3]**2*R*self.parsecToMeter/G/solarMass

			alphaVir=5*coreRow[Size3]**2*R*self.parsecToMeter/G/solarMass/coreSolar
 			
		else: #w345
			W345D=self.getDisW345(coreL )
 

			lengthA=2.35*np.tan(np.radians(degSize1))*W345D #self.W345distance #*self.parsecToMeter*100. #cm 
			lengthB=2.35*np.tan(np.radians(degSize2))* W345D #self.W345distance #*self.parsecToMeter*100.
			sigmaV=coreRow[Size3]/1000.
 
 			R=np.sqrt( lengthA*lengthB)
 			
 			Mvirial=k2*R*(2.35*sigmaV)**2
 			
 			Mv1=k1*3*coreRow[Size3]**2*R*self.parsecToMeter/G/solarMass
			alphaVir=5*coreRow[Size3]**2*R*self.parsecToMeter/G/solarMass/coreSolar
			#print alphaVir, Mvirial/coreSolar
			 
 		#print Mv1, Mvirial	
		
		#print  coreSolar/Mvirial

		return coreSolar,Tex, [Mvirial, alphaVir ]



	def getDisG216(self,coreB,coreL):
		"""
		pass
		"""

		#W5distance=2300. #pc
		#W3distance= 1950. # pc
		#W345SplitDL= 135.6  # distance split  


		if  coreB >self.G214B and coreL<self.G214L:
			return self.G214Distance
			
		return self.distance


	def getDisW345(self,coreL):
		"""
		pass
		"""

		#W5distance=2300. #pc
		#W3distance= 1950. # pc
		#W345SplitDL= 135.6  # distance split  


		if coreL>self.W345SplitDL:
			return self.W5distance
		return self.W3distance
		

	def getMassForCoreWithDistance(self,coreRow,distance ):

		#print distance,"???????????????????????"


		
		W345=False
		
		
		Cen1="Cen1"
		Cen2="Cen2"
		Cen3="Cen3"
		
		
		Size1="Size1"
		Size2="Size2"
		Size3="Size3"
		
		Sum="Sum"
		
		colnames=['PIDENT',
							'Peak1',
							'Peak2',
							'Peak3',
							'Cen1',
							'Cen2',
							'Cen3',
							'Size1',
							'Size2',
							'Size3',
							'Sum',
							'Peak',
							'Volume']
		
		
		
		if coreRow["Peak1"]< 160:
			W345=True
		
		
		degSize1=coreRow[Size1]/60./60.
		
		degSize2=coreRow[Size2]/60./60.
		
		pixSize=30./3600.
		

		#looks like useless
		core_lRange=[ coreRow[ Cen1]-degSize1 ,    coreRow[Cen1] + degSize1   ]
		
		core_bRange=[ coreRow[ Cen2]-  degSize2,    coreRow[Cen2] +  degSize2  ]

		core_vRange=[ coreRow[ Cen3]- coreRow[Size3] ,    coreRow[Cen3] + coreRow[Size3]    ]

		
		core_vRange[0]=core_vRange[0]/1000.
		core_vRange[1]=core_vRange[1]/1000.

		#crop 12COFITS, then average  
		
		coreL=coreRow[ Cen1]
		
		LVB=[ coreRow[ Cen1], coreRow[ Cen2],coreRow[ Cen3]/1000. ]
		#print W345,"????????????"
		# better use peak values
		if not W345:
			T12=myFITS.getVoxValue(self.CO12DATA,self.WCSCO12,LVB ) 
		else:
			T12=myFITS.getVoxValue(self.W345CO12DATA,self.WCSW345CO12,LVB ) 

		
		
		Tex=self.getTexFromT12(T12)
		#this is wrong
		#ex= np.mean(coreData)
		#TexSTD=np.std(coreData,ddof=1)


 
		S=np.pi*degSize1*degSize2
		
		PartTex=1./(  1-np.exp(-5.29/Tex)  )
		
		
		
		NH2=1.49e20*PartTex*coreRow[Sum]*self.vResolution13 #0.158737644553 # cm-2

		#distance=2200 # pc
		
		#parsecToMeter= 3.0857e16 #m 
		#length1=np.radians(degSize1)*self.distance*self.parsecToMeter*100. #cm 
		#length2=np.radians(degSize2)*self.distance*self.parsecToMeter*100.


		if not 			W345:
 			length1=np.radians(30./60./60.)*self.distance*self.parsecToMeter*100. #cm 
 			
 			
		else:
			W345D=self.getDisW345(coreL )
			
 			length1=np.radians(30./60./60.)*W345D*self.parsecToMeter*100. #cm 

 
 
		mu=1.36
		
		Mh2=3.35e-27 #kg
		solarMass=1.9891e30 #kg
		#s=np.pi*length1*length2
		s= length1*length1

		coreSolar= s*NH2*Mh2*mu/solarMass

		#calculate virial mass
		Mvirial=0
		G=6.673e-11 #N m2/kg2
		
		k2=210. #use uniform 
		k1=5./3.
		if not 			W345: #G216
			lengthA=np.tan(np.radians(degSize1))*distance #*self.parsecToMeter*100. #cm 
			lengthB=np.tan(np.radians(degSize2))*distance #*self.parsecToMeter*100.

			sigmaV=coreRow[Size3]/1000.
			
 			R=np.sqrt( lengthA*lengthB)
 			
 			Mvirial=k2*R*(2.35*sigmaV)**2
 			
 			Mv1=k1*3*coreRow[Size3]**2*R*self.parsecToMeter/G/solarMass


 			
 			
		else: #w345
			lengthA=np.tan(np.radians(degSize1))*distance #*self.parsecToMeter*100. #cm 
			lengthB=np.tan(np.radians(degSize2))*distance #*self.parsecToMeter*100.
			sigmaV=coreRow[Size3]/1000.
 
 			R=np.sqrt( lengthA*lengthB)
 			
 			Mvirial=k2*R*(2.35*sigmaV)**2
 			
 			Mv1=k1*3*coreRow[Size3]**2*R*self.parsecToMeter/G/solarMass

 			
 		#print Mv1, Mvirial	
		
		
		return coreSolar,Tex,Mvirial


 

	def getTexListFromTB(self,coreTB):
 	
		TexList=[]
		
		for coreRow in  coreTB:
 
 
 
			
			coreSolar,Tex,Mvirial=self.getMassForCore(coreRow )
			TexList.append( Tex ) #,Tex
			#print coreSolar, coreRow[Sum], Tex
			
		#print "The minimum Tex is ",np.min(TexList)
	#	print "The maximum Tex is ",np.max(TexList)
		#print "The mean Tex is ",np.mean(TexList)

		return np.array(TexList)
		
		
		


	def getMassFromGauSum(self,coreTB,saveName=''):
		"""
		calculate the mass of cores with 12CO with Sum, and 12CO is used to calculate the excitation temperature
		
		
		1. According to the postion of molecular cores, get the average tempearture of 12CO as the excitation tempearture 
		
		2. Then use sum to calculte the mass
		"""
		

		#should return the data row
		print "Calculating masses for {} cores".format(len(coreTB) )
		
		
		#cocmf=myCMF()

		
		sourceName=""
		
		if "G216" in saveName:
			sourceName="G216"
			
		else:
			sourceName="W345"
		
		docmf=myCMF(sourceName)
		
		newRow= docmf.getEmptyRow()
 
		
		#read the fits file 
		#coreData,coHead=myFITS.readFITS(COFITS)

  
		coreN=len(coreTB)
		#distance 
		massList=[]
		MvirialList=[]
		
		alphaVirialList=[]
		
		for coreRow in  coreTB:
  
			coreSolar,Tex,MvirialAndAlpha=self.getMassForCore(coreRow )
			massList.append( coreSolar ) #,Tex
			
			Mvirial,alphaVirial=MvirialAndAlpha
			
			MvirialList.append(Mvirial  )
			
			
			alphaVirialList.append( alphaVirial )
			
			
			
			#print coreSolar, coreRow[Sum], Tex
 
		totalMass= np.sum(  massList)
		print totalMass, "total mass"
		print np.sum( MvirialList ),"Total virial mass"
		if 1: #draw test
			massList=np.array( massList)
			#massList=ICO13
			#low masses seems wrong...
			#massList=massList[massList>1]  
			
			alpha1Arrays=[]
			alpha2Arrays=[]
			MtArrays=[]

			chains=6

			
			for i in range(chains):
			
				sampleArray=self.fitPowerLawWithMCMC3p(massList )
				
				alpha1Arrays.append(sampleArray[0])
				alpha2Arrays.append(sampleArray[1])
				MtArrays.append(sampleArray[2])

			
			alpha1Arrays=np.array(alpha1Arrays)
			alpha2Arrays=np.array(alpha2Arrays)
			MtArrays=np.array(MtArrays)

			sampleN=len(alpha1Arrays)


			
			
			alpha1Arrays=alpha1Arrays.reshape(-1)
			alpha2Arrays=alpha2Arrays.reshape(-1)
			MtArrays=MtArrays.reshape(-1)

			#draw cornermap
			sampleArray=np.array([alpha1Arrays, alpha2Arrays,  MtArrays ]		 )	
			sampleArray=np.transpose(sampleArray)
 
 
			 
 
			#figureCorner = corner.corner(sampleArray, quantiles=[0.16, 0.5, 0.84])

			#figureCorner.savefig(saveName+"corner.png")
			
			#sampleArray1=self.fitPowerLawWithMCMC3p(massList )

			#print len(sampleArray[0])

			alpha1Mean= np.mean(alpha1Arrays)
			alpha2Mean=np.mean(alpha2Arrays)
			MtMean= np.mean(MtArrays)
			
			alpha1_std= np.std(alpha1Arrays,ddof=1)
			alpha2_std=np.std(alpha2Arrays,ddof=1)
			Mt_std= np.std(MtArrays,ddof=1)	
			
			
			turnOverlog=np.log10(MtMean)
			print "===================MCMC mean + std ======================"
			print "alpha1 {} +/- {}".format( alpha1Mean, alpha1_std )
			print "alpha2 {} +/- {}".format( alpha2Mean, alpha2_std )
			print "Mt {} +/- {}".format( MtMean, Mt_std )
			print "========================================="

			
			#print  ""turnOverlog, alpha1Mean, alpha2Mean

			#fig, ax = plt.subplots()
			fig, axs = plt.subplots( figsize=(14,6) )  
			rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
			#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	 
			rc('text', usetex=True)
			#fig = plt.figure()
			nn=1
			for i in range(3):
				plt.subplot2grid((3*nn,6*nn), (nn*i, 0) ,rowspan=1,colspan=1)
				plt.subplot2grid((3*nn,6*nn), (nn*i, 1) ,rowspan=1,colspan=1)
				plt.subplot2grid((3*nn,6*nn), (nn*i, 2) ,rowspan=1,colspan=1)

		
		
		
			figureCorner = corner.corner(sampleArray,fig=fig,  show_titles=True )

			ndim=3
			# Extract the axes
			axes = np.array(figureCorner.axes).reshape((ndim, ndim))
			
			alpha1Mean= np.mean(alpha1Arrays)
			alpha2Mean=np.mean(alpha2Arrays)
			MtMean= np.mean(MtArrays)
			
			alpha1_std= np.std(alpha1Arrays,ddof=1)
			alpha2_std=np.std(alpha2Arrays,ddof=1)
			Mt_std= np.std(MtArrays,ddof=1)	


			if 1:

				meanValues = [alpha1Mean, alpha2Mean,MtMean, ] #np.mean(sampleArray, axis=0)
				lowerPHD90 = [ alpha1_std, alpha2_std,   Mt_std ] #np.mean(sampleArray, axis=0)
				#upperPHD90 = [alpha1Mean+alpha1_std,alpha2Mean+alpha2_std, MtMean+ Mt_std ] #np.mean(sampleArray, axis=0)
				upperPHD90 = [ alpha1_std, alpha2_std,   Mt_std ] #np.mean(sampleArray, axis=0)

				titleAlpha1=r'$\alpha_1 =  {:.2f}\pm{:.2f} $'.format(alpha1Mean, alpha1_std)
				titleAlpha2=r'$\alpha_2 =  {:.2f}\pm{:.2f} $'.format(alpha2Mean, alpha2_std)
				#agmu1Sigma,loweragmu1Sigma,upperagmu1Sigma
				titleMturn=r'M$_{{\rm turn}} =  {:.0f}\pm{:.0f} \ M_{{\odot}}$  '.format(MtMean, Mt_std)
 
				labels=[titleAlpha1,  titleAlpha2,titleMturn ]
				# Loop over the diagonal
				for i in range(ndim):
					axCorner = axes[i, i]
					#pass
					#ax.axvline(meanValues[i], color="black" )
					#axCorner.title("aaa")
	
					axCorner.set_title(labels[i],fontsize=12)
					#if i==0  :
						#axCorner.set_title(labels[i],fontsize=14)
	
					axCorner.axvline(meanValues[i], color="black" ,linestyle='-',linewidth=1.5)
					axCorner.axvline(meanValues[i]-lowerPHD90[i], color="black" ,linestyle='--',linewidth=1.5)
					axCorner.axvline(meanValues[i]+upperPHD90[i], color="black" ,linestyle='--',linewidth=1.5)
		
		
		
			#plt.subplot2grid((8*nn,8*nn), (i*nn*2, 0) ,rowspan=2*nn,colspan=nn)
			#plt.subplot2grid((8*nn,8*nn), (i*nn*2,  nn),rowspan=2*nn,colspan=nn)
			#plt.subplot2grid((8*nn,8*nn), (i*nn*2,  2*nn),rowspan=2*nn,colspan=nn)
			ax=plt.subplot2grid((5,12), (0,7),rowspan=5,colspan=5)

		
		
			rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
			#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
	 
			rc('text', usetex=True)
			logM=np.log10(massList)

			DN,edges=np.histogram(logM,bins=7)
			centerEdge=(edges[1:] + edges[:-1]  )/2.
			
			dlogM=edges[1:]- edges[:-1]
			
			Y=np.log10( DN/dlogM )
			
			X=centerEdge
			

			#draw the Error 
			lowerY=Y-np.log10( (DN-np.sqrt(DN))/dlogM )
			
			upperY=np.log10( (DN+np.sqrt(DN))/dlogM )-Y
			
			#better draw the error bar

		#plt.hist(coreMass,100,density=True)
			#ax.plot(X,Y)
			ax.errorbar(X,Y,yerr=[lowerY,upperY],c='b',marker='o',capsize=1.5,elinewidth=0.8,lw=1,label="{} clumps in total".format(len(coreTB)))
			#ax.plot([turnOverlog,turnOverlog],[min(Y-lowerY)  ,max(Y+upperY)  ]  ,c='black',label="turnover mass")
			ax.axvline(turnOverlog, linestyle='--',color="black", label=r"M$_{{\rm turn}} = {:.0f}\pm{:.0f}  M_{{\odot}}$".format(MtMean,Mt_std),lw=1)
			
			
			
			x0=turnOverlog

			#determine y0;
			
			#frist, find the closest point in the right sise of turnover mass
			
			# make sure the line corss this point,
			
			leftOneX=   X[X>x0][1]
			leftOneY=   Y[X>x0][1]

			#print leftOneX,leftOneY
			
			
			y0=leftOneY+(1-alpha2Mean)*(x0-leftOneX)
 
			
			
			
			x1=max(X)
			y1=y0+(1-alpha2Mean)*(x1-x0)
 
			ax.plot([x0,x1],[y0,y1],'r--',lw=1,label="Fitted power law")
			
			#draw alpha1
			
			ax.text( x0+0.2, 1.4, r"$\alpha_2={:.2f} \pm {:.2f}$".format(alpha2Mean,alpha2_std),fontsize=11)
			ax.text( x0-0.7, 1.4, r"$\alpha_1={:.2f} \pm {:.2f}$".format(alpha1Mean,alpha1_std) ,fontsize=11)

			#y0=max(Y)
			
			ax.legend()

			
			x1=min(X)
			y1=y0+(1-alpha1Mean)*(x1-x0)
 
			ax.plot([x0,x1],[y0,y1],'r--',lw=1)
			
			ax.set_xlabel(r"$\log(\frac{M}{M_{\odot}})$")
			ax.set_ylabel(r"$\log(\frac{dN}{ d\log( \frac{M}{M_{\odot}} )})$")




			#cornerMap.
				
			at = AnchoredText(sourceName, loc=3, frameon=False)
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			ax.add_artist(at)



			#ax.hist(massList,bins=50 )
			plt.savefig(saveName+'maddCoreCMF.pdf',bbox_inches="tight")



		#compare the fits larger than 5 sigma


		newRow[docmf.Nclump]=coreN
		
		newRow[docmf.totalMass]=totalMass
		newRow[docmf.alpha1]=alpha1Mean
		newRow[docmf.alpha2]=alpha2Mean

		newRow[docmf.alpha1Std]=alpha1_std
		newRow[docmf.alpha2Std]=alpha2_std

		newRow[docmf.Mturn]=MtMean
		newRow[docmf.MturnStd]=Mt_std



		massList=np.array( massList)
		MvirialList=np.array( MvirialList)
		alphaVirialList=np.array( alphaVirialList)
 
		
		
		
		Gbound=massList[alphaVirialList<2]

		
		newRow[docmf.Mvirial]=np.sum(Gbound)


		#newRow[docmf.Mvirial]=np.sum(Gbound)


		newRow[docmf.MvirialRatio]=  np.sum(Gbound)/totalMass

		print np.sum(Gbound)/totalMass, "ratio of mass in gravitinal bound"
		#print docmf.Mvirial," vitial mass "

		#newRow[docmf.peakMin]=np.min( )
 
 
		return newRow

		#return  [coreN,totalMass,  [alpha1Mean, alpha1_std ], [alpha2Mean, alpha2_std ],  [MtMean, Mt_std]  ]   # massList



	def examineCOREFITS(self,coreName):
		"""
		
		compare the core fits with the the for all values large than 5 sigma
		
		"""
		coreFITS=coreName+ "_1.fits"
		
		if "G216" in coreName:
 
			CO13FITS="G216_13.fits" 
		else:
			CO13FITS="./data/cropW345CO13.fits" 

			
		
		
		
		coreData,coreHead=myFITS.readFITS(coreFITS)

		
		
		CO13Data,CO13Head=myFITS.readFITS(CO13FITS)

		different=coreData-CO13Data

		differentRMS5=different[CO13Data>self.RMS13*5]

		#print differentRMS5.shape
		return  np.std(differentRMS5,ddof=1)




	def getMassFromSum(self,coreTB):
		"""
		calculate the mass of cores with 12CO with Sum, and 12CO is used to calculate the excitation temperature
		
		
		1. According to the postion of molecular cores, get the average tempearture of 12CO as the excitation tempearture 
		
		2. Then use sum to calculte the mass
		"""
		
		
		
		
 
		Xco=2.0e20

 
		ICO=coreTB["Sum"]
		NH2=ICO*Xco
		
		size1=coreTB["Size1"]
		size2=coreTB["Size2"]
		
		
		#distance
		d=2300# pc
		parsecToCM=3.085677581e18
 
		pixLength=0.5/60.*np.pi/180*d*parsecToCM
		Area=pixLength*pixLength
		mu=1.36
		
		Mh2=3.35e-27 #kg
		solarMass=1.9891e30 #kg
		
		
		
		return Area*NH2*mu*Mh2/solarMass

	def getEdgeCenter(self,edges):
		
		return (edges[1:]+edges[0:-1])/2.



	

	def drawCMF(self,coreTB):
		"""
		draw the core mass function
		"""

		
		print "Drawing CMF with file:"+ coreTB
		core=Table.read(coreTB)
		
		coreMass=core["Sum"].data
		
		
		fig, ax = plt.subplots()
		
		plt.hist(np.log(coreMass),30,density=True)
		
		coreSaveName=coreTB.replace(".fit","")
		print "Saving the histgram as:"+coreSaveName+'_hist.eps'
		plt.savefig(coreSaveName+'_hist.eps')
		
		
		plt.clf()
		fig, ax = plt.subplots()
		
		#draw log normal of the CMF
		
		logM=np.log(coreMass)
		
		#fitting coreMass
		if 0:
			sampleArray=self.fitPowerLawWithMCMC3p(coreMass )
	 
			lpha1Mean= np.mean(sampleArray[0])
			alpha2Mean=np.mean(sampleArray[1])
			MtMean= np.mean(sampleArray[2])
			turnOverlog=np.log(MtMean)
		if 1: #log normal
			
			sampleArray=self.fitPowerLawWithMCMC_lognormal(coreMass )

			muMean= np.mean(sampleArray[1])
			sigmaMean= np.mean(   np.reciprocal(sampleArray[0]) )

			
			alpha2Mean=np.mean(sampleArray[3])
			MtMean= np.mean(sampleArray[2])
			turnOverlog=MtMean
			print sigmaMean,muMean,MtMean,alpha2Mean
			
		


		#print "The slopes are:",alpha1Mean,alpha2Mean
		print "The turn over point is:",MtMean,turnOverlog

		
		#do histgram
		
		DN,edges=np.histogram(logM,bins=15)
		centerEdge=(edges[1:] + edges[:-1]  )/2.
		
		dlogM=edges[1:]- edges[:-1]
		
		Y=np.log( DN/ dlogM )
		
		X=centerEdge
		
		
		
		#
		
		#plt.hist(coreMass,100,density=True)
		ax.plot(X,Y)
		
 
		ax.plot( [ turnOverlog, turnOverlog],[ min(Y ), max(Y)  ],lw=1 )
		
		
		
		ax.plot( [ muMean, muMean],[ min(Y ), max(Y)  ],lw=1 )
		# draw a normal distribution
		
		
		
		
		xRange=np.linspace(min(X),turnOverlog,100)
		yRange= 1./np.sqrt(2*np.pi)/sigmaMean*np.exp(-0.5/sigmaMean**2*( xRange-muMean )**2   )
		ax.plot( xRange,np.log(yRange)+6,lw=1 )

		

		self.drawLine(  ax,[ turnOverlog, max(Y) ] , alpha2Mean, [turnOverlog, max( X) ]   )
		#self.drawLine(  ax,[ turnOverlog, max(Y) ] , alpha2Mean, [turnOverlog, max( X) ]   )

		
		
		
		#self.drawLine(  ax,[ turnOverlog, max(Y) ] , 2.35, [turnOverlog, max( X) ]   )

		#self.drawLine(  ax,[ turnOverlog, max(Y) ] , alpha1Mean, [turnOverlog, min( X) ]   )

		coreSaveName=coreTB.replace(".fit","")
		
		
		print "Saving the log log diagram as:"+coreSaveName+'_log_log.eps'
		
		
		
		
		
		plt.savefig(coreSaveName+'_log_log.eps')
		
		
		#

	def drawLine(self,ax,passPoint, alpha,xRange):
		
		"""
		draw a line 
		"""
		
		
		slope=1-alpha
		
		intercept=passPoint[1]-passPoint[0]*slope
		
		
		yRange=[0,0]
		
		yRange= np.array(xRange)*slope+intercept
		
		ax.plot(xRange, yRange )
		
		

		
	def maddCMF(self):
		"""
		Do the CMF for maddalena
		"""
		
		
		#coreFile
		
		core1=Table.read("./maddalenaCore/maddCore.fit")
		core2=Table.read("./maddalenaCore/maddCoreIt2.fit")
		
		
		
		
		
		gaussclumpCore=  Table.read("./maddalenaCore/GCCore.fit")
		gaussclumpCoreRes1=  Table.read("./maddalenaCore/GCCoreVRes1.fit")

		print len(gaussclumpCore),len(gaussclumpCoreRes1)

		
		gaussclumpCore=self.clearCores(gaussclumpCore)
		core1=self.clearCores(core1)

		
		
		#draw hist gram
		fig, ax = plt.subplots()
		
		#plt.hist(core2["Sum"],100,density=True)
		#plt.hist(core1["Sum"],200,density=True)

		#plt.hist(gaussclumpCore["Sum"],200,density=True)

		from scipy.stats import expon
		
		a=np.sort(gaussclumpCore["Sum"])


		solarMass=self.getMassFromSum(gaussclumpCoreRes1 )
		
		
		solarMassFellwalker=self.getMassFromSum(core1 )

 
		alpha1, alpha2, Mt=self.doPowerLaw.fitPowerLawWithMCMC5p(solarMass)
 
 
		#alpha1  =self.doPowerLaw.fitPowerLawWithMCMC1p(solarMass[solarMass<Mt])
		#alpha2  =self.doPowerLaw.fitPowerLawWithMCMC1p(solarMass[solarMass>3000])

		#alpha2  =self.doPowerLaw.fitPowerLawWithMCMC1p(solarMass[solarMass>2000])
		#alpha2  =self.doPowerLaw.fitPowerLawWithMCMC1p(solarMass[solarMass>1000])

		#alpha2  =self.doPowerLaw.fitPowerLawWithMCMC1p(solarMass[solarMass>800])

  
		hist,edges=np.histogram(np.log(solarMass),bins=50)
		histCenter=(edges[1:]+edges[0:-1])/2.
		dx=edges[1:]-edges[0:-1] 


		X=histCenter 
		
		Y= np.log(hist/dx)

		plt.plot( X,Y)
 
		
 
 
		#plt.hist(solarMass, bins=100)
		#plt.plot([7136,7136],[0,500] )

		#print np.sum(solarMass),"total mass?"
		#ax.plot(self.getEdgeCenter(edges), Number,'-o')
		#print   
	 
		
		#ax.plot(a,2*a**(-alpha),lw=0.5	)
		

		#print "Testing Fellwalker"

		#sampleArray=self.fitPowerLawWithMCMC3p(core1["Sum"]  )
		#print "Value of k:",np.mean(sampleArray[0]),np.std(sampleArray[0],ddof=1) 
		#print "Value of alpha:",np.mean(sampleArray[1]),np.std(sampleArray[1],ddof=1) 


		#self.madeCount(gaussclumpCore["Sum"]/100)
		#self.madeCount(core1["Sum"])

		#ax.plot(sortMass,countN)
		#a1,a2= expon.fit(core1["Sum"])
		
		#a=np.linspace(0,max(core1["Sum"]),1000)
		#print expon.pdf(a,a2)

		#ax.plot(a,expon.pdf(a,a2))
		
		plt.savefig('CMF_Madd.eps')


	def getEndPoints(self,centerXY,pathAngle,length,pixLength):
		"""
		aaaa
		"""
		
		
		pixs= length/pixLength 
		
		pathAngleRadians=  np.radians( pathAngle )	 
 

		endPointX=centerXY[0] -pixs/2.*np.cos(pathAngleRadians)
		endPointY=centerXY[1] -pixs/2.*np.sin(pathAngleRadians)

		beginPointX=centerXY[0] +pixs/2.*np.cos(pathAngleRadians)
		beginPointY=centerXY[1] +pixs/2.*np.sin(pathAngleRadians)
		
		return (endPointX, endPointY), (beginPointX,beginPointY)

		
	def drawOutFlow(self,centerLB,pathAngle,length,width,vRange,COFITS,survey="IRAS"):

		"""
		this function is used to draw the outflow or to test the outflow
		
		
		four parts, CO int,arrow, WISE, and  PC diagram, and spectral lines
		
		the velocity range is in km/s
		"""
		#read CO fits
		
		#
		
		#CO FITS
		
		coData,coHead=myFITS.readFITS(COFITS)
		
		coData,coHead=self.convert3DTo2D(coData,coHead)
		
		pixLengthCO= abs(coHead["CDELT1"])
		
		WCS_CO=WCS(coHead) 
		
		centerX,centerY = WCS_CO.wcs_world2pix(centerLB[0],centerLB[1],0)
 
		beginP,endP=self.getEndPoints([centerX,centerY],  pathAngle, length, pixLengthCO)
		
		widthPix=width/pixLengthCO
		
 


 
		#outflowSizeL=max( abs(starLB[0]-centerLB[0]), abs(endLB[0]-centerLB[0]))
		
		#outflowSizeB=max( abs(starLB[1]-centerLB[1]), abs(endLB[1]-centerLB[1]))
 
		size=length*1.3
		
 
		
		#double the sieze of outflow
		fig = plt.figure(figsize=(10, 6)) 
		
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		#WISE fits
		wiseData,wiseHead=self.getWISE(centerLB,size)

		#save wise file
		cropWISEName=self.WISEcropPath+"WISE_{}_{}_{}.fits".format(centerLB[0],centerLB[1],size)
		
		if not os.path.isfile(cropWISEName):
			
		
		
			fits.writeto(cropWISEName,wiseData,header=wiseHead)
		
		
		



		WCS_WISE=WCS(wiseHead)
		

		pixLengthWISE= abs(wiseHead["CDELT1"])
		

		centerXWISE,centerYWISE = WCS_WISE.wcs_world2pix(centerLB[0],centerLB[1],0)
		beginPWISE,endPWISE=self.getEndPoints([centerXWISE,centerYWISE],  pathAngle, length, pixLengthWISE)
		widthPixWISE=width/pixLengthWISE
		






		axWISE=pywcsgrid2.subplot(221, header=WCS_WISE)
		
		axWISE.set_ticklabel_type("absdeg", "absdeg")
		
		
		#axFUGIN.imshow(FUGINdata,origin="lower",cmap="bone",vmin=3,vmax=40,interpolation='none')
		axWISE.imshow(wiseData, origin="lower",cmap="bone" , norm=colors.LogNorm(vmin=wiseData.min() , vmax= 150), interpolation='none')
		
		
		
		atN = AnchoredText(r"WISE 22$\mu$m", loc=1, prop=dict(size=8), frameon=True, )
		axWISE.add_artist(atN)
		
		
		
		# the middle figure is  spectral, including, C18O, 13CO, #better show the map of C18O....
		
		#purple cross, position of the center star
		axWISE["gal"].scatter(centerLB[0],centerLB[1],s=15,color="purple", linewidths=2,marker="+", lw=0.8,alpha=0.5)
		axWISE.axis[:].major_ticks.set_color("w")

		# draw arrow
 
 
		axWISE.arrow(beginPWISE[0],beginPWISE[1],endPWISE[0]-beginPWISE[0],endPWISE[1]-beginPWISE[1] ,  head_width=35, fc='g',ec='g',lw=0.6, length_includes_head=True, capstyle="projecting",alpha=0.5  )
		
		
		#draw width
		
		
		thetaAngle=pathAngle-90
		
		thetaAngleRadians=np.radians(thetaAngle)
		vectorShiftX=widthPixWISE/2.*np.cos(thetaAngleRadians)
		vectorShiftY=widthPixWISE/2.*np.sin(thetaAngleRadians)

		axWISE.plot([beginPWISE[0]+vectorShiftX, endPWISE[0]+vectorShiftX ] ,  [beginPWISE[1]+vectorShiftY, endPWISE[1]+vectorShiftY ] , 'g--',lw=0.5) 

		axWISE.plot([beginPWISE[0]-vectorShiftX, endPWISE[0]-vectorShiftX ] ,  [beginPWISE[1]-vectorShiftY, endPWISE[1]-vectorShiftY ] , 'g--',lw=0.5) 

		axWISE.plot([beginPWISE[0]-vectorShiftX, beginPWISE[0]+vectorShiftX ] ,  [beginPWISE[1]-vectorShiftY, beginPWISE[1]+vectorShiftY ] , 'g--',lw=0.5) 
		axWISE.plot([endPWISE[0]-vectorShiftX, endPWISE[0]+vectorShiftX ] ,  [endPWISE[1]-vectorShiftY, endPWISE[1]+vectorShiftY ] , 'g--',lw=0.5) 




		axWISE.plot(    )

		#import matplotlib.patches as mpatch
		

		#arrow = mpatch.FancyArrowPatch((0, 0), (0.5, 0.5))
		#axWISE.add_artist(arrow)
		
	 
		
		
		
		### PV diagram

		
		from pvextractor import extract_pv_slice,Path

		
 
		endpoints = [beginP,endP]
		xy = Path(endpoints,width= widthPix )

		pv = extract_pv_slice(  self.CO12HDU, xy)

 
		pvCO13=extract_pv_slice(  self.CO13HDU, xy)

 
		outflowPV="my_slice.fits"
		
		os.system("rm "+outflowPV)
		pv.writeto(outflowPV)
 
		#read data and fits
		pvData,pvHead=pv.data,pv.header
		pvWCS=WCS(pvHead)
		
		a,vPixStart=pvWCS.wcs_world2pix(0,vRange[0]*1000,0)
		a,vPixEnd=pvWCS.wcs_world2pix(0,vRange[1]*1000,0)

		
		axPV=pywcsgrid2.subplot(122, header=WCS(pvHead))
		axPV.imshow(pvData, origin="lower",cmap="bone" ,vmin=0 ,vmax= 8  ,    interpolation='none')


		#axPV.set_ticklabel_type("manual", "manual")
		
		axPV.set_ticklabel1_type("manual",locs=[0.1,0.2,0.3,0.4,0.5],labels=["0.1","0.2","0.3","0.4","0.5"])
 
		
		vInterval=2 #km/s
		
		yLocs=np.arange(int(vRange[0])+vInterval, int(vRange[1])+vInterval,2)
		
		yLabels= map(int, yLocs)
		yLabels=map(str,yLabels)
 
		axPV.set_ticklabel2_type("manual",locs=yLocs*1000.,labels=yLabels)

		axPV.set_ylabel(r"Velocity ($\rm km \ s^{-1}$)")
		
		axPV.set_xlabel(r"Offset (degree)")
		axPV.set_ylim(vPixStart,vPixEnd)
		
 
		#replace the ticks
		axPV.axis[:].major_ticks.set_color("w")

		#axPV.set_xticks((0,0.1,0.2,0.3,0.4), ("a","b","c","d","e") )
		#axPV.set_yticks( np.array([36000,40000]))

		#draw spectral lines 
 
		axSpectral=plt.subplot(223)
		
		avgSpec=np.mean(pvData,axis=1)
		
		vIndex=np.arange( len(avgSpec) )
		
		a,velocities=pvWCS.wcs_pix2world(0,vIndex,0) 
 
		velocities=np.array(velocities)/1000.  
		
		WCSCO13=WCS(pvCO13.header)
		pv13Data=pvCO13.data
		avgSpec13=np.mean(pv13Data,axis=1)
		vIndex13=np.arange( len(avgSpec13) )

		a,velocities13=WCSCO13.wcs_pix2world(0,vIndex13,0) 
		velocities13=np.array(velocities13)/1000.  


		
		axSpectral.plot([0,70],[0,0],color="black")
		axSpectral.step(velocities,avgSpec,lw=0.8,color='dimgray',label=r"$^{12}\mathrm{CO}\left(J=1\rightarrow0\right)$")
		axSpectral.step(velocities13,avgSpec13,lw=0.8,color='darkviolet',label=r"$^{13}\mathrm{CO}\left(J=1\rightarrow0\right)$")


		
		axSpectral.legend()
		axSpectral.set_xlabel(r"Velocity ($\rm km \ s^{-1}$)")
		axSpectral.set_ylabel(r"T$_{\rm mb}$ (K)")


		outflowName='outflow_{}_{}_{}_{}.png'.format(centerLB[0],centerLB[1],pathAngle,survey)
		
		plt.savefig( self.outflowFigurePath+outflowName, bbox_inches="tight",dpi=300)
#






	def drawAllOutflow(self):
		
		pathAngles= [0.,45.,90.,135.] #degree
		
		length= 0.5 #degre
	 
		vRange=[0.,70.]
		width=0.05
	
		WISEcat=self.WISEHIICat.copy()
		
		colLWISE="GLong<br>(deg.)"
		colBWISE="GLat<br>(deg.)"
		
		maddaHII=self.filterTBByRange(WISEcat,colLWISE,[210,219.5])
		maddaHII=self.filterTBByRange(maddaHII,colBWISE,[-5.,5.])



		MaddaIRASHII=self.filterTBByRange(self.IRASHIICat,"glon",[210-360,219.5-360])
		MaddaIRASHII=self.filterTBByRange(MaddaIRASHII,"glat",[-5.,5.])

 
		colLIRAS="glon"
		colBIRAS="glat"
 



		for eachIRASHII in MaddaIRASHII:
			surveyIRAS='IRAS'
			#print eachIRASHII[colLIRAS],eachIRASHII[colBIRAS],surveyIRAS

			centerLB=[eachIRASHII[colLIRAS],eachIRASHII[colBIRAS]]
			
			for pathAngle in pathAngles:
				doMadd.drawOutFlow(centerLB,pathAngle,length,width,vRange,"CO13RGB_G.fits",survey=surveyIRAS)


		for eachWISEHII in maddaHII:
			surveyWISE='WISE'

			centerLB=[ eachWISEHII[colLWISE],eachWISEHII[colBWISE]   ]
			
			for pathAngle in pathAngles:
				doMadd.drawOutFlow(centerLB,pathAngle,length,width,vRange,"CO13RGB_G.fits",survey=surveyWISE)


	def drawWISEImage(self,l,b,size=0.5,survey="IRAS"):
		"""
		Draw WISE iamge of IRAS sources
		
		"""
		pass
		

	def filterMaddCore(self,coreTB):
		"""
		This function is used to select cores that belongs to the maddalena molecular cores,
		
		# we calculate the cores with the distance 2.2 kpc used by "rethinking a mysterious molecular cloud"
		"""
		
		#the l range of maddalena cloud 
		
		print "Filtering the maddalena cores"
		
		#lRange_Maddalena=[213.8,219.44] #deg 
		#bRange_Maddalena=[-4.9 ,-1.25 ]  #deg
		#vRange_Maddalena= [16604, 38687.4] #km/s
		
		
		maddCoreTB=coreTB.copy()
		
		maddCoreTB.add_index("Peak1")
		
		maddCoreTB=maddCoreTB.loc["Peak1",213.8:219.44]
		
		maddCoreTB.add_index("Peak2")
		
		maddCoreTB=maddCoreTB.loc["Peak2", -4.9:-1.25]
		
		maddCoreTB.add_index("Peak3")
		
		maddCoreTB=maddCoreTB.loc["Peak3", 16604:38687.4]
		#maddCoreTB=coreTB[0:1]
		#maddCoreTB.remove_row(0)
		
		
		#for eachCore in coreTB:
			
			
			#if 
			
			#maddCoreTB.add_row(eachCore)
		
		
		return   maddCoreTB


	def examineRing(self):
		
		"""
		Examine the elliptical ring
		"""
		
		#first do momment in velocity range [ , ]
		vRange=[25., 27.] #kms
		#vRange=[27., 29.] #kms

		
		#imageRms=
		
		myFITSDo=myFITS()
		
		Rdata,Rhead= myFITSDo.momentFITS(self.COFile13, vRange,0,outFITS="ellipticcalRing.fits")
		

	def outflowMannually(self,centerLB,length=0.5,width=0.05):
		
		
		"""
		
		
		"""
 
		
		pass
		
		

	def getRgcAlpha(self,d,l,b):
		
		"""
		
		
		
		"""
		#d=d*np.cos(np.radians(b) )
		
		
		R0=8.2 #34 the background would be changed
		Theta0=240.
		dThetadR=-0.2
		
		#Theta=theta0+dThetadR*(R-R0) #km/s
		

		d=d*np.cos( np.deg2rad(b) )

		angle0=np.radians(l-90)

		
		v1x=d*np.cos(angle0)
		v1y=d*np.sin(angle0)
		
		
		vx=v1x
		vy=v1y+R0
		
		r=np.sqrt( vx**2+vy**2)
		
		
		if vy>0:
			
			alpha=np.arccos(  vx/r)
		if vy<0 and vx>0:
			alpha=np.arcsin(  vy/r)

		if vy<0 and vx<0:
			alpha=np.arcsin(  abs(vy/r) )+np.pi 
 
		
		#with d cal R
		
		#smallL=np.radians(l-180.)

		#r2=(d*np.sin(smallL) )**2+( R0+d*np.cos(smallL) )**2

		#r=np.sqrt(r2)
		
 
 
 
		#alpha=np.arcsin(d*np.sin(smallL)/r)
		
		
		
		return r,alpha*180/np.pi #inradians 
		
		
		

	def ReidA5(self, l,b,v):
		
		"""
		return the distance with A5 model of Reid 2014
		
		#what about the errors?
		
		return distance 
		"""
		
		
		#first convert l,b to  the string use by the ma
		
		c = SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')
		
		
		cor=str(c.icrs.to_string( 'hmsdms'))
		
		raStr,decStr=cor.split()
		
		raStr=raStr.replace("h",'')
		raStr=raStr.replace("m",'')
		raStr=raStr.replace("s",'')
		
		decStr=decStr.replace("d",'')
		decStr=decStr.replace("m",'')
		decStr=decStr.replace("s",'')
		
		sourceStr="calsource {} {} {} 0".format(raStr,decStr,v)
		
		#our source
		
		outSourceStr= 'echo "{}" > ./2014revised_kd/source_file.dat'.format(sourceStr)
		
		
 
		
		
		os.system(outSourceStr)
		
		os.system("cd ./2014revised_kd; ./a.out > distanceInfo.txt; cd ..")

		with open("./2014revised_kd/distanceInfo.txt") as f:
			lines=f.readlines()
			
		disStr=""
		for eachL in lines:
			
			if eachL[0]=="!":
				continue
			
			disStr=eachL
			
		distance= float(disStr.split()[5])
				
		disErrorUp=  float(disStr.split()[6])
		disErrorLow=  float(disStr.split()[7])

		rGC,alpha=self.getRgcAlpha(distance,l,b)
				

		return [distance,disErrorLow,disErrorUp  ], rGC,alpha
		#print this to  ./2014revised_kd/source_file.data





	def drawArmsNew(self,coreTB):
		
		"""
		calculate the distances of cores, then draw, spiral arms with it
		"""
		
		xLim=[-10.05,10.05]
		yLim=[-3,17.55]
		
		mwData,mwHead=myFITS.readFITS(self.MWFITS)
		
		
		mwWCS=WCS(mwHead)
		
		fig=plt.figure(figsize=(6,6))
		
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
		
		rc('text', usetex=True)
		
		
		ax=pywcsgrid2.subplot(111,header=WCS(mwHead))
		
		
		
		#ax.imshow(mwData,origin='lower',cmap="gray",  norm=LogNorm(vmin=60 ,vmax=500),)
		
		
		#ax.imshow(mwData,origin='lower',cmap="Greys",  vmin=0.00001  ,vmax=0.03)

		#ax.imshow(mwData,origin='lower',cmap="Greys",  norm=LogNorm(vmin=0.0025 ,vmax=0.13),)


		if 1: #longtitude
			lRange=[-4,4]
			lInterval=1 #degree
			
			xLocs=  np.arange( (lRange[0]) ,  (lRange[1])+lInterval,lInterval)    #np.arange( (lRange[0]) ,  (lRange[1])+lInterval,lInterval)
			
			#xLabels= map(int, xLocs)
			xLabels=map(str,xLocs)
	
			locsX=xLocs  #abs(xLocs-219.5)
	
			#print locsX
	
			ax.set_ticklabel1_type("manual",locs=locsX,labels=xLabels)
	
			ax.set_xlabel(r"X (kpc)")

		
		if 1: #longtitude
			lRange=[4,12]
			lInterval=1 #degree
			
			xLocs=  np.arange( (lRange[0]) ,  (lRange[1])+lInterval,lInterval)    #np.arange( (lRange[0]) ,  (lRange[1])+lInterval,lInterval)
			
			#xLabels= map(int, xLocs)
			xLabels=map(str,xLocs)
	
			locsX=xLocs  #abs(xLocs-219.5)
	
			#print locsX
	
			ax.set_ticklabel2_type("manual",locs=locsX,labels=xLabels)
	
			ax.set_ylabel(r"Y (kpc)")
		
		
		#psun (0, 8.34)
		wcsMW=WCS(mwHead)
		
		lowerLeft=wcsMW.wcs_world2pix( xLim[0],yLim[0],0 )
		lowerRight=wcsMW.wcs_world2pix( xLim[1],yLim[1],0 )

		
		
		ax.set_xlim(  lowerLeft[0] ,  lowerRight[0]   )
		ax.set_ylim(  lowerLeft[1] ,  lowerRight[1]   )







		gaiaDR2OCTB=Table.read("/home/qzyan/WORK/projects/maddalena/data/GAIADR2OC.fit")
		MWSCOCTB=Table.read("/home/qzyan/WORK/projects/maddalena/data/MWSCOC.fit")




		OstarTB=Table.read("/home/qzyan/WORK/projects/maddalena/data/XuOstar.fit")



		gaiaDR2OCTB=gaiaDR2OCTB.copy()

		logt='logt'
		#maxAge=  np.log10(5.5e7) #7*np.log10(5.5) #8.5  #7*np.log10(1.2) #1.2e7 

		
		
		if 0:
			ageCol=gaiaDR2OCTB['X'].copy()
			ageCol=ageCol*0+1000.
			ageCol.name=logt
			gaiaDR2OCTB.add_column(ageCol)
 
			#youngMWSC=MWSCOCTB.loc[logt,:maxAge]
		
		#print "ploting OCs"
			for eachOC in gaiaDR2OCTB:
				
				OCName=eachOC['Cluster']
				
				mwscOC=  MWSCOCTB[MWSCOCTB['Name']==OCName  ] 
				#print len(mwscOC[logt])
				
				if len(mwscOC)==0:
					continue
				
				#print OCName, mwscOC[logt][0], mwscOC['Name'][0]
				  
				#if mwscOC[logt][0] <maxAge:
				 
				eachOC[logt]=mwscOC[logt] 
			#save fits
			
			
			
			gaiaDR2OCTB.write("GAISDR2COWithAge.fit")
				
		
		
		
		gaiaDR2OCTBWithAge=Table.read("GAISDR2COWithAge.fit")
		print "ploting OCs"
		maxAge=   np.log10(1.2e7) #7*np.log10(5.5) #8.5  #7*np.log10(1.2) #1.2e7 

		#maxAge=np.log10(5.5e7) 
		OCcolor='purple'
		for eachOC in gaiaDR2OCTBWithAge:
			break
			if eachOC[logt]<maxAge:
				
				l=eachOC["GLON"]
				b=eachOC["GLAT"]

				distance=eachOC["dmode"]/1000.
				
				
				error16= eachOC["d16"]/1000.
				error84= eachOC["d84"]/1000.
				
				
				Rgc,alpha=self.getRgcAlpha(distance,l,b)
				alpha=np.radians(alpha)
				#print distance,eachMaser["plx"] l,alpha
				
				ax["gal"].scatter( Rgc*np.cos( alpha ) , Rgc*np.sin( alpha ) ,s= 10,  facecolors='none', edgecolors=OCcolor,lw=0.3,marker="s")
				
				Rgc,alpha1=self.getRgcAlpha(error16,l,b)
				alpha1=np.radians(alpha1)
				
				xError16= Rgc*np.cos( alpha1 )
				yError16= Rgc*np.sin( alpha1 )
	
				
				#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")
	
				Rgc,alpha2=self.getRgcAlpha(error84,l,b)
				alpha2=np.radians(alpha2)
				#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")
				xError84= Rgc*np.cos( alpha2 )
				yError84= Rgc*np.sin( alpha2 )
				
				
				ax["gal"].plot([xError16,xError84], [yError16,yError84], lw=0.5,color=OCcolor )
					






		print "Plotting out arm in the second quandrant" #2016ApJS..224....7D
		 
		#better code the cloud size with mass 
		
		
		outArm2rdTB=Table.read("./data/outArm2rd.fit")
		outArmColor='blue'
		
 
		for eachMC2rd in outArm2rdTB:
			break
			s= eachMC2rd["Mass"]*10 
			

			l=eachMC2rd["_Glon"]
			b=eachMC2rd["_Glat"]
			distance=eachMC2rd["Dist"]
			Rgc,alpha=self.getRgcAlpha(distance,l,b)
			alpha=np.radians(alpha)
			#print distance,eachMaser["plx"] l,alpha
			ax["gal"].scatter( Rgc*np.cos( alpha ) , Rgc*np.sin( alpha ) ,s= s,  facecolors='none', edgecolors=outArmColor,lw=0.3,marker=".")
			
	
		

 

 

		#====================================
		
		print "ploting masers"

					
		#output those core to 
		
		maserTB=Table.read("maserTB.csv")
		
		
		
		
		for eachMaser in maserTB:
			l=eachMaser["_Glon"]
			
			#if l>210 and l<220:
				#print  eachMaser
			
			
			b=eachMaser["_Glat"]
			
			
			if l>210 and l<220:
				
				print l,b,"???????????????????????????????"
			
			
			if eachMaser["e_plx"]/eachMaser["plx"]>0.2:
				continue
				
				
			colorEdge='black'
			if eachMaser["Arm"]=='Loc':
				colorEdge='darkorange'
				#colorEdge='pi'

 
			if eachMaser["Arm"]=='Out':
				colorEdge='g'
				
			if eachMaser["Arm"]=='Per':
				colorEdge='b'
			if eachMaser["Arm"]=='Sgr':
				colorEdge='cyan'
				
			if 0:
				np.random.seed()
				distanceSample=np.random.normal( eachMaser["plx"], eachMaser["e_plx"],10000 ) 
				

				
				distanceSample=1./distanceSample  # kpc  
				
				distance=  np.mean(distanceSample)*np.cos( np.radians(b) )
				error16=np.percentile(distanceSample,16 )*np.cos( np.radians(b) )
				error84=np.percentile(distanceSample,84 )*np.cos( np.radians(b) )

			

			if 1:

				distance=1./eachMaser["plx"] 
				error16=1./(eachMaser["plx"]+ eachMaser["e_plx"])
				error84=1./(eachMaser["plx"]- eachMaser["e_plx"])

			#if distance>7:
 
				#print eachMaser["plx"], eachMaser["e_plx"], error16, distance,error84
			

			Rgc,alpha=self.getRgcAlpha(distance,l,b)
			alpha=np.radians(alpha)
			#print distance,eachMaser["plx"] l,alpha
			
			ax['gal'].scatter( Rgc*np.cos( alpha ) , Rgc*np.sin( alpha ) ,s= 15,  facecolors='none', edgecolors=colorEdge,lw=0.6,marker="^")
			
			Rgc,alpha1=self.getRgcAlpha(error16,l,b)
			alpha1=np.radians(alpha1)
			
			xError16= Rgc*np.cos( alpha1 )
			yError16= Rgc*np.sin( alpha1 )

			
			#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")

			Rgc,alpha2=self.getRgcAlpha(error84,l,b)
			alpha2=np.radians(alpha2)
			#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")
			xError84= Rgc*np.cos( alpha2 )
			yError84= Rgc*np.sin( alpha2 )
			
			
			ax['gal'].plot([xError16,xError84], [yError16,yError84], lw=0.6,color=colorEdge )
			
			#if distance>7:
			
				#print alpha1,alpha2,alpha
				
				#print xError16, xError84
			
			#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker="^")

		########################
		print "Draw Gaia distance clouds"
		
		
		gaiaMCTB=Table.read("gaiaMCDis.fit")
		#yan Et al . 
		for eachO in gaiaMCTB:
			
			#break #no need to draw, all local
			l=eachO["l"]
			
			#if l>210 and l<220:
				#print  eachMaser
 
			b=eachO["b"]
			
			
			distance= eachO["distance"]/1000.
			
			
			projRatio=np.cos( np.radians(b) )
			
			disStd=  (eachO["upperDis"] + eachO["lowerDis"] )/4./1000.
			
			sysError= disStd**2+ (distance*0.05 )**2
			
			sysError=np.sqrt( sysError ) 
			
			
			
			#error16=   #( eachO["distance"] -  eachO["lowerDis"])*np.cos( np.radians(b) ) /1000.
			#error84=  #( eachO["distance"] + eachO["upperDis"] ) *np.cos( np.radians(b) )/1000.

			
			error16=distance- sysError #5% systematic error
			error84=distance+ sysError


			
			distance=distance*projRatio
			error16=error16*projRatio
			error84=error84* projRatio

 

			Rgc,alpha=self.getRgcAlpha(distance,l,b)
			alpha=np.radians(alpha)
			#print distance,eachMaser["plx"] l,alpha
			
			resizeCloud=["Mon R2", "Rosette"]
			
			size= {"Mon R2":40*1.078594861373204 , "Rosette": 160*0.83265625 }
			
			if eachO["MBMID"] not in resizeCloud:
			
			
				ax['gal'].scatter( Rgc*np.cos( alpha ) , Rgc*np.sin( alpha ) ,s=3,  facecolors='none', edgecolors='purple',lw=0.3,marker="o")
			else:
				
				cise=size[ eachO["MBMID"]]
				
				ax['gal'].scatter( Rgc*np.cos( alpha ) , Rgc*np.sin( alpha ) ,s= cise/2.,  facecolors='none', edgecolors='purple',lw=0.3,marker="o")
			
			
			
			Rgc,alpha1=self.getRgcAlpha(error16,l,b)
			alpha1=np.radians(alpha1)
			
			xError16= Rgc*np.cos( alpha1 )
			yError16= Rgc*np.sin( alpha1 )

			
			#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")

			Rgc,alpha2=self.getRgcAlpha(error84,l,b)
			alpha2=np.radians(alpha2)
			#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")
			xError84= Rgc*np.cos( alpha2 )
			yError84= Rgc*np.sin( alpha2 )
			
			
			ax['gal'].plot([xError16,xError84], [yError16,yError84], lw=0.3,color='purple' )
		
 
		
		#draw G130150
		if 1:

			newDisCatG130150="/home/qzyan/WORK/projects/maddalena/dendroDisPath/G130150/G130150goodDisTB.fit"

			print "Draw clouds G130150"
			self.plotDistanceTB( newDisCatG130150,ax,isG210=False)


		#draw G2650
		if 1:

			newDisCatG2650="/home/qzyan/WORK/projects/maddalena/dendroDisPath/G2650/G2650goodDisTB.fit"

			print "Draw clouds G2650"
			self.plotDistanceTB( newDisCatG2650,ax,isG210=False)



		########################

		
		#print "Draw clouds G210
		#self.plotDistanceTB( "/home/qzyan/WORK/projects/maddalena/G130150Path/goodDistance.fit",ax)
 
		########################
		
		
		newDisCatG210="/home/qzyan/WORK/projects/maddalena/dendroDisPath/G210/disCatG210.fit"

		print "Draw clouds G210"
		self.plotDistanceTB( newDisCatG210,ax)
		
		########################





		########################
		print "Draw clouds G210 manuselected cloud region"
		
		#extraNameList=[ "G214.8+04.0", "G214.4-04.3", "G211.1-02.1" , "W5"   ]
		
		extraNameList=[ "G214.8+04.0", "G214.4-04.3", "G211.1-02.1"    ]

				
		gaiaTBCollect=Table.read( "/home/qzyan/WORK/projects/gaiaDistance/data/disCollection.fit" )  
		doG210Dis=disTB( "/home/qzyan/WORK/projects/gaiaDistance/data/disCollection.fit" )
		G210Color="red"




		for eachO in gaiaTBCollect:
			
			if eachO[disTB.sourceName] not in extraNameList:
				continue
			
			
			#break #no need to draw, all local
			
			
			l=eachO["l"]
			
 
			b=eachO["b"]

			if l==0.0:
				
				l=eachO["cloudBoxCenterL"]
				
		
				b=eachO["cloudBoxCenterB"]
			
			distance= eachO["distance"]/1000.
			
			disStd=eachO[disTB.disStd]/1000.
			sysError= (distance*0.05 )**2+disStd**2
			sysError=np.sqrt( sysError )

			
			error16= distance - sysError
			error84=  distance + sysError
			
			distance=distance*np.cos( np.radians(b) ) 
			error16=error16* np.cos( np.radians(b) ) 
			error84=error84* np.cos( np.radians(b) ) 



			#distance= eachO["distance"]*np.cos( np.radians(b) )/1000.
			
			#sysError= distance*0.05

			
			#error16= ( eachO["distance"] -  eachO[disTB.disHPDLow])*np.cos( np.radians(b)  ) /1000. #-sysError
			#error84= ( eachO["distance"] + eachO[disTB.disHPDUp] ) *np.cos( np.radians(b) )/1000. #+sysError


			Rgc,alpha=self.getRgcAlpha(distance,l,b)
			alpha=np.radians(alpha)
			#print distance,eachMaser["plx"] l,alpha
			mass=round(doG210Dis.calMassByRow( eachO)/1000.,1) #1000 solar

			ax['gal'].scatter( Rgc*np.cos( alpha ) , Rgc*np.sin( alpha ) ,s= mass/2.,  facecolors='none', edgecolors=G210Color,lw=0.3,marker="o")
			Rgc,alpha1=self.getRgcAlpha(error16,l,b)
			alpha1=np.radians(alpha1)
			
			xError16= Rgc*np.cos( alpha1 )
			yError16= Rgc*np.sin( alpha1 )
 
			Rgc,alpha2=self.getRgcAlpha(error84,l,b)
			alpha2=np.radians(alpha2)
			#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")
			xError84= Rgc*np.cos( alpha2 )
			yError84= Rgc*np.sin( alpha2 )
			
			
			ax['gal'].plot([xError16,xError84], [yError16,yError84], lw=0.3,color=G210Color )
		
 
		
	
		
		
		
		
		
		###############################
		
		print "Draw lines and Markers"
		
		xSun=0.
		ySun=8.2 #
		#draw sun symbol
		ax["gal"].scatter(xSun,ySun,s=1,facecolors="red")
		ax["gal"].scatter(xSun,ySun,s=20,facecolors='none', edgecolors='black',lw=0.3)

		#ax.text(xSun+0.2,ySun+0.2, 'Sun')

		
		pSun=np.array( [xSun,ySun  ]) 
		#draw anglues
		
		maxDrawR=8.5
		
		sunRadius=0.12 #
		
		#for drawL in [ 180, 190,200,210,220,230,240,250,260,270]:
		for drawL in [ 30, 60, 90,  120, 150,180,  210, 240,   270,300,330 ]:

		
			#true angles
			drawAngle=np.radians(drawL-90)
			
			unitvector=np.array( [ np.cos(drawAngle), np.sin(drawAngle) ]  )
			drawp1= pSun-unitvector*maxDrawR

			drawp1_end= pSun-unitvector*sunRadius

			 
			ax["gal"].plot( [drawp1_end[0],drawp1[0]],   [drawp1_end[1],drawp1[1]] ,'--',color='black',lw=0.5,alpha=0.5 )

			drawp2_start= pSun+unitvector*sunRadius

			drawp2= pSun+unitvector*maxDrawR
			ax["gal"].plot( [drawp2_start[0],drawp2[0]],   [drawp2_start[1],drawp2[1]] ,'--',color='black',lw=0.5,alpha=0.5 )

			#ax.plot( [drawp1[0],drawp2[0]],   [drawp1[1],drawp2[1]] ,'--',color='green',lw=0.8,alpha=0.7 )
			
						
			refPosition= pSun+unitvector*(4.05 )

			#ax.text(drawp2Text[0]-0.6,drawp2Text[1]+0.2, r"{}$^\circ$".format(drawL) )
			
			#ax.text(drawp2Text[0] ,drawp2Text[1]  , r"${}^\circ$".format(drawL),fontsize=6 )

			
			#refPosition= drawp2Text
			
			thetaShift=  np.radians( drawL  )
			
			shiftVector=  np.array( [ np.cos(thetaShift), np.sin(thetaShift) ]  )
			
			shiftVector=shiftVector *0.5
			
			drawp2Text=refPosition+shiftVector
			
			
			#ax["gal"].text(drawp2Text[0] ,drawp2Text[1]  , r"${}^\circ$".format(drawL),fontsize=6 ,rotation=  drawL -180)
			ax["gal"].text(drawp2Text[0] ,drawp2Text[1]  , r"${}^\circ$".format(drawL),fontsize=6,rotation=drawL-180, rotation_mode='anchor'  )

			#if drawL<120 and  drawL>0:
			
			
				#ax["gal"].text(drawp2Text[0] ,drawp2Text[1]-0.1 , r"${}^\circ$".format(drawL),fontsize=6 )
				
			
			#if drawL<270 and  drawL>=120:
			
			
				#ax["gal"].text(drawp2Text[0] -0.4,drawp2Text[1]+0.1 , r"{}$^\circ$".format(drawL),fontsize=6)
				
			#if  drawL>=270:
			
			
				#ax["gal"].text(drawp2Text[0]-0.6,drawp2Text[1] -0.2, r"{}$^\circ$".format(drawL) ,fontsize=6 )
				
			
				
			#else:
				#ax.text(drawp2Text[0]-0.6,drawp2Text[1]+0.2, r"{}$^\circ$".format(drawL) )

				

		for radiusCircle in [ 1,2,3,4]:
			"""
			"""
			#circle=plt.Circle((xSun,ySun),radius=3 )
			
			drawA=np.linspace(0,2*np.pi,200)
			
			x=radiusCircle*np.cos(drawA)
			y=radiusCircle*np.sin(drawA)
			ax["gal"].plot( x +xSun ,   y +ySun,'--',color='black',lw=0.5,alpha=0.5 )
			
			#if radiusCircle!=2:
			ax["gal"].text( radiusCircle*np.cos(np.pi/2.) +xSun+0.05,  radiusCircle*np.sin(np.pi/2.)+ySun+0.05, "{} kpc".format( radiusCircle) ,fontsize=7 )

			
			#ax.add_patch(circle)

		ax["gal"].scatter(0,0,s=15,facecolors="black")
		ax["gal"].text(-2.8,-0.8,"the Galactic Center")

		
		#ax.set_xlabel(r"kpc")	
		#ax.set_ylabel(r"kpc")	
		
		if 1: # set the limit
		
			extendN=20
			
			figureSize=4.5 # kpc
			
			minX,minY=mwWCS.wcs_world2pix( -figureSize, ySun- figureSize ,0  )
			maxX,maxY=mwWCS.wcs_world2pix(  figureSize, ySun+ figureSize ,0  )

			ax.set_xlim( minX, maxX )
			ax.set_ylim( minY, maxY   )



			#ax.set_xlim(2400-extendN,3200+extendN )
			#ax.set_ylim(3480-extendN,3480+800+extendN )


		# fix legend
		if 1:
			#ax.scatter( -100 , 0 ,s= 15,  facecolors='none', edgecolors='g',lw=0.6,marker="^",label='Perseus')
			#ax.scatter( -100 , 0 ,s= 15,  facecolors='none', edgecolors='r',lw=0.6,marker="^",label='Local')
			#ax.scatter( -100 , 0 ,s= 15,  facecolors='none', edgecolors='cyan',lw=0.6,marker="^",label='Sagittarius')
			#ax.scatter( -100 , 0 ,s= 15,  facecolors='none', edgecolors='b',lw=0.6,marker="^",label='Outer')
			
			#ax.scatter( -100 , 0 ,s= 20,  facecolors='none', edgecolors='purple',lw=0.6,marker="s",label=r'Open Clusters ($<$ 12 Myr)')
			ax.scatter( -100 , 0 ,s= 20,  facecolors='none', edgecolors='black',lw=0.6,marker="^",label='Masers (all colors)')
			#ax.scatter( -100 , 0 ,s= 20,  facecolors='none', edgecolors='black',lw=0.6,marker="*",label='O stars')
			ax.scatter( -100 , 0 ,s= 40,  facecolors='none', edgecolors='purple',lw=0.6,marker=".",label='Molecular Clouds (Yan et al. 2019)')
			ax.scatter( -100 , 0 ,s= 40,  facecolors='none', edgecolors='red',lw=0.6,marker=".",label='Molecular Clouds (this work)')

			#ax.scatter( -100 , 0 ,s= 15,  facecolors='none', edgecolors='b',lw=0.6,marker=".",label='Clouds')
	
			ax.legend(loc=3,fontsize='small',ncol=1)
		
		
		
		
		
		#replace 
		
		
		#ax["gal"].plot(0,8.34,'r.')
		#
		#ax["gal"].plot([-10,10],[0,10],'g')
		
		
		#ax["gal"].plot(-5,8.34,'r.')





		#ax.axis[:].major_ticks.set_color("b")
		
		
		#fig.savefig('MWBackWCS.pdf', bbox_inches="tight",dpi=300)
	
		fig.savefig( 'armMap.pdf', bbox_inches="tight")
		fig.savefig( 'armMap.png', bbox_inches="tight",dpi=600)


	def plotDistanceTB(self,disTBName,ax,colorDraw='red',isG210=True):
		
		"""
		plot distance TB
		"""



		processDis=disTB(disTBName)



		allDistances=processDis.getTB()

		if isG210:
			orderDict=self.getNDict()
		
		#cloudName= self.getCloudNameByID(ID) 		

		for eachO in allDistances:
			
			#break #no need to draw, all local
			l=eachO["l"]
			
 
			b=eachO["b"]
			
			
			cloudName=self.getCloudNameByLB(l,b)

			if 	 isG210:
				N=orderDict[cloudName]

				if N>11:
					continue

			distance= eachO["distance"]/1000.
			
			disStd=eachO[disTB.disStd]/1000.
			sysError= (distance*0.05 )**2+disStd**2
			sysError=np.sqrt( sysError )

			
			error16= distance - sysError
			error84=  distance + sysError
			
			
			#error16= ( eachO["distance"] -  eachO[disTB.disHPDLow])*np.cos( np.radians(b)  ) /1000.  #- distance*0.05*np.cos( np.radians(b)  )
			#error84= ( eachO["distance"] + eachO[disTB.disHPDUp] ) *np.cos( np.radians(b) )/1000.   #+distance*0.05*np.cos( np.radians(b)  )
			
 
			
			distance=distance*np.cos( np.radians(b) ) 
			error16=error16* np.cos( np.radians(b) ) 
			error84=error84* np.cos( np.radians(b) ) 


			#error16= ( eachO["distance"] -  eachO[disTB.disHPDLow])*np.cos( np.radians(b) ) /1000.  #- sysError
			#error84= ( eachO["distance"] + eachO[disTB.disHPDUp] ) *np.cos( np.radians(b) )/1000. #+ sysError


			Rgc,alpha=self.getRgcAlpha(distance,l,b)
			alpha=np.radians(alpha)
			#print distance,eachMaser["plx"] l,alpha
			
			#get mass
			
			mass=round(processDis.calMassByRow( eachO)/1000.,1) #1000 solar

			if np.isnan(mass):
				
				print mass,"??????????????"
				continue



			ax['gal'].scatter( Rgc*np.cos( alpha ) , Rgc*np.sin( alpha ) ,s= mass/2.,  facecolors='none', edgecolors=colorDraw,lw=0.3,marker="o")
	 
			
			Rgc,alpha1=self.getRgcAlpha(error16,l,b)
			alpha1=np.radians(alpha1)
			
			xError16= Rgc*np.cos( alpha1 )
			yError16= Rgc*np.sin( alpha1 )
 
			Rgc,alpha2=self.getRgcAlpha(error84,l,b)
			alpha2=np.radians(alpha2)
			#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")
			xError84= Rgc*np.cos( alpha2 )
			yError84= Rgc*np.sin( alpha2 )
			
			
			ax['gal'].plot([xError16,xError84], [yError16,yError84], lw=0.3,color=colorDraw )
		
 

	def aaaaaaadfasdf(self,coreTB):
		
		"""
		calculate the distances of cores, then draw, spiral arms with it
		"""
		

		
		#print coreTB.colnames
 
		return
		#cal core distances once for all
		
		
		#out put all source to 
		
		
		
		
		fig=plt.figure(figsize=(8,6))
		ax=fig.add_subplot(111)
		
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		gaiaDR2OCTB=Table.read("/home/qzyan/WORK/projects/maddalena/data/GAIADR2OC.fit")
		MWSCOCTB=Table.read("/home/qzyan/WORK/projects/maddalena/data/MWSCOC.fit")




		OstarTB=Table.read("/home/qzyan/WORK/projects/maddalena/data/XuOstar.fit")



		gaiaDR2OCTB=gaiaDR2OCTB.copy()

		logt='logt'
		#maxAge=  np.log10(5.5e7) #7*np.log10(5.5) #8.5  #7*np.log10(1.2) #1.2e7 

		
		
		if 0:
			ageCol=gaiaDR2OCTB['X'].copy()
			ageCol=ageCol*0+1000.
			ageCol.name=logt
			gaiaDR2OCTB.add_column(ageCol)
 
			#youngMWSC=MWSCOCTB.loc[logt,:maxAge]
		
		#print "ploting OCs"
			for eachOC in gaiaDR2OCTB:
				
				OCName=eachOC['Cluster']
				
				mwscOC=  MWSCOCTB[MWSCOCTB['Name']==OCName  ] 
				#print len(mwscOC[logt])
				
				if len(mwscOC)==0:
					continue
				
				#print OCName, mwscOC[logt][0], mwscOC['Name'][0]
				  
				#if mwscOC[logt][0] <maxAge:
				 
				eachOC[logt]=mwscOC[logt] 
			#save fits
			
			
			
			gaiaDR2OCTB.write("GAISDR2COWithAge.fit")
				
		
		
		
		gaiaDR2OCTBWithAge=Table.read("GAISDR2COWithAge.fit")
		print "ploting OCs"
		maxAge=   np.log10(1.2e7) #7*np.log10(5.5) #8.5  #7*np.log10(1.2) #1.2e7 

		#maxAge=np.log10(5.5e7) 
		OCcolor='black'
		for eachOC in gaiaDR2OCTBWithAge:
			if eachOC[logt]<maxAge:
				
				l=eachOC["GLON"]
				b=eachOC["GLAT"]

				distance=eachOC["dmode"]/1000.
				
				
				error16= eachOC["d16"]/1000.
				error84= eachOC["d84"]/1000.
				
				
				Rgc,alpha=self.getRgcAlpha(distance,l,b)
				alpha=np.radians(alpha)
				#print distance,eachMaser["plx"] l,alpha
				
				ax.scatter( Rgc*np.cos( alpha ) , Rgc*np.sin( alpha ) ,s= 10,  facecolors='none', edgecolors=OCcolor,lw=0.3,marker="s")
				
				Rgc,alpha1=self.getRgcAlpha(error16,l,b)
				alpha1=np.radians(alpha1)
				
				xError16= Rgc*np.cos( alpha1 )
				yError16= Rgc*np.sin( alpha1 )
	
				
				#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")
	
				Rgc,alpha2=self.getRgcAlpha(error84,l,b)
				alpha2=np.radians(alpha2)
				#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")
				xError84= Rgc*np.cos( alpha2 )
				yError84= Rgc*np.sin( alpha2 )
				
				
				ax.plot([xError16,xError84], [yError16,yError84], lw=0.5,color=OCcolor )
					
					
				#print eachOC['Cluster'], eachOC[logt]
		print "Plotting out arm in the second quandrant" #2016ApJS..224....7D
		 
		#better code the cloud size with mass 
		
		
		outArm2rdTB=Table.read("./data/outArm2rd.fit")
		outArmColor='blue'
		
 
		for eachMC2rd in outArm2rdTB:
			s= eachMC2rd["Mass"]*10 
			

			l=eachMC2rd["_Glon"]
			b=eachMC2rd["_Glat"]
			distance=eachMC2rd["Dist"]
			Rgc,alpha=self.getRgcAlpha(distance,l,b)
			alpha=np.radians(alpha)
			#print distance,eachMaser["plx"] l,alpha
			ax.scatter( Rgc*np.cos( alpha ) , Rgc*np.sin( alpha ) ,s= s,  facecolors='none', edgecolors=outArmColor,lw=0.3,marker=".")
			
	
		
		

		print "Plotting O stars"
		for eachO in OstarTB:
			l=eachO["_Glon"]
			
			#if l>210 and l<220:
				#print  eachMaser
			
			
			b=eachO["_Glat"]
			if eachO["e_plx"]/eachO["plx"]>0.2:
				continue
				
			distance=1./eachO["plx"] 
			error16=1./(eachO["plx"]+ eachO["e_plx"])
			error84=1./(eachO["plx"]- eachO["e_plx"])


			Rgc,alpha=self.getRgcAlpha(distance,l,b)
			alpha=np.radians(alpha)
			#print distance,eachMaser["plx"] l,alpha
			
			ax.scatter( Rgc*np.cos( alpha ) , Rgc*np.sin( alpha ) ,s= 10,  facecolors='none', edgecolors='black',lw=0.3,marker="*")
			Rgc,alpha1=self.getRgcAlpha(error16,l,b)
			alpha1=np.radians(alpha1)
			
			xError16= Rgc*np.cos( alpha1 )
			yError16= Rgc*np.sin( alpha1 )

			
			#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")

			Rgc,alpha2=self.getRgcAlpha(error84,l,b)
			alpha2=np.radians(alpha2)
			#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")
			xError84= Rgc*np.cos( alpha2 )
			yError84= Rgc*np.sin( alpha2 )
			
			
			ax.plot([xError16,xError84], [yError16,yError84], lw=0.5,color='black' )
		
 

			
		print "ploting masers"

					
		#output those core to 
		
		maserTB=Table.read("maserTB.csv")
		
		
		
		
		for eachMaser in maserTB:
			l=eachMaser["_Glon"]
			
			#if l>210 and l<220:
				#print  eachMaser
			
			
			b=eachMaser["_Glat"]
			
			
			if l>210 and l<220:
				
				print l,b,"???????????????????????????????"
			
			
			if eachMaser["e_plx"]/eachMaser["plx"]>0.2:
				continue
				
				
			colorEdge='black'
			if eachMaser["Arm"]=='Loc':
				colorEdge='r'
				
 
			if eachMaser["Arm"]=='Out':
				colorEdge='b'
				
			if eachMaser["Arm"]=='Per':
				colorEdge='g'
			if eachMaser["Arm"]=='Sgr':
				colorEdge='cyan'
				
			if 0:
				np.random.seed()
				distanceSample=np.random.normal( eachMaser["plx"], eachMaser["e_plx"],10000 ) 
				

				
				distanceSample=1./distanceSample  # kpc  
				
				distance=  np.mean(distanceSample)
				error16=np.percentile(distanceSample,16 )
				error84=np.percentile(distanceSample,84 )



			if 1:

				distance=1./eachMaser["plx"] 
				error16=1./(eachMaser["plx"]+ eachMaser["e_plx"])
				error84=1./(eachMaser["plx"]- eachMaser["e_plx"])

			#if distance>7:
 
				#print eachMaser["plx"], eachMaser["e_plx"], error16, distance,error84
			

			Rgc,alpha=self.getRgcAlpha(distance,l,b)
			alpha=np.radians(alpha)
			#print distance,eachMaser["plx"] l,alpha
			
			ax.scatter( Rgc*np.cos( alpha ) , Rgc*np.sin( alpha ) ,s= 10,  facecolors='none', edgecolors=colorEdge,lw=0.6,marker="^")
			
			Rgc,alpha1=self.getRgcAlpha(error16,l,b)
			alpha1=np.radians(alpha1)
			
			xError16= Rgc*np.cos( alpha1 )
			yError16= Rgc*np.sin( alpha1 )

			
			#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")

			Rgc,alpha2=self.getRgcAlpha(error84,l,b)
			alpha2=np.radians(alpha2)
			#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker=".")
			xError84= Rgc*np.cos( alpha2 )
			yError84= Rgc*np.sin( alpha2 )
			
			
			ax.plot([xError16,xError84], [yError16,yError84], lw=0.5,color=colorEdge )
			
			#if distance>7:
			
				#print alpha1,alpha2,alpha
				
				#print xError16, xError84
			
			#ax.scatter(-Rgc*np.sin( alpha ) , Rgc*np.cos( alpha ) ,s= 15,  facecolors='none', edgecolors='purple',lw=0.6,marker="^")

		
		#for eachCore in coreTB[0:40]:
		for eachCore in coreTB:
			
			l=eachCore["Peak1"]
			b=eachCore["Peak2"]
			v=eachCore["Peak3"]/1000.
			
			d,Rgc,alpha=self.ReidA5( l,b,v) 

 
			#calculate mass here
 
			clumpMass,clumpTex,clumpVirial =self.getMassForCoreWithDistance( eachCore,d*1000 )
			
			#print clumpMass,clumpTex,d*1000
			
			#print clumpMass
			
			
			size= clumpMass/10. #eachCore["Sum"]
			
			
			
			#print l,b,v,d,Rgc,alpha
		
			alpha=np.radians(alpha)
			
			colorEdge='b'
			
			if v< self.V_PerseusArm:
				colorEdge='r'
				
			if self.V_PerseusArm <=   v<= self.V_OutArm:
				colorEdge='g'	
				
			if v>self.V_OutArm:
				colorEdge='b'
			
			#print -Rgc*np.sin( alpha ) , Rgc*np.cos( alpha )
			ax.scatter( Rgc*np.cos( alpha ) , Rgc*np.sin( alpha ) ,s= size,  facecolors='none', edgecolors=colorEdge,lw=0.3)
		
		
		ax.set_aspect("equal")
		
		xSun=0.
		ySun=8.34
		#draw sun symbol
		ax.scatter(xSun,ySun,s=1,facecolors="red")
		ax.scatter(xSun,ySun,s=20,facecolors='none', edgecolors='red',lw=0.6)

		#ax.text(xSun+0.2,ySun+0.2, 'Sun')

		
		pSun=np.array( [xSun,ySun  ]) 
		#draw anglues
		
		maxDrawR=8.5
		
		sunRadius=0.12 #
		
		#for drawL in [ 180, 190,200,210,220,230,240,250,260,270]:
		for drawL in [ 30, 60, 90,  120, 150,180,  210, 240,   270,300,330 ]:

		
			#true angles
			drawAngle=np.radians(drawL-90)
			
			unitvector=np.array( [ np.cos(drawAngle), np.sin(drawAngle) ]  )
			drawp1= pSun-unitvector*maxDrawR

			drawp1_end= pSun-unitvector*sunRadius

			 
			ax.plot( [drawp1_end[0],drawp1[0]],   [drawp1_end[1],drawp1[1]] ,'--',color='black',lw=0.2,alpha=0.5 )

			drawp2_start= pSun+unitvector*sunRadius

			drawp2= pSun+unitvector*maxDrawR
			ax.plot( [drawp2_start[0],drawp2[0]],   [drawp2_start[1],drawp2[1]] ,'--',color='black',lw=0.2,alpha=0.5 )

			#ax.plot( [drawp1[0],drawp2[0]],   [drawp1[1],drawp2[1]] ,'--',color='green',lw=0.8,alpha=0.7 )
			
						
			drawp2Text= pSun+unitvector*(maxDrawR )

			#ax.text(drawp2Text[0]-0.6,drawp2Text[1]+0.2, r"{}$^\circ$".format(drawL) )
			
			#ax.text(drawp2Text[0] ,drawp2Text[1]  , r"${}^\circ$".format(drawL),fontsize=6 )

 
			
			if drawL<120 and  drawL>0:
			
			
				ax.text(drawp2Text[0] ,drawp2Text[1]-0.1 , r"${}^\circ$".format(drawL),fontsize=6 )
				
			
			if drawL<270 and  drawL>=120:
			
			
				ax.text(drawp2Text[0] -0.4,drawp2Text[1]+0.1 , r"{}$^\circ$".format(drawL),fontsize=6)
				
			if  drawL>=270:
			
			
				ax.text(drawp2Text[0]-0.6,drawp2Text[1] -0.2, r"{}$^\circ$".format(drawL) ,fontsize=6 )
				
			
				
			#else:
				#ax.text(drawp2Text[0]-0.6,drawp2Text[1]+0.2, r"{}$^\circ$".format(drawL) )

				

		for radiusCircle in [ 2, 4, 6,8]:
			"""
			"""
			#circle=plt.Circle((xSun,ySun),radius=3 )
			
			drawA=np.linspace(0,2*np.pi,200)
			
			x=radiusCircle*np.cos(drawA)
			y=radiusCircle*np.sin(drawA)
			ax.plot( x +xSun ,   y +ySun,'--',color='black',lw=0.2,alpha=0.5 )
			
			if radiusCircle!=2:
				ax.text( radiusCircle*np.cos(np.pi/2.) +xSun+0.1,  radiusCircle*np.sin(np.pi/2.)+ySun, "{} kpc".format( radiusCircle) ,fontsize=6  )

			
			#ax.add_patch(circle)

		ax.scatter(0,0,s=15,facecolors="black")
		ax.text(-2.8,-0.8,"the Galactic Center")

		
		ax.set_xlabel(r"kpc")	
		ax.set_ylabel(r"kpc")	
		



		# fix legend
		if 1:
			#ax.scatter( -100 , 0 ,s= 15,  facecolors='none', edgecolors='g',lw=0.6,marker="^",label='Perseus')
			#ax.scatter( -100 , 0 ,s= 15,  facecolors='none', edgecolors='r',lw=0.6,marker="^",label='Local')
			#ax.scatter( -100 , 0 ,s= 15,  facecolors='none', edgecolors='cyan',lw=0.6,marker="^",label='Sagittarius')
			#ax.scatter( -100 , 0 ,s= 15,  facecolors='none', edgecolors='b',lw=0.6,marker="^",label='Outer')
			
			ax.scatter( -100 , 0 ,s= 20,  facecolors='none', edgecolors='black',lw=0.6,marker="s",label=r'Open Clusters ($<$ 12 Myr)')
			ax.scatter( -100 , 0 ,s= 20,  facecolors='none', edgecolors='black',lw=0.6,marker="^",label='Masers (all colors)')
			ax.scatter( -100 , 0 ,s= 20,  facecolors='none', edgecolors='black',lw=0.6,marker="*",label='O stars')
			ax.scatter( -100 , 0 ,s= 20,  facecolors='none', edgecolors='black',lw=0.6,marker=".",label='Clouds (all colors)')

			#ax.scatter( -100 , 0 ,s= 15,  facecolors='none', edgecolors='b',lw=0.6,marker=".",label='Clouds')
	
			ax.legend(loc=3,fontsize='small',ncol=2)

		ax.set_xlim(-10,10)
		ax.set_ylim(-3,17.5)

		fig.savefig( 'armMap.pdf', bbox_inches="tight")

		
		

	def cropG216FITS(self):
		
		"""
		get fits of maddalena accoding to vRange, lRange, and bRange,
		"""
		
		#we do have to 
		#maddCoreTB=maddCoreTB.loc["Peak1",213.8:219.44]
		
		#maddCoreTB.add_index("Peak2")
		
		#maddCoreTB=maddCoreTB.loc["Peak2", -4.9:-1.25]
		
		#maddCoreTB.add_index("Peak3")
		
		#maddCoreTB=maddCoreTB.loc["Peak3", 16604:38687.4]
		#maddCoreTB=coreTB[0:1]
		#pass
		
		vRange=[14., 38.]
		
		lRange=[213.5,219.5]
		
		bRange=[-5,-1]
		
		myFITS.cropFITS(self.COFile13,outFITS="G216.fits",Vrange=vRange, Lrange=lRange, Brange=bRange)
		
		

	
		
	def getCoreG216(self):
		"""
		"""
		#pass
		
	def crops284FITS(self):
		"""
		get fits of s284 accoding to its vRange, lRange, and bRange,
		"""


		vRange=[35., 53.]
		
		lRange=[213,211.]
		
		bRange=[-2.,0.]
		
		#myFITS.cropFITS(self.COFile13,outFITS="S284.fits",Vrange=vRange, Lrange=lRange, Brange=bRange)

		myFITS.cropFITS(self.COFile12,outFITS="S284_CO12.fits",Vrange=vRange, Lrange=lRange, Brange=bRange)


	def getRMS(self,N,vResolution=0.166,singleRMS=0.5):
		return np.sqrt(N)*vResolution*singleRMS




	def integralS284(self,startV=38 ,endV=52):

		"""
		
		draw over lapping integral map over herschel
		
		
		"""


		saveFigure="herschelOverlap_S284.pdf"


 
		S284_Herschel="S284_Herschel.fits"

		dataCO12,headCO12=myFITS.readFITS(self.S284FITS12)
		vel12 = Velo(headCO12)
		
		
		dataCO13,headCO13=myFITS.readFITS(self.S284FITS13)
		vel13 = Velo(headCO13)
		
		
		
		
		dataHerschel,headHerschel=myFITS.readFITS(S284_Herschel)

		
		fig = plt.figure(1, figsize=(8,4.5) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		
		gh = pywcsgrid2.GridHelper(wcs=WCS(headHerschel))
		gh.locator_params(nbins=3)
		
		g = axes_grid.ImageGrid(fig, 111,
		                        nrows_ncols=(1, 2),
		                        ngrids=None,
		                        direction='row',
		                        axes_pad=0.02, add_all=True,
		                        share_all=True, aspect=True,
		                        label_mode='L', cbar_mode=None,
		                        axes_class=(pywcsgrid2.Axes, dict(grid_helper=gh)))


		ax12=g[0]
		ax13=g[1]
		
		ax12.imshow(dataHerschel,origin='lower',cmap="bone",vmin=25 ,vmax=250,interpolation='none')
		
		for axCO in [ax12,ax13]:
			axCO.imshow(dataHerschel,origin='lower',cmap="bone",vmin=25 ,vmax=250,interpolation='none')
	 
			axCO.set_ticklabel_type("absdeg", "absdeg")
			axCO.axis[:].major_ticks.set_color("w")
			
			
		at = AnchoredText(r"$^{12}\mathrm{CO}~(J=1\rightarrow0)$", loc=4, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		ax12.add_artist(at)

		#at = AnchoredText(r"Herschel 250 $\mu$m", loc=4, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		#ax12.add_artist(at)




		at = AnchoredText(r"$^{13}\mathrm{CO}~(J=1\rightarrow0)$", loc=4, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		ax13.add_artist(at)	
		

		at = AnchoredText(r"Herschel 250 $\mu m$", loc=1, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		ax13.add_artist(at)	
		


		at = AnchoredText(r"[{}, {}] km s$^{{-1}}$".format(startV,endV), loc=3, prop=dict(size=10 ,color="white"),  frameon=False )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		ax12.add_artist(at)	
		
		
		#draw 12
		v0=startV 
		v1=endV
		v0Index=vel12.to_pixel(v0*1000.)
		v1Index=vel12.to_pixel(v1*1000.)
		
		N=v1Index-v0Index
		
		noise12= np.sqrt(N)*self.vResolution12*self.RMS12
		 
		eachMap12=dataCO12[v0Index:v1Index]
		
		eachMap12=np.sum(eachMap12,axis=0)*self.vResolution12 
		
		#print np.sum(dataCO,axis=0).shape
		
		levels=np.linspace( noise12*5,  np.max(eachMap12),10 )

		ax12[WCS(headCO12)].contour(eachMap12,levels=levels, cmap='jet',linewidths=0.4)

		displayLrange=[211.5,212.7]
		displayBrange=[-1.8,-0.5]


		self.setLrange( ax12,WCS(headHerschel),displayLrange)
		self.setBrange( ax12,WCS(headHerschel),displayBrange)


		self.setLrange( ax13,WCS(headHerschel),displayLrange)
		self.setBrange( ax13,WCS(headHerschel),displayBrange)

		
		

		#draw 13
		v0=startV 
		v1=endV
		v0Index=vel13.to_pixel(v0*1000.)
		v1Index=vel13.to_pixel(v1*1000.)
		
		N=v1Index-v0Index
		
		noise13= np.sqrt(N)*self.vResolution12*self.RMS13
		 
		eachMap13=dataCO13[v0Index:v1Index]
		
		eachMap13=np.sum(eachMap13,axis=0)*self.vResolution13
		
		#print np.sum(dataCO,axis=0).shape
		
		levels=np.linspace( noise13*5,  np.max(eachMap13),10 )

		ax13[WCS(headCO13)].contour(eachMap13,levels=levels, cmap='jet',linewidths=0.4)




		#draw the circle by Anderson et al. 2014ApJS..212....1A
		
		wiseCat=Table.read("./data/WISEcat.fit")

		wiseCat.add_index("_Glon")
		
		G210WISE=wiseCat.loc["_Glon",210:215]
		
		G210WISE=G210WISE[G210WISE["Cl"]=="K"]

		for eachWISE in G210WISE:
			if "S284" in eachWISE["HIIName"]:
 
				#draw circle
				
				
			#	print eachWISE["GLON"] 
				#print  eachWISE["GLAT"]
				#print eachWISE["Rad"]
				


				centerL= eachWISE["GLON"] #212.021
				centerB= eachWISE["GLAT"] #-1.309
				radius=   eachWISE["Rad"]/60./60. # deg #1075
				
				
				
				
				angles=np.linspace(0,2*np.pi,50) 
					
				xx=[]
				yy=[]
				for eachAng in angles:
					
					xx.append(centerL+radius*np.cos(eachAng) )
					yy.append(centerB+ radius*np.sin(eachAng)) 

 
 
				ax13["gal"].plot(xx,yy,'r--',lw=0.6)
				ax12["gal"].plot(xx,yy,'r--',lw=0.6)
				break


		#set s y lim

		
		
		#at = AnchoredText(r"{}-{} km s$^{{-1}}$".format(v0,v1), loc=3, prop=dict(size=7 ,color="white"),  frameon=False )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		#g[i].add_artist(at)

		plt.savefig(saveFigure, bbox_inches="tight")




	def S284expand(self,startV=38 ,endV=50,vInterVal=1.,use12=True):
		
		"""
		draw channel maps of  S284
		"""
		saveFigure="S284expand.pdf"
		#using 12CO
		if use12:			
			vReslution=self.vResolution12 
			singleRMS=self.singleChannelRMS12
			#saveFigure="co12_channel_S284Contour.pdf"
			saveFigure="S284expandCO12.pdf"

			
			S284_COFITS="S284_CO12.fits"
 
		else:

			vReslution=self.vResolution13 
			singleRMS=self.singleChannelRMS13
			#saveFigure="co13_channel_S284Contour.pdf"
			S284_COFITS="S284_CO13.fits"
			saveFigure="S284expandCO13.pdf"

		#G216CO12FITS="G216_CO12.fits"
		
		#
		#first 
		
		# read fits
		
		
		#read Halpha=
		HalphaData,HalphaHead=myFITS.readFITS("./data/S284HalphaNew.fits")

		WCSHalpha=WCS( HalphaHead)

		S284_Herschel="S284_Herschel.fits"

		dataCO,headCO=myFITS.readFITS(S284_COFITS)
		vel = Velo(headCO)
		
		dataHerschel,headHerschel=myFITS.readFITS(S284_Herschel)

		
		fig = plt.figure(1, figsize=(9, 7) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		gs = gridspec.GridSpec(4, 2)
		
		#ax=pywcsgrid2.subplot(211, header=WCSHalpha)
		
		ax=pywcsgrid2.subplot(gs[0:2,0], header=WCSHalpha)

		
		# make colorbar
		#ax = g[-1]
		#cax = inset_axes(ax,
		                # width="8%", # width = 10% of parent_bbox width
		                # height="302%", # height : 50%
		                 #loc=3,
		                 #bbox_to_anchor=(1.01, 0, 1, 1),
		                 #bbox_transform=ax.transAxes,
		                # borderpad=0.
		                # )
		
 
		
		vMax=np.max( dataCO )

		vMax=0

		totalN=1


		#BlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlueBlue
		blueVrange=[40,44] #
		v0Index=vel.to_pixel(blueVrange[0]*1000.)
		v1Index=vel.to_pixel(blueVrange[1]*1000.)
		
 	
			
		blueMap=dataCO[v0Index:v1Index]
		blueMap=np.sum(blueMap,axis=0)*vReslution

 

		vInterVal=abs(blueVrange[0]-blueVrange[1]   )
		rmsNoise= self.getRMS(vInterVal/vReslution,singleRMS=singleRMS )

		#levelsB=np.arange(2,40,1)*2*rmsNoise

		levelsB= np.array( [ 4,8,16,32,64,128  ]  ) *rmsNoise #np.arange(2,40,1)*2*rmsNoise


		# RedRedRedRedRedRedRedRedRedRedRedRedRedRedRedRedRedRedRedRedRedRedRedRedRed
		redVrange=[44,52] #
		v0Index=vel.to_pixel(redVrange[0]*1000.)
		v1Index=vel.to_pixel(redVrange[1]*1000.)
		
 	
			
		redMap=dataCO[v0Index:v1Index]
		redMap=np.sum(redMap,axis=0)*vReslution
 
		vInterVal=abs(redVrange[0]-redVrange[1]   )
		rmsNoise= self.getRMS(vInterVal/vReslution,singleRMS=singleRMS )

		#levelsR=np.arange(2,20,1)*2*rmsNoise
		levelsR= np.array( [ 4,8,16,32,64,128  ]  ) *rmsNoise #np.arange(2,40,1)*2*rmsNoise
 
		
		im=ax.imshow(HalphaData,origin='lower',cmap="bone",norm=LogNorm(vmin=2500 ,vmax=6000),interpolation='none')

		#im=g[i].imshow(eachMap,origin='lower',cmap="bone",vmin=rmsNoise ,vmax=rmsNoise*10,interpolation='none')
	#	im=g[i][WCS(headCO)].contour(eachMap,origin='lower',cmap="bone", levels=levels,)
		#im=g[i][WCS(headCO)].contour(eachMap, )

		#im=g[i].imshow(dataHerschel,origin='lower',cmap="bone",vmin=25 ,vmax=250,interpolation='none')
		
		#g[i][WCS(headCO)].contour(eachMap,levels=levels,cmap='jet',linewidths=0.4)
		
		#im=g[i].imshow(eachMap, origin='lower',cmap='bone',vmin=rmsNoise*0.5 ,vmax=vMax*0.5,interpolation='none')
		#im=g[i].contour(eachMap, origin='lower',cmap='bone',vmin=rmsNoise*3 ,vmax=vMax*0.5 )
		im=ax[WCS(headCO)].contour(blueMap, origin='lower', levels=levelsB,  colors='b' ,linewidths=0.5,alpha=0.8)
		
		im=ax[WCS(headCO)].contour(redMap, origin='lower', levels=levelsR,  colors='r' ,linewidths=0.5,alpha=0.6)
		#im=g[i][WCS(headCO)].contour(greenMap, origin='lower', levels=levelsG,  colors='g' ,linewidths=0.5,alpha=0.6)

		#g[i].imshow(eachMap, origin='lower',cmap='bone', interpolation='none')


		
		#at = AnchoredText(r"{}-{} km s$^{{-1}}$".format(v0,v1), loc=3, prop=dict(size=7 ,color="white"),  frameon=False )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		#ax.add_artist(at)
		ax.set_ticklabel_type("absdeg", "absdeg")
		ax.axis[:].major_ticks.set_color("w")
		
		at = AnchoredText(r"$\rm H_{\alpha} \left(SHS\right)$", loc=2, prop=dict(size=12 ,color="limegreen"),  frameon=False )
 
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		ax.add_artist(at) 
 



		#draw the circle by Anderson et al. 2014ApJS..212....1A
		
		wiseCat=Table.read("./data/WISEcat.fit")

		wiseCat.add_index("_Glon")
		
		G210WISE=wiseCat.loc["_Glon",210:215]
		
		G210WISE=G210WISE[G210WISE["Cl"]=="K"]

		S284WISERow=None

		for eachWISE in G210WISE:
			if "S284" in eachWISE["HIIName"]:
 
				#draw circle
				S284WISERow=eachWISE
				break
 
		centerLS284= eachWISE["GLON"] #212.021
		centerBS284= eachWISE["GLAT"] #-1.309
		radiusS284=   eachWISE["Rad"]/60./60. # deg #1075

 
		angles=np.linspace(0,2*np.pi,50) 
			
		xx=[]
		yy=[]
		for eachAng in angles:
			
			xx.append(centerLS284+radiusS284*np.cos(eachAng) )
			yy.append(centerBS284+ radiusS284*np.sin(eachAng)) 

		ax["gal"].plot(xx,yy,'--',lw=1.,color='limegreen')

		drawRadius=radiusS284+0.05
		
		self.setAxLimX(ax,   WCSHalpha , [centerLS284-drawRadius ,centerLS284+drawRadius] )
		self.setAxLimY(ax,   WCSHalpha , [centerBS284-drawRadius, centerBS284+drawRadius ] )




		#draw spectral AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa
		
		axSpec1=plt.subplot(gs[0,1]  )
		
		PositionA=[211.9844825, -1.1292447]
		self.drawSpecByLB( ax,axSpec1,PositionA,"A",dataCO,headCO, blueVrange, redVrange) 


		#draw spectral BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

		axSpec2=plt.subplot(gs[1,1] ,sharex=axSpec1,sharey=axSpec1 )
		
		PositionB=[ 212.0655977,-1.2434648] 
		
		self.drawSpecByLB( ax,axSpec2,PositionB,"B",dataCO,headCO, blueVrange, redVrange) 

		
 

 
		#pistal   CCCCCCCCCCCCCC
		axSpec3=plt.subplot(gs[2,1] ,sharex=axSpec1,sharey=axSpec1 )

		PositionC=[212.1915751, -1.3425466]
		self.drawSpecByLB( ax,axSpec3,PositionC,"C",dataCO,headCO, blueVrange, redVrange) 

		
		#pistal   DDDDDDDDDDDD
		axSpec4=plt.subplot(gs[3,1] ,sharex=axSpec1,sharey=axSpec1 )

		PositionD=[ 211.9451520,-1.4422921]
		self.drawSpecByLB( ax,axSpec4,PositionD,"D",dataCO,headCO, blueVrange, redVrange) 

		#212.0655977
		axSpec4.set_ylim(-1.5,6)



          
		#draw PV diagram
		pvFITSS284="S284PV.fits"
		
		if 1:
			

			centerLS284= eachWISE["GLON"] #212.021
			centerBS284= eachWISE["GLAT"] #-1.309
			radiusS284=   eachWISE["Rad"]/60./60. # deg #1075
				
			#centerLB= [212.0046000, -1.3174597]
			centerLB=[centerLS284,centerBS284 ]
			
			angle=44.
			#theta=np.radians( angle )  
			width= radiusS284*2*0.9 #0.5 
			length=radiusS284*2*1.2 #0.7
			
		if 0:
			centerLB =[ 212.0655977,-1.2434648] 
		
			angle=135.
			#theta=np.radians( angle )  
			width=0.2  
			length=0.7
		
		 
		
		if 1:
		
			doMadd.getPVFitsByLB(self.COFile12,centerLB ,width,length,angle,pvFITSS284) 
			#doMadd.getPVFitsByLB( "./data/S284MaskFITSCO12.fits",centerLB ,width,length,angle,pvFITSS284) 

			
		#draw pvRange
          
		XY,HW,angleRect,arrowStart,arrowEnd=self.getRectangleByWCS(HalphaHead  ,centerLB,length,width, angle) # the input angle is the pixel angle
		rect=Rectangle((XY[0], XY[1]), HW[1],HW[0] , angle=angleRect, facecolor='None', edgecolor = 'cyan', linestyle='--', linewidth=0.8)
		ax.add_patch( rect)
		
		ax.arrow(arrowStart[0],arrowStart[1],arrowEnd[0]-arrowStart[0],arrowEnd[1]-arrowStart[1] ,   head_width=60, fc='cyan',ec='cyan',lw=1., length_includes_head=True, capstyle="projecting"  )

		pvData,pvHead=myFITS.readFITS(pvFITSS284)

		
		axPV=pywcsgrid2.subplot(gs[2:,0], header=WCS(pvHead))
 
		
		disPlayVRange=[39,50] #km /s
		
		#draw AxPv
		
		axPV.axis[:].major_ticks.set_color("w") 
		
		pvRMS=0.08
		maxSlice=np.max( pvData )
		

		contLevels= np.arange( pvRMS*3,maxSlice ,pvRMS*2 )
		
		contLevels=list(contLevels)
		
		contLevels.insert(0,np.min( pvData) )
		
		
		#axPV.contourf(  pvData ,levels=contLevels,cmap='bone' )
		axPV.imshow(  pvData ,origin="lower", cmap='jet' )
		axPV.contour (  pvData ,levels=contLevels,cmap='jet',linewidths=0.5 )

		pvNx,pvNy=pvData.shape


			# modify the labels
		minV,maxV=disPlayVRange
		
		pvWCS=WCS(pvHead)
		
		a,minY =pvWCS.wcs_world2pix(0,minV*1000,  0)
		a,maxY =pvWCS.wcs_world2pix(0,maxV*1000,  0)


			
		if 1:
		
			vInterval=2 #km/s
			
			yLocs=np.arange(int(minV) , int(maxV)+vInterval,vInterval)
			
			yLabels= map(int, yLocs)
			yLabels=map(str,yLabels)
			
			#print yLocs*1000.
			#print yLabels
			
			
					 
			axPV.set_ticklabel2_type("manual",locs=yLocs*1000.,labels=yLabels)
			#axPV.set_ticklabel2_type("manual",locs=[35000.],labels=["35000"])
	
		axPV.set_ylim( minY, maxY)			

 


		axPV.set_ylabel(r"$\rm V_{LSR} \left(km \ s^{-1}\right)$")
		
		axPV.set_xlabel(r"$\rm Offset\ \left(degree\right)$")


		
		
		
                                                       
		#axAspec.set_xlabel(r"Velocity ($\rm km \ s^{-1}$)")
		#axAspec.set_ylabel(r"T$_{\rm mb}$ (K)")            

		fig.tight_layout()
		plt.savefig(saveFigure, bbox_inches="tight" )
		plt.savefig("S284expandCO12.png", bbox_inches="tight",dpi=300 )

 

	def drawSpecByLB(self,axImage,axSpec,LB,pName,dataCO,headCO,blueV,redV):
		
		l1,b1=LB
		spec1,velo1=myFITS.getSpectraByLB(dataCO,headCO,l1,b1)
		axSpec.step(velo1,spec1,lw=0.8, color='black' )
		
		
		hatchLW=0.2
		
		plt.rcParams["hatch.linewidth"]=hatchLW

		axSpec.fill_between(blueV, [7,7] , -2 , facecolor='none' ,hatch='/',edgecolor='b', linewidth=hatchLW)
		axSpec.fill_between(redV, [7,7] , -2 , facecolor='none' ,hatch='\\',edgecolor='r', linewidth=hatchLW)

		
		axImage["gal"].plot(l1,b1,'.',color='limegreen')
		axImage["gal"].text(l1-0.01,b1-0.01,pName,color='limegreen')

		#at = AnchoredText(pName, loc=1, prop=dict(size=12 ,color="black"),  frameon=False )
		at = AnchoredText(pName, loc=1,  frameon=False )

			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axSpec.add_artist(at)
		axSpec.axhline(  linestyle='--',linewidth=0.5,color='black' )
		axSpec.set_ylabel(r"$\rm T_{\rm mb} \left(K\right)$")            
		#axSpec.set_ylabel(r"$\rm T_{\rm mb} \left(K\right)$")            
		axSpec.set_xlabel(r"$\rm V_{LSR} \left(km \ s^{-1}\right)$")
		#axSpec3.set_xlabel(r"$\rm V_{LSR} \left(km \ s^{-1}\right)$")
          

	def 		 setAxLimX(self,ax,wcs,lRange):
		"""
		"""


		leftL= max(lRange)
		rightL= min(lRange)
 
		leftN,aa=wcs.wcs_world2pix(leftL,0,0)
		rightN,aa=wcs.wcs_world2pix(rightL,0,0)

		
		ax.set_xlim(leftN,rightN)


	def 		 setAxLimY(self,ax,wcs,bRange):
		"""
		"""


		bottomB= min(bRange)
		topB= max(bRange)
 
		aa,lowerN=wcs.wcs_world2pix(212., bottomB, 0)
		aa,upN=wcs.wcs_world2pix(212., topB, 0)

		
		ax.set_ylim(lowerN,upN )




	def channelS284(self,startV=38 ,endV=50,vInterVal=1.,use12=True):
		
		"""
		draw channel maps of  S284
		"""
		saveFigure="co12_channel_S284.pdf"
		#using 12CO
		if use12:			
			vReslution=self.vResolution12 
			singleRMS=self.singleChannelRMS12
			saveFigure="co12_channel_S284.pdf"
			S284_COFITS="S284_CO12.fits"
 
		else:

			vReslution=self.vResolution13 
			singleRMS=self.singleChannelRMS13
			saveFigure="co13_channel_S284.pdf"
			S284_COFITS="S284_CO13.fits"
	
		#G216CO12FITS="G216_CO12.fits"
		
		#
		#first 
		
		# read fits
		
		

		S284_Herschel="S284_Herschel.fits"

		dataCO,headCO=myFITS.readFITS(S284_COFITS)
		vel = Velo(headCO)
		
		dataHerschel,headHerschel=myFITS.readFITS(S284_Herschel)

		
		fig = plt.figure(1, figsize=(9, 7) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		
		#gh = pywcsgrid2.GridHelper(wcs=WCS(headHerschel))
		
		gh = pywcsgrid2.GridHelper(wcs=WCS(headCO))

		
		gh.locator_params(nbins=3)
		
		g = axes_grid.ImageGrid(fig, 111,
		                        nrows_ncols=(3, 4),
		                        ngrids=None,
		                        direction='row',
		                        axes_pad=0.02, add_all=True,
		                        share_all=True, aspect=True,
		                        label_mode='L', cbar_mode=None,
		                        axes_class=(pywcsgrid2.Axes, dict(grid_helper=gh)))
		
 
		# make colorbar
		ax = g[-1]
		cax = inset_axes(ax,
		                 width="8%", # width = 10% of parent_bbox width
		                 height="302%", # height : 50%
		                 loc=3,
		                 bbox_to_anchor=(1.01, 0, 1, 1),
		                 bbox_transform=ax.transAxes,
		                 borderpad=0.
		                 )
		
		rmsNoise= self.getRMS(vInterVal/vReslution,singleRMS=singleRMS )
		
		ims=[]
		
		vMax=np.max( dataCO )

		vMax=0

		totalN=12

		for i in range(totalN):
			
			v0=startV+i*vInterVal
			v1=v0+vInterVal
			v0Index=vel.to_pixel(v0*1000.)
			v1Index=vel.to_pixel(v1*1000.)
 
			eachMap=dataCO[v0Index:v1Index]
			
			eachMap=np.sum(eachMap,axis=0)*vReslution
			if  np.max( eachMap )>vMax:
				vMax=np.max( eachMap )
		
		levels=np.linspace( rmsNoise*5,  vMax,8 )
		
		print levels
		for i in range(totalN):
			
			v0=startV+i*vInterVal
			v1=v0+vInterVal
			v0Index=vel.to_pixel(v0*1000.)
			v1Index=vel.to_pixel(v1*1000.)
 
			eachMap=dataCO[v0Index:v1Index]
			
			eachMap=np.sum(eachMap,axis=0)*vReslution
 			
			#im=g[i].imshow(eachMap,origin='lower',cmap="bone",vmin=rmsNoise ,vmax=rmsNoise*10,interpolation='none')
		#	im=g[i][WCS(headCO)].contour(eachMap,origin='lower',cmap="bone", levels=levels,)
			#im=g[i][WCS(headCO)].contour(eachMap, )

			#im=g[i].imshow(dataHerschel,origin='lower',cmap="bone",vmin=25 ,vmax=250,interpolation='none')
			
			#g[i][WCS(headCO)].contour(eachMap,levels=levels,cmap='jet',linewidths=0.4)
			
			im=g[i].imshow(eachMap, origin='lower',cmap='bone',vmin=0 ,vmax=vMax*0.3,interpolation='none') #good


			#im=g[i].imshow(eachMap, origin='lower',cmap='bone',norm=LogNorm(vmin=rmsNoise ,vmax=vMax*0.4),interpolation='none')


			#g[i].imshow(eachMap, origin='lower',cmap='bone', interpolation='none')


			ims.append(im)

		
			at = AnchoredText(r"{}-{} km s$^{{-1}}$".format(v0,v1), loc=4, prop=dict(size=7 ,color="white"),  frameon=False )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			g[i].add_artist(at)
			g[i].set_ticklabel_type("absdeg", "absdeg")
			g[i].axis[:].major_ticks.set_color("w")
		
		
			#draw circle 
			
			if 1:
				centerL=  212.021
				centerB=  -1.309
				radius=  1075./60./60. # deg #
 
				angles=np.linspace(0,2*np.pi,50) 
					
				xx=[]
				yy=[]
				for eachAng in angles:
					
					xx.append(centerL+radius*np.cos(eachAng) )
					yy.append(centerB+ radius*np.sin(eachAng)) 
	
 
				g[i]["gal"].plot(xx,yy,'r--',lw=0.2)
				g[i]["gal"].plot(xx,yy,'r--',lw=0.2)
	
		
		
		if use12:
	 
			at = AnchoredText(r"$^{12}\mathrm{CO}~(J=1\rightarrow0)$", loc=2, frameon=False,prop={"color":"w"} )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			g[0].add_artist(at)
		
		else:
			at = AnchoredText(r"$^{13}\mathrm{CO}~(J=1\rightarrow0)$", loc=2, frameon=False,prop={"color":"w"} )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			g[0].add_artist(at)
			
		
		
		cb = plt.colorbar(ims[6], cax=cax )
		cb.set_label("K km s$^{-1}$")
		#cb.set_ticks([0, 1, 2, 3])

		#set range

		displayLrange=[211.5,212.7]
		displayBrange=[-1.8,-0.5]


		self.setLrange( g[0],WCS(headCO),displayLrange)
		self.setBrange( g[0],WCS(headCO),displayBrange)



		
		plt.savefig(saveFigure, bbox_inches="tight")



	def channelG216(self,startV,endV,vInterVal=2.,use12=True):
		
		"""
		draw channel maps of 
		"""

		#using 12CO
		
		#G216CO12FITS="G216_CO12.fits"
		
		use12=False
		
		if use12:			
			vReslution=self.vResolution12 
			singleRMS=self.singleChannelRMS12
			G216CO12FITS= "G216_CO12.fits"

 
		else:

			vReslution=self.vResolution13 
			singleRMS=self.singleChannelRMS13
			G216CO12FITS="G216_13.fits"
	

	 
		
		#noiseRMS=np.sqrt(vInterVal/vReslution)*


		data,head=myFITS.readFITS(G216CO12FITS)
		vel = Velo(head)
		 
		fig = plt.figure(1, figsize=(9, 7), dpi=70)
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		


		gh = pywcsgrid2.GridHelper(wcs=WCS(head))
		gh.locator_params(nbins=3)
		
		g = axes_grid.ImageGrid(fig, 111,
		                        nrows_ncols=(4, 3),
		                        ngrids=None,
		                        direction='row',
		                        axes_pad=0.02, add_all=True,
		                        share_all=True, aspect=True,
		                        label_mode='L', cbar_mode=None,
		                        axes_class=(pywcsgrid2.Axes, dict(grid_helper=gh)))
		
		# make colorbar
		ax = g[-1]
		cax = inset_axes(ax,
		                 width="8%", # width = 10% of parent_bbox width
		                 height="403%", # height : 50%
		                 loc=3,
		                 bbox_to_anchor=(1.01, 0, 1, 1),
		                 bbox_transform=ax.transAxes,
		                 borderpad=0.
		                 )

		rmsNoise= self.getRMS(vInterVal/vReslution,singleRMS=singleRMS )
		print rmsNoise,"?????????????????"
		ims=[]
		
		for i in range(12):
			
			v0=startV+i*vInterVal
			v1=v0+vInterVal
			v0Index=vel.to_pixel(v0*1000.)
			v1Index=vel.to_pixel(v1*1000.)
 
			eachMap=data[v0Index:v1Index]
			
			eachMap=np.sum(eachMap,axis=0)*vReslution
 			
 			
 			#use LogNormal
			#cmap = plt.cm.gray_r
			#import matplotlib.colors as mcolors
			#norm = mcolors.Normalize()
		
			#im=g[i].imshow(eachMap,origin='lower',  norm=norm, cmap="jet")

			im=g[i].imshow(eachMap,origin='lower',cmap="bone",vmin=rmsNoise*0.5 ,vmax=rmsNoise*8,interpolation='none')
			#sim=g[i].imshow(eachMap,origin='lower',cmap="bone",norm=LogNorm(vmin= rmsNoise  , vmax=rmsNoise*10),  interpolation='none')

 
			
			
			ims.append(im)
			at = AnchoredText(r"{}-{} km s$^{{-1}}$".format(v0,v1), loc=3, prop=dict(size=8 ,color="white"),  frameon=False )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			g[i].add_artist(at)
			g[i].set_ticklabel_type("absdeg", "absdeg")
			g[i].axis[:].major_ticks.set_color("w")
		# make colorbar
		cb = plt.colorbar(ims[4], cax=cax)
		cb.set_label("K km s$^{-1}$")
		#cb.set_ticks([0, 1, 2, 3])
		

		at = AnchoredText(r"$^{13}\mathrm{CO}~(J=1\rightarrow0)$", loc=1, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		g[2].add_artist(at)



		plt.savefig("co13channelG216.pdf", bbox_inches="tight")

		#plt.savefig("co13_channel_G216.pdf", bbox_inches="tight")
		
		
	def outPutDiskStars(self,TB,outFile="GAIAFILE.txt"):
		f = open(outFile,"w") 
 
		
		for eachYSO in  TB: #self.YSOTB:
			
			l= eachYSO["_Glon"]
			b=eachYSO["_Glat"]
 		
			c = SkyCoord(l=l*u.degree, b=b*u.degree, frame='galactic')
			cor=str(c.icrs.to_string( 'decimal'))
			if c.icrs.dec<0:
			
				outStr=  cor.replace(" ","")
			else:
				outStr=  cor.replace(" ","+")

 
			f.write(outStr+"\n")
		 
		f.close() 
		
		
		
		
	def calRMS(self,FITSName):
		
		"""
		
		calculate the RMS with truncated Gaussion
		"""




	def examineOneCore(self,coreName):
		
		"""
		should return the row
		"""
 
		coreTBName = coreName + "Cat_1.fit"
		coreFITSName = coreName + "_1.fits"
 
		coreTB=Table.read(coreTBName)
  
		coreCMF=self.getMassFromGauSum(coreTB, coreName) 

		coreCMF[myCMF.peakMin]=np.min(coreTB["Peak"])


		return coreCMF



	def checkAllCores(self,testRandom=True,source="G216",saveTBName="",startI=0,endI=50):
		
		"""
		"""
		recoreCMF=myCMF(source)

		returnTB=recoreCMF.getEmptyTB()

		
		if testRandom: #if test random cores
			#coreTests=[]
			
			#for i in range(1):
				
			for i in np.arange(startI,endI ):
				caseName="./gaussClumpTest/"+"RANDOM_{}_{}GaussClumpscore".format(source,i)
				#coreTests.append(saveName)
		
			#for   coreTest in coreTests:
				
			#for   coreTest in coreTests[0:5]: #for te
				
				print "Proceeding ",caseName,"..."
				coreInfo=self.examineOneCore(caseName)
				#		return  [coreN,totalMass,  [alpha1Mean, alpha1_std ], [alpha2Mean, alpha2_std ],  [MtMean, Mt_std]  ]   # massList
	
				residualRMS=self.examineCOREFITS(caseName) 
				#print "residual RMS",residualRMS
				
				coreInfo[myCMF.peakMin]= coreInfo[myCMF.peakMin]/self.RMS13
				coreInfo[myCMF.residual]=  residualRMS/self.RMS13
				coreInfo[myCMF.case]= i

				returnTB.add_row(coreInfo)
			
			
			#write returnTB
			
			
			os.system("rm "+saveTBName)
			returnTB.write(saveTBName,format="fits")
			#get mean and std with a list of mean and std
	 
			#self.getMeanAndStd(alpha2Mean,alpha2Std) 
 
		#coreTest4="TEST_THRESH6_0GaussClumpscore"
 
		#calculate excitation temperatures
		if 0:
			for TBName in coreTests:
				
				coreTBTest=Table.read( TBName+"Cat_1.fit" )
				TexList=doMadd.getTexListFromTB(coreTBTest)
				
			
		
 
		

			
			


	def getMeanAndStd(self,meanList,stdList):
		
		"""
		"""
		
		pass
		
		
		
		
		
		

	def mathCTBWithLB(self,TB,LB):
		
		"""
		"""
		
		disMin=10.
		
		pl,pb=LB
		rowMatch=TB[0]
		for eachRow in TB:
			l=eachRow["Peak1"]	
			b=eachRow["Peak2"]		
			distance=np.sqrt( (l-pl)**2 + (b- pb)**2  )
			
			if distance<disMin:
				
				rowMatch=eachRow
				disMin=distance
				
		return 	rowMatch
			


	def drawRingPV(self,fitsName,centerLB,r1,r2,saveFITS="ring.fits"):
		"""
		all in degree
		
		"""
		
		#ring the PPV diagram along the circle
		
		data,head=myFITS.readFITS(fitsName)
		
		resolutionPix = abs(head["CDELT2"] )
		
		ring_wcs=WCS(head)

		widthDeg=r2-r1

		r0=np.mean([r1,r2]) 

		R = np.arange(0, 2*np.pi, 0.1)
 
		xx = centerLB[0] + r0*np.cos(R)  #latitude
		yy = centerLB[1] + r0*np.sin(R)  #longtitude

		ring_x,ring_y,dd= ring_wcs.all_world2pix(xx,yy,45000.,0)


		pathPoints=[]
		
		
		for i in range( len(ring_x)):
			
			pathPoints.append( [ ring_x[i], ring_y[i] ] )

 
 
		widthPix= widthDeg/resolutionPix
		
		
		
		print xx
		from pvextractor import extract_pv_slice,Path
 
		xy = Path(pathPoints,width= widthPix )

		#xy=Path( [ [315.15652, 315.15652], [ 388.05878, 320.25191] ], widthPix )
		
		ringHDU= fits.open(fitsName)[0]
		
		
		#ringHDU= fits.open( "G216_CO12.fits")[0]

		pv = extract_pv_slice(   ringHDU, xy)
		
 

		os.system("rm "+saveFITS)
		pv.writeto(saveFITS)
		




	def fittingEliptical(self,vRange):
		
		"""
		moment the vRange, find the cores by hand, and fitting the eliptical ring
		
		this ring confirms the existing of shells in G216, which is not an important discovery, but what kind of thing do we have to display in the paper?
		"""

		#doing moment
		outRingFITS="ringMoment.fits"
		outRingFITSM1="ringMoment_1.fits"

		
		
		ringCoreTBName="ringCoresGaussClumpscoreCat_1.fit"
		
		ringTB=Table.read(ringCoreTBName)
		
		if 0:
		
			#ring_Data,ring_head=doMadd.doFITS.momentFITS( self.G216FITS, vRange,0,outFITS=outRingFITS)
	
			RMS=0.20
			
			#
			from starlinkTest import myCore
			doCore=myCore()
			doCore.testCoreMethod("ringMoment.sdf",doCore.GaussClumps, RMS,useConfig=True,saveName="ringCores")


		if 0: #moment 1
			doMadd.doFITS.momentFITS( self.G216FITS, vRange,1,outFITS=outRingFITSM1)

			
			
		#get the cores
		
		
		# Open the file with read only permit
		f = open('ringCorePosition.txt', "r")
		# use readlines to read all lines in the file
		# The variable "lines" is a list containing all lines in the file
		lines = f.readlines()
		# close the file after reading the lines.
		f.close()
				
		
		peaks_L=[]
		peaks_B=[]
		
		
		
		onRingCores=ringTB[0:1]
		onRingCores.remove_row(0)
		
		
		#
		
		for eachLine in lines:
			
			lStr,bStr= eachLine.split(',')
			LB=[ float(lStr),float(bStr)]
			matchedCore=self.mathCTBWithLB(  ringTB,LB )
			
			
			onRingCores.add_row(matchedCore)
			
			#peaks_L.append(matchedCore["Peak1"]   )
			#speaks_B.append(matchedCore["Peak2"]   )



		#peaks_L=onRingCores["Peak1"].data  #np.array(peaks_L)
		#peaks_B= onRingCores["Peak2"].data  #np.array(peaks_B)

		peaks_L=onRingCores["Cen1"].data  #np.array(peaks_L)
		peaks_B= onRingCores["Cen2"].data  #np.array(peaks_B)



		#print onRingCores["Peak3"]

		#get the velocity
		fig, ax = plt.subplots()
		
		maxVList=[]
		
		for eachRingCore in onRingCores:
			#print eachRingCore
			
			#according to size average spectrum
			#seachRingCore["Size1"],eachRingCore["Size2"]
			
			l=eachRingCore["Peak1"]
			b=eachRingCore["Peak2"]


			#l=eachRingCore["Cen1"]
			#b=eachRingCore["Cen2"]

			data,header=myFITS.readFITS(self.G216FITS) #
			
			spec,velo=myFITS.getSpectraByLB(data,header,l,b)
			

 
			
			#ax.plot([vRange[0],vRange[0]],[0,1],'r--')
			#ax.plot([vRange[1],vRange[1]],[0,1],'r--')
			specP=spec[velo>vRange[0]]
			veloP=velo[velo>vRange[0]]



 
			
			specPP=specP[veloP<vRange[1]]
			veloPP=veloP[veloP<vRange[1]]
			
			
			if 1: #usepeak
			
				maxCoreV=veloPP[specPP.argmax() ]
			if 0: #use average
				maxCoreV= np.sum(veloP*specP)/np.sum(specP)
 
				
				
			# it seems the sensitivity cannot detect the velocity correctely
			
			maxVList.append(maxCoreV)
			#print maxCoreV,"????????????"
			ax.scatter(l,b,c='r')
			ax.text(l,b,"{:.2f}".format(maxCoreV) )

		RMSV= max(maxVList)-min(maxVList) #np.std(maxVList,ddof=1  )
		print RMSV,"Vdispersion"

		
		#plt.show() 
		
		
		from ellipse import myEllipse

		#from ellipse import myEllipse
		doE=myEllipse()
		doE.doFitAll(peaks_L,peaks_B)
				
		#fig = plt.figure(figsize=(6,6))
		#ax = fig.add_subplot(111)
		#ax.axis('equal')
		
		ringData,ringHead=myFITS.readFITS(outRingFITS)
		
		ax=pywcsgrid2.subplot(111,header=WCS(ringHead))
		
		rms=0.16
		
		ax.imshow(ringData[0],origin='lower',cmap="bone",vmin=rms ,vmax=rms*10,interpolation='none')
		
		#ax.imshow(bData, origin="lower")
		#from astropy.visualization import make_lupton_rgb
		#image = make_lupton_rgb(rData, gData, bData, stretch=0.5)
		#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		ax.set_ticklabel_type("absdeg", "absdeg")
		
		
		
		
		#ax.plot(peaks_L, peaks_B, 'bo', label='test data', zorder=1)




		centerE =doE.centerE
		#phi = ellipse_angle_of_rotation(a)
		phi =  doE.angle 
		axes =doE.axes
		print phi*180/np.pi,"Angle?"

		print "Center coordinates",centerE
		
		print "axes in degree: ",axes
		
		arc = 4
		R = np.arange(0,arc*np.pi, 0.1)
		b,a = axes
		xx = centerE[0] + a*np.cos(R)*np.cos(phi) - b*np.sin(R)*np.sin(phi)
		yy = centerE[1] + a*np.cos(R)*np.sin(phi) + b*np.sin(R)*np.cos(phi)

		ring_wcs=WCS(ringHead)


		ring_x,ring_y,dd= ring_wcs.all_world2pix(xx,yy,0,0)


		pathPoints=[]
		
		
		for i in range( len(ring_x)):
			
			pathPoints.append( [ ring_x[i], ring_y[i] ] )

 
		#get pv diagrame 
		widthPix=4
		from pvextractor import extract_pv_slice,Path
 
		xy = Path(pathPoints,width= widthPix )

		#xy=Path( [ [315.15652, 315.15652], [ 388.05878, 320.25191] ], widthPix )
		
		ringHDU= fits.open(self.G216FITS)[0]
		
		
		#ringHDU= fits.open( "G216_CO12.fits")[0]

		pv = extract_pv_slice(   ringHDU, xy)
		
		ringPV="ringCircle.fits"

		os.system("rm "+ringPV)
		pv.writeto(ringPV)
		
		



		print "inclination angle",np.arccos(b/a)*180/np.pi
		print  RMSV*a/b,"The true velocity?"

		expandV=RMSV*a/b
		
		T=  np.sin(np.radians(2*a))*2200*3.0857e16/1000./expandV
		
 
		T=T/60./60./24/365 #year
		
		print  T/1e6,"The age Myr?"

		ax["gal"].plot(xx,yy,'r--',lw=0.5,alpha=0.6)
		
		ax["gal"].plot(centerE[0] ,centerE[1],'ro',markersize=1,alpha=0.6)
		
		
	
		at = AnchoredText(r"$^{13}\mathrm{CO}~(J=1\rightarrow0)$ " , loc=2, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		ax.add_artist(at)
		at = AnchoredText(r"$ [{},\ {}] \ \rm km \ s^{{-1}}$".format(vRange[0],vRange[1]), loc=3, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		ax.add_artist(at)
		#at = AnchoredText(r"$^{13}\mathrm{CO}~(J=1\rightarrow0)$", loc=1, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		#ax.add_artist(at)
		ax.axis[:].major_ticks.set_color("w")
		
		
		ax.set_xlim(100,630)
		ax.set_ylim(84,444)	
		
		#dx=a
		#dy=dx*np.tan(phi)
		
	#	X1=center[0]+dx
		
		#Y1=center[1]+dy
		
		#ax["gal"].plot([center[0], X1 ] , [ center[1],Y1 ] )


		#display    vRange



		if 0:
		
			ax.plot(xx,yy, color = 'red')
			
			ax.plot(center[0] ,center[1],'ro')
			
			dx=a
			dy=dx*np.tan(phi)
			
			X1=center[0]+dx
			
			Y1=center[1]+dy
			
			ax["gal"].plot([center[0], X1 ] , [ center[1],Y1 ] )
	
	
			dx=-a
			dy=dx*np.tan(phi)
			
			X1=center[0]+dx
			
			Y1=center[1]+dy
			
			ax.plot([center[0], X1 ] , [ center[1],Y1 ] )
			
	
			dx=b
			dy=dx*np.tan(phi-np.pi/2)
			
			X1=center[0]+dx
			
			Y1=center[1]+dy
			
			ax.plot([center[0], X1 ] , [ center[1],Y1 ] )
	
	
			dx=-b
			dy=dx*np.tan(phi-np.pi/2)
			
			X1=center[0]+dx
			
			Y1=center[1]+dy
			
			ax.plot([center[0], X1 ] , [ center[1],Y1 ] )
			
	
			
			plt.legend()
		
		
		
		
		#plt.show()
		plt.savefig( 'ringG216.pdf', bbox_inches="tight")		
		
	def calDisByXfactor(self,coInt,distance,xFactor=2.0e20):
		
		"""
		The unit of coInt must be K km/s, 
		distance pc
		not calDis, calMass
		"""
		
		NH2= coInt *xFactor # cm-2 #
		

		#distance=2200 # pc
		
		#parsecToMeter= 3.0857e16 #m 
		#length1=np.radians(degSize1)*self.distance*self.parsecToMeter*100. #cm 
		#length2=np.radians(degSize2)*self.distance*self.parsecToMeter*100.

		#should use single pix for 12CO
		
		length1=np.radians(30./60./60.)*distance*self.parsecToMeter*100. #cm 
		#length2=np.radians(degSize2)*self.distance*self.parsecToMeter*100.

		mu=1.36
		
		Mh2=3.35e-27 #kg
		solarMass=1.9891e30 #kg
		#s=np.pi*length1*length1
		s= length1*length1

		coreSolar= s*NH2*Mh2*mu/solarMass
		return coreSolar
		
	def massOfCo12Core(self,coreRow,xFactor=2.0e20,distance=2350):
		"""
		The default Xfactor about 2.2e20 cm-2, according to 2015ApJ...803...38I , 
		
		For eachCore, we calculate the mass. better use fellwalker
		
		"""
		
		Cen1="Cen1"
		Cen2="Cen2"
		Cen3="Cen3"
		
		
		Size1="Size1"
		Size2="Size2"
		Size3="Size3"
		
		Sum="Sum"
		
		colnames=['PIDENT',
							'Peak1',
							'Peak2',
							'Peak3',
							'Cen1',
							'Cen2',
							'Cen3',
							'Size1',
							'Size2',
							'Size3',
							'Sum',
							'Peak',
							'Volume']
		
		
		degSize1=coreRow[Size1]/60./60.
		
		degSize2=coreRow[Size2]/60./60. 
		
		pixSize=30./3600.
 
		return self.calDisByXfactor( coreRow[Sum]*0.1587376445530, distance, xFactor=xFactor)
		


		return coreSolar
			

	def calMassWithCO12(self,TBName):
		"""
		"""
		
		#calculate the mass
		
		#calculate the area
		
		TB=Table.read(TBName )
		
		massList=[]
		for eachRow in TB:
			coreMassSingle= self.massOfCo12Core(eachRow)
			massList.append(coreMassSingle)
			
			
		print "total Mass",np.sum(massList)
		
	
	def compareLTEVirial(self,TBName):
		
		coreTB=Table.read(TBName )
		massList=[]
		MvirialList=[]
		alphaVirialList=[]
		
		for coreRow in  coreTB:
  
			coreSolar,Tex,MvirialAndAlpha=self.getMassForCore(coreRow )
			massList.append( coreSolar ) #,Tex
			Mvirial,alphaVirial=MvirialAndAlpha
			MvirialList.append( Mvirial )
			alphaVirialList.append(alphaVirial)
			#if coreSolar>Mvirial:
				
		massList=np.array(massList)
		MvirialList=np.array(MvirialList)
		alphaVirialList=np.array(alphaVirialList)

		a=massList[alphaVirialList<2]
		print len(a),len(coreTB),len(a)*1./len(coreTB),np.sum(a)/np.sum(massList) 
		
		print "total mass",  np.sum(massList)  , "mass gravitation bound", np.sum(a)
		
		
		
	def setLrange(self,ax,wcs,lRange):
		
		"""
		2D
		
		"""
		
		try: #2D
		
			
			index1,dddd=wcs.wcs_world2pix( lRange[0],0, 0)	
			
			index2,dddd=wcs.wcs_world2pix( lRange[1],0, 0)		
		except: #3D
			index1,dddd,oo=wcs.wcs_world2pix( lRange[0],0,0, 0)	 
			
			index2,dddd,oo=wcs.wcs_world2pix( lRange[1],0,0, 0)		
			


		ax.set_xlim(min(index1,index2),max(index1,index2))

		
	def setBrange(self,ax,wcs,bRange):
		
		"""
		2D
		
		"""
		try: #2D
		
			dddd,index1 =wcs.wcs_world2pix( 212., bRange[0],  0)	
			
			dddd,index2=wcs.wcs_world2pix(212., bRange[1],  0)		
		except:
			dddd,index1,oo=wcs.wcs_world2pix( 212., bRange[0],0,  0)	
			
			dddd,index2,oo=wcs.wcs_world2pix(212., bRange[1], 0, 0)	
			

		ax.set_ylim(min(index1,index2),max(index1,index2))

		
	def G214Multiwave(self):
		"""
		
		draw a over of S284 with multi band, including optical, WISE12, WSIE 22, Herschel 250, others?
		
		"""
 
		dataPath='./data/'

		saveFigure="multiwaveG214.pdf"

		lRange=[215.2,214.1]
		
		bRange=[-2.2, -1.6]

 
		G214_Herschel=dataPath+"G214Herschel250.fits"
		#S284_DSS2Red=dataPath+"S284DSS2Red.fits"
		

		G214_WISE12=dataPath+"G214WISE12.fits"
		G214_WISE22=dataPath+"G214WISE22.fits"
		G214_CO13=dataPath+"G214CO13.fits"
		G214_CO12=dataPath+"G214CO12.fits"

		
		
		dataHerschel,headHerschel=myFITS.readFITS(G214_Herschel)
		#dataDSS2Red,headDSS2Red=myFITS.readFITS(S284_DSS2Red)

		dataWISE12,headWISE12=myFITS.readFITS(G214_WISE12)
		dataWISE22,headWISE22=myFITS.readFITS(G214_WISE22)


		#13CO
		dataCO13,headCO13=myFITS.readFITS(G214_CO13)
		
		dataCO12,headCO12=myFITS.readFITS(G214_CO12)



		fig = plt.figure(1, figsize=(8,4.5) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)

		if 1:
			#herschel
			axHerschel=  pywcsgrid2.subplot(221, header=  WCS( headHerschel) )
			axHerschel.imshow(dataHerschel,origin='lower',cmap="bone",vmin=25 ,vmax=160,interpolation='none')
			
			self.setLrange(axHerschel,WCS( headHerschel), lRange)
			self.setBrange(axHerschel,WCS( headHerschel), bRange)
 
			at = AnchoredText(r"Herchel 250 $\mu$m" , loc=4, frameon=False,prop={"color":"w"} )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			axHerschel.add_artist(at)
			axHerschel.set_ticklabel_type("absdeg", "absdeg")
			axHerschel.axis[:].major_ticks.set_color("w")

		if 0:

			axWISE12= pywcsgrid2.subplot(222, header=  WCS( headWISE12) )
			
			axWISE12.imshow(dataWISE12,origin='lower',cmap="bone" , norm=LogNorm(vmin=480, vmax=500),   interpolation='none')
			self.setLrange(axWISE12,WCS( headWISE12), lRange)
			self.setBrange(axWISE12,WCS( headWISE12), bRange) 


		if 1:

			axWISE22= pywcsgrid2.subplot(223, header=  WCS( headWISE22) )
			
			axWISE22.imshow(dataWISE22,origin='lower',cmap="bone" ,  norm=LogNorm(vmin=139.5,  vmax=142),   interpolation='none')
			
			#axWISE22.imshow(dataWISE22,origin='lower',cmap="bone" ,  vmin=150,  vmax=170 ,   interpolation='none')

			#draw WISE YSO
			
			WISEGLON="GLON"
			WISEGLAT="GLAT"
			print "Total number of YSO in WISE",len (self.WISEYSOTB)

			WISEYSOG214=self.filterTBByRange(self.WISEYSOTB,"GLON",lRange)
			WISEYSOG214=self.filterTBByRange(WISEYSOG214,"GLAT",bRange)

			#for eachYSO in WISEYSOG214:
				
				#print eachYSO[WISEGLON],eachYSO[WISEGLAT] ,eachYSO["SDist"]
			print "Total number of YSO in G214",len (WISEYSOG214)
			axWISE22["gal"].scatter( WISEYSOG214[WISEGLON],WISEYSOG214[WISEGLAT],s=10,facecolor='none' ,edgecolor="limegreen",lw=0.5)


			 




			self.setLrange(axWISE22,WCS( headWISE22), lRange)
			self.setBrange(axWISE22,WCS( headWISE22), bRange)
 
			
			
			at = AnchoredText(r"WISE 22 $\mu$m" , loc=4, frameon=False,prop={"color":"w"} )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			axWISE22.add_artist(at)

	
		if 2: #CO12

			axCO12= pywcsgrid2.subplot(222, header=  WCS( headCO12) )
			
			#axWISE22.imshow(dataWISE22,origin='lower',cmap="bone" ,  norm=LogNorm(vmin=139.5,  vmax=142),   interpolation='none')
			
			#axWISE22.imshow(dataWISE22,origin='lower',cmap="bone" ,  vmin=150,  vmax=170 ,   interpolation='none')
			
			int12CO=np.sum(dataCO12,axis=0)*self.vResolution12
			
			axCO12.imshow(int12CO,origin='lower',cmap="bone" ,  vmin=0, vmax=30,   interpolation='none')

			
			#axWISE22[WCS(headCO13)].contour(int13CO,levels=[3,6,9,12,15] ,cmap='jet',linewidths=0.3)
			
			
			self.setLrange(axCO12,WCS( headCO12), lRange)
			self.setBrange(axCO12,WCS( headCO12), bRange)
	
			at = AnchoredText(r"$^{12}\mathrm{CO}~(J=1\rightarrow0)$", loc=4, frameon=False,prop={"color":"w"} )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			axCO12.add_artist(at)
	
		if 1: #CO13

			axCO13= pywcsgrid2.subplot(224, header=  WCS( headCO13) )
			
			#axWISE22.imshow(dataWISE22,origin='lower',cmap="bone" ,  norm=LogNorm(vmin=139.5,  vmax=142),   interpolation='none')
			
			#axWISE22.imshow(dataWISE22,origin='lower',cmap="bone" ,  vmin=150,  vmax=170 ,   interpolation='none')
			
			int13CO=np.sum(dataCO13,axis=0)*self.vResolution13
			
			axCO13.imshow(int13CO,origin='lower',cmap="bone" ,  vmin=0, vmax=10,   interpolation='none')

			
			#axWISE22[WCS(headCO13)].contour(int13CO,levels=[3,6,9,12,15] ,cmap='jet',linewidths=0.3)
		 
			
			self.setLrange(axCO13,WCS( headCO13), lRange)
			self.setBrange(axCO13,WCS( headCO13), bRange)
	

					#axCO13,headCO13=myFITS.readFITS(G214_CO13)

	
			at = AnchoredText(r"$^{13}\mathrm{CO}~(J=1\rightarrow0)$", loc=4, frameon=False,prop={"color":"w"} )
			axCO13.add_artist(at)

 
			
		for eachAx in [ axCO12,axWISE22,axCO13]:
		
 
			eachAx[WCS(headHerschel)].contour(dataHerschel,levels=[75] ,colors='r',linewidths=0.3)
			eachAx.set_ticklabel_type("absdeg", "absdeg")
			eachAx.axis[:].major_ticks.set_color("w")
					
					

		fig.tight_layout()
		plt.savefig(saveFigure, bbox_inches="tight")



		
		
	def S284Multiwave(self):
		"""
		
		draw a over of S284 with multi band, including optical, WISE12, WSIE 22, Herschel 250, others?
		
		"""
 
		dataPath='./data/'

		saveFigure="multiwave_S284.pdf"


 
		S284_Herschel=dataPath+"S284_Herschel250.fits"
		S284_DSS2Red=dataPath+"S284DSS2Red.fits"
		
 
		S284_WISE12=dataPath+"S284WISE12.fits"
		S284_WISE22=dataPath+"S284WISE22.fits"
		
		
		
		dataHerschel,headHerschel=myFITS.readFITS(S284_Herschel)
		dataDSS2Red,headDSS2Red=myFITS.readFITS(S284_DSS2Red)

		dataWISE12,headWISE12=myFITS.readFITS(S284_WISE12)
		dataWISE22,headWISE22=myFITS.readFITS(S284_WISE22)


		
		fig = plt.figure(1, figsize=(8,4.5) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
 
 
		axWISE12= pywcsgrid2.subplot(121, header=  WCS( headWISE12) )
		
		axWISE12.imshow(dataWISE12,origin='lower',cmap="bone" , norm=LogNorm(vmin=500, vmax=700),   interpolation='none')
		self.setLrange(axWISE12,WCS( headWISE12), [211.6,212.4])
		self.setBrange(axWISE12,WCS( headWISE12), [-1.7, -0.9 ]) 
		
		
		
		
		axWISE12.set_ticklabel_type("absdeg", "absdeg")
		axWISE12.axis[:].major_ticks.set_color("w")


		at = AnchoredText(r"WISE 12 $\mu$m", loc=3, frameon=True)
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axWISE12.add_artist(at)





		if 0:# WISE2
			axWISE22= pywcsgrid2.subplot(222, header=  WCS( headWISE22) )
			
			axWISE22.imshow(dataWISE22,origin='lower',cmap="bone" ,  norm=LogNorm(vmin=150,  vmax=165),   interpolation='none')
			
			#axWISE22.imshow(dataWISE22,origin='lower',cmap="bone" ,  vmin=150,  vmax=170 ,   interpolation='none')
	
			
			self.setLrange(axWISE22,WCS( headWISE22), [211.6,212.4])
			self.setBrange(axWISE22,WCS( headWISE22), [-1.7, -0.9 ])
	
			axWISE22.set_ticklabel_type("absdeg", "absdeg")
			axWISE22.axis[:].major_ticks.set_color("w")
	
		 
			
		
		#axDSS2Red["gal"].set_xlim(211.5,212)
		if 0:
		
			axHerschel=  pywcsgrid2.subplot(223, header=  WCS( headHerschel) )
			axHerschel.imshow(dataHerschel,origin='lower',cmap="bone",vmin=25 ,vmax=250,interpolation='none')
			
			self.setLrange(axHerschel,WCS( headHerschel), [211.6,212.4])
			self.setBrange(axHerschel,WCS( headHerschel), [-1.7, -0.9 ])
			axHerschel.set_ticklabel_type("absdeg", "absdeg")
			axHerschel.axis[:].major_ticks.set_color("w")
		
		#ax12.imshow(dataHerschel,origin='lower',cmap="bone",vmin=25 ,vmax=250,interpolation='none')
		
 
		axDSS2Red= pywcsgrid2.subplot(122, header=  WCS( headDSS2Red) )
		
		axDSS2Red.imshow(dataDSS2Red,origin='lower',cmap="bone" , norm=LogNorm(vmin=4600,  vmax=8000),   interpolation='none')
		self.setLrange(axDSS2Red,WCS( headDSS2Red), [211.6,212.4])
		self.setBrange(axDSS2Red,WCS( headDSS2Red), [-1.7, -0.9 ])

		axDSS2Red.set_ticklabel_type("absdeg", "absdeg")
		axDSS2Red.axis[:].major_ticks.set_color("w")

		
		at = AnchoredText(r"DSS2 Red", loc=3, frameon=True)
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axDSS2Red.add_artist(at)

		
		 


		fig.tight_layout()
		plt.savefig(saveFigure, bbox_inches="tight")




 




	def drawG214Outflow(self):
		"""
		A possible high mass star forming region in a fillament
		
		"""
		
		dataPath='./data/'

		saveFigure="G214Outflow.pdf"

		#lRange=[215.3,214.1]
		#bRange=[-2.3, -1.5]
		centerL=214.493 #def
		centerB=-1.811 #deg
		
		
		bluePL=214.46667
		bluePB=-1.8166667
		
		
		sizeL=0.6 
		sizeB=0.4
		
		sizeGlimpseL=0.065
		sizeGlimpseB=0.045
		
		centerLB=[centerL,centerB]
		
		lRange=[ centerL-sizeL/2,  centerL+ sizeL/2 ]
		
		
		bRange=[ centerB-sizeB/2,  centerB +sizeB/2 ]
 
		doFITS=myFITS()


		radias=120./60./60. #deg

		

 
		G214_Herschel=dataPath+"G214Herschel250.fits"
		#S284_DSS2Red=dataPath+"S284DSS2Red.fits"
		

		G214_WISE12=dataPath+"G214WISE12.fits"
		G214_WISE22=dataPath+"G214WISE22.fits"
		G214_CO13=dataPath+"G214CO13.fits"
		G214_CO12=dataPath+"G214CO12.fits"

		G214_GLM45=dataPath+"G214GLM45.fits"
		
		dataGLM45,headGLM45=myFITS.readFITS(G214_GLM45)
 
		dataHerschel,headHerschel=myFITS.readFITS(G214_Herschel)
		#dataDSS2Red,headDSS2Red=myFITS.readFITS(S284_DSS2Red)

		dataWISE12,headWISE12=myFITS.readFITS(G214_WISE12)
		dataWISE22,headWISE22=myFITS.readFITS(G214_WISE22)




		#13CO
		
		dataCO13,headCO13=myFITS.readFITS(G214_CO13)
		dataCO12,headCO12=myFITS.readFITS(G214_CO12)

		
		data12all,head12all=self.doFITS.readFITS(self.COFile12)
		data13all,head13all=self.doFITS.readFITS(self.COFile13)


		#get the spectral of CO 13
		specCO12C,veloCO12C=myFITS.getSpectraByLB(data12all,head12all,centerL,centerB) 
		specCO13C,veloCO13C=myFITS.getSpectraByLB(data13all,head13all,centerL,centerB) 

		
		blueWingeRange=[20,24]

		
		momentFITS="outflowBlueMoment.fits"
		doFITS=myFITS()
		blueData,blueHead=doFITS.momentFITS( self.COFile12,blueWingeRange,0,outFITS=momentFITS)
		cropBlue="cropblue.fits"
		#100%
 
		
		
		myFITS.cropFITS2D(momentFITS,outFITS=cropBlue,Lrange= lRange, Brange=bRange,overWrite=True )
		blueData,blueHead=myFITS.readFITS(cropBlue)
 
		
		
		



		blueL=214.4833300
		blueB= -1.8166667
		
		
		specCO12B,veloCO12B = myFITS.getSpectraByLB(data12all,head12all,blueL,blueB) 
		specCO13B,veloCO13B = myFITS.getSpectraByLB(data13all,head13all,blueL,blueB) 


		fig = plt.figure(1, figsize=(8,6) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)

 

		if 1:
			#herschel
			axHerschel=  pywcsgrid2.subplot(221, header=  WCS( headHerschel) )
			axHerschel.imshow(dataHerschel,origin='lower',cmap="bone",vmin=30 ,vmax=160,interpolation='none')
			
			self.setLrange(axHerschel,WCS( headHerschel), lRange)
			self.setBrange(axHerschel,WCS( headHerschel), bRange)


			axHerschel.set_ticklabel_type("absdeg", "absdeg")
			axHerschel.axis[:].major_ticks.set_color("w") 

			#draw contour
			
			
			noise12=self.getNoise12(blueWingeRange) 
			#print noise12,"noise????????"
			axHerschel[WCS(blueHead)].contour( blueData, levels=np.arange(2,8)*2*noise12,colors='b' ,linewidths=0.5,alpha=0.6)

			#axHerschel["gal"].plot(bluePL,bluePB,'g.')

			#axHerschel["gal"].plot( centerL,centerB,'g.')
 			#axHerschel["gal"].plot( blueL,blueB,'b.')
 			
 			pixBeam=52./3600./headHerschel["CDELT1"]
 			pixBeam=abs(pixBeam)
			axHerschel.add_beam_size(pixBeam, pixBeam, 0,  3 , patch_props={"facecolor":'None','edgecolor':'w',"linewidth":0.8}  )
			at = AnchoredText(r"Herchel 250 $\mu$m" , loc=4, frameon=False,prop={"color":"w"} )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			axHerschel.add_artist(at)

		if 2:# GLIMPSE, I2, 4.5 micron, I 1, 3.6 Micron
			axGLM45=  pywcsgrid2.subplot(222, header=  WCS( headGLM45) )
			#axGLM45.imshow(dataGLM45,origin='lower',cmap="bone", interpolation='None',norm=colors.LogNorm(vmin=0.3 , vmax= 200))
			#axGLM45.imshow(dataGLM45,origin='lower',cmap="bone", interpolation='none', vmin=0.2 , vmax= 10 )

			
		#cmap.set_bad('white',1.)
		
		
			cmap = plt.cm.bone
			cmap.set_bad('black',1.)
			axGLM45.imshow(dataGLM45,origin='lower',cmap=cmap,   vmin=0.2 , vmax= 200 , interpolation='None',norm=colors.LogNorm())

		
		
			self.setLrange(axGLM45,WCS( headGLM45), [ centerL-sizeGlimpseL/2., centerL+sizeGlimpseL/2. ])
			self.setBrange(axGLM45,WCS( headGLM45), [ centerB-sizeGlimpseB/2., centerB+sizeGlimpseB/2. ])

			axGLM45.set_ticklabel_type("absdeg", "absdeg")
			axGLM45.axis[:].major_ticks.set_color("w") 

			at = AnchoredText(r"GLIMPSE360 4.5 $\mu$m" , loc=4, frameon=False,prop={"color":"w"} )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			axGLM45.add_artist(at)
		if 0: #draw spectral
			
			axSpecCenter=plt.subplot(223)
			
			#axSpecCenter.step(veloCO13C,specCO13C,color='black',lw=0.5)
			#axSpecCenter.step(veloCO13B,specCO13B,color='blue',lw=0.5)


			axSpecCenter.step(veloCO12C,specCO12C,color='green',lw=0.5)
			axSpecCenter.step(veloCO12B,specCO12B,color='blue',lw=0.5)


			axSpecCenter.step(veloCO13C,specCO13C,color='green',lw=0.5,linestyle='--')
			axSpecCenter.step(veloCO13B,specCO13B,color='blue',lw=0.5,linestyle='--')



			minV=10
			maxV=45
 
			axSpecCenter.plot([minV,maxV],[0,0],color='black',lw=0.8)
			axSpecCenter.set_xlim(minV,maxV)






		if 1: #draw line on Herschel

			minV=15
			maxV=40
			wcsHerschel=WCS(headHerschel )
			axSpecAll= pywcsgrid2.subplot(223, header=  wcsHerschel )

			axSpecAll.imshow(dataHerschel,origin='lower',cmap="bone",vmin=30 ,vmax=400,interpolation='None',  )
			axSpecAll[WCS(blueHead)].contour( blueData, levels=np.arange(2,8)*2*noise12,colors='b' ,linewidths=0.5,alpha=0.6)


			axSpecAll.set_ticklabel_type("absdeg", "absdeg")
			axSpecAll.axis[:].major_ticks.set_color("w") 


			#axSpecAll["gal"].plot(bluePL,bluePB,'g.')



			if  0:#draw beam size

				pixBeam=52./3600./headHerschel["CDELT1"]
 				pixBeam=abs(pixBeam)
				axSpecAll.add_beam_size(pixBeam, pixBeam, 0,  3 , patch_props={"facecolor":'None','edgecolor':'w',"linewidth":0.3}  )
			#at = AnchoredText(r"Herchel (250 $\mu$m)" , loc=4, frameon=False,prop={"color":"w"} )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			#axSpecAll.add_artist(at)




			#draw 12CO spectral
			
			#start from centerLB
			#find nearstest  point around centerLB
			
			pixCenterL,pixCenterB,a=self.WCSCO12.wcs_world2pix(centerL,centerB, 27,0 )


			#pixCenterL,pixCenterB,a=self.WCSCO12.wcs_world2pix(centerL,centerB, 27,0 )

			a,b,minVindex=self.WCSCO12.wcs_world2pix(centerL,centerB, minV*1000,0 )
			a,b,maxVindex=self.WCSCO12.wcs_world2pix(centerL,centerB, maxV*1000,0 )


			#minvIndex=
			
			
			
			pixCenterL =  int(round(pixCenterL) )
			pixCenterB =  int(round(pixCenterB) )
			
			
			
			minVindex =  int(round(minVindex) )
			maxVindex =  int(round(maxVindex) )
			
			
			
			
			drawPixs=[  (pixCenterL,pixCenterB)    ]
			
			#creast pixels
			
			xN=3
			yN=2
			
			
			
			for i in np.arange(-xN,xN+1):
				for j in np.arange(-yN,yN+1):
 
					drawPixs.append(( pixCenterL+i, pixCenterB +j ) )
 
 
			
			
			
			for position in drawPixs:
			
				spL,spB= position
				
 
				centerLCO12,centerBCO12,a= self.WCSCO12.wcs_pix2world(spL,spB, 27,0 )
	 
	
				sp1,vo1=myFITS.getSpectraByIndex(self.CO12DATA,self.COHEADER,spL,spB) 
				
				sp1=sp1[minVindex:maxVindex+1]
				vo1=vo1[minVindex:maxVindex+1]

				
				#better Cur spectral
				
				
				#sp11,vo11=myFITS.getSpectraByLB( self.CO12DATA,self.COHEADER, centerLCO12,centerBCO12 ) 
	

				#convert spectral to pixels on herschel image
				#left
				LB=[ centerLCO12, centerBCO12 ]
 
 				
				herchelSpX,herschelSpY=self.getHerschelSpXY(wcsHerschel,sp1,vo1,LB)
				
  
				axSpecAll['gal'].plot(herchelSpX,herschelSpY,'r-',lw=0.1)
				
			
			
			
			self.setLrange(axSpecAll,WCS( headHerschel), [ centerL-sizeGlimpseL/2., centerL+sizeGlimpseL/2. ])
			self.setBrange(axSpecAll,WCS( headHerschel), [ centerB-sizeGlimpseB/2., centerB+sizeGlimpseB/2. ])



		if 1: #draw PV diagram
			
			width=4./60. #deg
			
			widthPix=width/(30./3600.)
			
			
			length=20./60. #deg
			
			
			
			startPL=[centerL-length/2.  ]
			
			startPB=[centerB-length/2. ]

			slope=(centerB-blueB)/(centerL-blueL)
			
			theta=np.arctan(slope)
			
						
			startPL=  centerL+length/2.*np.cos( theta)   
			startPB=  centerB+length/2.*np.sin( theta)   

			endPL=  centerL-length/2.*np.cos( theta)   
			endPB=  centerB-length/2.*np.sin( theta)   

			#self.CO12DATA,self.COHEADE
			SPX1,SPY1,a=self.WCSCO12.wcs_world2pix( startPL,startPB,27000, 0 )
			
			SPX2,SPY2,a=self.WCSCO12.wcs_world2pix( endPL,endPB,27000, 0 )

			pvFITS="G214CenterPV.fits"

 

			#os.system("rm "+pvFITS)

			self.drawPV( self.COFile12,[SPX1,SPY1],[SPX2,SPY2],widthPix=widthPix,saveName= pvFITS ) 


			# plot the staring point of path
			#axHerschel["gal"].plot( [startPL,endPL], [startPB, endPB]   )			
			 
						
			XY,HW,angleRect,arrowStart,arrowEnd=self.getRectangleByWCS( headHerschel,centerLB,length,width,-theta*180/np.pi) # the input angle is the pixel angle
			rect=Rectangle((XY[0], XY[1]), HW[1],HW[0] , angle=angleRect, facecolor='None', edgecolor = 'cyan', linestyle='--', linewidth=0.8)
			axHerschel.add_patch( rect)
			
			axHerschel.arrow(arrowStart[0],arrowStart[1],arrowEnd[0]-arrowStart[0],arrowEnd[1]-arrowStart[1] ,   head_width=10, fc='cyan',ec='cyan',lw=1., length_includes_head=True, capstyle="projecting"  )
			
			
			XY,HW,angleRect,arrowStart,arrowEnd=self.getRectangleByWCS( headHerschel,centerLB,sizeGlimpseL,sizeGlimpseB,0) # the input angle is the pixel angle
			rect=Rectangle((XY[0], XY[1]), HW[1],HW[0] , angle=angleRect, facecolor='None', edgecolor = 'purple', linewidth=0.8,linestyle='--')
			axHerschel.add_patch( rect)
			
			
			#read pvFITS
			pvData,pvHead=myFITS.readFITS( pvFITS )
			
			
			axPV=  pywcsgrid2.subplot(224, header=  WCS( pvHead) )

			#axPV.imshow( pvData,origin='lower',aspect='auto')
			axPV.axis[:].major_ticks.set_color("w") 
			
			pvRMS=0.2
			maxSlice=np.max( pvData )
			
			

			contLevels= np.arange( pvRMS*3,maxSlice+pvRMS*3,pvRMS*3)

			#contLevels= np.arange( pvRMS*3,maxSlice+pvRMS*5,pvRMS*5)


			contLevels=list(contLevels)
			
			contLevels.insert(0,np.min( pvData) )
 




			
			axPV.contourf(  pvData ,levels=contLevels,cmap='bone' )

			pvNx,pvNy=pvData.shape

			
			# modify the labels
			minV=18
			
			pvWCS=WCS(pvHead)
			
			a,minY =pvWCS.wcs_world2pix(0,minV*1000,  0)
			a,maxY =pvWCS.wcs_world2pix(0,maxV*1000,  0)

		
			axPV.set_ylim( minY, maxY)			

			a,blueIndex1=pvWCS.wcs_world2pix(0,blueWingeRange[0]*1000,  0)
			a,blueIndex2=pvWCS.wcs_world2pix(0,blueWingeRange[1]*1000,  0)


			#axPV.plot([0,pvNx],[ blueIndex1,blueIndex1] ,'b--',lw=0.5  )

			axPV.axhline(y=blueIndex1, linestyle='--',linewidth=0.5,color='b' )
			axPV.axhline(y=blueIndex2, linestyle='--',linewidth=0.5,color='b' )

			
			
			#axPV['ga'].plot( [0,0.32] , [blueWingeRange[0], blueWingeRange[0] ],'w--',lw=0.5  )
	
			#axPV.set_ticklabel2_type("manual",locs=[35000.],labels=["35000"])

			
			if 1:
			
				vInterval=2 #km/s
				
				yLocs=np.arange(int(minV) , int(maxV)+vInterval,vInterval)
				
				yLabels= map(int, yLocs)
				yLabels=map(str,yLabels)
				
				#print yLocs*1000.
				#print yLabels
				
				
						 
				axPV.set_ticklabel2_type("manual",locs=yLocs*1000.,labels=yLabels)
				#axPV.set_ticklabel2_type("manual",locs=[35000.],labels=["35000"])



			axPV.set_ylabel(r"Velocity ($\rm km \ s^{-1}$)")
			
			axPV.set_xlabel(r"Offset (degree)")


 

		fig.tight_layout()
		plt.savefig(saveFigure, bbox_inches="tight")



	def getPVFitsByLB(self,FITSName,centerLB,width,length,angle,pvFITS):
		
		
		"""
		
		Better write this in to myPYTHON
		
		widthSize,length in degree 
		
		"""
		#
		
		#resolution
		
 
		fitsData,fitsHead=myFITS.readFITS( FITSName)
		
		wcs=WCS(fitsHead)
		
		centerL,centerB=centerLB
		
		resolution=abs(fitsHead["CDELT1"] ) #in degree
		
 
		
		widthPix=width/resolution
		
 
		#slope=(centerB-blueB)/(centerL-blueL)
		
		theta=   np.radians(angle )
		
					
		startPL=  centerL+length/2.*np.cos( theta)   
		startPB=  centerB-length/2.*np.sin( theta)   

		endPL=  centerL-length/2.*np.cos( theta)   
		endPB=  centerB+length/2.*np.sin( theta)   


		print startPL,startPB
		
		print endPL,endPB

		#self.CO12DATA,self.COHEADE
		SPX1,SPY1,a=wcs.wcs_world2pix( startPL,startPB,27000, 0 )
		
		SPX2,SPY2,a=wcs.wcs_world2pix( endPL,endPB,27000, 0 )

		self.drawPV(  FITSName, [SPX1,SPY1],[SPX2,SPY2],widthPix=widthPix,saveName= pvFITS ) 

		#os.system("rm "+pvFITS)

		# /home/qzyan/WORK/projects/maddalena/data/S284MaskFITSCO12.fits, test S284
		
		
		
		#self.drawPV(  "./data/S284MaskFITSCO12.fits", [SPX1,SPY1],[SPX2,SPY2],widthPix=widthPix,saveName= pvFITS ) 

		#self.drawPV(  "./data/S284MaskFITSCO12.fits", [SPX1,SPY1],[SPX2,SPY2],widthPix=widthPix,saveName= pvFITS ) 

		#self.drawPV( self.COFile12, [SPX1,SPY1],[SPX2,SPY2],widthPix=widthPix,saveName= pvFITS ) 

	def getHerschelSpXY(self,herschelWCS,sp,velo,LB):
		"""
		"""
		
		dPix=30./60./60. #deg
		
		
		leftL=LB[0]+ dPix/2.
		rightL=LB[0]-dPix/2.

		lowerB=LB[1]- dPix/2.
		upperB=LB[1]+ dPix/2.

 
		
		X=velo/(np.max(velo)-np.min(velo) )*( rightL- leftL )+ leftL- np.min(velo)*( rightL- leftL)/(np.max(velo)-np.min(velo) )
		#Y=velo/(np.max(velo)-np.min(velo) )*( rightL- leftL )+ leftL- np.min(velo)*( rightL- leftL)/(np.max(velo)-np.min(velo) )

		Tmax=15.
		Tmin=-2.
		
		
 
		
		Y=sp/(Tmax-Tmin)*( upperB- lowerB )+ lowerB- Tmin*( upperB- lowerB)/(Tmax-Tmin)

 

		return X,Y



	def getRectangleByWCS(self,wiseHead,centerLB,length,width,angle):
		
		"""
		angle, in degree
		"""
		WCS_WISE=WCS(wiseHead)
		

		pixLengthWISE= abs(wiseHead["CDELT1"])
 

		centerXWISE,centerYWISE = WCS_WISE.wcs_world2pix(centerLB[0],centerLB[1],0)
		#beginPWISE,endPWISE=doMadd.getEndPoints([centerXWISE,centerYWISE],  pathAngle, length, pixLengthWISE)
		widthPixWISE=width/pixLengthWISE
		
		lengthPixWISE= length/pixLengthWISE
		
		
		heightAngle=angle 
		
		heightAngleRadians=np.radians(heightAngle)
		
		vectorHeightX=lengthPixWISE*np.cos(heightAngleRadians)
		vectorHeightY=lengthPixWISE*np.sin(heightAngleRadians)
		
		widthAngle=    angle+90
		
		widthAngleeRadians=np.radians(widthAngle)
		
		vectorWidthX=widthPixWISE*np.cos(widthAngleeRadians)
		vectorWidthY=widthPixWISE*np.sin(widthAngleeRadians)


		XY=[ centerXWISE -vectorHeightX/2.+vectorWidthX/2.  ,  centerYWISE  -vectorHeightY/2.+vectorWidthY/2.  ]
		
		HW=[lengthPixWISE,widthPixWISE ]




		arrowStart=[ centerXWISE-vectorHeightX/2., centerYWISE  -vectorHeightY/2.  ]
		arrowEnd=[ centerXWISE+vectorHeightX/2., centerYWISE  +vectorHeightY/2.  ]


		return XY,HW,270+angle,arrowStart,arrowEnd
		


	def plotNegative(self,moment=False):
		
		"""
		plot three negative velocity components
		
		"""

		#integal the data ,
		#box(210.2094116,4.3903029,272.541",329.741",0)
		#box(216.3728503,2.4680128,298.565",206.075",0)
		#box(215.1582116,1.6176686,298.610",217.478",0)
		#A(-6.8, -5)
		
		saveFigure="negativeV.pdf"
		negativeCO12FITS="./data/mosaic_U_minus70_0.fits"


		##########################draw AAAAAAAAAAAAAAAAAAAAAAA

		vRangeA=[-6.8,-5]
		
		lRangeA,bRangeA,centerLBA,sizeLBA=box(216.3728503,2.4680128,298.565,206.075,0)
 
 
 
		fitsA="momentFITSA.fits"
		RMSA=0.33
		
		if moment:
			Adata,Ahead=doMadd.doFITS.momentFITS(negativeCO12FITS, vRangeA,0,outFITS=fitsA)
		else:
			Adata,Ahead=myFITS.readFITS(fitsA)
			
		Adata,Ahead=self. convert3DTo2D(Adata,Ahead)
		#draw three parts
		fig = plt.figure(1, figsize=(8,6) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		wcsA=WCS( Ahead)
 
		
		grid = ImageGrid(fig, (3, 2, 3), nrows_ncols = (1, 1),
		                 cbar_mode="single", cbar_pad="2%",
		                 cbar_location="right",
		                 axes_class=(pywcsgrid2.Axes, dict(header=wcsA)))
		
		axAimg = grid[0]

		cb_axes = grid.cbar_axes[0] # colorbar axes
		im=axAimg.imshow(Adata ,origin='lower',cmap="bone",vmin=RMSA*0.5,vmax=RMSA*10,interpolation='none') 
		#im = main_axes.imshow(d, origin="lower", cmap=plt.cm.gray_r,
		                     # vmin=4.e-05, vmax=0.00018,
		                      #interpolation="nearest")
		#cb_axes.set_ylabel("aa")
		cb_axes.colorbar(im)
		#cb_axes.set_xlabel("cc")
		cb_axes.set_ylabel(r"$\rm K \ km \ s^{-1}$")

		
		
		
		 
		
		
		axAimg.set_ticklabel_type("absdeg", "absdeg")
		axAimg.axis[:].major_ticks.set_color("w") 
		axAimg.locator_params(axis="x", nbins=3)

		XY,HW,angleRect,arrowStart,arrowEnd=self.getRectangleByWCS( Ahead,centerLBA,sizeLBA[0],sizeLBA[1],0) # the input angle is the pixel angle
		rect=Rectangle((XY[0], XY[1]), HW[1],HW[0] , angle=angleRect, facecolor='None', edgecolor = 'green', linewidth=0.5,linestyle='--')
		axAimg.add_patch( rect)
		
	
		disLrange=[lRangeA[0]-0.1, lRangeA[1]+0.1 ]
		disBrange=[bRangeA[0]-0.1, bRangeA[1]+0.1]

		self.setLrange( axAimg,wcsA,disLrange)
		self.setBrange( axAimg,wcsA,disBrange)




		axAspec=plt.subplot(324)	

		# get average spectral

		#crop fits
		
		#averageSpecFITS="cropSpectral.fits"

		#myFITS.cropFITS(negativeCO12FITS, outFITS=	averageSpecFITS  ,  Lrange=lRangeA, Brange=bRangeA, overWrite=True)

		avgSpec,avgVs=doMadd.doFITS.getAverageSpecByLBrange( negativeCO12FITS,lRangeA,bRangeA )
		axAspec.step(avgVs,avgSpec,color='blue',lw=0.5)
		axAspec.axhline(y=0,color='black',lw=0.4)
		axAspec.set_xlim(min(vRangeA)-10, 0)


		axAspec.set_xlabel(r"$\rm V_{LSR} \left(km \ s^{-1}\right)$")
		axAspec.set_ylabel(r"$\rm T_{\rm mb} \left(K\right)$")   

 

		at = AnchoredText(r"B", loc=1, frameon=False)
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axAspec.add_artist(at)


		at = AnchoredText(r"B" , loc=2, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axAimg.add_artist(at)



		##########################draw AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa

		vRangeA=[-7.7,-6.2]   

		
		lRangeA,bRangeA,centerLBA,sizeLBA= box(215.1582116,1.6176686,298.610,217.478,0)
 

 
		fitsA="momentFITSB.fits"
		RMSA=0.39
		
		if moment:
			Adata,Ahead=doMadd.doFITS.momentFITS(negativeCO12FITS, vRangeA,0,outFITS=fitsA)
		else:
			Adata,Ahead=myFITS.readFITS(fitsA)
			
		Adata,Ahead=self. convert3DTo2D(Adata,Ahead)
		#draw three parts
		fig = plt.figure(1, figsize=(8,6) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		wcsA=WCS( Ahead)
 


		grid = ImageGrid(fig, (3, 2, 1), nrows_ncols = (1, 1),
		                 cbar_mode="single", cbar_pad="2%",
		                 cbar_location="right",
		                 axes_class=(pywcsgrid2.Axes, dict(header=wcsA)))
		
		axAimg = grid[0]
		 
		
		cb_axes = grid.cbar_axes[0] # colorbar axes
		im=axAimg.imshow(Adata ,origin='lower',cmap="bone",vmin=RMSA*0.5,vmax=RMSA*6,interpolation='none') 
		#im = main_axes.imshow(d, origin="lower", cmap=plt.cm.gray_r,
		                     # vmin=4.e-05, vmax=0.00018,
		                      #interpolation="nearest")
		#cb_axes.set_ylabel("aa")
		cb_axes.colorbar(im)
		#cb_axes.set_xlabel("cc")
		cb_axes.set_ylabel(r"$\rm K \ km \ s^{-1}$")





		
		
		axAimg.set_ticklabel_type("absdeg", "absdeg")
		axAimg.axis[:].major_ticks.set_color("w") 
		axAimg.locator_params(axis="x", nbins=3)

		XY,HW,angleRect,arrowStart,arrowEnd=self.getRectangleByWCS( Ahead,centerLBA,sizeLBA[0],sizeLBA[1],0) # the input angle is the pixel angle
		rect=Rectangle((XY[0], XY[1]), HW[1],HW[0] , angle=angleRect, facecolor='None', edgecolor = 'green', linewidth=0.5,linestyle='--')
		axAimg.add_patch( rect)
		
	
		disLrange=[lRangeA[0]-0.1, lRangeA[1]+0.1 ]
		disBrange=[bRangeA[0]-0.1, bRangeA[1]+0.1]


		#crop crops
		cropedWISE22='cropedWISE22.fits'
		
		myFITS.cropFITS2D(self.WISEFITS, outFITS=cropedWISE22,Lrange=disLrange,Brange=disBrange,overWrite=True)

		rmsWISE=0.0476529
		
		startWISE=144.84
		
		levels=np.arange(10)*rmsWISE*5+startWISE
		
		cropWISEData,cropWISEHead= myFITS.readFITS( cropedWISE22 )
		axAimg[WCS(cropWISEHead)].contour( cropWISEData,levels=levels, colors='r',linewidths=0.3 )
		

		self.setLrange( axAimg,wcsA,disLrange)
		self.setBrange( axAimg,wcsA,disBrange)




		axAspec=plt.subplot(322)	

		# get average spectral

		#crop fits
		
		#averageSpecFITS="cropSpectral.fits"

		#myFITS.cropFITS(negativeCO12FITS, outFITS=	averageSpecFITS  ,  Lrange=lRangeA, Brange=bRangeA, overWrite=True)

		avgSpec,avgVs=doMadd.doFITS.getAverageSpecByLBrange( negativeCO12FITS,lRangeA,bRangeA )
		axAspec.step(avgVs,avgSpec,color='blue',lw=0.5)
		axAspec.axhline(y=0,color='black',lw=0.4)
		axAspec.set_xlim(min(vRangeA)-10, 0)


		axAspec.set_xlabel(r"$\rm V_{LSR} \left(km \ s^{-1}\right)$")
		axAspec.set_ylabel(r"$\rm T_{\rm mb} \left(K\right)$")   

 
		at = AnchoredText(r"A", loc=1, frameon=False)
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axAspec.add_artist(at)



		at = AnchoredText(r"A" , loc=2, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axAimg.add_artist(at)






		##########################draw CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

		vRangeA=[-3.8, -2.6]   

		
		lRangeA,bRangeA,centerLBA,sizeLBA=   box(217.7170544,1.5355732,178.329,168.235,0)




 

 
		fitsA="momentFITSC.fits"
		RMSA=0.27
		
		if moment:
			Adata,Ahead=doMadd.doFITS.momentFITS(negativeCO12FITS, vRangeA,0,outFITS=fitsA)
		else:
			Adata,Ahead=myFITS.readFITS(fitsA)
			
		Adata,Ahead=self. convert3DTo2D(Adata,Ahead)
		#draw three parts
		fig = plt.figure(1, figsize=(8,6) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		wcsA=WCS( Ahead)
		
		
		
		grid = ImageGrid(fig, (3, 2, 5), nrows_ncols = (1, 1),
		                 cbar_mode="single", cbar_pad="2%",
		                 cbar_location="right",
		                 axes_class=(pywcsgrid2.Axes, dict(header=wcsA)))
		
		axAimg = grid[0]
		 
		
		cb_axes = grid.cbar_axes[0] # colorbar axes
		im=axAimg.imshow(Adata ,origin='lower',cmap="bone",vmin=RMSA*0.5,vmax=RMSA*7,interpolation='none') 
		#im = main_axes.imshow(d, origin="lower", cmap=plt.cm.gray_r,
		                     # vmin=4.e-05, vmax=0.00018,
		                      #interpolation="nearest")
		#cb_axes.set_ylabel("aa")
		cb_axes.colorbar(im)
		#cb_axes.set_xlabel("cc")
		cb_axes.set_ylabel(r"$\rm K \ km \ s^{-1}$")
		#cb_axes.set_xlabel("K km s$^{-1}$")
		#axAimg=  pywcsgrid2.subplot(325, header= wcsA)
		
		
		
		axAimg.set_ticklabel_type("absdeg", "absdeg")
		axAimg.axis[:].major_ticks.set_color("w") 
		axAimg.locator_params(axis="x", nbins=3)

		XY,HW,angleRect,arrowStart,arrowEnd=self.getRectangleByWCS( Ahead,centerLBA,sizeLBA[0],sizeLBA[1],0) # the input angle is the pixel angle
		rect=Rectangle((XY[0], XY[1]), HW[1],HW[0] , angle=angleRect, facecolor='None', edgecolor = 'green', linewidth=0.5,linestyle='--')
		axAimg.add_patch( rect)
		
	
		disLrange=[lRangeA[0]-0.15, lRangeA[1]+0.15 ]
		disBrange=[bRangeA[0]-0.15, bRangeA[1]+0.15]

		self.setLrange( axAimg,wcsA,disLrange)
		self.setBrange( axAimg,wcsA,disBrange)




		axAspec=plt.subplot(326)	

		# get average spectral

		#crop fits
		
		#averageSpecFITS="cropSpectral.fits"

		#myFITS.cropFITS(negativeCO12FITS, outFITS=	averageSpecFITS  ,  Lrange=lRangeA, Brange=bRangeA, overWrite=True)

		avgSpec,avgVs=doMadd.doFITS.getAverageSpecByLBrange( negativeCO12FITS,lRangeA,bRangeA )
		axAspec.step(avgVs,avgSpec,color='blue',lw=0.5)
		axAspec.axhline(y=0,color='black',lw=0.4)
		axAspec.set_xlim(min(vRangeA)-10, 0)


		axAspec.set_xlabel(r"$\rm V_{LSR} \left(km \ s^{-1}\right)$")
		axAspec.set_ylabel(r"$\rm T_{\rm mb} \left(K\right)$")   

 
        
        
		at = AnchoredText(r"C", loc=1, frameon=False)
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axAspec.add_artist(at)



		at = AnchoredText(r"C" , loc=2, frameon=False,prop={"color":"w"} )
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axAimg.add_artist(at)



 


		fig.tight_layout()
		plt.savefig(saveFigure, bbox_inches="tight")





	def calMassWithCO12(self,FitsFile,lRange,bRange,vRange,distance,xFactor=2.0e20):
		"""
		distance in pc
		"""
		
		#by default ,the xFactor is 2.2e20 
		
		
		RMS=0.5
		
		massFITS="cropForMass.fits"
		vResolution=0.158737644553 #km s


		
		
		
		#crop fits
		
		myFITS.cropFITS(FitsFile,Lrange=lRange,Brange=bRange,Vrange=vRange,outFITS=massFITS,overWrite=True)
		
		fitsData,fitsHead=myFITS.readFITS(massFITS)
		
		#use duchamp to mask
		
		fitsData[fitsData<3*RMS]=0
		
		
		allCo=np.sum(fitsData)*vResolution #K km/s
		
		lenPix= np.tan(  np.radians(30./60./60.) )*distance*self.parsecToMeter*100 #cm
		
		
		
		#pixArea=
		
		NH2=allCo*xFactor
		
		mu=1.36
		
		Mh2=3.35e-27 #kg
		solarMass=1.9891e30 #kg
		#s=np.pi*length1*length2
		s=  lenPix*lenPix

		coreSolar= s*NH2*Mh2*mu/solarMass
		
		
		return coreSolar



	def removeIndex(self,TB):
		
		TestTB=Table() 
		
		for eachCol in TB.colnames:
			
			
			TB[eachCol].mask=False #do not mask
			
			aaa=Column( TB[eachCol]  ,name=eachCol )	
			
 
			TestTB.add_column(aaa )
		return TestTB

	def examineHIIRegionBySp(self,source=None):
		
		"""
		"""
		pass


		#Bdata,Bhead=doMadd.doFITS.readFITS( backGroundFITS)
		
		#if len(Bdata.shape)==3:
			#Bdata=Bdata[0]
		
		colNameWISE="WISE Name"
		colVLSRWISE="VLSR<br>(km/s)" 
		colVLSRMolWISE="VLSR (Mol.)<br>(km/s)" 
		colMembershipWISE="Membership"
		
		colRadiusWISE="Radius<br>(arcsec.)"

		WISEcatalog="Catalog"
 
		
		CO12Mask="./data/G210MaskCO13.fits"
		
		colL="GLong<br>(deg.)"
		colB="GLat<br>(deg.)"


		WISEcat=self.WISEHIICat.copy()

		#WISEcat.sort([ colL,colB ])

		lMin,lMax=min( self.G210LRange), max( self.G210LRange), 
		
		
		maddaHII=self.filterTBByRange(WISEcat,colL, [lMin,lMax ])

		maddaHII=self.filterTBByRange(maddaHII,colB,[-5.,5.])

		
		MaddaIRASHII=self.filterTBByRange(self.IRASHIICat,"glon",[lMin-360,lMax-360])

		MaddaIRASHII=self.filterTBByRange(MaddaIRASHII,"glat",[-5.,5.])


		WC89TB=Table.read("IRASWC89.fits")
		
		WC89G210TB=self.filterTBByRange(WC89TB,"glon",[lMin,lMax])

		WC89G210TB=self.filterTBByRange(WC89G210TB,"glat",[-5.,5.])

		
		#examine WISE HII

		
		#print maddaHII.colnames
		

		print len(maddaHII  ), len(MaddaIRASHII  )


		#eliminate Duplicated Sources
		
		#print MaddaIRASHII.colnames
		
		searchRadius=1./60. #degree 
		
		dupN=0
		
		#MaddaIRASHII.remove_indices("glon")
		#MaddaIRASHII.remove_indices("glat")

		MaddaIRASHII=self.removeIndex(MaddaIRASHII)
		
		newIRASHII=MaddaIRASHII[0:1]
 
		newIRASHII.remove_row(0)
		
 
		for eachIRASHII in MaddaIRASHII:
			"""
			"""
			
			IRASl= eachIRASHII["glon"] 
			IRASb= eachIRASHII["glat"] 
			IRASl=IRASl+360
			
			dup=False
			
			for eachWiseHII in maddaHII:
				
				radiusDegree= eachWiseHII[colRadiusWISE]/60./60. 
				HIIl = eachWiseHII[colL] 
				HIIb = eachWiseHII[colB] 
				
				middleB=np.mean( [ HIIb,IRASb ] )
				
				IRAS_WISE_dis=np.sqrt(  (HIIl- IRASl )**2*( np.cos( np.radians( middleB ) ) )**2 + (HIIb - IRASb )**2    )
		 
				if IRAS_WISE_dis<radiusDegree:

				#if IRAS_WISE_dis<searchRadius:

					dup=True
					dupN=dupN+1
					break
			###
			
			if not dup:
				newIRASHII.add_row(eachIRASHII)
				
				
				
		print "In total duplicated IRAS and WISE ",dupN

		print "New",len( newIRASHII )
		
		
		
		 

		#get WC89 
		if 1:
			
			newIRASWC89=WC89TB[0:1]
	 
			newIRASWC89.remove_row(0)
			
			
			print WC89G210TB.colnames
			
			for eachIRASHII in WC89G210TB:
				"""
				"""
				
				IRASl= eachIRASHII["glon"] 
				IRASb= eachIRASHII["glat"] 
				#IRASl=IRASl+360
				
				dup=False
				
				for eachWiseHII in maddaHII:
					
					radiusDegree= eachWiseHII[colRadiusWISE]/60./60. 
					HIIl = eachWiseHII[colL] 
					HIIb = eachWiseHII[colB] 
					
					middleB=np.mean(  [ HIIb, IRASb  ] )
					
					IRAS_WISE_dis=np.sqrt(  (HIIl- IRASl )**2 *( np.cos( np.radians( middleB ) ) )**2+ (HIIb - IRASb )**2    )
					#print IRAS_WISE_dis,IRASl
					if IRAS_WISE_dis<radiusDegree:
	
					#if IRAS_WISE_dis<searchRadius:
	
						dup=True
						dupN=dupN+1
						break
				###
				
				
				for eachIRASHIIA in newIRASHII:
					
					radiusDegree= 1./60.  
					HIIl = eachIRASHIIA["glon"]+360 #eachWiseHII[colL] 
					HIIb =eachIRASHIIA["glat"]  #eachWiseHII[colB] 
					
					#middleB=np.mean(  HIIb-)
					middleB=np.mean(  [ HIIb, IRASb  ] )

					
					IRAS_WISE_dis=np.sqrt(  (HIIl- IRASl )**2*( np.cos( np.radians( middleB ) ) )**2 + (HIIb - IRASb )**2    )
					#print IRAS_WISE_dis,IRASl
					if IRAS_WISE_dis<radiusDegree:
	
					#if IRAS_WISE_dis<searchRadius:
	
						dup=True
						dupN=dupN+1
						break
				
				
				
				if not dup:
					
					print eachIRASHII["pscname"]
					newIRASWC89.add_row(eachIRASHII)

		
		#if 1:

			#print len(newIRASWC89),"new IRAS WC89?"
		
		print newIRASWC89
		
 
		aaaaa
		cubeG210Raw12 = SpectralCube.read(self.COFile12)  
		cubeG21012  = cubeG210Raw12.with_spectral_unit(u.km / u.s)


		cubeMaskCO12Raw = SpectralCube.read(CO12Mask)  
		cubeMaskCO12  = cubeMaskCO12Raw.with_spectral_unit(u.km / u.s)

 


		cubeG210Raw13 = SpectralCube.read(self.COFile13)  
		cubeG21013  = cubeG210Raw13.with_spectral_unit(u.km / u.s)



		velo13= cubeG21013.spectral_axis
		velo12= cubeG21012.spectral_axis
		veloMask=cubeMaskCO12.spectral_axis

		d=DS9("drawWiseHII")

		doHII=myHII()

		tempCO13Region="tempCO13Region.reg"
		
 
		
		maddaHII=self.removeIndex(maddaHII)
		
		maddaHII.sort( [colL,  colB  ] )
		
		

		if 1:# use DS9 WISE
			useCO12=False
			
			if useCO12:
				d.set("file "+self.COFile12)

			else:
				d.set("file "+self.COFile13)
			for eachWiseHII in maddaHII:
	 			sourceName=  eachWiseHII[colNameWISE] 
	 			
	 			
	 			if source is not None:
	 				
	 				if sourceName!=source:
	 					continue
	 			
	 			
	 			wiseRow=doHII.getRowByName(sourceName)
 
				wiseRow[myHII.sourceName]= sourceName
	 			
				radiusDegree= eachWiseHII[colRadiusWISE]/60./60. 
				HIIl = eachWiseHII[colL] 
				HIIb = eachWiseHII[colB] 
				
				
				
				
				if wiseRow[ myHII.CO13L]  ==0:  #
					
					region_str = "galactic; circle({},{}, {})".format(  HIIl, HIIb,radiusDegree) 
					
				else: #show 13 if it has 12
 
 
					region_str = "galactic; circle({},{}, {})".format(  wiseRow[myHII.CO13L], wiseRow[myHII.CO13B],wiseRow[myHII.CO13Radius]) 

					
					
				print HIIl,HIIb ,radiusDegree
				
				print "Previous VLSR:",		eachWiseHII[colVLSRWISE] , "#############",eachWiseHII[colVLSRMolWISE]  
				
				
				wiseRow[ myHII.sourceName]= sourceName

				wiseRow[ myHII.l]= HIIl
				wiseRow[ myHII.b]= HIIb
				wiseRow[ myHII.radius]= radiusDegree
				wiseRow[ myHII.catalog]=   eachWiseHII[WISEcatalog]  


				d.set('regions delete all' )

				d.set('regions',region_str)

				#plot
				
				#sub_cube13 = cubeG21013.subcube_from_ds9region(region_str)
				#sp13=sub_cube13.mean(axis=(1,2))
				
			#	d.set('plot new LINE  stdin   ')
				
				#d.set('plot new' )
				defaultV=wiseRow[ myHII.velocity]


				vlsr=raw_input("Molecular velocity of  WISE {} ({} km/s)  : ".format(sourceName ,defaultV ) )
				vlsr=str.strip(vlsr)
				
				
 


				if wiseRow[ myHII.CO13L]  ==0:
					wiseRow[ myHII.velocity]= -999.

				
				try :
					vlsr=float(vlsr)
					wiseRow[ myHII.velocity]= float(vlsr) 
					if useCO12:
						wiseRow[ myHII.note]= "12CO"
					else:
						wiseRow[ myHII.note]= ""
				except :
					pass
					#if 
					
					#vlsr= -999.
					
					

				#get the region
				
				d.set('regions system wcs')
				d.set('regions sky galactic')
		
				d.set('regions save '+tempCO13Region)
				#d.set('quit')
 
				regions = pyregion.open(tempCO13Region)

				eachRegion= regions[0]
				
				reginCor= eachRegion.coord_list
				
				wiseRow[myHII.CO13L]= reginCor[0]
				
				wiseRow[myHII.CO13B]= reginCor[1]

				wiseRow[myHII.CO13Radius]= reginCor[2]
 
				
				doHII.addRowToTB(wiseRow)   #=doHII.getRowByName(sourceName)

				raw_input("Next region?")
 
				
			
		#print newIRASHII.colnames
		useCO12=False
		if 0:# use DS9 IRAS
			
			if useCO12:
				d.set("file "+self.COFile12)
			else:
				d.set("file "+self.COFile13)

			
			#draw 
			
			#draw all wise 
			
 
			
			for eachIRASHII in newIRASHII: #
				
				d.set('regions delete all' )

  
	 			sourceName=  eachIRASHII["pscname"] 
	 			
	 			wiseRow=doHII.getRowByName(sourceName)
 
				wiseRow[myHII.sourceName]= sourceName
	 			
				radiusDegree= 3/60. 
				HIIl = eachIRASHII["glon"] +360.
				HIIb = eachIRASHII["glat"] 
				
 
				if wiseRow[ myHII.CO13L]  ==0:  #
				
					region_str = "galactic; circle({},{}, {}) # color=red".format(  HIIl, HIIb,radiusDegree) 
					
				else: #show 13 if it has 12
 
 
					region_str = "galactic; circle({},{}, {}) # color=red".format(  wiseRow[myHII.CO13L], wiseRow[myHII.CO13B],wiseRow[myHII.CO13Radius]) 

					
					
				print HIIl,HIIb ,radiusDegree
				
				#print "Previous VLSR:",		eachWiseHII[colVLSRWISE] , "#############",eachWiseHII[colVLSRMolWISE]  
				
				wiseRow[ myHII.sourceName]= sourceName

				wiseRow[ myHII.l]= HIIl
				wiseRow[ myHII.b]= HIIb
				wiseRow[ myHII.radius]= radiusDegree
				#wiseRow[ myHII.catalog]=   eachWiseHII[WISEcatalog]  



				d.set('regions',region_str)

				#plot
				
				#sub_cube13 = cubeG21013.subcube_from_ds9region(region_str)
				#sp13=sub_cube13.mean(axis=(1,2))
				
			#	d.set('plot new LINE  stdin   ')
				
				#d.set('plot new' )
				defaultV=wiseRow[ myHII.velocity]


				vlsr=raw_input("Molecular velocity of  IRAS{} ({} km/s)  : ".format(sourceName ,defaultV ) )
				vlsr=str.strip(vlsr)
				
		 
				if wiseRow[ myHII.CO13L]  ==0:
					wiseRow[ myHII.velocity]= -999.

				
				try :
					vlsr=float(vlsr)
					wiseRow[ myHII.velocity]= float(vlsr) 
					if useCO12:
						wiseRow[ myHII.note]= "12CO"
					else:
						wiseRow[ myHII.note]= ""

						

				except :
					pass
					#if 
					
					#vlsr= -999.
					
					

				#get the region
				
				d.set('regions system wcs')
				d.set('regions sky galactic')
		
				d.set('regions save '+tempCO13Region)
				#d.set('quit')
 
				regions = pyregion.open(tempCO13Region)

				eachRegion= regions[0]
				
				reginCor= eachRegion.coord_list
				
				wiseRow[myHII.CO13L]= reginCor[0]
				
				wiseRow[myHII.CO13B]= reginCor[1]

				wiseRow[myHII.CO13Radius]= reginCor[2]
 
				
				doHII.addRowToTB(wiseRow)   #=doHII.getRowByName(sourceName)

				raw_input("Next region?")


 
			
			
		#print WC89.colnames
		useCO12=False
		if 0 :# use DS9 WC89
			
			if useCO12:
				d.set("file "+self.COFile12)
			else:
				d.set("file "+self.COFile13)

			
			#draw 
			
			#draw all wise 
			
 
			
			for eachIRASHII in newIRASWC89: #
				
				d.set('regions delete all' )

  
	 			sourceName=  eachIRASHII["pscname"] 
	 			
	 			
	 			print sourceName,"?????????????????"
	 			
	 			
	 			wiseRow=doHII.getRowByName(sourceName)
 
				wiseRow[myHII.sourceName]= sourceName
	 			
				radiusDegree= 3/60. 
				HIIl = eachIRASHII["glon"]  
				HIIb = eachIRASHII["glat"] 
				
 
				if wiseRow[ myHII.CO13L]  ==0:  #
				
					region_str = "galactic; circle({},{}, {}) # color=red".format(  HIIl, HIIb,radiusDegree) 
					
				else: #show 13 if it has 12
 
 
					region_str = "galactic; circle({},{}, {}) # color=red".format(  wiseRow[myHII.CO13L], wiseRow[myHII.CO13B],wiseRow[myHII.CO13Radius]) 

					
					
				print HIIl,HIIb ,radiusDegree
				
				#print "Previous VLSR:",		eachWiseHII[colVLSRWISE] , "#############",eachWiseHII[colVLSRMolWISE]  
				
				wiseRow[ myHII.sourceName]= sourceName

				wiseRow[ myHII.l]= HIIl
				wiseRow[ myHII.b]= HIIb
				wiseRow[ myHII.radius]= radiusDegree
				#wiseRow[ myHII.catalog]=   eachWiseHII[WISEcatalog]  



				d.set('regions',region_str)

				#plot
				
				#sub_cube13 = cubeG21013.subcube_from_ds9region(region_str)
				#sp13=sub_cube13.mean(axis=(1,2))
				
			#	d.set('plot new LINE  stdin   ')
				
				#d.set('plot new' )
				defaultV=wiseRow[ myHII.velocity]


				vlsr=raw_input("Molecular velocity of  IRAS{} ({} km/s)  : ".format(sourceName ,defaultV ) )
				vlsr=str.strip(vlsr)
				
		 
				if wiseRow[ myHII.CO13L]  ==0:
					wiseRow[ myHII.velocity]= -999.

				
				try :
					vlsr=float(vlsr)
					wiseRow[ myHII.velocity]= float(vlsr) 
					if useCO12:
						wiseRow[ myHII.note]= "12CO"
					else:
						wiseRow[ myHII.note]= ""

						

				except :
					pass
					#if 
					
					#vlsr= -999.
					
					

				#get the region
				
				d.set('regions system wcs')
				d.set('regions sky galactic')
		
				d.set('regions save '+tempCO13Region)
				#d.set('quit')
 
				regions = pyregion.open(tempCO13Region)

				eachRegion= regions[0]
				
				reginCor= eachRegion.coord_list
				
				wiseRow[myHII.CO13L]= reginCor[0]
				
				wiseRow[myHII.CO13B]= reginCor[1]

				wiseRow[myHII.CO13Radius]= reginCor[2]
 
				
				doHII.addRowToTB(wiseRow)   #=doHII.getRowByName(sourceName)

				raw_input("Next region?")


			
			
			
			
			
			#print eachWiseHII["WISE Name"], eachWiseHII["VLSR<br>(km/s)"], "#####",eachWiseHII["VLSR (Mol.)<br>(km/s)"], eachWiseHII["Membership"]



		print "Total number of HII reions?",len( MaddaIRASHII)+len(MaddaIRASHII)
		

	def fittingCO13ForHII(self,draw=True):
		"""
		
		fitting Gaussian profile for spectrals
		"""

		doHII=myHII()
		
		allHII= doHII.getAllTB()

		cubeG210Raw13 = SpectralCube.read(self.COFile13)  
		cubeG21013  = cubeG210Raw13.with_spectral_unit(u.km / u.s)

		velo13=cubeG21013.spectral_axis 



		#f#ig, ax = plt.subplots()

		from pylab import rcParams
		rcParams['figure.figsize'] = 20, 12
		
		fig = plt.figure(   )

 		for eachHII in allHII:
			
			HIIName=eachHII[myHII.sourceName]
			
			
			wiseRow=doHII.getRowByName(HIIName)
 
			if eachHII[myHII.velocity]<=0:
				
				print "Skipping "+ HIIName

				continue
			
			
			l,b,r=eachHII[myHII.CO13L] , eachHII[myHII.CO13B], eachHII[myHII.CO13Radius]
 	

			region_str13 = "galactic; circle({},{}, {})  ".format( l,b,r) 

			sub_cube = cubeG21013.subcube_from_ds9region(region_str13) 

			sp13=sub_cube.mean(axis=(1,2))
			
			
			
			print "Processing "+ HIIName
			
			#print "at {} km/s".format( eachHII[myHII.velocity])
			
			
			
			
		#for 

		#region_str = "galactic; circle({},{}, {})".format(  HIIl, HIIb,radiusDegree) 


			g_init = models.Gaussian1D(amplitude=1., mean=eachHII[myHII.velocity] , stddev=1.)
			fit_g = fitting.LevMarLSQFitter()
			g = fit_g(g_init, velo13, sp13)


			#print dir(g )
			#print g.parameters
			itPeak=round(  g.parameters[0],2)

			fitV=round(  g.parameters[1],2)
			fitStd=round(  g.parameters[2],2)

			
			print "Equal?", fitV, '-----------',eachHII[myHII.velocity]
			 
			wiseRow[myHII.CO13VFit]= fitV
			wiseRow[myHII.vStd]= fitStd
			wiseRow[myHII.CO13Peak]= itPeak

			#print wiseRow
			#plot
			
			ax=plt.subplot(111 )
 
 			ax.plot(velo13,  sp13,lw=2,color='b')
 			
 			
 			ax.plot(velo13,g(velo13),'r--',lw=1.5,label=HIIName )
 			
 			ax.legend(fontsize=14)
 			
 			if draw:
				plt.show()
			#plt.show(figsize=(20,12))			
		
		
		
		
		
		
		
		
		
		
		
		
		
			doHII.addRowToTB(wiseRow)   #=doHII.getRowByName(sourceName)




	def printHII(self):
		
		"""
		"""

		doHII=myHII()
		doHII.printHIILatex()


 

	def suYang(self):
		
		"""
		Draw a figure for Su yang
		"""
		trimMix=1
		trimMax=20
		
		CO12RFile= "maskMomentLocal.fits"
		CO12GFile= "maskMomentPerseus.fits"
		CO12BFile= "maskMomentOuter.fits"
		
		
 
		Rdata,Rhead=doMadd.doFITS.readFITS( CO12RFile)
		Gdata,Ghead=doMadd.doFITS.readFITS( CO12GFile)
		Bdata,Bhead=doMadd.doFITS.readFITS( CO12BFile)

		rData=Rdata[0]
		gData=Gdata[0]
		bData=Bdata[0]
		
		
		#draw RGB$^{12}CO ()
			
		img = np.zeros((bData.shape[0], bData.shape[1], 3), dtype=float)

		img[:,:,0]=img_scale.sqrt(rData, scale_min=trimMix, scale_max=15)
		img[:,:,1]= img_scale.sqrt(gData, scale_min=trimMix, scale_max=trimMax)
		img[:,:,2]=img_scale.sqrt(bData, scale_min=trimMix, scale_max=8)
			
		plt.clf()
		#print dir(pywcsgrid2)
		#plt.rcParams['font.sans-serif'] = ['Helvetica']
		fig = plt.figure()
		ax=pywcsgrid2.subplot(111,header=WCS(Bhead))
		#ax.imshow(bData, origin="lower")
		#from astropy.visualization import make_lupton_rgb
		#image = make_lupton_rgb(rData, gData, bData, stretch=0.5)
		#rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		ax.set_ticklabel_type("absdeg", "absdeg")
		
		#ax.set_xlabel(r'$l$')
		
		#rc('text.latex', preamble='\usepackage{xcolor}')
		#plt.rcParams['font.sans-serif'] = ['Helvetica']
		#ax.imshow(bData,origin="lower",cmap="jet")zz
		#trimMax=10
		#trimMix=0
		#i = imshow_rgb(ax, self.doFITS.mytrim(rData,trimMix,trimMax) ,self.doFITS.mytrim(gData,trimMix,trimMax) , self.doFITS.mytrim(bData,trimMix,trimMax) , origin="lower", interpolation="none",vmin=trimMix,vmax=trimMax )
		i = imshow_rgb(ax,  img[:,:,0]  , img[:,:,1]  ,  img[:,:,2]  , origin="lower", interpolation="none" )

		#ax.imshow(img[:,:,0],interpolation="none",vmin=0.3,vmax=3 )
 
		#i = imshow_rgb(ax, img[:,:,0] ,img[:,:,1] ,img[:,:,2] , origin="lower", interpolation="none",vmin=2,vmax=15 )
		#axins = zoomed_inset_axes(ax,   loc=3)
		#axins=inset_axes(ax, width="32%", height="25%", loc=1)

		#if "12" in outFileName:
		patchLine=mpatches.Patch(color='None',label=r"$^{12}$CO ($J=1\rightarrow 0$)")
		
 
		patchRed=mpatches.Patch(color='red',label=r"0-14 km s$^{-1}$ (the Local Arm)")
		patchGreen=mpatches.Patch(color='green',label=r"14-38 km s$^{-1}$ (the Perseus Arm)")
		patchBlue=mpatches.Patch(color='blue',label=r"38-70 km s$^{-1}$ (the Outer Arm)")
		
		patches=[patchLine,patchRed,patchGreen,patchBlue]
		
		
		l=ax.legend(handles=patches,loc=2,handlelength=0.5,handleheight=0.3,fontsize=7,framealpha=0.6)
		
		cs=['black','darkred','green','blue' ]
		for countN,text in enumerate(l.get_texts()):
			text.set_color(cs[countN])
		
 
		#plt.subplots_adjust(left=0.02, right=1.05, top=0.99, bottom=0.05)

		#draw HII regions
		
		WISEcat=self.WISEHIICat.copy()
		
		#print WISEcat.colnames
		
		colL="GLong<br>(deg.)"
		colB="GLat<br>(deg.)"
		
		colRadius='Radius<br>(arcsec.)'
		
		lMin,lMax=min( self.G210LRange), max( self.G210LRange), 
		
		
		maddaHII=self.filterTBByRange(WISEcat,colL, [lMin,lMax ])

		maddaHII=self.filterTBByRange(maddaHII,colB,[-5.,5.])
		
		
		angles=np.linspace(0,2*np.pi,50) 

		if 0:
			for eachWISEHII in  maddaHII:
				
				
				#radius
				radius=eachWISEHII[colRadius]/3600.
				
				if radius<5/60.:
					continue
				
				centerL,centerB=eachWISEHII[colL],eachWISEHII[colB]
					
				xx=[]
				yy=[]
				for eachAng in angles:
					
					xx.append(centerL+radius*np.cos(eachAng) )
					yy.append(centerB+ radius*np.sin(eachAng)) 
		
				#ax["gal"].plot(xx,yy,'-',lw=0.5 ,color='white', dashes=(1, 0.5) )			
				ax["gal"].plot(xx,yy, lw=0.4 ,color='white'  )			
	
				#ax["gal"].scatter(eachWISEHII[colL],eachWISEHII[colB],color= "white",s=10,facecolors='None',lw=0.4)
	
		if 1:
			
			circles=[circle(214.3203262,3.4995615,1173.056),
			circle(217.6312840,-0.4818238,1242.016),
			circle(217.1470939,-2.3896891,824.634),
			circle(218.4923009,-2.7871671,693.873),
			circle(216.6762937,-0.2353826,1938.413),
			circle(210.4175060,-3.2235031,1138.628),
			circle(219.1551333,-4.0490156,1524.884),
			circle(213.8966435,1.5981049,1378.063),
			circle(215.0374080,4.5926157,1015.330),
			circle(212.0215331,-1.2978158,597.795),
			circle(211.1602321,-2.2231002,2051.810),
			circle(216.8013072,0.4425537,1011.627),
			circle(217.1891190,0.1542092,742.830) 
			]
			
			
			for eachCircle in circles:
				
				centerL,centerB,radius=eachCircle
					
				xx=[]
				yy=[]
				for eachAng in angles:
					
					xx.append(centerL+radius*np.cos(eachAng) )
					yy.append(centerB+ radius*np.sin(eachAng)) 
		
				#ax["gal"].plot(xx,yy,'-',lw=0.5 ,color='white', dashes=(1, 0.5) )			
				ax["gal"].plot(xx,yy, lw=0.4 ,color='white' ,alpha=0.8 )			
	

			

		ax.axis[:].major_ticks.set_color("w")
			
 
			
		#set lRange
		
		ax.set_xlim(0,1231)
			
			
		plt.tight_layout(pad=0)
 



		fig.savefig( 'suyang.pdf', bbox_inches="tight")
		fig.savefig( 'suyang.png', bbox_inches="tight")

		#fig.savefig(outFileName+'.png', bbox_inches="tight",dpi=300)

	def drawFITSCompare(self,parameter,FITSList):
		
		"""
		the FITSList is generated by starlink
		"""
		
		saveName="./figures/"+parameter+"fitsCompare.png"
		fig, axs = plt.subplots( figsize=(15,11) )  
		dv=		0.166040420532 

 
		AX = gridspec.GridSpec(4,3)
		AX.update(wspace = 0.25, hspace = 0.25)
		
		observedFITS="G216_13.fits"
		
		FITSList.insert(0,observedFITS )
		
 
		for i,FITSName in enumerate( FITSList):
			
			
			rowIndex= (i/3)
			colIndex= i%3
			
			#print i,rowIndex, colIndex
			
			data,header=myFITS.readFITS(FITSName)
			axMy=pywcsgrid2.subplot(AX[rowIndex,colIndex],header=WCS(header))  
			inTdata=np.sum(data,axis=0)*dv
			axMy.imshow(inTdata,origin='lower', vmin=0.01 ,vmax=10 ,cmap='jet')

  
		
		
		
		
		
		plt.savefig( saveName,bbox_inches='tight',dpi=600)

		
		
		

	def  examineTestPara(self,parameter):
		"""
		
		drawAll figures by para
		"""
		doCore=myCore()
		doCMF=myCMF()
		TBList,FITSList= doCore.getTBFITSNameListByPara(parameter)


		self.drawFITSCompare(parameter,FITSList)

		  
		allCMF=doCMF.getEmptyTB()
 		
		saveIndex=0
		
		
		for eachTBName in TBList:
			
			saveName="./figures/G216"+parameter+str(saveIndex)
			eachTB=Table.read(eachTBName )
			#CMFRows.append( self.getMassFromGauSum(eachTB)  )
			allCMF.add_row(self.getMassFromGauSum(eachTB,saveName) )
			saveIndex=saveIndex+1
		
		#print allCMF
		
 
		#draw founr 
		
		fig, axs = plt.subplots( figsize=(15,11) )  
		
 
		AX = gridspec.GridSpec(2,3)
		AX.update(wspace = 0.25, hspace = 0.25)
		
		axNclump  = plt.subplot(AX[0,0])
		axTotal = plt.subplot(AX[0,1])
		axRatio = plt.subplot(AX[0,2])
		
		
		axAlpha1  = plt.subplot(AX[1,1])

		
		axAlpha2 = plt.subplot(AX[1,2])
		axMturn = plt.subplot(AX[1,0])


		
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		######################
		#axNclump= #plt.subplot2grid((2,2), (0, 0) ,rowspan=1,colspan=1)
		NclumpG216=allCMF[myCMF.Nclump]
		axNclump.plot(NclumpG216 , 'g.'  )
		#axNclump.set_aspect("equal")

		#axNclump.set_xlabel("Clump number of G216")
		axNclump.set_ylabel("Clump number of G216")

		
		######################
		#axAlpha1=plt.subplot2grid((2,2), (0, 1) ,rowspan=1,colspan=1)
		#axAlpha1.errorbar(alpha1ListG216,alpha1ListW345,xerr=alpha1ErrorListG216, capsize=1.0, yerr= alpha1ErrorListW345,marker='o', ms=2.5, ls = 'none',color='blue',lw=0.4)
		#axAlpha1.errorbar( allCMF[myCMF.Nclump],  capsize=1.0, yerr= alpha1ErrorListW345,marker='o', ms=2.5, ls = 'none',color='blue',lw=0.4)

		alpha1ListG216=allCMF[myCMF.alpha1]
		axAlpha1.plot(alpha1ListG216 , 'b.' )
 
		#axAlpha1.set_aspect("equal")
		#axAlpha1.set_xlabel(r"$\alpha_1$ of G216")
		axAlpha1.set_ylabel(r"$\alpha_1$ of G216")




		#axAlpha2=plt.subplot2grid((2,2), (1, 0) ,rowspan=1,colspan=1)
		#axAlpha2.errorbar(alpha2ListG216,alpha2ListW345,xerr=alpha2ErrorListG216, capsize=1.0, yerr= alpha2ErrorListW345,marker='o', ms=2.5, ls = 'none',color='blue',lw=0.4,zorder=1)
	#	axAlpha2.set_aspect("equal")
		#axAlpha2.errorbar([G216Alpha2mean] ,[W345Alpha2mean],xerr=[G216Alpha2error], capsize=1.0, yerr= [W345Alpha2error],marker='o', ms=2.5, ls = 'none',color='red',lw=0.4)

		alpha2ListG216=allCMF[myCMF.alpha2]
		axAlpha2.plot(alpha2ListG216 , 'b.' )
 
		#axAlpha2.set_xlabel(r"$\alpha_2$ of G216")
		axAlpha2.set_ylabel(r"$\alpha_2$ of G216")



		#axMturn=plt.subplot2grid((2,2), (1, 1) ,rowspan=1,colspan=1)
		#axMturn.errorbar(MturnListG216,MturnListW345,xerr=MturnErrorListG216, capsize=1.0, yerr= MturnErrorListW345,marker='o', ms=2.5, ls = 'none',color='blue',lw=0.4)
		#axMturn.set_aspect("equal")
		
		MturnListG216=allCMF[myCMF.Mturn]
		axMturn.plot(MturnListG216 , 'b.' )
 
		
		
		
		#axMturn.set_xlabel(r"$M_{\rm turn}$ of G216 ($M_\odot$)")
		axMturn.set_ylabel(r"$M_{\rm turn}$ of G216 ($M_\odot$)")
		

		#totalMassG216= fitTBG216[self.totalMass]
		#totalMassW345= fitTBW345[self.totalMass]

		#gBoundG216= fitTBG216[self.MvirialRatio]
		#gBoundW345= fitTBW345[self.MvirialRatio]
		
		#axTotal.scatter( )
		totalMassG216=allCMF[myCMF.totalMass]
		axTotal.plot(totalMassG216/1e4,  'g.' )
		
		#axNclump.set_aspect("equal")
		
		#axTotal.set_xticks([1.0,1.5,2.0])
		#axTotal.set_yticks([4,6,8,10,12])

		
		#locs, labels = plt.xticks()            # Get locations and labels
		#plt.xticks(ticks, [labels], **kwargs)  # Set locations and labels

		#xticks(np.arange(0, 1, step=0.2))


		#axTotal.set_xlabel(r" Total $M_{\rm LTE}$ of G216 ($10^4\ M_\odot$)")
		axTotal.set_ylabel(r" Total $M_{\rm LTE}$ of G216 ($10^4\ M_\odot$)")

		gBoundG216=allCMF[myCMF.MvirialRatio]

		axRatio.plot(gBoundG216  , 'g.', label='Ratio of mass in gravitationally bound')
		#axRatio.set_xticks([0.1,0.2,0.3,0.4])
		#axRatio.set_yticks([0.6,0.7,0.8,0.9,1.0])
		
		#axRatio.set_xlabel(r"$M_{\rm LTE}$ ratio in gravitationally bound (G216)")
		axRatio.set_ylabel(r"$M_{\rm LTE}$ ratio in gravitationally bound (G216)")
		#axRatio.legend( loc=4)
		#plt.tight_layout()
		#plt.savefig(parameter+"TEST.pdf",bbox_inches='tight')

		plt.savefig("./figures/"+parameter+"TEST.png",bbox_inches='tight',dpi=600)


	def creatOutAndPeseusMaskFITS(self,doCO12=True):
		"""
		"""
		
		#create mask fits according to 
		#v=a*l+b a= 1.65580467802 b= -318.196012159, out arm line
		
		#perseus

		#read 

		allMaskData=None

		if doCO12:

			doFile=self.COFile12
			allMaskData,aaa=myFITS.readFITS("/home/qzyan/WORK/projects/maddalena/data/duchampG210/mosaic_U.MASK.fits") 

		else:
			doFile=self.COFile13
			allMaskData,aaa=myFITS.readFITS("/home/qzyan/WORK/projects/maddalena/data/duchamp13CO/mosaic_L.MASK.fits") 

		fitsData,fitsHead=myFITS.readFITS(doFile)
		
		coWCS=WCS(fitsHead)
		
		a= 1.65580467802
		b= -318.196012159
		
		maskFITSDataOuter= fitsData*0
		maskFITSDataPerseus= fitsData*0+1

		#v > a*l+b , should be outarm
		
		#self.V_PerseusArm=15.
		Nz,Ny,Nx=fitsData.shape
		
		
		vIndex= np.arange(Nz)
		
		vIndex=vIndex/1000.
		
		aa,aa,vValue=coWCS.wcs_pix2world(   0,0,vIndex,0)
		
		#for i in range(Nz):
		correspondL= (vValue-b)/a 
		
		correspondLIndex,aa,aa=coWCS.wcs_world2pix(  correspondL,0,vValue*1000.,0)

		correspondLIndex=map(round, correspondLIndex )
		correspondLIndex=map(int, correspondLIndex )

		
		#print correspondLIndex
		for k in range(Nz):
			cutIndexL=min( [correspondLIndex[k], Nx-1])
			cutIndexL=max( [ cutIndexL, 0])

			#print cutIndexL
			if cutIndexL==Nx-1:
				continue
			
			maskFITSDataOuter[k,:,cutIndexL:]= 1 
			maskFITSDataPerseus[k,:,cutIndexL:]= 0

		aa,aa,perseusIndex=coWCS.wcs_world2pix(  215. , 0 , self.V_PerseusArm*1000 ,0)

		cutVIndex= int(round( perseusIndex ))
		
		#maskFITSDataOuter[0:cutVIndex,:,:]=0
		
		maskFITSDataPerseus[0:cutVIndex,:,:]=0

		#save 

		if doCO12:
			outFile=self.outMask12
			perFile=self.perMask12
		else:
			outFile=self.outMask13
			perFile=self.perMask13
			
			
		try:
			os.remove( self.dataPath+outFile )
		except:
			pass
	
		try:
			os.remove( self.dataPath+ perFile )
		except:
			pass


		

		fits.writeto(self.dataPath+outFile ,maskFITSDataOuter*fitsData*allMaskData,header=fitsHead)
		fits.writeto(self.dataPath+ perFile ,maskFITSDataPerseus*fitsData*allMaskData,header=fitsHead)

			#print    #(vValue[k]-b )/a
		
		#for eachV in maskFITSDataOuter:
			#print 
		
		#print fitsData.shape
		

	def intOutAndPer(self):
		
		"""
		integrate outer and perseus mar
		"""
		# int CO12
		
		CO12RGB="CO12RGB"
		#CO13RGB="CO13RGB"
		outFileName= CO12RGB
		CO12RFile=outFileName+"_R.fits"
		CO12GFile=outFileName+"_G.fits"
		CO12BFile=outFileName+"_B.fits"
		
 		
		
		Gdata,Ghead=doMadd.doFITS.momentFITS(self.dataPath+self.perMask12, [0, 70],0,outFITS=CO12GFile)
		Bdata,Bhead=doMadd.doFITS.momentFITS(self.dataPath+self.outMask12, [0, 70],0,outFITS=CO12BFile)

		CO12RGB="CO13RGB"
		#CO13RGB="CO13RGB"
		outFileName= CO12RGB
		CO12RFile=outFileName+"_R.fits"
		CO12GFile=outFileName+"_G.fits"
		CO12BFile=outFileName+"_B.fits"
		
 		
		
		Gdata,Ghead=doMadd.doFITS.momentFITS(self.dataPath+self.perMask13, [0, 70],0,outFITS=CO12GFile)
		Bdata,Bhead=doMadd.doFITS.momentFITS(self.dataPath+self.outMask13, [0, 70],0,outFITS=CO12BFile)


	def getCloudNameByLB(self, l,b):
	
		#if b>=0:
	
			lStr= str(l)
			
			bStr="{:+f}".format(b)
			
			
			if '.' in lStr:
			
				lStrPart1,lStrPart2=lStr.split('.')
			
			else:
				lStrPart1 =lStr
				lStrPart2='0'
	
	
			if '.' in bStr:
	
				bStrPart1,bStrPart2=bStr.split('.')
			else:
				bStrPart1 =bStr
				bStrPart2='0'
				
				
			lStr=lStrPart1+'.'+lStrPart2[0:1]
	
			
			bStr=bStrPart1+'.'+bStrPart2[0:1]
	
			
			lStr=lStr.zfill(5)
			
			
			#bStr="{:+.1f}".format(b)
			bStrNumberPart=bStr[1:]
			bStr=bStr[0:1]+  bStrNumberPart.zfill(4)
			
			cName="G{}{}".format(lStr,bStr)
			
			return cName

	def getLBVByID(self,ID):
		dendroCat="/home/qzyan/WORK/projects/maddalena/markedCatalog.fit"
		
		TB=Table.read(dendroCat)
		
		foundRow=None
		
		for eachRow in TB:
			
			if eachRow["_idx"]==ID:
				
				foundRow=eachRow
				break
		
		l=eachRow["x_cen"]
		b=eachRow["y_cen"]
		v=eachRow["v_cen"]

		return l,b,v
	def getCloudNameByID(self, ID):

		"""
		getCloudName (generated by l and b) by Identities
		"""
		l,b,v=self.getLBVByID(ID)

		return self.getCloudNameByLB( l,b)
 

	def getNDict(self):
		"""
		get the dict for the number of clouds
		
		Then number is not the ID given by astrodendro, but the order number to display the clouds easily
		
		"""
		
		returnDic={}
		
		G210OrderFile="G210Order.txt"
		
		with open(G210OrderFile) as f:
			lines = f.readlines()
		
		for eachLine in lines:
			splitL= eachLine.split("&")
		
		
			if len(splitL) != 2:
				continue
				
			
			splitL=map( str.strip, splitL   )
			N,cloudName= splitL
			
			cloudName=cloudName.replace("$","")
			N=int(N)
			returnDic[cloudName ]=N
			
		return returnDic

 	def ZZZ(self):
		
		pass




#/home/qzyan/WORK/projects/maddalena/
		

		
doMadd=myG210()



if  0: #draw face one
	"""
	
	"""

	#coreTB=Table.read("fellwalkercoreCat_1.fit")

	coreTB=Table.read("allCoresGaussClumpscoreCat.fit")

	#coreTB=Table.read("GaussClumpscoreCat_1_limit_1.fit")

	print len(coreTB)

	doMadd.drawArmsNew(coreTB)




if  0:  #draw three color # one by one ,the second could be affected by the frist one
	#run two times to keep 
 
	doMadd.drawThreeColor(doMadd.COFile13,"CO13RGB",doMoment=False,trimMix=-0,trimMax=6)
	doMadd.drawThreeColor(doMadd.COFile13,"CO13RGB",doMoment=False,trimMix=0,trimMax=6)

	doMadd.drawThreeColor(doMadd.COFile12,"CO12RGB",doMoment=False,trimMix=0,trimMax=20)
	doMadd.drawThreeColor(doMadd.COFile12,"CO12RGB",doMoment=False,trimMix=0,trimMax=20)





if 0:

	print cloudName
	print orderDict[cloudName]

if 0:
	doMadd.creatOutAndPeseusMaskFITS()

	doMadd.creatOutAndPeseusMaskFITS(doCO12=False)





if 0: #draw compare
	doMadd.drawCompareCOFUGINDame()





if 0:
	doMadd.intOutAndPer()



if 0:

	#doMadd.checkAllCores(testRandom=True,source="G216",saveTBName="G216CMFnoUpper.fit",startI=0,endI=200)
	
	doMadd.checkAllCores(testRandom=True,source="W345",saveTBName="W345CMFnoUpper.fit",startI=0,endI=200)


if 0:
	#examine G214

	doMadd.drawG214Outflow()




if 0:
	for k,v in myCore.gaussPara.iteritems():
		doMadd.examineTestPara(k)

if 0:
	
	doMadd.examineTestPara('S0')

if 0:
	doMadd.plotNegative(moment=False)




if 0:
	doMadd.suYang()
	#doMadd.suYang()

if 0: #draw all regions
	#doMadd.drawHIIRegions("mosaic_U_M0.fits")
	doMadd.drawHIIRegions("outArmDuchamp.fits")

	#doMadd.drawHIIRegions("WISE22.fits")


if 0:
	 
	doMadd.examineHIIRegionBySp( )


if 0:
	doMadd.G214Multiwave()






if 0: #mass of arms

	#self.setLrange(axDSS2Red,WCS( headDSS2Red), [211.6,212.4])
	#self.setBrange(axDSS2Red,WCS( headDSS2Red), [-1.7, -0.9 ])
	lRange_Maddalena=   doMadd.G210LRange
	bRange_Maddalena=  doMadd.G210BRange
 

 

	#lRange_G214,bRange_G214,aaaa,aaaa = box(214.6240124,-1.8017578,4638.337,3909.934,0)
 
 

	vRange=[38,70]
	
	massOut=doMadd.calMassWithCO12( "./data/mosaic_U.fits",lRange_Maddalena,bRange_Maddalena,vRange, 5300.)
	
	print "Outer Arm {:.2f}e10^6".format(massOut/1e6)
 
	vRange=[14,38]
	massPerseus=doMadd.calMassWithCO12( "./data/mosaic_U.fits",lRange_Maddalena,bRange_Maddalena,vRange, 2500.)
 
	print "Perseus Arm {:.2f}e10^6".format(massPerseus/1e6)

	vRange=[0,14]
	massLocal=doMadd.calMassWithCO12( "./data/mosaic_U.fits",lRange_Maddalena,bRange_Maddalena,vRange, 2500.)
 
	print "Local Arm {:.2f}e10^6".format(massLocal/1e6)




if 0:
	doMadd.printHII()


if 0:
	##
	doMadd.fittingCO13ForHII(draw=False)






if 0:

	doMadd.checkAllCores(testRandom=True,source="G216",saveTBName="G216CMFnoUpperDefault.fit",startI=210,endI=211)
	
	#doMadd.checkAllCores(testRandom=True,source="W345",saveTBName="W345CMFnoUpperDefault.fit",startI=210,endI=211)
 


if 0:

	
	#doMadd.integralS284()
 	
	#
	
	doMadd.S284expand(use12=True)





if 0: #mass 214

	#self.setLrange(axDSS2Red,WCS( headDSS2Red), [211.6,212.4])
	#self.setBrange(axDSS2Red,WCS( headDSS2Red), [-1.7, -0.9 ])
	



	#lRange_G214,bRange_G214,aaaa,aaaa = box(214.6240124,-1.8017578,4638.337,3909.934,0)
	
	
	lRange_Maddalena=   self.lRange_G214
	bRange_Maddalena=  self.bRange_G214
 
	vRange=[22,33]
	
	mass214=doMadd.calMassWithCO12( "./data/mosaic_U.fits",lRange_Maddalena,bRange_Maddalena,vRange,2450.)
	
	print "Mass 214",mass214
 


	#mass 216
	
 
	
	lRange_Maddalena= self.lRange_Maddalena
	bRange_Maddalena= self.bRange_Maddalena
 
	vRange= self.vRange_Maddalena
	
	mass216=doMadd.calMassWithCO12( "./data/mosaic_U.fits",lRange_Maddalena,bRange_Maddalena,vRange,2450.)
	
	print "Mass 216",mass216
 


	#mass 284
 
	
	displayLrange=[211.5,212.7]
	displayBrange=[-1.8,-0.5]
	vRange=[38,52]
	
	mass284=doMadd.calMassWithCO12( "./data/S284_CO12.fits",displayLrange,displayBrange,vRange,5000.)
	
	print "Mass S284",mass284
	
	
	
	
	
	





	#doMadd.S284Multiwave()



if 0:
	#doMadd.integralS284()
	
	doMadd.channelS284(use12=True)


	#doMadd.S284Multiwave()




if 0: #de faultValues
	coreTest1="./gaussClumpTest/RANDOM_G216_1GaussClumpscore"
	doMadd.examineOneCore(coreTest1)

	#coreTest1="./gaussClumpTest/RANDOM_W345_210GaussClumpscore"
	#doMadd.examineOneCore(coreTest1)







 




if 0:
 	print "S287?"
	#doMadd.ReidA5(217.794,-0.023,   28.5 )


	print "G216?"
	#doMadd.ReidA5(216.596,-2.947,   23.5 )

	#print "G214?"
	#doMadd.ReidA5(214.624,-1.802,   27 )

	#print "G211?"
	#doMadd.ReidA5(211.784,2.451,   7.1 )
	print "S284?"
	#doMadd.ReidA5(212.1, -1.1,  45 )

	print "V838 Mon"
	#doMadd.ReidA5(217.797, 1.052, 49 )
	
	print "Tesing"
 	#doMadd.ReidA5(203, 0.1, 40 )



	
	print "W5"
	doMadd.ReidA5(137.222, 0.920, -36.6 )





if 0:
	
	doMadd.channelG216(14,38)




 


if 0: #draw ppvL
	
	#doMadd=myG210()
	#doMadd.drawPVL(doMadd.COFile13) 
	doMadd.drawPVL(doMadd.COFile12) 

	



if 0:
 	print "S287?"
	#doMadd.ReidA5(217.794,-0.023,   28.5 )


	print "G216?"
	#doMadd.ReidA5(216.596,-2.947,   23.5 )

	#print "G214?"
	#doMadd.ReidA5(214.624,-1.802,   27 )

	#print "G211?"
	#doMadd.ReidA5(211.784,2.451,   7.1 )
	print "S284?"
	#doMadd.ReidA5(212.1, -1.1,  45 )

	print "V838 Mon"
	#doMadd.ReidA5(217.797, 1.052, 49 )
	
	print "Tesing"
 	#doMadd.ReidA5(203, 0.1, 40 )






 















if  0: #examine virial mass, to see, if the virial mass is less than the LTE mass, then the core will be gravitational bound
	
	
	coreTest1="./gaussClumpTest/RANDOM_W345_210GaussClumpscoreCat_1.fit"
	doMadd.compareLTEVirial(coreTest1)

	coreTest1="./gaussClumpTest/RANDOM_G216_210GaussClumpscoreCat_1.fit"
	doMadd.compareLTEVirial(coreTest1)

 








	#doMadd.checkAllCores(testRandom=True,source="G216",saveTBName="G216CMF.fit",startI=0,endI=200)

	#doMadd.checkAllCores(testRandom=True,source="W345",saveTBName="W345CMF.fit",startI=0,endI=200)

	#doMadd.checkAllCores(testRandom=True,source="W345",saveTBName="W345CMF.fit")

	
	#coreTest1="./gaussClumpTest/RANDOM_G216_0GaussClumpscore"
	#doMadd.examineOneCore(coreTest1)
	
	#coreTest1="./gaussClumpTest/RANDOM_W345_6GaussClumpscore"
	#doMadd.examineOneCore(coreTest1)
	
	
	











	#source="W345"
	
	#i=7
	
	#coreTest1="./gaussClumpTest/RANDOM_{}_{}GaussClumpscore".format(source,i)
	#doMadd.examineOneCore(coreTest1)
	
	#coreTest1="./gaussClumpTest/RANDOM_G216_6GaussClumpscore"
	#doMadd.examineOneCore(coreTest1)
	#coreTBTest=Table.read(coreTest1+"Cat_1.fit") 
	#doMadd.getTexListFromTB(coreTBTest)
	















if 0:#examine the core mass of a possible G214
	
	massList=[]
	
	for i in range(200):
	
		examineTB="./gaussClumpTest/RANDOM_G216_{}GaussClumpscoreCat_1.fit".format(i)
		
		examineTB=Table.read(examineTB)
		
		examineTB.add_index("Peak1")
		examineTB=examineTB.loc["Peak1",214.45:214.51 ]

		#print len(examineTB),"????????????"
		if len(examineTB)==13: #Only one selected core
			
			Mass, Tex,Mvirial = doMadd.getMassForCore(examineTB)
			
			massList.append(Mass)
			continue


		
		examineTB.add_index("Peak2")
		examineTB=examineTB.loc["Peak2",-1.81:-1.8]
	
	
		#print len(examineTB),"????????????"
		if len(examineTB)==13: #Only one selected core
			
			Mass, Tex,Mvirial = doMadd.getMassForCore(examineTB)
			
			massList.append(Mass)
			continue
	
		examineTB.add_index("Peak3")
		examineTB=examineTB.loc["Peak3", 28000: ]
		#print len(examineTB),"????????????"
		if len(examineTB)==13: #Only one selected core
			
			Mass, Tex,Mvirial = doMadd.getMassForCore(examineTB)
			
			massList.append(Mass)
			continue
		#print examineTB
		#print len(examineTB),"What the hell"
 
	print len(massList)
	print np.mean(massList)




if 0:
	
	doMadd.drawPV('./data/G214CO13.fits',beginP=[201.22345,90.423686],endP=[149.63139,54.4214],widthPix=10) 










if 0:
	#doMadd.channelS284()

	#doMadd.integralS284()
	centerLB=[212.0432736, -1.2724404  ]
	
	r1=0.1498044
	
	r2=0.2347144
	
	doMadd.drawRingPV( doMadd.S284FITS12 ,centerLB,r1,r2,saveFITS="ringS284.fits") 
 
if 0:
	
	doMadd.calMassWithCO12( "FellwalerG216CO12fellwalkercoreCat_1.fit" )
	#doMadd.calMassWithCO12( "TEST_MODELLIM0_5GaussClumpscoreCat_1.fit" )

if 0:
	
	doMadd.fittingEliptical( [25.,27])








if 0:
	doMadd=myG210()

	#testing new core mass function
	#coreTB=Table.read("GaussClumpscoreCat_1_NoFig.fit") 
	
	coreTB=Table.read("GaussClumpscoreCat_1_limit_1.fit") 
	maddCoreTB=doMadd.filterMaddCore(coreTB)
	print len(maddCoreTB)

	#maddCoreTB=Table.read("G216GaussClumpscoreCat_1.fit") 

	print len(maddCoreTB)
	doMadd.getMassFromGauSum(maddCoreTB)


if 0:
	
	
	TB=Table.read('stars.fit')
	print len(TB)
	doMadd.outPutDiskStars(TB)
	#matching YSO with Gaia Failed.





if 0:
	pass
	#doMadd.cropG216FITS()



	#doMadd.getPPVFits( doMadd.COFile13)



if 0:
	pass
	#doMadd.ReidA5(217.278,-0.237,15.92)
	#doMadd.ReidA5(217.752,-0.148,27.17)
	#calsource 064807.2691 -035209.0678 25 0
	
	#local arm 
	#doMadd.ReidA5(217.752,-0.148,8)
	#calsource 064807.2691 -035209.0678 25 0
	
	
	#print doMadd.ReidA5(211.8493593, 1.8963028,8)

	#km distance of maddalens
	#print doMadd.ReidA5(216., -2.5   ,25.)
	
	
	#G216
	print "Reid?"
	doMadd.ReidA5(214.8,-1.9,    26.6)

	#! Source     Gal Long  Gal Lat    V_lsr     V_rev    Rev. D_k     +/-  
	#!              (deg)    (deg)    (km/s)    (km/s)     (kpc)      (kpc)
	# calsource    214.800  -1.900      26.6      27.9      2.49   0.82  -0.73
	
	





	#calsource 064807.2691 -035209.0678 25 0
	# calsource    216.000  -2.500      25.0      26.2      2.25   0.76  -0.69
	
	#doMadd.ReidA5(211.6, 2.5, 8 )
	#calsource 065754.8606 +021927.833 8 0
	
	#calsource    211.600   2.500       8.0       9.3      0.77   0.66  -0.60



	#alsource 
	
	#calculate S284, 
	#doMadd.ReidA5(  212.1 , -1.1,  45 )
	#calsource    212.5  -1.100      45.0      46.7      5.34   1.35  -1.14

	print "what?"










if 0:# test an eliptical ring seen by 13CO
	pass
	#
	doMadd=myG210()
	doMadd.examineRing()

if 0: #select cores belong to maddalena cores
	
	doMadd=myG210()
	coreTB=Table.read("GaussClumpscoreCat_1_NoFig.fit") 
	doMadd.filterMaddCore(coreTB)







if 0: #test cores
	#doMadd.drawCMF("./maddalenaCore/fellwalker1.fit")

	doMadd.drawCMF("./maddalenaCore/GaussClumps1.fit")


 
if 0:
	centerLB= [ 213.0972, -3.5569 ]
	
	pathAngle= 0. #degree
	
	length= 0.5
 
	vRange=[45.,56.]

	doMadd.drawOutFlow(centerLB,pathAngle,length,0.05,vRange,"CO13RGB_G.fits",survey="TestOne")


if 0:
	
	doMadd.drawAllOutflow()

if 0:
	doMadd.maddCMF()




if 0:
	doMadd.drawRegions()

if 0: #draw distance
	doMadd.drawDistance() 


	

if 0:
	
	pass
	
	#doMadd.getGaiaDR2Stars("mosaic_U_M0.fits", 3,10) 
	
	#doMadd.getGaiaDR2Stars("CO12RGB_R.fits", 3,10)   
	
	#doMadd.getGaiaDR2Stars("CO12RGB_B.fits", 3,10)   
	
	
	

if 0 : # examine the sitribution with repsect X, Y, better check V later
	doMadd.examineXYDis("mosaic_U_M0.fits")


if 0: #test pvslice
	#doMadd.drawPVslice([210.2262188,-1.333511053],[1,1],45,0.1)
	#doMadd.drawPVslice([210.465194,-2.334800181],[1,1],45,0.1)

	doMadd.drawPVslice([-146.9885,0.9862],[0.5,0.5],45,0.1)


if 0: #test download skyview survey
	
	doMadd.downLoadSkyView()
	
 
 


#/home/qzyan/WORK/projects/maddalena/mosaic_U_M0.fits

if 0:

	doMadd.drawCMFByVrange([0,70],nameSuffix="All")
	doMadd.drawCMFByVrange([15,35],nameSuffix="Maddalena")

	

if 0:
	doMadd.drawYSOSpectraAll(drawClass="Proto",smooth=True)
	doMadd.drawYSOSpectraAll(drawClass="Disk",smooth=True)
if 0:
	#examine the spectra of proto stars given by 
	doMadd.extractProtoStarSpectra("to2D_CO12RGB_R.fits")
 



#doMadd.moment(doMadd.COFile12,doMadd.vRange)
#doMadd.moment(doMadd.COFile13,doMadd.vRange)
#doMadd.moment(doMadd.COFile18,doMadd.vRange)
