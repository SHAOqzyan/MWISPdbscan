from astropy.table import   Table
import numpy as np
import os
#os.environ['QT_QPA_PLATFORM']='offscreen'
#gaiaDIS.py
from myPYTHON import *
import matplotlib as mpl
from pyds9 import DS9
from spectral_cube import SpectralCube

from progressbar import * 
import pymc3 as pm
from astropy.table import Column
from scipy.stats import norm,expon
from scipy import special
import multiprocessing
import scipy
import corner
from astropy.wcs import WCS
import pywcsgrid2
from matplotlib import rc
import os
from myTable import myTB
from astropy.io import fits
from scipy.stats import truncnorm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import interpolate
import pyregion

import matplotlib.gridspec as gridspec

from  distanceTB import disTB
#from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText
from matplotlib.offsetbox import AnchoredText
#from distanceMCMC5p import getDisAndErrorMCMCTrucatedGau5p
#from sklearn.isotonic import IsotonicRegression
#from distanceMCMCLinearGauss import getDisAndErrorMCMCTrucatedGau6p

#from distanceMCMCTwoGauss  import getDisAndErrorMCMCG2
#from termcolor import colored
import readline
#this  module is used to calculate the distance with Gaia ag extinction

#from gaiaAVTB import GAIAAV

#from gaiaAVFlagTB import GAIAAVFlag


from gaiaTB import GAIATB
#from gaiaTBCor import GAIATB

doAG= GAIATB()

doAV= GAIAAVFlag()  #avTB

class GAIADIS:
	#name="GSF_Gaia2_all"

	
	GAIA_parallax=doAG.GAIA_parallax
	GAIA_a_g_val=doAG.GAIA_a_g_val
	GAIA_parallax_err=doAG.GAIA_parallax_err
	
	GAIA_parallax_err=doAG.GAIA_parallax_err

	GAIA_a_g_percentile_lower = doAG.GAIA_a_g_percentile_lower   
	
	GAIA_a_g_percentile_upper = doAG.GAIA_a_g_percentile_upper
	
	agError =doAG.agError
	
	
	GAIA_distance=doAG.GAIA_distance
	GAIA_l=doAG.GAIA_l
	GAIA_b=doAG.GAIA_b
	
	GAIA_distanceError=doAG.GAIA_distanceError
	
	
	
	
	coint="coint"
	
	gaiaDo=GAIATB()
	
	nChains=10
	
	chainSample= 1000 #1000 #1000
	burn_in=50
	
	thinning=15 #get on sample every XX stars.
	
	maxDisCut=1000.
	


	parametersSave="para.txt"
	
	regionName='disRegions.reg'

	smScale="smScale"
	
	maskCol="mask"
	foreCOCol="foreCO"
	
	
	
	foregroundCOCut=3.  #K km/s

	#foregroundCOCut=10.  #K km/s

	#the lRang of G210 is 209.75-219.75
	
	extraDis=300. # select extra 300 pc stars to make the tail of off cloud star normal
	agThresh=0.05 #mag, lower
 

	tmpPath="/home/qzyan/WORK/myDownloads/GaiaTmpFiles/"

	offColor="dodgerblue"

	def __init__(self,sourceName="test",useAV=False, maskFITS=None,useBaseLine=False,disTBPath=None,redoPara=False,foregroundFITS=None): 


		self.useAV=useAV #otherwise use AG


		if self.useAV:
			
			GAIA_parallax=doAV.GAIA_parallax
			GAIA_a_g_val=doAV.GAIA_a_g_val
			GAIA_parallax_err=doAV.GAIA_parallax_err
			
			GAIA_parallax_err=doAV.GAIA_parallax_err
		

			
		  
			GAIA_distance=doAV.dist50
			
			GAIA_l=doAG.GAIA_l
			GAIA_b=doAG.GAIA_b
			
 
					
					
			

		
		self.redoPara=redoPara 


		self.maskFITS=maskFITS
		self.useBaseLine= useBaseLine

		self.disTBPath=disTBPath
		
		

		self.sourceRow=None
		
		if self.disTBPath !=None:
		
			self.disRecord=disTB(self.disTBPath )
		else:
			self.disRecord=disTB( )


		self.sourceName=sourceName
		self.sourcename=sourceName

		self.foregroundFITS=foregroundFITS


	def getAgBase(self,dataDis,baseLine):
		
		
		#it is possible that the nearst star is closer than off sttars
 
		baseAG= baseLine.predict(  dataDis )
		
		
		#print baseLine.predict([  baseLine.X_min_  ])[0], baseLine.predict([  baseLine.X_max_  ])[0]
		
		baseAG[ dataDis<  baseLine.X_min_  ]  = baseLine.predict([  baseLine.X_min_  ])[0]
		
		baseAG[ dataDis>  baseLine.X_max_  ]  =baseLine.predict([  baseLine.X_max_  ])[0]

		return baseAG
 
		#starInBin[GAIATB.GAIA_a_g_val]= starInBin[GAIATB.GAIA_a_g_val]- baseAG    

 

	def calDisWithG2(self,oncloudStars,offcloudStars,baseLine,smScale=50,correctBaseline=True  ):
		"""
		calculate the distance with oncloud stars, no figure will be provided, and only sample will be returned
		
		5 parameers are used 
		two gaussian model
		
		subtracted the base line
		"""
		#pass
 

		dataDis, dataAG,disError, AGError=self.getDisAndAGFromTB(oncloudStars)

		dataDisOff, dataAGOff,disErrorOff, AGErrorOff=self.getDisAndAGFromTB(offcloudStars)

		
 


		rawDataOn=[dataDis.copy(), dataAG.copy(),disError.copy(), AGError.copy()] 
		rawDataOff=[dataDisOff.copy(), dataAGOff.copy(),disErrorOff.copy(), AGErrorOff.copy()] 

		#rawDataBeforeBase=[ dataDis, dataAG,disError, AGError ]

		#smDisOff, smAGOff,a, a=self.getDisAndAGFromTB(offcloudStars)
		#dataDisRawOff, dataAGRawOff,disErrorRawOff, AGErrorRawOff=self.getDisAndAGFromTB(offcloudStars)
		#smDisOff,smAGOff,aa,aa=self.getSmoothArrayWithTB(offcloudStars,smScale)

		if baseLine ==None:
			baseLine = IsotonicRegression()

			offcloudStars=offcloudStars[ AGErrorOff> 0.05 ] 
			dataDisOff, dataAGOff,disErrorOff, AGErrorOff=self.getDisAndAGFromTB(offcloudStars)
			
			rawDataOff=[dataDisOff.copy(), dataAGOff.copy(),disErrorOff.copy(), AGErrorOff.copy()] 

			
			agWeight=1./AGErrorOff**2
		
			#agWeight=  1./offcloudStars[self.agError] 

 
			baseLine.fit(dataDisOff,dataAGOff,agWeight)

		
		#baseLine.fit_transform(smDisOff,smAGOff)
		
		#rawDataOff=[smDisOff,smAGOff,aa,aa]
		

		if correctBaseline: #subtract base line
			#axRaw.plot(drawNOCODis,y_ ,color='red',label=r"A$_{\rm G}$ baseline")	
 
 
			agBase=self.getAgBase(dataDis, baseLine)

			#dataAG=dataAG- agBase
			
			dataAG=dataAG - agBase


		#dataAG=dataAG-baseLine.predict(dataDis)
		
		#print baseLine.predict(dataDis)
		

		rawDataBase=[dataDis, dataAG,disError, AGError]

		nChain=10 # =================you could use 8 chain, if your computer allows it
		print "Calculating {} chains and each chain has {} samples, thinned by {}.".format(self.nChains,self.chainSample,self.thinning)  

		procs = []
	
		manager = multiprocessing.Manager()
		returnSampleDic = manager.dict()

		# instantiating process with arguments
		for name in range(nChain):
			# print(name)
			print "Starting process ",name
			
			proc = multiprocessing.Process(target=getDisAndErrorMCMCG2, args=(dataDis,dataAG,disError,AGError,name, returnSampleDic,self.chainSample,self.burn_in,self.thinning, self.maxDisCut))
			procs.append(proc)
			proc.start()



		for proc in procs:
			proc.join()
 
		sampleArray=[]
		

		for j in range(len(returnSampleDic[0])):
			
			sampleArrayj=[]
			for i in  range(len( returnSampleDic) ):
				
				sampleArrayj.append(returnSampleDic[i][j])
			
			sampleArray.append( np.concatenate(sampleArrayj) )

		return sampleArray,rawDataOn,rawDataOff,rawDataBase,baseLine 
		
		
		

	def calDisWith6Ps(self,oncloudStars ):
		"""
		calculate the distance with oncloud stars, no figure will be provided, and only sample will be returned
		
		5 parameers are used 
		mu1 is assgned to be zero
		"""
		#pass
 
		
		dataDis, dataAG,disError, AGError=self.getDisAndAGFromTB(oncloudStars)
 
		rawData=[dataDis, dataAG,disError, AGError]

		nChain=4 # =================you could use 8 chain, if your computer allows it
		print "Calculating {} chains and each chain has {} samples, thinned by {}.".format(self.nChains,self.chainSample,self.thinning)  

		procs = []
	
		manager = multiprocessing.Manager()
		returnSampleDic = manager.dict()

		# instantiating process with arguments
		for name in range(nChain):
			# print(name)
			print "Starting process ",name
			
			proc = multiprocessing.Process(target=getDisAndErrorMCMCTrucatedGau6p, args=(dataDis,dataAG,disError,AGError,name, returnSampleDic,self.chainSample,self.burn_in,self.thinning, self.maxDisCut))
			procs.append(proc)
			proc.start()



		for proc in procs:
			proc.join()
 
		sampleArray=[]
		

		for j in range(len(returnSampleDic[0])):
			
			sampleArrayj=[]
			for i in  range(len( returnSampleDic) ):
				
				sampleArrayj.append(returnSampleDic[i][j])
			
			sampleArray.append( np.concatenate(sampleArrayj) )

		return sampleArray,rawData
 


	

	def calDisWith5Ps(self,oncloudStars ):
		"""
		calculate the distance with oncloud stars, no figure will be provided, and only sample will be returned
		
		5 parameers are used 
		mu1 is assgned to be zero
		"""
		#pass
		print "Using truncated gauss"
		
		dataDis, dataAG,disError, AGError=self.getDisAndAGFromTB(oncloudStars)

		#print len(dataDis )


		rawData=[dataDis, dataAG,disError, AGError]

		nChain=self.nChains # =================you could use 8 chain, if your computer allows it
		print "Calculating {} chains and each chain has {} samples, thinned by {}.".format(self.nChains,self.chainSample,self.thinning)  

		procs = []
	
		manager = multiprocessing.Manager()
		returnSampleDic = manager.dict()

		# instantiating process with arguments
		for name in range(nChain):
			# print(name)
			print "Starting process ",name
			
			proc = multiprocessing.Process(target=getDisAndErrorMCMCTrucatedGau5p, args=(dataDis,dataAG,disError,AGError,name, returnSampleDic,self.chainSample,self.burn_in,self.thinning, self.maxDisCut))
			procs.append(proc)
			proc.start()



		for proc in procs:
			proc.join()
 
		sampleArray=[]
		

		for j in range(len(returnSampleDic[0])):
			
			sampleArrayj=[]
			for i in  range(len( returnSampleDic) ):
				
				sampleArrayj.append(returnSampleDic[i][j])
			
			sampleArray.append( np.concatenate(sampleArrayj) )


		print len(sampleArray[0]),"total samples????"

		return sampleArray,rawData



	def removeIndex(self,TB):
		
		TestTB=Table() 
		
		for eachCol in TB.colnames:
			
			
			#TB[eachCol].mask=False #do not mask
			
			aaa=Column( TB[eachCol]  ,name=eachCol )	
			
 
			TestTB.add_column(aaa )
		return TestTB
		
	def getSmoothDis(self,TB,smScale=10):
		"""
		the unnie of smoothScale is pc
		"""
		TB=TB.copy()
		#TB.sort(self.GAIA_dist50)
		TB=self.removeIndex(TB)

		newDisList=[]
		newExtictionList=[]
		
		error1=[]
		error2=[]

		distances, av,  disError,avError =self. getDisAndAGFromTB( TB)

		minDis=int( min( distances) )
		maxDis=int( max(distances) ) + 1
		
		newDisRange=np.arange(minDis,maxDis+smScale,smScale)

 
		TB.sort('distance') #need to clean index
		
 
		for countN,beginD in enumerate(newDisRange[0:-1]):
 
 
			endD=newDisRange[countN+1]
 
			selectCriteria=   (distances>=beginD) &  (distances<=endD)


			cut1=distances[ selectCriteria ]
 
			cut2=av[  selectCriteria ]

			if len(cut1)==0:
				continue


 

			cut1Err=disError[  selectCriteria]  #this is the error of distanes
			cut2Err=avError[  selectCriteria ]  #this is the error of distanes

			weight1=1./cut1Err**2/(np.sum(1./cut1Err**2))
			weight2=1./cut2Err**2/(np.sum(1./cut2Err**2))


			#print np.average(cut1,weights=weight1)
			
			
			newDisList.append( np.average(cut1,weights=weight1) ) #distance in pc
			
			newExtictionList.append( np.average(cut2,weights=weight2) ) # ag average
			
		#no errors are needed
			
  
		return np.array(newDisList),np.array(newExtictionList),None,None



	def getSmoothArrayWithTB(self,gaiaTB,smScale=10):
		"""
		the unnie of smoothScale is pc
		"""
		
		# do this consider if this is AV of AG..... 
		
 
 

		#better do this with distance directely
		
		if self.GAIA_distance in gaiaTB.colnames   :
			
			return self.getSmoothDis(gaiaTB,smScale=smScale)
 
		#print "aaaaaaaaaaaaaaaaaaaaaa"
		
		p_coParallax= gaiaTB[self.GAIA_parallax] 
		p_coExtention=gaiaTB[self.GAIA_a_g_val]
		p_coParallaxErr=gaiaTB[self.GAIA_parallax_err] 
		p_coExtentionErr= gaiaTB[self.agError]*p_coExtention
		
 
		#get somoothed on     cloud star
		return self.getSmoothArrayWithParallaxByMean(p_coParallax.data,p_coExtention.data,list1Error=p_coParallaxErr.data,list2Error=p_coExtentionErr.data,smoothNumber=smScale)
		#smsortCOPara,smsortCOExtention,smweigtCOPara,smweightCOExtention=self.getSmoothArrayWithParallax(p_coParallax.data,p_coExtention.data,list1Error=p_coParallaxErr.data,list2Error=p_coExtentionErr.data,smoothNumber=smoothNumber)
		#return 
		# off cloud stars		


	def getSmoothArrayWithParallaxByMean(self, list1,list2,list1Error=None,list2Error=None,smoothNumber=10):
		"""
		This function is used to sort list1 and list1Error, is parallax
		list1, parallax
		list2, extinction
		
		#perhaphs use 
		
		"""
		
		#radius=smoothNumber/2. #not a number ,but pc
		
		list1Copy=list1.copy()
 
		NewList2=[]
		
		minPara=list1.min()
		maxPara=list1.max()
		
		maxDis=1./minPara*1000
		minDis=1./maxPara*1000
		
		minDis=int(minDis)
		maxDis=int(maxDis) 
		 
		newDisRange=np.arange(minDis,maxDis+smoothNumber,smoothNumber)
		newParaRange=1000./newDisRange 
		
 
		newDisList=[]
		newExtictionList=[]
		
		error1=[]
		error2=[]
		
		
		
		for countN,eachD in enumerate(newParaRange[0:-1]):
 
			endPara=eachD
			
 
			beginPara=newParaRange[countN+1]
			 
			cut1=list1[ (list1Copy>=beginPara) &  (list1Copy<=endPara)]
 
			cut2=list2[ (list1Copy>=beginPara) &  (list1Copy<=endPara)]

			if len(cut1)==0:
				continue

			cut1Err=list1Error[ (list1Copy>=beginPara) &  (list1Copy<=endPara)]  #this is the error of distanes
			cut2Err=list2Error[ (list1Copy>=beginPara) &  (list1Copy<=endPara)]  #this is the error of distanes

			weight1=1./cut1Err**2/(np.sum(1./cut1Err**2))
			weight2=1./cut2Err**2/(np.sum(1./cut2Err**2))


			if np.nan in weight1:
				print "?????????????????????????????????"


			newDisList.append(1./np.average(cut1,weights=weight1)*1000 ) #distance in pc
			
			newExtictionList.append( np.average(cut2,weights=weight2) ) # ag average
			
		#no errors are needed
			
  
		return np.array(newDisList),np.array(newExtictionList),None,None






	def getDisArray(self,lRange,bRange,fitsName, signalLevel,noiseLevel, saveTBName, lowerDisCut=1,upperDisCut=2000,paraErrorCut=0.2, saveFigureName=None, newQuery=False, useGAIAFITS=False,disRow=None,maskFITS=None,foregroundCOFITS=None):
		"""
		
		
		calculate the distance with MCMC
		
		Parameters:
		
		lRange: the galactic longtitude range of the molecular cloud
		lRange: the galactic latitude range of the molecular cloud
		
		fitsName: the background fits file, which is used to classify Gaia stars
		signalLevel,noiseLevel: the signal and noise level cut
		
		saveTBName: the name of TB file, which is saved in the current path.
		
		lwoerDisCut,upperDisCut: the distance range we care about #unit pc
		
		paraErrorCut: the error limit of parallax, default value is 20%
		
		newQuery: upate the saveTBName, otherwise we would use the saved file to avoiding querying the database every time.
		
		
		
		#this function selection on and off cloud stars accroding to the signalLevel and noiseLevel, but there are many way to selct on and off cloud stars,
		
		#so, isolate the 

		"""
		
		
		#display the uncat
		
		#read fits file

		print "Calculating distances normally..."
		lowerPara=1./upperDisCut*1000
		upperPara= 1./lowerDisCut*1000
		
		if not useGAIAFITS:
		
			if newQuery: #save the table
				gaiaOnCloudStars=self.gaiaDo.getByLBRange(lRange,bRange,lowerPara=lowerPara,upperPara=1000., paraError=paraErrorCut)
				
				
				
				os.system("rm "+saveTBName)
				gaiaOnCloudStars.write(saveTBName)
				
	
			gaiaAllStars=Table.read(saveTBName)
		else:
			#use gaiaFITS
			gaiaAllStars=Table.read(saveTBName)
			
			gaiaAllStars=myTB.filterByRange(gaiaAllStars,self.GAIA_parallax,[lowerPara,None] )
			gaiaAllStars=myTB.filterByRange(gaiaAllStars,self.GAIA_l,lRange)
			gaiaAllStars=myTB.filterByRange(gaiaAllStars,self.GAIA_b,bRange)
			
			gaiaAllStars=myTB.filterByRange(gaiaAllStars,self.relative_error, [None,paraErrorCut] )
 
		#gaiaOnCloudStars=gaiaAllStars.loc[self.coint,signalLevel:]
		#gaiaOffCloudStars=gaiaAllStars.loc[self.coint,  :noiseLevel]

		gaiaOnCloudStars,gaiaOffCloudStars=self.getOnAndOffStars( gaiaAllStars,NL=noiseLevel,SL=signalLevel,maskFITS=self.maskFITS,intFITS=fitsName,foregroundFITS= foregroundCOFITS) 
 
		
		
		
		
		selectRows=  1- np.isnan(gaiaOnCloudStars[self.coint] )
		selectRows=selectRows==1
		
		gaiaOnCloudStars=gaiaOnCloudStars[selectRows]

		gaiaOnsource=gaiaOnCloudStars
 
		selectRows= 1- np.isnan(gaiaOffCloudStars[self.coint] )
		selectRows=selectRows==1

		gaiaOffCloudStars=gaiaOffCloudStars[selectRows]
		gaiaOffsource=gaiaOffCloudStars

		##############################################save new rows
		#disCatTB=disTB( )
		
		#newDisRow=disCatTB.getEmptyRow()
 
		
		
			
		print "On source star number: {}".format(len( gaiaOnCloudStars)  )


		#print gaiaOnCloudStars
		#print gaiaOffCloudStars


		TestTBOn=Table() 
		
		
		
		
		
		for eachCol in gaiaOnCloudStars.colnames:
			aaa=Column( gaiaOnCloudStars[eachCol],name=eachCol )	
			
 
			TestTBOn.add_column(aaa )
		#gaiaOffCloudStars.remove_column(self.coint)
		gaiaOnCloudStars=TestTBOn


		#gaiaOnCloudStars.remove_indices(self.coint)


		gaiaOnCloudStars.sort( self.GAIA_parallax )

		gaiaOnCloudStars.reverse() #
		

 
		#print gaiaOffCloudStars.colnames

		TestTB=Table() 
		
		for eachCol in gaiaOffCloudStars.colnames:
			aaa=Column( gaiaOffCloudStars[eachCol],name=eachCol )	
			
 
			TestTB.add_column(aaa )
		#gaiaOffCloudStars.remove_column(self.coint)
		gaiaOffCloudStars=TestTB
		#gaiaOffCloudStars.remove_indices(self.coint)
		
		#gaiaOffCloudStars.remove_indices(self.coint)

		gaiaOffCloudStars.sort( self.GAIA_parallax )
		gaiaOffCloudStars.reverse()

		print "Off source star number: {}".format(len( gaiaOffCloudStars)  )


		#p_coParallax= gaiaOnCloudStars[self.GAIA_parallax] 
		#p_coExtention=gaiaOnCloudStars[self.GAIA_a_g_val]
		#p_coParallaxErr=gaiaOnCloudStars[self.GAIA_parallax_err] 
		#p_coExtentionErr= gaiaOnCloudStars[self.agError]*p_coExtention
		
		#p_nocoParallax= gaiaOffCloudStars[self.GAIA_parallax]  #pc
		#p_nocoExtention=gaiaOffCloudStars[self.GAIA_a_g_val]
		#p_nocoParallaxErr=gaiaOffCloudStars[self.GAIA_parallax_err]  #pc
		#p_nocoExtentionErr= gaiaOffCloudStars[self.agError]*p_nocoExtention


		
		smoothNumber=10 #pc
		#get somoothed on     cloud star
		#smsortCOPara,smsortCOExtention,smweigtCOPara,smweightCOExtention=self.getSmoothArrayWithParallaxByMean(p_coParallax.data,p_coExtention.data,list1Error=p_coParallaxErr.data,list2Error=p_coExtentionErr.data,smoothNumber=smoothNumber)
		#smsortCOPara,smsortCOExtention,smweigtCOPara,smweightCOExtention=self.getSmoothArrayWithParallax(p_coParallax.data,p_coExtention.data,list1Error=p_coParallaxErr.data,list2Error=p_coExtentionErr.data,smoothNumber=smoothNumber)

		# off cloud stars
		#drawNOCOPara,drawNOCOAG,smweightNOCOPara,smweightNOCOExtention=self.getSmoothArrayWithParallaxByMean(p_nocoParallax.data,\
		#p_nocoExtention.data,list1Error=p_nocoParallaxErr.data,list2Error=p_nocoExtentionErr.data,smoothNumber=smoothNumber)
		

		#
		
		truegaiaOOnCloudStars=gaiaOnCloudStars.copy()
		truegaiaOOnCloudStars.add_index(self.GAIA_parallax)
 
		
		truegaiaOnCloudStars_1=truegaiaOOnCloudStars.loc[self.GAIA_parallax,  :upperPara]

		foreGroundStars=None

		if lowerDisCut>1:
			foreGroundStars=truegaiaOOnCloudStars.loc[self.GAIA_parallax,  upperPara:]


			dataDisFore, dataAGFore,disErrorFore, AGErrorFore=self.getDisAndAGFromTB(foreGroundStars)


		
		
		sampleArray,rawData=self.calDisWith5Ps(truegaiaOnCloudStars_1 )

		

		dataDis, dataAG,disError, AGError=rawData 
		

		self.saveParaToRow(  disRow , sampleArray)

		newDisRow= disRow

		distance=newDisRow[ disTB.distance ]  
		disStd=newDisRow[ disTB.disStd ]  
		lowerDis=newDisRow[ disTB.disHPDLow ] 
		upperDis= newDisRow[ disTB.disHPDUp ] 


		
		disStd= np.std(sampleArray[0],ddof=1) 
		disStd= round(disStd)
		disStd=int(disStd)
		distanceStd=disStd
		print "The distance of {} is: {}+/-{} pc. ".format(self.sourceName,distance,disStd),
		
		totalStd=( 0.05*distance )**2 +disStd**2
		
		totalStd=totalStd**0.5

		totalStd= round(totalStd)
		totalStd=int(totalStd)
		
		
		#print "Including the systematic error is:{}+/-{} pc".format(distance,totalStd)

 
		disStr="{}+/-{} pc".format(distance,totalStd)
		
		print "Including the 5% systematic error is:", colored(disStr,"red") #.format(distance,totalStd)
		
  
		#if not draw:
			
			#return distance,lowerDis,upperDis
			#return distance, distanceStd #lowerDis,upperDis
 
		agmu1= newDisRow[ disTB.AG1 ] 
		#newDisRow[ disTB.AG1Std ]= round( np.std( sampleArray[1],ddof=1),3) #distanceStd
		loweragmu1=newDisRow[ disTB.AG1HPDLow ] 
		upperagmu1=newDisRow[ disTB.AG1HPDUp ] 


		agmu1Sigma=newDisRow[ disTB.sigma1 ] 
		#newDisRow[ disTB.sigma1Std ]=  round(np.std( sampleArray[2],ddof=1),3 )  #distanceStd
		loweragmu1Sigma=newDisRow[ disTB.Sigma1HPDLow ] 
		upperagmu1Sigma=newDisRow[ disTB.Sigma1HPDUp ] 
		
		agmu2Sigma= newDisRow[ disTB.sigma2 ] 
		#newDisRow[ disTB.sigma2Std ]= round( np.std( sampleArray[4],ddof=1) , 3) #distanceStd
		loweragmu2Sigma= newDisRow[ disTB.Sigma2HPDLow ] 
		upperagmu2Sigma= newDisRow[ disTB.Sigma2HPDUp ] 


		agmu2= newDisRow[ disTB.AG2 ] 
		#newDisRow[ disTB.AG2Std ]= round( np.std( sampleArray[3],ddof=1), 3) #distanceStd
		loweragmu2=newDisRow[ disTB.AG2HPDLow ] 
		upperagmu2=newDisRow[ disTB.AG2HPDUp ] 

 

		#draw figuers

		fig, axs = plt.subplots(ncols=1, nrows=2,figsize=(16.5,8) )  
		rc('text', usetex=True )
		rc('text.latex',  preamble=r'\usepackage{upgreek}')

		for ax in  axs  :
			ax.remove()

		##########draw corner maps
		nn=1
		
		
		self.drawCorner(fig,sampleArray,newDisRow,nn)
  
 
	 

		#draw gia Rawstars
		#draw gaia raw
		if 1:
			ax=plt.subplot2grid((10*nn,10*nn), (0,4),colspan=6,rowspan=4)

			
			#ax=plt.subplot(gs[6:,6: ])
			#axRaw=plt.subplot(gs[0:5,6: ] )
			#axRaw=plt.subplot2grid((10*nn,10*nn), (0,4),colspan=6,rowspan=4)
			self.drawRawGaia( ax,gaiaOnsource,gaiaOffsource,gaiaForeCut=foreGroundStars,row=newDisRow,figNameMark=disRow[disTB.sourceName],markLoc=4  ) 
			
		
 
				
		#draw CO FITS

		#draw fits
		import pywcsgrid2
		from mpl_toolkits.axes_grid1.inset_locator import inset_axes
		import mpl_toolkits.axes_grid1.inset_locator as inset_locator
		
		

		backFISTHDR=fits.open(fitsName)

		backData=backFISTHDR[0].data
		backHead=backFISTHDR[0].header


		if len( backData.shape)==3:
			temp3D='temp3d.fits'
			myFITS.downTo2D(fitsName,outPUT=temp3D,overwrite=True) 

			backFISTHDR=fits.open(temp3D)
	
			backData=backFISTHDR[0].data
			backHead=backFISTHDR[0].header
				
		if len( backData.shape)==4:
			backData=backData[0]
			backData=backData[0]

		
 

		axBack=pywcsgrid2.subplot(224,  header= WCS(backHead))



		self.drawCloudFITS( fitsName, axBack,backData, backHead,lRange,bRange,noiseLevel,signalLevel , gaiaOnsource,gaiaOffsource, maskFITS=maskFITS) 


		if saveFigureName==None:
 
			AGAVText="Ag"
			
			if self.useAV:
				AGAVText="Av"
	
	
			saveFigureName="{}_{}".format( self.sourceRow[disTB.sourceName]  ,'Normal_{}.pdf'.format(AGAVText) )
	
			#plt.savefig(  self.sourceName+'_extinctionGaiaAg.pdf', bbox_inches="tight")
			plt.savefig( saveFigureName, bbox_inches="tight")
			saveFigureName="{}_{}".format( self.sourceRow[disTB.sourceName]  ,'Normal_{}.png'.format(AGAVText) )
	
			print "saving as: " ,saveFigureName
	
			#plt.savefig(  self.sourceName+'_extinctionGaiaAg.pdf', bbox_inches="tight")
			plt.savefig( saveFigureName, bbox_inches="tight",dpi=600)
			return distance,lowerDis,upperDis
		else:
			
			plt.savefig( saveFigureName, bbox_inches="tight")

			if ".pdf" in  saveFigureName:
				saveFigureName=saveFigureName.replace( ".pdf", ".png" )
				plt.savefig( saveFigureName, bbox_inches="tight",dpi=600)


			print "saving as: " ,saveFigureName
		
			return distance,lowerDis,upperDis

	
	def getDisAndAGFromTB(self,oncloudStars):
		#print "Calculating distances and their errors...."
		
		
		
		
		
		if "distance"  in oncloudStars.colnames and doAV.dist50 in oncloudStars.colnames:
			# has  
			#print "distance column found!!!!!!!!!!!!"
			
			#oncloudStars.sort("distance") 
			#offcloudStars.sort("distance") 
			dataDis= oncloudStars["distance"] 
			disError= oncloudStars["distance_err"]  
			
			dataAG=oncloudStars[doAV.av50]
			AGError= oncloudStars["agError"]* dataAG

			return np.array(dataDis),np.array( dataAG),np.array(disError), np.array(AGError)
 
		if "distance"  in oncloudStars.colnames  :
			# has  
			#print "distance column found!!!!!!!!!!!!"
			
			#oncloudStars.sort("distance") 
			#offcloudStars.sort("distance") 
			dataDis= oncloudStars["distance"] 
			disError= oncloudStars["distance_err"]  
			
			dataAG=oncloudStars["a_g_val"]
			AGError= oncloudStars["agError"]* dataAG

			return np.array(dataDis),np.array( dataAG),np.array(disError), np.array(AGError)
 



		if doAV.dist50 in oncloudStars.colnames:
			# has  
			#print "distance column found!!!!!!!!!!!!"
			
			#oncloudStars.sort("distance") 
			#offcloudStars.sort("distance") 
			dataDis= oncloudStars[ doAV.dist50 ] 
			disError= (oncloudStars[ doAV.dist84  ]  - oncloudStars[ doAV.dist16 ]  )/2.
			
			dataAG=oncloudStars[doAV.av50]			
			AGError=  (oncloudStars[ doAV.av84  ]  - oncloudStars[ doAV.av16 ]  )/2.

			


			return np.array(dataDis)*1000.,np.array( dataAG),np.array(disError)*1000., np.array(AGError)

		
		
		
		dataDis= oncloudStars["parallax"]*0
		disError=dataDis*0
		dataAG=oncloudStars["a_g_val"]
		AGError= oncloudStars["agError"]* dataAG

		#if 1 :#:#calculate parallax with mcmc
		for i in range(len(dataDis)):
			para=oncloudStars[i]["parallax"]
			paraError=oncloudStars[i][doAG.GAIA_parallax_err] #+0.1 #   #
			
			dA=1./np.random.normal(para,paraError,20000)*1000
			
			dataDis[i]=np.mean(dA)
			disError[i]=np.std(dA,ddof=1)
		#no bining
		
		return np.array(dataDis),np.array( dataAG),np.array(disError), np.array(AGError)
 


	def getDisAndHPD(self,disSample):
		
		

		disSample=np.array(disSample)
		
		meanDis=np.mean(disSample)
		
		#disHPD=pm.stats.hpd(disSample,alpha=0.1)
		disHPD=pm.stats.hpd(disSample,alpha=0.05)

		
		print disHPD,"===95%HPD===="
		#print disSample[-10:]
		lowerDis=meanDis-min(disHPD)
		upperDis= max(disHPD)-meanDis
		
		meanDis,lowerDis,upperDis=map(round,[meanDis,lowerDis,upperDis])
		meanDis,lowerDis,upperDis=map(int,[meanDis,lowerDis,upperDis])
		
		return meanDis, lowerDis,upperDis 


	def getmuAndHPD(self,muSample):
		muSample=np.array(muSample)
		
		meanMu=np.mean(muSample)
		
		muHPD=pm.stats.hpd(muSample,alpha=0.1)
		lowerMu=meanMu-min(muHPD)
		upperMu= max(muHPD)-meanMu
 

		return  round(meanMu,3),  round(lowerMu,3), round(upperMu,3) 


	def assignOnsourceGaia(self,GaiaTB,bgFITS,colName=None):


		GaiaTB= GaiaTB.copy()


		bgData,bgHead= fits.open(bgFITS)[0].data,fits.open(bgFITS)[0].header    #myFITS.readFITS(bgFITS)

		if colName ==None:
			colName=self.coint
			#addCol=Column(name=self.coint,data=np.zeros(len(GaiaTB)))


		addCol=Column(name=colName,data=np.zeros(len(GaiaTB)))

		GaiaTB.add_column(addCol)

		bgWCS=WCS(bgHead,naxis=2)




		Ls=GaiaTB["l"]
		Bs=GaiaTB["b"]

		Xs,Ys=bgWCS.all_world2pix(Ls,Bs,0)

		if Xs[0]<0:
			Xs,Ys=bgWCS.all_world2pix(Ls-360,Bs,0)


		Xs=map(round,Xs)
		Ys=map(round,Ys)
		Xs=map(int,Xs)
		Ys=map(int,Ys)

		#[y,x]



		sizeData=bgData.shape

		if len(sizeData)==3:
			bgData=bgData[0]
		if len(sizeData)==4:
			bgData=bgData[0]
			bgData=bgData[0]

		maxY,maxX=bgData.shape


		for i in range(len(Xs)):

			if Ys[i] >maxY-1 or  Xs[i]>maxX-1:
				GaiaTB[i][colName]= -1000 #self.getCointByLB(bgData, WCS(bgHead), eachStar["l"],eachStar["b"])

			else:

				GaiaTB[i][colName]= bgData[Ys[i]][Xs[i]]  #self.getCointByLB(bgData, WCS(bgHead), eachStar["l"],eachStar["b"])



		return  GaiaTB


	def calWithFITS(self,fitsName,Lrange,Brange,noise,signalLevel,noiseLevel):
		
		"""
		the noiseLevel is only used to compare, not for calculation
		"""
		#selectColName

		#first select Gaia by  LB range
		
		gaiaAllSourceStars=self.gaiaDo.getByLBRange(Lrange,Brange,upperPara= 1/2.2,paraError=0.2) # error<20%
		gaiaAllSourceStars=self.assignOnsourceGaia(gaiaAllSourceStars,fitsName)

		print len(gaiaAllSourceStars)
		gaiaAllSourceStars.add_index(self.coint)
		gaiaOnSourceStars=gaiaAllSourceStars.loc[self.coint,noise*signalLevel:]
			
		print "The total number of on-cloud gaia stars: ", len(gaiaOnSourceStars)
		gaiaOffSourceStars=gaiaAllSourceStars.loc[self.coint,:noise*noiseLevel]
	
		dataDis, dataAG,disError, AGError=self.getDisAndAGFromTB(gaiaOnSourceStars)
		dataDisOff, dataAGOff,disErrorOff, AGErrorOff=self.getDisAndAGFromTB(gaiaOffSourceStars)

 
		sampleArray=getDisAndErrorMCMCTrucatedGau6p(dataDis, dataAG,disError, AGError)
		#	return  [disArray,disSigmaArray,mu1Array,mu1SigmaArray,mu2Array,mu2SigmaArray] #disEqual,errorEqual,round()
		

		
		distance,lowerDis,upperDis=self.getDisAndHPD(sampleArray[0])
		
		disError,lowerDisError,upperDisError=self.getDisAndHPD(sampleArray[1])
		
		agmu1,loweragmu1,upperagmu1=self.getmuAndHPD(sampleArray[2])

		
		
		agmu1Sigma,loweragmu1Sigma,upperagmu1Sigma=self.getmuAndHPD(sampleArray[3])
		
		
		agmu2,loweragmu2,upperagmu2=self.getmuAndHPD(sampleArray[4])
		agmu2Sigma,loweragmu2Sigma,upperagmu2Sigma=self.getmuAndHPD(sampleArray[5])
		

		if 1: #draw
			
			#fig, ax = plt.subplots(figsize=(8, 6))
			#fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(8, 6))
			bgData,bgHead=myFITS.readFITS(fitsName)

			axFITS=pywcsgrid2.subplot(122,header=WCS(bgHead) )
			
			axStar=plt.subplot(121 )

 
			
			axStar.scatter(dataDis, dataAG,s=3)
			#ax.scatter(dataDisOff, dataAGOff)
			axFITS.imshow(bgData[0],origin='lower',interpolation=None)
			
			
			axFITS["gal"].plot([Lrange[0],Lrange[0]],Brange,'b--',lw=0.8)
			axFITS["gal"].plot([Lrange[1],Lrange[1]],Brange,'b--',lw=0.8)
			
			axFITS["gal"].plot(Lrange,[Brange[0],Brange[0]],'b--',lw=0.8)
			axFITS["gal"].plot(Lrange,[Brange[1],Brange[1]],'b--',lw=0.8)



			#crop fits
			
			outCropFITS= " FTSIDiscrop.fits"
			myFITS.cropFITS2D(fitsName,outFITS=outCropFITS,Lrange=Lrange,Brange=Brange,overWrite=True)
			
			cropData,cropHead=myFITS.readFITS(outCropFITS)
			cmap = plt.cm.winter
			cmap.set_bad('white',1.)
			axFITS[WCS(cropHead)].contour(cropData,  [noise*noiseLevel,noise*signalLevel],cmap=cmap,vmin=noise*noiseLevel*0.7,vmax=noise*signalLevel*1.1, linewidths=0.1)
 
			
			plt.savefig(  'distance.png', bbox_inches="tight",dpi=300)



	def box(self,centerL,centerB,lSize,bSize,dummy=0):

		"""
		return lRange and B Range
		"""

		lSize=lSize/3600.
		bSize=bSize/3600.
		
		
		return [centerL-lSize/2.,  centerL+lSize/2.   ], [centerB-bSize/2., centerB+bSize/2. ] 
		


	#def self.calDisWithRow( newRow)
	def calDisNormal(self,sName,getPara=True,maskFITS=None,foregroundCOFITS=None):
		
		
		calRow=self.disRecord.getRowByName(sName)

		if getPara or calRow==None:
			self.getParaWithDS9(sName)
		
		#self.calDisWithRow( newRow)
	
		calRow=self.disRecord.getRowByName(sName)

		
		
		self.sourceRow=calRow
		self.calDisWithRowNormal( calRow,maskFITS=maskFITS, foregroundCOFITS=foregroundCOFITS )
		
		
		
	
		
	def getParaWithDS9(self,sName):


	#def calDisWithDS9(self,sName):
		
		"""
		calculate FITS absed on the FITS File
		
		
		"""
		
		
		#first examine the fits file 
  
 
		if 		sName=="":
			print "You have to provide a name!Existing..."
			return 
		self.sourceRow=self.disRecord.getRowByName(sName)
		
		if self.sourceRow==None:
			print "Source does not exist, creating a new one...."
 
			self.sourceRow=self.disRecord.getEmptyRow()
			
			self.sourceRow[disTB.sourceName]=sName
		
		

		#print self.sourceRow
				
		#self.disRecord.addRowToTB( self.sourceRow )
		
		#get lowerDisCut
		
		disCutLow=self.sourceRow[disTB.cutDistanceLower] 
		
		if disCutLow<1:
			disCutLow=1
			
		
		
		NewdisCutLow = str.strip( raw_input("Lower discut ({} pc): ".format(disCutLow) ) )
		if NewdisCutLow!="":
			NewdisCutLow=float(NewdisCutLow)
			if NewdisCutLow>=1:
				disCutLow=NewdisCutLow
		self.sourceRow[disTB.cutDistanceLower] =disCutLow

 
 
 
		disCutUpper=self.sourceRow[disTB.cutDistanceUpper] 
		
		if disCutUpper<self.sourceRow[disTB.cutDistanceLower]:
			self.sourceRow[disTB.cutDistanceUpper]=self.sourceRow[disTB.cutDistanceUpper]+1000
			disCutUpper=self.sourceRow[disTB.cutDistanceUpper] 
 

		NewdisCutUpper= str.strip( raw_input("Upper discut ({} pc): ".format(disCutUpper) ) )
 
		
		if NewdisCutUpper!="":
			NewdisCutUpper=float(NewdisCutUpper)
			if NewdisCutUpper> self.sourceRow[disTB.cutDistanceLower] :
				disCutUpper=NewdisCutUpper
		self.sourceRow[disTB.cutDistanceUpper] =disCutUpper

		#print self.sourceRow[disTB.cutDistanceLower], self.sourceRow[disTB.cutDistanceUpper]
 
 
		#get fits Name
		
		fitsFile=self.sourceRow[disTB.fitsFile] 
 


		readline.parse_and_bind("tab: complete")



		newfitFile= str.strip( raw_input("FITS image ({}): ".format(fitsFile) ) )

		if ".fits" in newfitFile and  os.path.exists(newfitFile):
			fitsFile=   os.path.abspath(newfitFile)
 
		self.sourceRow[disTB.fitsFile] =fitsFile

		if not os.path.exists(self.sourceRow[disTB.fitsFile]):
			
			print "No usable FITS is found! Existing..."
	 
			return
			
		
		# signal cut
		
		
		noiseLevel=self.sourceRow[disTB.noiseLevel] 

		NEWnoiseLevel= str.strip( raw_input("Noise level ({} K): ".format(noiseLevel) ) )

		if NEWnoiseLevel!="":
			NEWnoiseLevel=float(NEWnoiseLevel)
			noiseLevel=NEWnoiseLevel
		self.sourceRow[disTB.noiseLevel] =noiseLevel



		signalLevel=self.sourceRow[disTB.signalLevel] 

		NEWsignalLevel= str.strip( raw_input("Sigma level ({} K): ".format(signalLevel) ) )

		if NEWsignalLevel!="":
			NEWsignalLevel=float(NEWsignalLevel)
			
			if NEWsignalLevel> self.sourceRow[disTB.noiseLevel]:
				signalLevel=NEWsignalLevel
 
			
		self.sourceRow[disTB.signalLevel] =signalLevel

		
 

		paraErrorCut=self.sourceRow[disTB.paraErrorCut] 

		NEWparaErrorCut= str.strip( raw_input("Maximum parallax error ({}): ".format(paraErrorCut) ) )

		if NEWparaErrorCut!="":
			NEWparaErrorCut=float(NEWparaErrorCut)
			
			if NEWparaErrorCut<1 and NEWparaErrorCut>0:
			
				paraErrorCut=NEWparaErrorCut
				
		self.sourceRow[disTB.paraErrorCut] =paraErrorCut


		

		#save Parameters
		



		print "Now openning ds9 to draw a region"
		
		
		
		 
		d=DS9("drawBoxDis")
		
		
		#fitsFile="local0.fits"
		
		d.set('file '+self.sourceRow[disTB.fitsFile])
		
		
		if self.sourceRow[disTB.boxSizeL]>0.:
			regionStr='galactic; box({},{},{}",{}",0)'.format(  self.sourceRow[disTB.boxCenterL],   self.sourceRow[disTB.boxCenterB],   self.sourceRow[disTB.boxSizeL]*3600.,   self.sourceRow[disTB.boxSizeB]*3600. )
 
			d.set('regions',regionStr)
			
		else:
			#put an   rectangular in the center of the fits
			
			#read fits
			data,head=myFITS.readFITS(self.sourceRow[disTB.fitsFile])
			
			dataShape= data.shape
			
			Ny,Nx=dataShape[-2],dataShape[-1]
 
			wcsBackground=WCS(head,naxis=(1,2))
			
			initialL,intialB=wcsBackground.wcs_pix2world(Nx/2,Ny/2,0)
 
			
			
			regionStr='galactic; box({},{},{}",{}",0)'.format( initialL,  intialB,   0.5*3600.,   0.5*3600. )
 
			d.set('regions',regionStr)
			
			#open the maskfits
		
		
		newfitFileMask= str.strip( raw_input("FITS MASK ({}): ".format("") ) )

		if ".fits" in newfitFileMask and  os.path.exists(newfitFileMask):
			fitsFileMask=   os.path.abspath(newfitFile)
 
			print fitsFileMask,"???????"
 
 
		#self.sourceRow[disTB.fitsFile] =fitsFile
			d.set("frame new") 
			d.set("frame lock  wcs") 
			d.set('file '+ newfitFileMask)

			
			d.set("contour")
 
			#	d.set("contour nlevels 2")
			d.set("contour color red")
			

 
			d.set("contour save maskTmp.ctr")
			
			#d.set('regions',regionStr)
			d.set(" frame move back") 
			d.set(" contour load maskTmp.ctr") 

			
			
			
		else :
			
			print "No usable mask FITS is found! "
	 
 



		raw_input("When you have drawn regions, press Enter to save regions and calculate distances...do not close ds9...")
		

		
		#d.set('regions','Galactic;')
		d.set('regions system wcs')
		d.set('regions sky galactic')

		d.set('regions save '+self.regionName )
		d.set('quit')


		regions = pyregion.open(self.regionName)

		print len(regions), " regions found!"
		
 
		for i,eachRegion in  enumerate(regions):
			
			newRow=self.sourceRow.table
			newRow=newRow.copy()
			newRow=newRow[0]
 
			newNameForEachRegion=newRow[disTB.sourceName] 
			
			

			if i>0:
				newRow[ disTB.sourceName ]=newNameForEachRegion+ "_" +str(i)
			
			
			boxPosition=eachRegion.coord_list
				
			
			lRange,bRange= self.boxRegion(boxPosition )
			
 
			newRow[ disTB.boxCenterL  ]= boxPosition[0]
			newRow[ disTB.boxCenterB  ]= boxPosition[1]

			newRow[ disTB.boxSizeL  ]= round( boxPosition[2],7)
			newRow[ disTB.boxSizeB  ]= round(boxPosition[3],7)


			
			
			#print lRange, newRow[ disTB.sourceName ]
			self.disRecord.addRowToTB(newRow) #save first in case discalculate are wrong..
			
			
			
			# stop here
			
			
			
			
			
			
			#self.getDisArray(lRange,bRange, newRow[disTB.fitsFile] ,newRow[disTB.signalLevel],newRow[disTB.noiseLevel],\
			  #newRow[disTB.sourceName] +"TB.fit",lowerDisCut=newRow[disTB.cutDistanceLower],upperDisCut=newRow[disTB.cutDistanceUpper],paraErrorCut=newRow[disTB.paraErrorCut],newQuery=True,disRow=newRow ) #taurus
 
			#self.disRecord.addRowToTB(newRow) #save second

 
		#cloudDis.getDisArray(lRange,bRange, "W3.fits", 4,1 ,"W3.fit",lowerDisCut=1000,upperDisCut=4000,paraErrorCut=0.2,newQuery=True ) #taur



	def findBestParameters(self,testSource ):
		"""
		test wich is best parameters
		
		sigma 1.1 K for G211, use fits  CO12RGB_R.fits 
		"""
		levelInterval=0.5
		
		#signaLevels=np.arange(3,8+levelInterval,levelInterval)
		#noiselevels= np.arange(1,3+levelInterval,levelInterval)
		#smScales=np.arange(5,55, 5)
		
		
		#array([ 3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13])
		
		RMS=1.0 #k km/s
		
		signaLevels=np.arange(3,10.2,0.2)  
		#noiselevels=np.arange(0,3,0.5)
		noiselevels=np.arange(0.,3.2,0.2)
		smScales=np.arange(5,55,5)
		

		initialValues={ disTB.signalLevel:6.0, disTB.noiseLevel:1., self.smScale:10.  }

		#get the distance of initial values
		#testRow=self.disRecord.getRowByName(sourceName)
		#eachScale, eachNoise, eachSignal=initialValues[self.smScale], initialValues[disTB.noiseLevel], initialValues[disTB.signalLevel]  



 
		self.testOneP( testSource,RMS,disTB.noiseLevel,noiselevels, initialValues )

		self.testOneP( testSource,RMS,disTB.signalLevel,signaLevels, initialValues )
		
		self.testOneP( testSource,RMS,self.smScale,smScales, initialValues )

	def testOneP(self,sourceName,RMS,paraName,testValues,initialValues):
		"""
		"""
		testRow=self.disRecord.getRowByName(sourceName)

		parameters=[]
		distanceAndStd=[]
		#name of saving files
		parameterFile= "{}_{}_testPara.txt".format(sourceName, paraName  )
		distanceFile= "{}_{}_testDis.txt".format( sourceName,paraName  )
 
		
	 
		initialValues=initialValues.copy()

		for eachV in testValues:

			initialValues[paraName]=eachV
			print "Tesing ..",paraName,"---: ",eachV
			
			
			eachScale, eachNoise, eachSignal=initialValues[self.smScale], initialValues[disTB.noiseLevel], initialValues[disTB.signalLevel]  
 
			parameters.append(    [ eachScale, eachNoise, eachSignal]   )
			
			testRow[disTB.noiseLevel]= eachNoise*RMS
			testRow[disTB.signalLevel]=eachSignal*RMS
			
			distance,std= self.calDisByRowBase(testRow, smScale= eachScale, draw=False) 

			distanceAndStd.append( [ distance,std  ] )
  
		np.savetxt( parameterFile, parameters )
		
		np.savetxt( distanceFile, distanceAndStd )
		
 

	def gaiaDisBaseline( self,sourceName, getNewPara=True,maskFITS=None,foregroundCOFITS=None,legendLoc=4):
		
		if getNewPara:
			
			self.getParaWithDS9(sourceName)		
		
		self.sourceRow=self.disRecord.getRowByName(sourceName)
		
		
		
		self.calDisByRowBase(self.sourceRow,maskFITS=maskFITS ,foregroundCOFITS=foregroundCOFITS,legendLoc=legendLoc)

	def getOnAndOffStars(self,TB,NL=1,SL=5,maskFITS=None,intFITS=None,foregroundFITS=None,useMask=True):
		
		"""
		"""
		
		#print NL,SL
		
		
		
 
		
		#print NL,SL,maskFITS,intFITS,foregroundFITS
		
		
		
		TB=TB.copy()
		#1  first, remove stars with foregroundFITS
		if foregroundFITS != None: # assuming this is foreground CO  fits  
			
			TB= self.assignOnsourceGaia(TB,foregroundCOFITS,self.foreCOCol )

			TB= TB[ TB[self.foreCOCol]<self.foregroundCOCut  ]	
 
		
		if maskFITS!=None and useMask:
			#add self.maskCol and further clean the data
			TB= self.assignOnsourceGaia(TB,maskFITS,self.maskCol )
			
		#
		
		#int fits should not be none 
		TB=self.assignOnsourceGaia(TB,intFITS)
		
		TB.add_index(self.coint)
		gaiaOnCloudStars=TB.loc[self.coint,SL:]
				
				
		gaiaOffCloudStars=TB.loc[self.coint,  :NL]
		
		#
		
		
		
 
		
		if maskFITS!=None and useMask:
			
			gaiaOffCloudStars=Table(gaiaOffCloudStars)
			
			gaiaOnCloudStars= 	gaiaOnCloudStars[ gaiaOnCloudStars[self.maskCol]>0.5   ]	
			
			gaiaOffCloudStars=	gaiaOffCloudStars[ gaiaOffCloudStars[self.maskCol]<0.5   ]	
			gaiaOffCloudStars=	gaiaOffCloudStars[ gaiaOffCloudStars[self.maskCol]>-1   ]	
				
				
		
		
			
			
		return gaiaOnCloudStars, gaiaOffCloudStars
		
			
 
	def saveParaToRow(self,newDisRow,sampleArray):
		
		"""
		"""
		
		#By default, there are five parameters
		
		
		distance,lowerDis,upperDis=self.getDisAndHPD(sampleArray[0])
		agmu1,loweragmu1,upperagmu1=self.getmuAndHPD(sampleArray[1])
		agmu1Sigma,loweragmu1Sigma,upperagmu1Sigma=self.getmuAndHPD(sampleArray[2])
		agmu2,loweragmu2,upperagmu2=self.getmuAndHPD(sampleArray[3])
		agmu2Sigma,loweragmu2Sigma,upperagmu2Sigma=self.getmuAndHPD(sampleArray[4])


		#save data 
		
		newDisRow[ disTB.distance ]=  distance
		newDisRow[ disTB.disStd ]=   	 round(np.std(sampleArray[0],ddof=1),0 ) 
		newDisRow[ disTB.disHPDLow ]=  lowerDis
		newDisRow[ disTB.disHPDUp ]=  upperDis


		newDisRow[ disTB.AG1 ]=  agmu1
		newDisRow[ disTB.AG1Std ]= round( np.std( sampleArray[1],ddof=1),3) #distanceStd
		newDisRow[ disTB.AG1HPDLow ]=  loweragmu1
		newDisRow[ disTB.AG1HPDUp ]=  upperagmu1


		newDisRow[ disTB.sigma1 ]=  agmu1Sigma
		newDisRow[ disTB.sigma1Std ]=  round(np.std( sampleArray[2],ddof=1),3 )  #distanceStd
		newDisRow[ disTB.Sigma1HPDLow ]=  loweragmu1Sigma
		newDisRow[ disTB.Sigma1HPDUp ]=  upperagmu1Sigma
		
		newDisRow[ disTB.sigma2 ]=  agmu2Sigma
		newDisRow[ disTB.sigma2Std ]= round( np.std( sampleArray[4],ddof=1) , 3) #distanceStd
		newDisRow[ disTB.Sigma2HPDLow ]=  loweragmu2Sigma
		newDisRow[ disTB.Sigma2HPDUp ]=  upperagmu2Sigma


		newDisRow[ disTB.AG2 ]=  agmu2
		newDisRow[ disTB.AG2Std ]= round( np.std( sampleArray[3],ddof=1), 3) #distanceStd
		newDisRow[ disTB.AG2HPDLow ]=  loweragmu2
		newDisRow[ disTB.AG2HPDUp ]=  upperagmu2



	def drawCorner(self,fig,sampleArray,newDisRow,nn):


		distance=newDisRow[ disTB.distance ]  
		lowerDis=newDisRow[ disTB.disHPDLow ]
		upperDis= newDisRow[ disTB.disHPDUp ] 
		distance,lowerDis,upperDis=map(int, [distance,lowerDis,upperDis] )

		agmu1= newDisRow[ disTB.AG1 ] 
		loweragmu1=newDisRow[ disTB.AG1HPDLow ] 
		upperagmu1=newDisRow[ disTB.AG1HPDUp ] 


		agmu1Sigma=newDisRow[ disTB.sigma1 ] 
		loweragmu1Sigma=newDisRow[ disTB.Sigma1HPDLow ] 
		upperagmu1Sigma=newDisRow[ disTB.Sigma1HPDUp ] 
		
		agmu2Sigma= newDisRow[ disTB.sigma2 ] 
		loweragmu2Sigma= newDisRow[ disTB.Sigma2HPDLow ] 
		upperagmu2Sigma= newDisRow[ disTB.Sigma2HPDUp ] 


		agmu2= newDisRow[ disTB.AG2 ] 
		loweragmu2=newDisRow[ disTB.AG2HPDLow ] 
		upperagmu2=newDisRow[ disTB.AG2HPDUp ] 

		

 

 
		meanValues = [distance,agmu1, agmu1Sigma,agmu2,agmu2Sigma] #np.mean(sampleArray, axis=0)
		lowerPHD90 = [lowerDis,loweragmu1, loweragmu1Sigma,loweragmu2,loweragmu2Sigma] #np.mean(sampleArray, axis=0)
		upperPHD90 = [upperDis, upperagmu1,upperagmu1Sigma,upperagmu2,upperagmu2Sigma] #np.mean(sampleArray, axis=0)

		#agmu1Sigma,loweragmu1Sigma,upperagmu1Sigma

		titleDis=r'$ D \rm  = {}_{{-{}}}^{{+{}}}$ pc'.format(distance,lowerDis,upperDis)
		
		titleMu1=r'$\mu_1 =  {:.3f}_{{-{:.3f}}}^{{+{:.3f}}}$ mag'.format(agmu1,loweragmu1,upperagmu1)

		titleMu1Sigma=r'$\sigma_1 =  {:.3f}_{{-{:.3f}}}^{{+{:.3f}}}$ mag'.format(agmu1Sigma,loweragmu1Sigma,upperagmu1Sigma)

		titleMu2=r'$\mu_2 =  {:.3f}_{{-{:.3f}}}^{{+{:.3f}}}$ mag'.format(agmu2,loweragmu2,upperagmu2)
		
		titleMu2Sigma=r'$ \sigma_2 =  {:.3f}_{{-{:.3f}}}^{{+{:.3f}}}$ mag'.format(agmu2Sigma,loweragmu2Sigma,upperagmu2Sigma)

		labels=[titleDis,titleMu1, titleMu1Sigma,titleMu2,titleMu2Sigma]
		for i in [0,1,2,3,4 ]:
			plt.subplot2grid((10*nn,10*nn), (i*nn*2, 0) ,rowspan=2*nn,colspan=nn)
			
			plt.subplot2grid((10*nn,10*nn), (i*nn*2,  nn),rowspan=2*nn,colspan=nn)
			plt.subplot2grid((10*nn,10*nn), (i*nn*2,  2*nn),rowspan=2*nn,colspan=nn)
			
			plt.subplot2grid((10*nn,10*nn), (i*nn*2,  3*nn),rowspan=2*nn,colspan=nn)
			plt.subplot2grid((10*nn,10*nn), (i*nn*2,  4*nn),rowspan=2*nn,colspan=nn)
		
		sampleArray=np.array(sampleArray)

		nx,ny=sampleArray.shape
		
		sampleArray=sampleArray.transpose()
		sampleArray=sampleArray.reshape(ny,nx)

		axLabels=["$D$ (pc)","$\mu_1$ (mag)", "$\sigma_1$ (mag)", "$\mu_2$ (mag)", "$\sigma_2$ (mag)", ]

		figureCorner = corner.corner(sampleArray, labels=axLabels, fig=fig,  show_titles=True )


		ndim=5
		axes = np.array(figureCorner.axes).reshape((ndim, ndim))
 
		

		# Loop over the diagonal
		for i in range(ndim):
			axCorner = axes[i, i]
			#pass
			#ax.axvline(meanValues[i], color="black" )
			#axCorner.title("aaa")

			axCorner.set_title(labels[i],fontsize=11)
			if i==0  :
				axCorner.set_title(labels[i],fontsize=12)

			axCorner.axvline(meanValues[i], color="black" ,linestyle='-',linewidth=1.5)
			axCorner.axvline(meanValues[i]-lowerPHD90[i], color="black" ,linestyle='--',linewidth=1.5)
			axCorner.axvline(meanValues[i]+upperPHD90[i], color="black" ,linestyle='--',linewidth=1.5)
			
		#Add labels

		#axDisLabelX=axes[0,0]
		#ax1.set_xlabel("minPts", fontsize=overAllFontSize)
		#axDisLabelX.set_xlabel("D (pc)", fontsize=11)

	def drawRawGaia(self,axRaw,gaiaOnsource,gaiaOffsource,gaiaForeCut=None,row=None,figNameMark=None,markLoc=4,baseline=None,drawOff=True ):
		#draw smooth reulst of on and off ccloud stars,
		
		#if distance results are provided, draw them



		dataDisOff,dataAGOff,aa,aa=self.getSmoothArrayWithTB(gaiaOffsource )
		
		smsortCODis,smsortCOExtention,aa,aa=self.getSmoothArrayWithTB( gaiaOnsource)

		#ax=plt.subplot2grid((10*nn,10*nn), (0,6),colspan=6,rowspan=4)
		# smooth the data 
		
		dataDis, dataAG,disError, AGError=self.getDisAndAGFromTB(gaiaOnsource)




		maxY= max(smsortCOExtention)+ (np.max(dataAG)-max(smsortCOExtention))*0.8



		minY= min(dataAG )-0.5


		minY=max( [minY, -1]  )
		
		
		
		#axRaw.scatter(dataDis, dataAG, edgecolors='darkgray', facecolors='darkgray',s=2,label="Raw on-cloud stars") # raw data

		if gaiaForeCut!=None:
			dataDisFore, dataAGFore,disError, AGError=self.getDisAndAGFromTB(gaiaForeCut)
			
			#should draw smooth dis fore
			
			axRaw.scatter(dataDisFore, dataAGFore, edgecolors='black', facecolors='black',s=2,label="Removed on-cloud foreground stars") # raw data

		sourceName=""
		
		if row!=None: #mark distances
			sourceName=row[disTB.sourceName]
			distance=row[ disTB.distance ]  
			lowerDis=row[ disTB.disHPDLow ]
			upperDis= row[ disTB.disHPDUp ] 
			
			sourceName = row[ disTB.sourceName ] 
			distance,lowerDis,upperDis=map(int, [distance,lowerDis,upperDis] )


			#change maxY
			backGroundStars = dataAG[dataDis>distance]

			if len(backGroundStars)==1:

				maxY=np.mean(backGroundStars)+1.9*0.5
			else:
				maxY = np.mean(backGroundStars) + 1.9 * np.std(backGroundStars, ddof=1)

			agmu1= row[ disTB.AG1 ]
			agmu2= row[ disTB.AG2 ] 
			axRaw.fill_betweenx(y=[minY, maxY*0.9], x1=distance-lowerDis, x2=upperDis+distance, color='moccasin',lw=0.1 );
			
			if baseline==None:
			
				axRaw.plot([min(dataDis),distance],[agmu1,agmu1], 'r--',lw=1.5  ,dashes=(4, 2))
				axRaw.plot([distance,max(dataDis)],[agmu2,agmu2], 'r--',lw=1.5  ,dashes=(4, 2))
				
			else:

				disB=np.arange(min(dataDisOff),max(dataDisOff),1)
				
				extLabel=r"$A_{G}$ baseline"
				if doAV.dist50 in gaiaOnsource.colnames: 
					extLabel=  r"$A_{V}$ baseline"
				
				axRaw.plot(disB,baseline.predict(disB), '--',lw=1.5,color='red',label=extLabel,   dashes=(3, 1))
			#ax.scatter(sortNOCOPara,sortNOCOExtention,lw=0.3,facecolors='b',s=5, edgecolors='b',label="Off-cloud stars" )

			axRaw.plot([distance,distance],[minY,maxY*0.9],lw=1,color="black")


			axRaw.text(distance-dataDis.max()/100.*6,maxY*0.92, r'${}_{{-{}}}^{{+{}}}$ pc'.format(distance,lowerDis,upperDis),fontsize=13)
			
		
		axRaw.scatter(dataDis, dataAG, edgecolors='darkgray', facecolors='darkgray',s=2,label="Raw on-cloud stars") # raw data

 
		if drawOff:
			axRaw.scatter(dataDisOff,dataAGOff, facecolors=self.offColor,s=5, edgecolors=self.offColor,label="Binned off-cloud stars" )
		
 
		axRaw.scatter(smsortCODis, smsortCOExtention, edgecolors="g" ,facecolors='g',s=5,label="Binned on-cloud stars")


		#sourceName=""


		if figNameMark==None:
		
			at = AnchoredText(sourceName, loc=markLoc, frameon=False,prop={"color":"black"} )
		else:
			at = AnchoredText(figNameMark, loc=markLoc, frameon=False,prop={"color":"black"} )

			
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axRaw.add_artist(at)
			

		
		
		#axRaw.plot([distance,distance],[-0.5,maxY*0.9],lw=1,color="black")
		#axRaw.text(distance-dataDis.max()/100.*5,maxY*0.91, r'${}_{{-{}}}^{{+{}}}$ pc'.format(distance,lowerDis,upperDis),fontsize=15)
		
		#axRaw.plot([min(dataDis),distance],[agmu1,agmu1], 'r--',lw=1.5  ,dashes=(4, 2))
		#axRaw.plot([distance,max(dataDis)],[agmu2,agmu2], 'r--',lw=1.5  ,dashes=(4, 2))

		#ax.set_xlabel(r"parallax (mas)",fontsize=17)
		axRaw.set_xlabel(r"distance (pc)",fontsize=13)
		#ax.set_ylabel(r"$\delta$ (Distance to nearest higher density pixels)",fontsize=17)
		showAV=False
		
		if  "dist50" in gaiaOnsource.colnames:
			showAV=True
		
		if showAV:
		
			axRaw.set_ylabel(r"$A_{V}$ (mag)",fontsize=13)
		else:
			
			axRaw.set_ylabel(r"$A_{G}$ (mag)",fontsize=13)

		
		if markLoc==4:
			axRaw.legend(loc=2,fontsize=13)
			#print "??????22222222222222?????????"

			
		else:
			axRaw.legend(loc=4,fontsize=13)
			#print "??????4444444???????????"
		axRaw.set_ylim(minY-0.5, max( max(smsortCOExtention)+1.0,minY+0.4)  )


	def drawBaselinePanel(self,ax,dataDisOnBase,dataAGOnBase, gaiaForeCut=None,row=None, baseline=None ):
		"""
		"""
		
		if self.useAV:
			maxY=  np.max(dataAGOnBase)*0.8
		
		else:
			maxY=  np.max(dataAGOnBase)*0.9

			
		
		minY=  np.min(dataAGOnBase)-0.1

		minY=max( [minY, -1]  )

		
		if row!=None: #mark distances
			sourceName=row[disTB.sourceName]
			distance=row[ disTB.distance ]  
			lowerDis=row[ disTB.disHPDLow ]
			upperDis= row[ disTB.disHPDUp ] 
			
			distance,lowerDis,upperDis=map(int, [distance,lowerDis,upperDis] )
			
			
			agmu1= row[ disTB.AG1 ] 
			agmu2= row[ disTB.AG2 ] 
		#ax.scatter(sortNOCOPara,sortNOCOExtention,lw=0.3,facecolors='b',s=5, edgecolors='b',label="Off-cloud stars" )
			ax.fill_betweenx(y=[minY, maxY*0.9], x1=distance-lowerDis, x2=upperDis+distance, color='moccasin',lw=0.1 );
			ax.plot([distance,distance],[minY,maxY*0.9],lw=1,color="black")
		
			ax.text(distance-dataDisOnBase.max()/100.*5,maxY*0.92, r'${}_{{-{}}}^{{+{}}}$ pc'.format(distance,lowerDis,upperDis),fontsize=13)
		
		ax.scatter(dataDisOnBase,dataAGOnBase, edgecolors='darkgray', facecolors='darkgray',s=3,label="Baseline-subtracted raw on-cloud stars ") # raw data

		


		if gaiaForeCut!=None:
		
			ax.scatter(dataDisFore, dataAGFore-baseline.predict(dataDisFore), edgecolors='black', facecolors='black',s=2,label="Removed on-cloud foreground stars") # raw data


		#ax.scatter(smsortCOPara, smsortCOExtention, edgecolors="g" ,facecolors='g',s=8,label="Binned on-cloud stars")



		ax.plot([min(dataDisOnBase),distance],[agmu1,agmu1], '--',lw=1.5  ,dashes=(4, 2),color="purple")
		ax.plot([distance,max(dataDisOnBase)],[agmu2,agmu2], '--',lw=1.5 ,dashes=(4, 2),color='purple')

		#ax.set_xlabel(r"parallax (mas)",fontsize=17)
		ax.set_xlabel(r"distance (pc)",fontsize=13)
		#ax.set_ylabel(r"$\delta$ (Distance to nearest higher density pixels)",fontsize=17)
		if self.useAV:
		
			ax.set_ylabel(r"$A_{V}$ (mag)",fontsize=13)
		
		else:
			ax.set_ylabel(r"$A_{G}$ (mag)",fontsize=13)

			
		ax.legend(loc=2,fontsize=13)

		ax.set_ylim(minY-0.5, max(dataAGOnBase)+1.0)


	def drawCloudFITS(self,fitsName, axBack,backData, backHead,lRange,bRange,noiseLevel,signalLevel,gaiaOnsource,gaiaOffsource, maskFITS=None,drawOffStar=True):
		"""
		"""
		if len( backData.shape)==3:
			temp3D='temp3d.fits'
			myFITS.downTo2D(fitsName,outPUT=temp3D,overwrite=True) 

			backFISTHDR=fits.open(temp3D)
	
			backData=backFISTHDR[0].data
			backHead=backFISTHDR[0].header
				
		if len( backData.shape)==4:
			backData=backData[0]
			backData=backData[0]
			
		
			
		#cropFITS

		outCropFITS= "CODrawCrop.fits"
		myFITS.cropFITS2D(fitsName,outFITS=outCropFITS,Lrange=lRange,Brange=bRange,overWrite=True)

 
		cropData,cropHead=myFITS.readFITS(outCropFITS)

		cropDataMask=cropData*0+1
	

		if backData.shape[0]==1:
			backData=backData[0]
 
		#axBack=	pywcsgrid2.plt.subplot2grid( (10*nn,10*nn), (0,3),colspan=3,rowspan=3, header= WCS(backHead) ) 
		#axBack=pywcsgrid2.subplot(224,  header= WCS(backHead))
		
		#noiseLevel=1.
		#signalLevel=4.
		

		np.mean(cropData )

		vmin= noiseLevel/2.

		aa=cropData[cropData>signalLevel]






		vmax=   np.mean(aa )+3*np.std(aa,ddof=1)  #signalLevel*3
		
		
		#print vmin, vmax,"?????"

		cmap = plt.cm.bone
		cmap.set_bad('black',1.)

		axBack.imshow(backData, origin="lower",cmap=cmap,norm=LogNorm(vmin=vmin, vmax=vmax),interpolation='none')
		#axBack.imshow(backData, origin="lower",cmap="bone", vmin=vmin, vmax=vmax ,interpolation='none')

		axBack.axis[:].major_ticks.set_color("w")
			
			
		cmap = plt.cm.hsv
		cmap.set_bad('white',1.)
		
		
		#axBack[WCS(cropHead)].contour(cropData,  [noiseLevel,signalLevel],cmap=cmap,vmin=noiseLevel*0.7,vmax=signalLevel*1.1, linewidths=0.1)

		#contourMaskFITS
		
		#axBack[WCS(cropHeadMask)].contour(cropDataMask,  [1],cmap=cmap,vmin=noiseLevel*0.7,vmax=signalLevel*1.1, linewidths=0.1)
		
		if maskFITS!=None:
			outMaskFITS= "CODrawCropMask.fits"
			myFITS.cropFITS2D(maskFITS,outFITS=outMaskFITS,Lrange=lRange,Brange=bRange,overWrite=True)
	
	 
			cropDataMask,cropHeadMask=myFITS.readFITS(outMaskFITS)
		
 
			axBack[WCS(cropHeadMask)].contour(cropDataMask,  [1],cmap=cmap,vmin= 0.99,vmax=signalLevel*1.1, linewidths=0.5)
		
		
		
		cmap2 = plt.cm.winter
		cmap2.set_bad('white',1.)


		print "..............."


		print np.max(  cropData)

		print np.min(  cropData),signalLevel




		axBack[WCS(cropHead)].contour(cropData*cropDataMask,  [  signalLevel],cmap=cmap2,vmin=noiseLevel*0.7,vmax=signalLevel*1.1, linewidths=0.3)

		#axBack[WCS(cropHead)].contour(cropData,  [noiseLevel ],cmap=cmap,vmin=noiseLevel*0.7,vmax=signalLevel*1.1, linewidths=0.1)
		#axBack[WCS(cropHeadMask)].contour
		drawOn=self.getRandomRows(gaiaOnsource,N=3000)
		
		
		axBack["gal"].scatter( drawOn[self.GAIA_l],drawOn[self.GAIA_b],s=0.5,color='green' ,marker='o' ,lw=0.3  )
		
		

		drawOff=self.getRandomRows(gaiaOffsource,N=3000)
		
		if drawOffStar:
			axBack["gal"].scatter( drawOff[self.GAIA_l],drawOff[self.GAIA_b],s=0.5,color=self.offColor ,marker='o'   ,lw=0.3)
		
		#cmap = plt.cm.winter
		#cmap.set_bad('white',1.)
		#axBack[WCS(backHead)].contour(cropData,  [IRASNoise,IRASSignal],cmap=cmap,vmin=IRASNoise*0.7,vmax=IRASSignal*1.1, linewidths=0.7)
		
		
		axBack["gal"].plot([lRange[0],lRange[0]],bRange,'b--',lw=0.8)
		axBack["gal"].plot([lRange[1],lRange[1]],bRange,'b--',lw=0.8)
		
		axBack["gal"].plot(lRange,[bRange[0],bRange[0]],'b--',lw=0.8)
		axBack["gal"].plot(lRange,[bRange[1],bRange[1]],'b--',lw=0.8)

		ExtendDeg  = 0  #1.5
		
		cutlSmall=min(lRange)-ExtendDeg
		cutlLarge=max(lRange)+ExtendDeg
		
		cutbSmall=min(bRange)-ExtendDeg
		cutbLarge=max(bRange)+ExtendDeg
		
		backWCS=WCS(backHead,naxis=2) #WCS(allHead)
		x0,y0=backWCS.all_world2pix(cutlLarge,cutbSmall,0)
		#
		x1,y1=backWCS.all_world2pix(cutlSmall,cutbLarge,0)	
		
		xRange=[x0,x1];
		yRange=[y0,y1]
		axBack.set_xlim(min(xRange),max(xRange)),  
		axBack.set_ylim(min(yRange),max(yRange)),


		axBack.set_ticklabel_type("absdeg", "absdeg")

		axBack.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
		axBack.set_ylabel(r"Galactic Latitude ($^{\circ}$)")

	def calDisAndrawWithOnAndOffCloudStars(self,sourceName,fitsName,maskFITS, gaiaOnsource,gaiaOffsource,baseLine,lRange,bRange,saveFigureName,noiseLevel=1.,signalLevel=5.,smScale=50,draw=True,lowerDisCut=1.,upperDisCut=2000,figNameMark=None,gaiaForeCut=None,inputRow=None,correctBaseline=True):
		
		"""
		
		smScale is 50 pc,
		
		"""
		
		if inputRow==None:
			disCatTB=disTB( )

			newDisRow=disCatTB.getEmptyRow()

		else:
			newDisRow=inputRow
			
		newDisRow[disTB.onStarN]=len(gaiaOnsource)
		newDisRow[disTB.offStarN]=len( gaiaOffsource)

		newDisRow[disTB.sourceName]= sourceName
		newDisRow[disTB.Note]= figNameMark
 

		newDisRow[disTB.boxCenterL]=  np.mean(lRange)
		newDisRow[disTB.boxCenterB]=  np.mean(bRange)

		newDisRow[disTB.boxSizeL]=  abs( lRange[0] -lRange[1]  )
		newDisRow[disTB.boxSizeB]=  abs( bRange[0] -bRange[1]  )


		#noiseLevel=1.,signalLevel
		
		#if baseLine==None:
			
			#self.getBaseLine(gaiaOffsource
			
		
		
		
		sampleArray, rawDataOn,rawDataOff,rawDataBase,baseline=self.calDisWithG2( gaiaOnsource,gaiaOffsource,baseLine,smScale,correctBaseline=correctBaseline )
		
		
		
		dataDis, dataAG,disError, AGError=rawDataOn 

		dataDisOff,dataAGOff,aa,aa=self.getSmoothArrayWithTB(gaiaOffsource )
		
		
		
		smsortCODis,smsortCOExtention,aa,aa=self.getSmoothArrayWithTB( gaiaOnsource)

		dataDisOnBase=dataDis

		if correctBaseline:
			dataAGOnBase=dataAG- self.getAgBase(dataDis, baseline)
		else:
			dataAGOnBase=dataAG #- self.getAgBase(dataDis, baseline)


		#dataDisOff, dataAGOff,disErrorOff, AGErrorOff=  rawDataOff  ##rawDataBase  rawDataOff
		#dataDisOnBase, dataAGOnBase,disErrorOnBase, AGErrorOnBase=  rawDataBase  ##rawDataBase  rawDataOff

		#print baseline.predict(dataDisOff)


		dataAGBase=dataAG-baseline.predict(dataDis)
		
		self.saveParaToRow( newDisRow , sampleArray)

 

		distance=newDisRow[ disTB.distance ]  
		disStd=newDisRow[ disTB.disStd ]  
		lowerDis=newDisRow[ disTB.disHPDLow ] 
		upperDis= newDisRow[ disTB.disHPDUp ] 


		
		disStd= np.std(sampleArray[0],ddof=1) 
		disStd= round(disStd)
		disStd=int(disStd)
		distanceStd=disStd
		print "The distance of {} is: {}+/-{} pc. ".format(sourceName,distance,disStd),
		
		totalStd=( 0.05*distance )**2 +disStd**2
		
		totalStd=totalStd**0.5

		totalStd= round(totalStd)
		totalStd=int(totalStd)
		
		
		#print "Including the systematic error is:{}+/-{} pc".format(distance,totalStd)


		disStr="{}+/-{} pc".format(distance,totalStd)
		
		print "Including the 5% systematic error is:", colored(disStr,"red") #.format(distance,totalStd)
		
  
		if not draw:
			
			return distance,lowerDis,upperDis
			#return distance, distanceStd #lowerDis,upperDis
 
		agmu1= newDisRow[ disTB.AG1 ] 
		#newDisRow[ disTB.AG1Std ]= round( np.std( sampleArray[1],ddof=1),3) #distanceStd
		loweragmu1=newDisRow[ disTB.AG1HPDLow ] 
		upperagmu1=newDisRow[ disTB.AG1HPDUp ] 


		agmu1Sigma=newDisRow[ disTB.sigma1 ] 
		#newDisRow[ disTB.sigma1Std ]=  round(np.std( sampleArray[2],ddof=1),3 )  #distanceStd
		loweragmu1Sigma=newDisRow[ disTB.Sigma1HPDLow ] 
		upperagmu1Sigma=newDisRow[ disTB.Sigma1HPDUp ] 
		
		agmu2Sigma= newDisRow[ disTB.sigma2 ] 
		#newDisRow[ disTB.sigma2Std ]= round( np.std( sampleArray[4],ddof=1) , 3) #distanceStd
		loweragmu2Sigma= newDisRow[ disTB.Sigma2HPDLow ] 
		upperagmu2Sigma= newDisRow[ disTB.Sigma2HPDUp ] 


		agmu2= newDisRow[ disTB.AG2 ] 
		#newDisRow[ disTB.AG2Std ]= round( np.std( sampleArray[3],ddof=1), 3) #distanceStd
		loweragmu2=newDisRow[ disTB.AG2HPDLow ] 
		upperagmu2=newDisRow[ disTB.AG2HPDUp ] 




		fig, axs = plt.subplots(ncols=1, nrows=2,figsize=(16.5,8) )  
		rc('text', usetex=True )
		#rc('text.latex',  preamble=r'\usepackage{upgreek}')
		rc('font', **{'family': 'sans-serif',  'size'   : 11,  'serif': ['Helvetica'] })
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]




		for ax in  axs  :
			ax.remove()
		#gs = gridspec.GridSpec(3, 3)
		gs = gridspec.GridSpec(10, 10)
		nn=2
		

		
		self.drawCorner(fig,sampleArray,newDisRow,nn)
  
		#plt.savefig( "test.png", bbox_inches="tight")

		#draw gaia raw
		if 1:
			#ax=plt.subplot(gs[6:,6: ])
			axRaw=plt.subplot(gs[0:5,6: ] )
			#axRaw=plt.subplot2grid((10*nn,10*nn), (0,4),colspan=6,rowspan=4)
			self.drawRawGaia( axRaw,gaiaOnsource,gaiaOffsource,gaiaForeCut=gaiaForeCut,row=newDisRow,figNameMark=figNameMark,markLoc=4,baseline=baseLine ) 
			

  
		#draw gia baselines 
		if 1:
			ax=plt.subplot(gs[6:,6: ] ,sharex=axRaw)
			#ax=plt.subplot(gs[0:5,6: ])

			#ax=plt.subplot2grid((10*nn,10*nn), (0,6),colspan=6,rowspan=4)
	 
			#maxY= max(smsortCOExtention)+ (np.max(dataAGBase)-max(smsortCOExtention)*0.8)
			
			
			self.drawBaselinePanel( ax,dataDisOnBase,dataAGOnBase,  row=newDisRow  )
				
  

		#draw CO FITS

		#draw fits
		import pywcsgrid2
		from mpl_toolkits.axes_grid1.inset_locator import inset_axes
		import mpl_toolkits.axes_grid1.inset_locator as inset_locator
		


		backFISTHDR=fits.open(fitsName)

		backData=backFISTHDR[0].data
		backHead=backFISTHDR[0].header

		axBack=pywcsgrid2.subplot(gs[0:4,3:5], header=WCS(backHead) )

		

		self.drawCloudFITS( fitsName, axBack,backData, backHead,lRange,bRange,noiseLevel,signalLevel , gaiaOnsource,gaiaOffsource, maskFITS=maskFITS) 

		#plt.show()

		#saveFigureName="{}_{}_{}_{}{}".format( self.sourceRow[disTB.sourceName] ,lowerDisCut,upperDisCut,paraErrorCut,'_extinctionGaiaAgBaseline.pdf')
		#fig.tight_layout()
		plt.savefig( saveFigureName, bbox_inches="tight")
		#plt.savefig( saveFigureName )


		#plt.savefig(  self.sourceName+'_extinctionGaiaAg.pdf', bbox_inches="tight")
		#plt.show()
		if "pdf" in saveFigureName:
			saveFigureName=saveFigureName.replace(".pdf",".png")
			plt.savefig( saveFigureName, bbox_inches="tight",dpi=600)
 
		#saveFigureName=sourceName #+  #"{}_{}_{}_{}{}".format(self.sourceRow[disTB.sourceName],lowerDisCut,upperDisCut,paraErrorCut,'_extinctionGaiaAgBaseline.png')

		print "saving as: " ,saveFigureName

		#plt.savefig(  self.sourceName+'_extinctionGaiaAg.pdf', bbox_inches="tight")
		#plt.savefig( saveFigureName, bbox_inches="tight")
		
		return newDisRow
		#return distance, distanceStd #lowerDis,upperDis



	def getRandomRows(self,TB,N=3000):
		"""
		"""

		if len(TB)<=N:
			return TB.copy()
		N=int(N)
 
		selectRowOff=np.zeros(len( TB ))
		selectRowOff[-N:]=1
		np.random.shuffle(selectRowOff)
		return TB[ selectRowOff==1 ]
		
	

	def addDisToTB(self,gaiaTB ):
		
		"""
		
		calculate distance, and add them to gaiaTB,
		
		usually only do once
		"""
		#copy file
		newTB=gaiaTB.copy()

 
		if doAV.dist50 in newTB.colnames and  doAG.GAIA_distance in newTB.colnames :

			return newTB
		
		if doAV.dist50 in newTB.colnames:
			
			
			
			disCol= ( gaiaTB[doAV.dist50])*1000 
			disCol.name= doAG.GAIA_distance
			
	 
			disErrCol= ( gaiaTB[ doAV.dist84] -  gaiaTB[ doAV.dist16] )/2.*1000.
			disErrCol.name= doAG.GAIA_distanceError
	

			
			agErrorCol= ( gaiaTB[ doAV.av84] -  gaiaTB[ doAV.av16] )/2./gaiaTB[ doAV.av50]
			agErrorCol.name= doAG.agError
			
			
 
			newTB.add_columns( [disCol,  disErrCol ,agErrorCol ])

			return newTB
		
		
		
 
		
		#print "?????????????????????????????????????"
		disCol= gaiaTB["parallax"]*0
		disCol.name= self.GAIA_distance
		
 
		disErrCol= gaiaTB["parallax"]*0
		disErrCol.name= self.GAIA_distanceError


		parallaxErrCol="parallax_err"

		if "parallax_error" in newTB.colnames: # for new gaia TB
			parallaxErrCol="parallax_error"


		
		newTB.add_columns( [disCol,  disErrCol  ])
		
		
		for eachRow in newTB: 
			para=eachRow["parallax"]
			paraError= eachRow[parallaxErrCol] #+0.1 #   #
			
			dA=1./np.random.normal(para,paraError,20000)*1000
			
			eachRow[self.GAIA_distance]=  round(np.mean(dA), 3) 
			
			eachRow[self.GAIA_distanceError]=   round(np.std(dA,ddof=1), 3) 
		

		return newTB

	def cleanOffTBByAgThresh(self,gaiaOffsource,agErrorThresh=0.05):
		
 
		
		AGErrorOff=  gaiaOffsource[self.agError]*gaiaOffsource[self.GAIA_a_g_val]  

			
 
		
		gaiaOffsource=gaiaOffsource[ AGErrorOff> agErrorThresh] 
		
		return gaiaOffsource



	def getBaseLine(self,offcloudStars):
		
		

		
		baseLine = IsotonicRegression()

		#AGErrorOff= offcloudStars[self.agError]*offcloudStars[self.GAIA_a_g_val]  


		#if self.useAV:
			#AGErrorOff= ( offcloudStars[doAV.av84] - offcloudStars[ doAV.av16]    )/2.  #   gaiaOffsource[doAV.agError]*gaiaOffsource[self.GAIA_a_g_val]  
		
		#else:
		AGErrorOff= offcloudStars[self.agError]*offcloudStars[self.GAIA_a_g_val]  

		
		offcloudStars=offcloudStars[ AGErrorOff>self.agThresh ] 
		
		dataDisOff, dataAGOff,disErrorOff, AGErrorOff=self.getDisAndAGFromTB(offcloudStars)
		
		#rawDataOff=[dataDisOff.copy(), dataAGOff.copy(),disErrorOff.copy(), AGErrorOff.copy()] 
		
		
		agWeight=1./AGErrorOff**2

		#agWeight=None #1./AGErrorOff**2


		#agWeight=  1./offcloudStars[self.agError] 

		baseLine.fit(dataDisOff,dataAGOff,agWeight)
		
		return baseLine

	
	def cutByUpperDis(self,TB,maxDis=3000 ):
		
		TB=TB.copy()
		
		newTB=TB[ TB[disTB.distance] < maxDis   ]
		return newTB

	def calDisByRowBase(self,testRow,smScale=10,draw=True,maskFITS=None,foregroundCOFITS=None,legendLoc=4):
		"""
		"""
		#reast source name
		print "Calculating distances with baseline subtracting"
		
		lRange,bRange=self.getLBRangeByRow( testRow) #used to crop fits
		
		paraErrorCut=testRow[disTB.paraErrorCut]
  
		#bgFITS="S287M0.fits"  #testRow[disTB.fitsFile] 
		fitsName= testRow[disTB.fitsFile] 

		sourceName=	 testRow[disTB.sourceName] 


		gaiaAll=self.getAllGaiaByRow(testRow,extraCut=self.extraDis) # cut extradis to fit the off cloud stars




		if self.useAV:
			gaiaAll=self.addDisToTB(gaiaAll)
 
 
 
 
		lowerDisCut=testRow[disTB.cutDistanceLower]
		upperDisCut=testRow[disTB.cutDistanceUpper]

		lowerPara=1./upperDisCut*1000
		upperPara= 1./lowerDisCut*1000

		
		
		if foregroundCOFITS != None: # assuming this is foreground CO  fits  
			
			gaiaAll= self.assignOnsourceGaia(gaiaAll,foregroundCOFITS,self.foreCOCol )

			#gaiaAll= gaiaAll[ gaiaAll[self.foreCOCol]<self.foregroundCOCut  ]
			
			
		gaiaAll=self.assignOnsourceGaia( gaiaAll,fitsName)

		
		if maskFITS!=None:
			#add self.maskCol and further clean the data
			gaiaAll= self.assignOnsourceGaia(gaiaAll,maskFITS,self.maskCol )
			
		


		#gaiaAll.write("testGaiaS287.fit")

		noiseLevel=testRow[disTB.noiseLevel]
		signalLevel=testRow[disTB.signalLevel]

		#upperPara
		#cut the distance according the distance range
		gaiaAllDisCut= myTB.filterByRange(gaiaAll,self.GAIA_parallax,[lowerPara, None ] ) #all stars

		
 
		#fore ground star cut 
		#gaiaForeCut= myTB.filterByRange(gaiaAll,self.GAIA_parallax,[upperPara, None] )
		gaiaOffsource=myTB.filterByRange(gaiaAllDisCut,self.coint,[ None,noiseLevel ] )

		if foregroundCOFITS != None: # assuming this is foreground CO  fits
			gaiaOffsource= gaiaOffsource[ gaiaOffsource[self.foreCOCol]<self.foregroundCOCut  ] #only do this on off cloud star
		
		
		
		onCloudStars= myTB.filterByRange(gaiaAllDisCut,self.coint,[ signalLevel,None  ] )
 
		gaiaOnsource= myTB.filterByRange(onCloudStars,self.GAIA_parallax,[None, upperPara] ) #cutSources

		gaiaOnsource=self.addDisToTB(gaiaOnsource)
		gaiaOffsource=self.addDisToTB(gaiaOffsource)
		
		gaiaOffsource=self.cleanOffTBByAgThresh(gaiaOffsource)
		

		#fit baseline 
		

		baseline=self.getBaseLine(gaiaOffsource)
		gaiaOnsource=self.cutByUpperDis(gaiaOnsource, upperDisCut )
		gaiaOffsource=self.cutByUpperDis(gaiaOffsource,upperDisCut) 


		if maskFITS!=None:
			gaiaOnsource= 	gaiaOnsource[ gaiaOnsource[self.maskCol]>0.5   ]	
			
			gaiaOffsource=	gaiaOffsource[ gaiaOffsource[self.maskCol]<0.5   ]	
			gaiaOffsource=	gaiaOffsource[ gaiaOffsource[self.maskCol]>-1   ]	




		gaiaForeCut=None
		if lowerDisCut>1: # usually not with baseline
		
			gaiaForeCut=myTB.filterByRange(onCloudStars,self.GAIA_parallax,[upperPara, None] )



		#adding distancecolumn ////////////////////////
		
		gaiaOnsource=gaiaOnsource[gaiaOnsource[disTB.distance]>0]
		gaiaOffsource=gaiaOffsource[gaiaOffsource[disTB.distance]>0] # one distances is nan, noidea what 

		
		print "Total on source stars",len(gaiaOnsource)


		#seems like useless anymore
		drawNOCODis,drawNOCOAG,aa,aa=self.getSmoothArrayWithTB(gaiaOffsource,smScale) # the para here is actually distance
		smsortCODis,smsortCOExtention,aa,aa=self.getSmoothArrayWithTB(gaiaOnsource,smScale) #
		
 
		
		if lowerDisCut>1: # usually not with baseline
			#gaiaForeCut= myTB.filterByRange(gaiaAll,self.GAIA_parallax,[upperPara, None] )
			#gaiaForeCut =myTB.filterByRange(gaiaForeCut,self.coint,[ signalLevel,None  ] )

			#dataDisForeRawSM,dataAGForeRawSM,aa,aa=self.getSmoothArrayWithTB(gaiaForeCut, smScale)

			#drawNOCODis,drawNOCOAG,aa,aa=self.getSmoothArrayWithTB(gaiaOffsource,smScale) # the para here is actually distance
			dataDisFore, dataAGFore,disError, AGError=self.getDisAndAGFromTB(gaiaForeCut)
		#ax.scatter(drawNOCOPara,drawNOCOAG,lw=0.3,facecolors='b',s=8, edgecolors='b',label="Binned off-cloud stars" )

 
 

		#save on and off cloud stars
		#

		if self.useAV:
			saveTBOnName=self.tmpPath+self.sourceName+"AV.fit"

		else:
			saveTBOnName=self.tmpPath+self.sourceName+"AG.fit"

		gaiaOnsource.write(saveTBOnName,overwrite=True)

		sampleArray, rawDataOn,rawDataOff,rawDataBase,baseline=self.calDisWithG2( gaiaOnsource,gaiaOffsource,baseline,smScale )
		dataDis, dataAG,disError, AGError=rawDataOn 


		
		
		


		dataDisOff, dataAGOff,disErrorOff, AGErrorOff=  rawDataOff  ##rawDataBase  rawDataOff
		
		dataDisOnBase, dataAGOnBase,disErrorOnBase, AGErrorOnBase=  rawDataBase  ##rawDataBase  rawDataOff
		
		
		dataDisOff,dataAGOff,aa,aa=self.getSmoothArrayWithTB(gaiaOffsource )

		#print baseline.predict(dataDisOff)


		dataAGBase=dataAG-baseline.predict(dataDis)
		
		distance,lowerDis,upperDis=self.getDisAndHPD(sampleArray[0])
 
		
		disStd= np.std(sampleArray[0],ddof=1) 
		disStd= round(disStd)
		disStd=int(disStd)

		distanceStd=disStd
		print "The distance of {} is: {}+/-{} pc. ".format(sourceName,distance,disStd),
		
		totalStd=( 0.05*distance )**2 +disStd**2
		
		totalStd=totalStd**0.5

		totalStd= round(totalStd)
		totalStd=int(totalStd)
		
		
		#print "Including the systematic error is:{}+/-{} pc".format(distance,totalStd)

 
		disStr="{}+/-{} pc".format(distance,totalStd)
		
		print "Including the 5% systematic error is:", colored(disStr,"red") #.format(distance,totalStd)
		
		
		
		

		if not draw:
			return distance, distanceStd #lowerDis,upperDis
			


		
		agmu1,loweragmu1,upperagmu1=self.getmuAndHPD(sampleArray[1])
		agmu1Sigma,loweragmu1Sigma,upperagmu1Sigma=self.getmuAndHPD(sampleArray[2])
		agmu2,loweragmu2,upperagmu2=self.getmuAndHPD(sampleArray[3])
		agmu2Sigma,loweragmu2Sigma,upperagmu2Sigma=self.getmuAndHPD(sampleArray[4])



		#save distance results
		#save data 
		if 1:
			testRow[ disTB.distance ]=  distance
			testRow[ disTB.disStd ]=  distanceStd
			testRow[ disTB.disHPDLow ]=  lowerDis
			testRow[ disTB.disHPDUp ]=  upperDis
	
	
			testRow[ disTB.sigma1 ]=  agmu1Sigma
			testRow[ disTB.sigma1Std ]=  round(np.std( sampleArray[2],ddof=1),3 )  #distanceStd
			testRow[ disTB.Sigma1HPDLow ]=  loweragmu1Sigma
			testRow[ disTB.Sigma1HPDUp ]=  upperagmu1Sigma
			
			testRow[ disTB.sigma2 ]=  agmu2Sigma
			testRow[ disTB.sigma2Std ]= round( np.std( sampleArray[4],ddof=1) , 3) #distanceStd
			testRow[ disTB.Sigma2HPDLow ]=  loweragmu2Sigma
			testRow[ disTB.Sigma2HPDUp ]=  upperagmu2Sigma
	
	
			testRow[ disTB.AG1 ]=  agmu1
			testRow[ disTB.AG1Std ]= round( np.std( sampleArray[1],ddof=1),3) #distanceStd
			testRow[ disTB.AG1HPDLow ]=  loweragmu1
			testRow[ disTB.AG1HPDUp ]=  upperagmu1
	
	
	
	
			testRow[ disTB.AG2 ]=  agmu2
			testRow[ disTB.AG2Std ]= round( np.std( sampleArray[3],ddof=1), 3) #distanceStd
			testRow[ disTB.AG2HPDLow ]=  loweragmu2
			testRow[ disTB.AG2HPDUp ]=  upperagmu2



		fig, axs = plt.subplots(ncols=1, nrows=2,figsize=(16.5,8) )  
		rc('text', usetex=True )
		rc('text.latex',  preamble=r'\usepackage{upgreek}')

		for ax in  axs  :
			ax.remove()
		#gs = gridspec.GridSpec(3, 3)
		gs = gridspec.GridSpec(10, 10)
 
		nn=2
		self.drawCorner(fig,sampleArray,testRow,nn)

 
 
 
			
		#draw gia baselines #baseline plot
		if 1:
			ax=plt.subplot(gs[6:,6: ])
			#ax=plt.subplot(gs[0:5,6: ])

			#ax=plt.subplot2grid((10*nn,10*nn), (0,6),colspan=6,rowspan=4)
			self.drawBaselinePanel( ax,dataDisOnBase,dataAGOnBase,  row=testRow ) 
				
	 

		#draw gaia raw
		if 1:
			#ax=plt.subplot(gs[6:,6: ])
			axRaw=plt.subplot(gs[0:5,6: ],sharex=ax)
			self.drawRawGaia( axRaw,gaiaOnsource,gaiaOffsource,gaiaForeCut=gaiaForeCut,row=testRow,figNameMark=None,markLoc=legendLoc,baseline=baseline ) 
			

	   

		#draw CO FITS

		#draw fits
		import pywcsgrid2
		from mpl_toolkits.axes_grid1.inset_locator import inset_axes
		import mpl_toolkits.axes_grid1.inset_locator as inset_locator
		


		backFISTHDR=fits.open(fitsName)

		backData=backFISTHDR[0].data
		backHead=backFISTHDR[0].header


		if len( backData.shape)==3:
			temp3D='temp3d.fits'
			myFITS.downTo2D(fitsName,outPUT=temp3D,overwrite=True) 

			backFISTHDR=fits.open(temp3D)
	
			backData=backFISTHDR[0].data
			backHead=backFISTHDR[0].header
				
		if len( backData.shape)==4:
			backData=backData[0]
			backData=backData[0]
			
 


 
		if backData.shape[0]==1:
			backData=backData[0]


		axBack=pywcsgrid2.subplot(gs[0:4,3:5], header=WCS(backHead) )
		
  
		self.drawCloudFITS( fitsName, axBack,backData, backHead,lRange,bRange,noiseLevel,signalLevel , gaiaOnsource,gaiaOffsource, maskFITS=maskFITS) 

		AGAVText="Ag"
		
		if self.useAV:
			AGAVText="Av"


		saveFigureName="{}_{}".format( self.sourceRow[disTB.sourceName] , '{}Baseline.pdf'.format(AGAVText) )

		#plt.savefig(  self.sourceName+'_extinctionGaiaAg.pdf', bbox_inches="tight")
		plt.savefig( saveFigureName, bbox_inches="tight")
		saveFigureName="{}_{}".format( self.sourceRow[disTB.sourceName] , '{}Baseline.png'.format(AGAVText) )

		print "saving as: " ,saveFigureName

		#plt.savefig(  self.sourceName+'_extinctionGaiaAg.pdf', bbox_inches="tight")
		plt.savefig( saveFigureName, bbox_inches="tight",dpi=600)
		return distance, distanceStd #lowerDis,upperDis

	def getArea(self, cropData,cutV=0.):
		
		"""
		Find the area of 3 sigma 
		
		The sigma is fitting with values less than 0.5
		
		"""
		
		if len(cropData.shape)==3 and cropData.shape[2]==1:
			cropData=cropData[0]

		
		
		shapeY,shapeX= cropData.shape
		
		oneD=cropData.reshape(shapeY*shapeX )
		
		fitGaussArr=oneD[oneD<cutV]
		fitGaussArr=fitGaussArr[fitGaussArr>-100]

		#truncaFit= norm.fit(fitGaussArr)  #truncnorm.fit(fitGaussArr)
		
		noiseBiase=  np.sqrt( np.mean( fitGaussArr**2) )
		 
		
		noise=noiseBiase +noiseBiase/4./len( fitGaussArr) #truncaFit[-1]
		

		
		
		
		signalPix=oneD[ oneD> noise*3]
		
		COsum=np.sum(signalPix)
		
		pixN= len(signalPix)
		areaCrop=round(pixN*30.*30/3600./3600.,2)
		print "Noise mean: {:.2f}------Noise Std: {:.2f}  area: {} ".format( 0  , noise,areaCrop)
		
		
		return areaCrop, COsum
		
		
		#print "?????????????????????"
		#print truncnorm.fit(fitGaussArr)
		
		
 
 

	def getAllGaiaByRow(self,newRow,extraCut=0. ):
		"""
		"""
		lowerDisCut=newRow[disTB.cutDistanceLower]
		upperDisCut=newRow[disTB.cutDistanceUpper]

		lowerPara=1./(upperDisCut+extraCut)*1000
		upperPara= 1./lowerDisCut*1000
		lRange,bRange=self.getLBRangeByRow(newRow)
		paraErrorCut=newRow[disTB.paraErrorCut]
		
		if self.useAV:
			
			#gaiaAllStars=doAV.getByLBRange(lRange,bRange,lowerDis=lowerDisCut,upperDis=upperDisCut+ extraCut)

			#print lRange,bRange,"Is this the right range fro Ursa Major?"


			gaiaAllStars=doAV.getByLBRange(lRange,bRange  )

			return gaiaAllStars 

		
		else:
			gaiaAllStars=self.gaiaDo.getByLBRange(lRange,bRange,lowerPara=lowerPara,upperPara=1000 , paraError=paraErrorCut)
			
			return gaiaAllStars 
		
		
		
	def getBoxPosition(self,newRow):
		return [ newRow[ disTB.boxCenterL  ] , newRow[ disTB.boxCenterB  ] ,newRow[ disTB.boxSizeL  ] ,newRow[ disTB.boxSizeB  ],0 ]
		
		
	def getLBRangeByRow(self, newRow ):
		boxPosition=self.getBoxPosition(newRow)

		return self.boxRegion(boxPosition )

	
	def calDisWithRowNormal(self,newRow,maskFITS=None,foregroundCOFITS=None,saveFigureName=None):
		#self.disRecord.addRowToTB(newRow) #save first in case discalculate are wrong..
		lRange,bRange=self.getLBRangeByRow(newRow)
		

		self.getDisArray(lRange,bRange, newRow[disTB.fitsFile] ,newRow[disTB.signalLevel],newRow[disTB.noiseLevel],\
		  newRow[disTB.sourceName] +"TB.fit",lowerDisCut=newRow[disTB.cutDistanceLower],upperDisCut=newRow[disTB.cutDistanceUpper],\
		  paraErrorCut=newRow[disTB.paraErrorCut],newQuery=True,disRow=newRow,maskFITS=maskFITS,foregroundCOFITS=foregroundCOFITS,saveFigureName=saveFigureName) #taurus
		
		self.disRecord.addRowToTB(newRow) #save second
		
		
		
		#cutDistanceLower="cutDistanceLower"
		#="cutDistanceUpper"
	def boxRegion(self,regionCoord):

		"""
		return lRange and B Range
		"""
		centerL,centerB,lSize,bSize,dummy=regionCoord
 
 
		return [centerL-lSize/2.,  centerL+lSize/2.   ], [centerB-bSize/2., centerB+bSize/2. ] 


	def drawTest(self):
		
		"""
		
		draw figure for eye confirm
		"""


	def calDisMain (self):
		#the interface of calculation distances


		


		if self.useBaseLine:
			
			self.gaiaDisBaseline(self.sourceName, getNewPara=self.redoPara, maskFITS=self.maskFITS, foregroundCOFITS=self.foregroundFITS) 
			
		else:
			
			self.calDisNormal(  self.sourcename,getPara=self.redoPara,maskFITS=self.maskFITS, foregroundCOFITS=self.foregroundFITS ) 
			

		
		
		

	def ZZZ(self):
		pass

	#def testDis



	