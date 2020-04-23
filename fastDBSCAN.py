
import os
import numpy as np
from astropy.table import Table,vstack
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
from skimage import measure
#fast way to perform DBSCAN
#
from progressbar import *
import math
from myPYTHON import *
from skimage.morphology import watershed
import sys
from skimage.morphology import erosion, dilation
from scipy.ndimage import label, generate_binary_structure,binary_erosion,binary_dilation
from sklearn.cluster import DBSCAN
from madda import  myG210
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
import seaborn as sns

from spectral_cube import SpectralCube
from astropy import modeling

from astropy import units as u
from  myGAIA import GAIADIS
import glob



gaiaDis=GAIADIS()

doFITS=myFITS()

doG210 = myG210()

def weighted_avg_and_std(values, weights):
	"""
	Return the weighted average and standard deviation.

	values, weights -- Numpy ndarrays with the same shape.
	"""
	average = np.average(values, weights=weights)
	# Fast and numerically precise:
	variance = np.average((values-average)**2, weights=weights)


	#if variance<0:
		#print weights

	return (average, math.sqrt(variance))



class myDBSCAN(object):

	rms = 0.5
	TBModel="minV3minP16_dendroCatTrunk.fit"

	tmpPath="./tmpFiles/"


	rawCOFITS="G2650Local30.fits"

	#
	pixelArea=0.25 # arcmins

	parsecToMeter= 3.0857e16 #m

	def __init__(self):
		pass

	def sumEdgeByCon1(self,extendMask): #7 in total
		raw=extendMask[1:-1,1:-1,1:-1]

		leftShiftZ=extendMask[0:-2, 1:-1, 1:-1]
		rightShiftZ=extendMask[2:, 1:-1 ,1:-1]

		leftShiftY=extendMask[1:-1, 0 : -2, 1:-1]
		rightShiftY=extendMask[1:-1, 2 : ,1:-1]

		leftShiftX=extendMask[1:-1, 1:-1,  0:-2]
		rightShiftX=extendMask[1:-1, 1:-1 , 2: ]


		sumAll=raw+leftShiftZ+rightShiftZ+leftShiftY+rightShiftY+leftShiftX+rightShiftX

		return  sumAll



	def sumEdgeByCon2(self,extendMask): #27 in total
		sumAll=extendMask[1:-1,1:-1,1:-1]*0
		Nz,Ny,Nx= sumAll.shape
		for i in [-1,0,1]:
			for j in [-1,0,1]:
				for k in [-1,0,1]:

					if np.sqrt( abs(i)+abs(j)+abs(k))>1.5:
						continue

					sumAll=sumAll+  extendMask[ 1+i:Nz+1+i , j+1:Ny+1+j , k+1: Nx+1+k  ]

		return  sumAll



	def sumEdgeByCon3(self,extendMask): #27 in total
		raw=extendMask[1:-1,1:-1,1:-1]
		Nz,Ny,Nx= raw.shape
		sumAll=raw*0
		for i in [-1,0,1]:
			for j in [-1,0,1]:
				for k in [-1,0,1]:
					sumAll=sumAll+  extendMask[ 1+i:Nz+1+i , j+1:Ny+1+j , k+1: Nx+1+k  ]


		return  sumAll




	def slowDBSCAN(self,COdata,COHead, min_sigma=2, min_pix=16, connectivity=2 ,region="" ,saveFITS=None ):
		"""
		Use the sklearn DBSCAN to calculate, just for comparison, to test the computeDBSCAN is right
		:param COdata:
		:param COHead:
		:param min_sigma:
		:param min_pix:
		:param connectivity:
		:param region:
		:return:
		"""
		###

		goodIndices=np.where(COdata>= min_sigma*self.rms   )

		coordinates= zip( goodIndices[0] , goodIndices[1] ,goodIndices[2]  )

		#eps=1.5 form connectivity 2,
		if connectivity==2:
			db = DBSCAN(eps=1.5  , min_samples= min_pix      ).fit(coordinates)

		if connectivity==1:
			db = DBSCAN(eps=1.1  , min_samples= min_pix      ).fit(coordinates)
		if connectivity==3:
			db = DBSCAN(eps=1.8  , min_samples= min_pix      ).fit(coordinates)

		labels = db.labels_
		print min(labels),"minimumLabel?"
		#u,c= np.unique(labels,return_counts=True)

		#print len(u)

		mask=np.zeros_like(COdata)-1

		mask[goodIndices]= labels
		if saveFITS==None:
			fits.writeto("dbscanMask1Sigma.fits",mask,header=COHead,overwrite=True)

		else:
			fits.writeto(saveFITS ,mask,header=COHead,overwrite=True)




	def computeDBSCAN(self,COdata,COHead, min_sigma=2, min_pix=16, connectivity=2 ,region="" , getMask=False,savePath="" ,mimicDendro=False, rmsFITS=None,inputRMS=None):
		"""
		min_pix the the minimum adjacen number for are core point, minPts
		:param COdata:
		:param min_sigma:
		:param min_pix:
		:param connectivity:
		:return:
		"""
		#pass

		if inputRMS==None:
			minValue = min_sigma*self.rms

		else:
			minValue = min_sigma* inputRMS



		Nz,Ny,Nx  = COdata.shape
		extendMask = np.zeros([Nz+2,Ny+2,Nx+2] ,dtype=int)


		if rmsFITS==None:

			goodValues= COdata>=minValue

		else:
			#pass,need to devide
			rmsData,rmsHead = myFITS.readFITS(rmsFITS)
			rmsCOData=COdata/rmsData
			goodValues= rmsCOData>=min_sigma

			#fits.writeto( "tttt.fits" , rmsCOData ,  header= COHead )

			#aaaa

		extendMask[1:-1,1:-1,1:-1] =  goodValues  #[COdata>=minValue]=1

		s=generate_binary_structure(3,connectivity)


		if connectivity==1:
			coreArray=self.sumEdgeByCon1(extendMask)

		if connectivity==2:
			coreArray=self.sumEdgeByCon2(extendMask)

		if connectivity==3:
			coreArray=self.sumEdgeByCon3(extendMask)

		coreArray = coreArray>=min_pix
		coreArray[  ~goodValues ]=False  # nan could be, #remove falsely, there is a possibility that, a bad value may have lots of pixels around and clould be
		coreArray=coreArray+0

		labeled_core, num_features=label(coreArray,structure=s) #first label core, then expand, otherwise, the expanding would wrongly connected

		selectExpand= np.logical_and(labeled_core==0,  goodValues   )
		#expand labeled_core
		#coreLabelCopy=labeled_core.copy()

		expandTry = dilation(labeled_core , s  ) # first try to expand, then only keep those region that are not occupied previously
		#it is possible  that a molecular cloud may have less pixN than 8, because of the closeness of two

		labeled_core[  selectExpand  ] =  expandTry[ selectExpand  ]

		labeled_array = labeled_core

		if mimicDendro:
			print "Mimicing dendrogram.."
			extendedArray=labeled_array>0
			extendedArray=extendedArray+0
			labeled_array, num_features = label(extendedArray, structure=s)

		saveName="{}dbscanS{}P{}Con{}.fits".format( region,min_sigma,min_pix,connectivity )

		if getMask:

			return labeled_array>0 #actually return mask


		print num_features,"features found!"

		fits.writeto(savePath+saveName, labeled_array, header=COHead, overwrite=True)
		return savePath+saveName



	def maskByGrow(self,COFITS,peakSigma=3,minV=1.):

		COData,COHead=myFITS.readFITS( COFITS )
		markers=np.zeros_like(COData )

		COData[COData<minV* self.rms]=0

		markers[COData>peakSigma*self.rms] = 1

		labels=watershed(COData,markers)
		fits.writeto("growMaskPeak3Min1.fits",labels,header=COHead,overwrite=True)



	def myDilation(self,scimesFITS,rawCOFITS,startSigma=20,endSigma=2, saveName="", maskCOFITS=None,savePath="" ):
		"""
		#because SCIMES removes weak emissions in the envelop of clouds, we need to add them back
		#one possible way is to use svm to split the trunk, test this con the  /home/qzyan/WORK/myDownloads/MWISPcloud/ClusterAsgn_ComplicateVe.fits

		:return:
		"""

		#cloudData,cloudHead = myFITS.readFITS("/home/qzyan/WORK/myDownloads/MWISPcloud/ClusterAsgn_ComplicateVe.fits")

		cloudData,cloudHead = myFITS.readFITS(scimesFITS)

		#rawFITS= rawCOFITS #"/home/qzyan/WORK/myDownloads/testScimes/complicatedTest.fits"


		if rawCOFITS!=None:
			rawCO,rawHead=   myFITS.readFITS( rawCOFITS )

		if maskCOFITS!=None:
			maskData,maskHead=myFITS.readFITS(maskCOFITS)

		#the expansion should stars from high coValue, to low CO values, to avoid cloud cross wak bounarires
		#sCon=generate_binary_structure(3,2)
		print "Expanding clous..."

		sigmaSteps= np.arange(startSigma,endSigma-1,-1)

		if endSigma not in sigmaSteps:
			sigmaSteps=list(sigmaSteps)
			sigmaSteps.append(endSigma)
		print "Step of sigmas, ", sigmaSteps

		cloudData = cloudData + 1  # to keep noise reagion  as 0
		dilationPath="./dilationMiddle/"

		for sigmas in sigmaSteps:

			#produceMask withDBSCAN
			#if sigmas>2:
			if maskCOFITS==None:
				COMask = self.computeDBSCAN( rawCO,rawHead, min_sigma=sigmas, min_pix=8, connectivity=2 ,region="" , getMask=True ) #use all of them
			else:
				COMask=maskData>sigmas*self.rms
			#else:
				#COMask = self.computeDBSCAN(  rawCO,rawHead, min_sigma=sigmas, min_pix=16, connectivity=2 ,region="" , getMask=True )


			#save middle FITS
			saveName=dilationPath+"dilationTo_{}sigma.fits".format(sigmas)
			fits.writeto( saveName  ,cloudData ,header=cloudHead,overwrite=True)

			for i in range(2000):

				rawAssign=cloudData.copy()
				d1Try=dilation(cloudData  ) #expand with connectivity 1, by defaults
				assignRegion=  np.logical_and(cloudData==0 , COMask )
				cloudData[ assignRegion ] = d1Try[ assignRegion ]


				diff=cloudData -rawAssign
				sumAll=np.sum(diff )
				if sumAll==0:
					print  "Sigmas: {}, Loop:{},  all difference:{}, beak".format(sigmas, i, sumAll)
					break

				else:
					print  "Sigmas: {}, Loop:{}, all difference:{} continue".format(sigmas, i, sumAll)
					continue




		cloudData=cloudData-1

		fits.writeto( savePath+saveName+"_extend.fits",cloudData ,header=cloudHead,overwrite=True)
		return savePath+saveName+"_extend.fits"



	def directLabel(self,COFITS,DBMaskFITS,min_sigma=3,min_pix=8,calCat=True ,useMask=True, peakSigma=3. ):

		saveMarker=""
		COData,COHead=myFITS.readFITS( CO12FITS )

		if useMask:
			DBMaskData,_=  myFITS.readFITS(  DBMaskFITS )

			maskData=np.zeros_like( DBMaskData )

			maskData[COData>min_sigma*self.rms]=1
			maskData[DBMaskData==0]=0
			saveLabel= "LabelSigma_{}_P{}.fits".format( min_sigma,min_pix )

		else:

			#use peak sigma to grow a mask

			maskData=np.zeros_like( COData )
			maskData[COData>min_sigma*self.rms]=1


			saveLabel= "NoMaskLabelSigma_{}_P{}.fits".format( min_sigma,min_pix )

			saveMarker="growMask"

		labels=measure.label(maskData,connectivity=1)
		fits.writeto(  saveLabel, labels,header=COHead, overwrite=True)


		if calCat  :

			self.getCatFromLabelArray(COFITS,saveLabel,self.TBModel,saveMarker=saveMarker,  minPix=min_pix,rms= min_sigma  )



	def getCatFromLabelArray(self,  CO12FITS,labelFITS,TBModel,minPix=8,rms=2 ,   saveMarker="", peakSigma=3,region="",pureDBSCAN=False):
		"""
		Extract catalog from label fits, the minPix and rms is only used for saving
		:param labelArray:
		:param head:
		:return:
		"""

		#do not make any selection here, because of the closeness of many clouds, some cloud may have pixels less than 8, we should keep them
		#they are usefull to mask edge sources..., and to clean fits

		if saveMarker=="":
			saveName= region+"DBSCAN{}_P{}Cat.fit".format(rms,minPix)

		else:
			saveName=saveMarker+".fit"

		clusterTBOld=Table.read( TBModel )

		###
		dataCO, headCO = myFITS.readFITS( CO12FITS )

		#dataCO=np.nan_to_num(dataCO), should not have Nan values


		dataCluster , headCluster=myFITS.readFITS( labelFITS )

		#create a data cube to find points that satisfy to have three consecutive points, can we find them ?
		minV = np.nanmin(dataCluster[0])


		Nz, Ny, Nx = dataCluster.shape

		#consecutiveData= np.zeros( ( Nz+2,Ny,Nx   )   )


		#b= dataCluster> minV


		#consecutiveData[1:-1]= dataCluster> minV
		#doSum= consecutiveData[0:-2]+consecutiveData[1:-1]+ consecutiveData[2: ]
		#consecutivepoints
		#P3 = doSum>=3

		wcsCloud=WCS( headCluster )

		clusterIndex1D= np.where( dataCluster>minV )
		clusterValue1D=  dataCluster[clusterIndex1D ]

		Z0,Y0,X0 = clusterIndex1D

		newTB= Table( clusterTBOld[0])
		newTB["sum"]=newTB["flux"]

		newTB["l_rms"]=newTB["v_rms"]
		newTB["b_rms"]=newTB["v_rms"]

		newTB["pixN"]=newTB["v_rms"]
		newTB["peak"]=newTB["v_rms"]

		newTB["peakL"]=newTB["v_rms"]
		newTB["peakB"]=newTB["v_rms"]
		newTB["peakV"]=newTB["v_rms"]
		newTB["area_accurate"]=newTB["v_rms"] # the column of area_accurate need to cos(b) facor

		#newTB["Nchannel"]=newTB["v_rms"] #number of channels, whis is used to sellect

		newTB["allChannel"]=newTB["v_rms"] # number channel involved
		newTB["has22"]=newTB["v_rms"] # number channel involved


		zeroProjection =  np.zeros( ( Ny , Nx  ) ) # one zero channel, used to get the projection area and
		zeroProjectionExtend = np.zeros( ( Ny+1, Nx+1 ) )



		idCol="_idx"


		#count all clusters

		#ids,count=np.unique(dataCluster,return_counts=True )
		ids,count=np.unique(  clusterValue1D,  return_counts=True  )
		GoodIDs=ids
		GoodCount=count
		print "Total number of turnks? ",len(ids)
		#print "Total number of Good Trunks? ",len(GoodIDs)

		#dataCO,headCO=doFITS.readFITS( CO12FITS )
		widgets = ['Recalculating cloud parameters: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),  ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options

		pbar = ProgressBar(widgets=widgets, maxval=len(GoodIDs))
		pbar.start()

		catTB=newTB.copy()
		catTB.remove_row(0)

		#zeroP

		for i in  range(len(GoodIDs)) :

			#i would be the newID
			newID= GoodIDs[i]
			pixN=GoodCount[i]

			newRow=newTB[0]

			newRow[idCol] = newID

			cloudIndex=self.getIndices(Z0,Y0,X0,clusterValue1D,newID)

			coValues=  dataCO[ cloudIndex ]

			#P3Value= P3[cloudIndex] # used to find

			#sortedCO=np.sort(coValues)
			#peak = sortedCO[-1] #np.max( coValues)
			#peak2=sortedCO[-2]
			peak= np.max(coValues)
			cloudV=cloudIndex[0]
			cloudB=cloudIndex[1]
			cloudL=cloudIndex[2]

			peakIndex=coValues.argmax()

			peakV = cloudV[peakIndex]
			peakB = cloudB[peakIndex]
			peakL = cloudL[peakIndex]

			#get the exact peak position, which would be used to

			projectIndex= tuple( [cloudB, cloudL ] )

			zeroProjection[projectIndex] =1

			#calculate the accurate

			indexB2D,indexL2D=np.where(zeroProjection==1 )



			_,BS2D, LS2D = wcsCloud.wcs_pix2world(indexL2D,indexB2D,  0, 0)




			area_accurate=np.sum( np.cos( np.deg2rad(BS2D) )    )*0.25

			newRow["area_accurate"]= area_accurate



			zeroProjectionExtend[0:-1,0:-1]=zeroProjection

			sumCO=np.sum( coValues )




			Vcen,Vrms= weighted_avg_and_std(cloudV, coValues )
			Bcen,Brms= weighted_avg_and_std(cloudB, coValues )
			Lcen,Lrms= weighted_avg_and_std(cloudL, coValues )

			#calculate the exact area

			#LBcore = zip(cloudB, cloudL)
			#pixelsN= {}.fromkeys(LBcore).keys() #len( set(LBcore) )
			#area_exact=len(pixelsN)*0.25 #arc mins square
			area_exact= np.sum( zeroProjection )*0.25
			#print area_exact,area_accurate
			#find the 2*2 patter

			sum22= zeroProjectionExtend[0:-1,0:-1] +  zeroProjectionExtend[0:-1,1: ]+ zeroProjectionExtend[1: , 0 :-1]+zeroProjectionExtend[1: , 1: ]

			#if any pixel>4:

			if 4 in sum22:
				newRow["has22"] = 1
			else:
				newRow["has22"] = 0

			diffVs= np.unique(cloudV)

			#dataClusterNew[cloudIndex] =newID

			#save values
			newRow["pixN"]= pixN
			newRow["peak"]= peak

			newRow["peakV"]= peakV
			newRow["peakB"]= peakB
			newRow["peakL"]= peakL

			#newRow["peak2"]= peak2

			newRow["sum"]= sumCO
			newRow["area_exact"]= area_exact

			newRow["x_cen"],  newRow["y_cen"], newRow["v_cen"]= wcsCloud.wcs_pix2world( Lcen, Bcen,Vcen ,0)
			newRow["v_cen"]= newRow["v_cen"]/1000.
			dv=headCluster["CDELT3"]/1000. #km/s

			dl= abs( headCluster["CDELT1"] ) #deg

			newRow["v_rms"] = Vrms*dv

			newRow["l_rms"] = Lrms*dl
			newRow["b_rms"] = Brms*dl

			#_, Nchan=np.unique( cloudV, return_counts=True)

			#newRow["Nchannel"] =    np.max(P3Value)# if there is a three consecutive spectra in the cloud
			newRow["allChannel"] =   len( diffVs )



			catTB.add_row(newRow)

			zeroProjection[projectIndex] = 0
			zeroProjectionExtend[ 0:-1,0:-1 ]=zeroProjection

			pbar.update(i)



		pbar.finish()
		#save the clouds

		#fits.writeto(self.regionName+"NewCloud.fits", dataClusterNew,header=headCluster,overwrite=True   )
		catTB.write( saveName ,overwrite=True)

		return saveName

	def getSumToFluxFactor(self):

		theta =  np.deg2rad(0.5/60)
		omega = theta * theta
		f=115.271202000
		waveLength =299792458/(f*1e9)
		k= 1.38064852e3 #has converted to jansky
		factorSumToFlux=  2*k*omega/waveLength/waveLength
		
		return factorSumToFlux
	def converSumToFlux(self,sumRow):

		factorSumToFlux=self.getSumToFluxFactor(factorSumToFlux)

		return sumRow* factorSumToFlux #jansky

	def converFluxToSum(self, fluxRow):
		factorSumToFlux=self.getSumToFluxFactor(factorSumToFlux)

		return fluxRow/factorSumToFlux



	def getIndices(self,Z0,Y0,X0,values1D,choseID):


		cloudIndices = np.where(values1D==choseID )

		cX0=X0[cloudIndices ]
		cY0=Y0[cloudIndices ]
		cZ0=Z0[cloudIndices ]

		return tuple( [ cZ0, cY0, cX0 ]  )



	def getIndices2D(self, Y0,X0,values1D,choseID):



		cloudIndices = np.where(values1D==choseID )

		cX0=X0[cloudIndices ]
		cY0=Y0[cloudIndices ]


		return tuple( [   cY0, cX0 ]  )




	def draw(self ):
		"""
		#draw compare of
		:return:
		"""



		fig=plt.figure(figsize=(12,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })

		drawTB=Table.read( "Sigma1_P25FastDendro.fit" )


		axNumber=fig.add_subplot(1,2,1)




		axArea= fig.add_subplot(1,2,2)

		areaEdges=np.linspace(0,6,1000)
		areaCenter=self.getEdgeCenter( areaEdges )


		totalTB=  [drawTB] #TBList1+TBList2

		for i in range( len(totalTB) ):

			eachTB = totalTB[i]

			binN,binEdges=np.histogram(eachTB["area_exact"]/3600., bins=areaEdges  )


			axArea.plot( areaCenter[binN>0],binN[binN>0], 'o-'  , markersize=1, lw=0.8  ,alpha= 0.5 )


		axArea.set_yscale('log')
		axArea.set_xscale('log')


		axArea.legend()



		axArea.set_xlabel(r"Area (deg$^2$)")
		axArea.set_ylabel(r"Number of trunks")


		plt.savefig( "compareDendroParaDBMask.pdf" ,  bbox_inches='tight')
		plt.savefig( "compareDendroParaDBMask.png" ,  bbox_inches='tight',dpi=300)


	def getEdgeCenter(self,edges):

		areaCenters= ( edges[1:] + edges[0:-1] )/2.

		return  areaCenters

	def cleanDBTB(self,dbTB,pixN=8,minV=3,minDelta=3):

		"""
		The minimum Peak, should be relative to the minValue
		:param dbTB:
		:param pixN:
		:param peak:
		:return:
		"""
		peakV=(minV + minDelta)*self.rms



		if type(dbTB)==list:
			newList=[]
			for eachT in dbTB:

				goodT=eachT.copy()
				goodT=goodT[ goodT["pixN"] >= pixN ]
				goodT=goodT[ goodT["peak"] >= peakV ]

				newList.append(goodT)

			return newList



		else:

			goodT=dbTB.copy()
			goodT=goodT[ goodT["pixN"] >= pixN ]
			goodT=goodT[ goodT["peak"] >= peakV ]

			return goodT





	def drawDBSCANArea(self):

		TB2_16= "G2650CO12DBCatS2P16Con2.fit"
		#TB2_16= "DBSCAN2_9Sigma1_P1FastDBSCAN.fit"
		TB25_9="G2650CO12DBCatS2.5P9Con2.fit"
		TB35_9="G2650CO12DBCatS3.5P9Con2.fit"
		TB45_9="G2650CO12DBCatS4.5P9Con2.fit"
		TB55_9="G2650CO12DBCatS5.5P9Con2.fit"
		TB65_9="G2650CO12DBCatS6.5P9Con2.fit"
		TB75_9="G2650CO12DBCatS7.5P9Con2.fit"

		TB3_9= "G2650CO12DBCatS3.0P9Con2.fit"
		TB4_9= "G2650CO12DBCatS4.0P9Con2.fit"
		TB5_9= "G2650CO12DBCatS5.0P9Con2.fit"
		TB6_9= "G2650CO12DBCatS6.0P9Con2.fit"
		TB7_9= "G2650CO12DBCatS7.0P9Con2.fit"



		TBFiles=[TB2_16,TB25_9,TB3_9, TB35_9, TB4_9, TB45_9,TB5_9, TB55_9, TB6_9 , TB65_9, TB7_9, TB75_9   ]



		sigmas=[2,2.5,3,3.5,4,4.5, 5, 5.5, 6, 6.5,7,7.5]

		labelStr=[  r"2$\sigma$, P16" ,   r"2.5$\sigma$, P16" ,  r"3$\sigma$, P16" ,  r"3.5$\sigma$, P16" ,   r"4$\sigma$, P16" ,  r"4.5$\sigma$, P16" , \
		            r"5$\sigma$, P16" ,  r"5.5$\sigma$, P16" , r"6$\sigma$, P16"  , r"6.5$\sigma$, P16"  , r"7$\sigma$, P16", r"7.5 $\sigma$, P16"   ]


		TBList=[]



		areaEdges=np.linspace(0,6,1000)
		areaCenter=self.getEdgeCenter( areaEdges )
		fig=plt.figure(figsize=(12,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })



		axArea=fig.add_subplot(1,2,1)



		for eachTBF,eachLab in zip(TBFiles,labelStr):
			tb=Table.read(eachTBF)

			tb=self.removeWrongEdges(tb)

			TBList.append( tb )


			goodT=tb


			goodT=goodT[ goodT["pixN"]>=16 ]

			goodT=goodT[ goodT["peak"]>=1.5 ]

			#
			binN,binEdges=np.histogram(goodT["area_exact"]/3600., bins=areaEdges  )




			axArea.plot( areaCenter[binN>0],binN[binN>0], 'o-'  , markersize=1, lw=0.8,label=eachLab ,alpha= 0.5 )


		axArea.set_yscale('log')
		axArea.set_xscale('log')

		axArea.legend()
		axArea.set_title("Plot of Area distribution with DBSCAN")

		plt.savefig( "dbscanArea.png" ,  bbox_inches='tight',dpi=300)


	def drawDBSCANNumber(self):

		minPix=8

		TB2_16="G2650CO12DBCatS2.0P{}Con2.fit".format(minPix)
		TB25_9="G2650CO12DBCatS2.5P{}Con2.fit".format(minPix)
		TB35_9="G2650CO12DBCatS3.5P{}Con2.fit".format(minPix)
		TB45_9="G2650CO12DBCatS4.5P{}Con2.fit".format(minPix)
		TB55_9="G2650CO12DBCatS5.5P{}Con2.fit".format(minPix)
		TB65_9="G2650CO12DBCatS6.5P{}Con2.fit".format(minPix)
		TB75_9="G2650CO12DBCatS7.5P{}Con2.fit".format(minPix)

		TB3_9= "G2650CO12DBCatS3.0P{}Con2.fit".format(minPix)
		TB4_9= "G2650CO12DBCatS4.0P{}Con2.fit".format(minPix)
		TB5_9= "G2650CO12DBCatS5.0P{}Con2.fit".format(minPix)
		TB6_9= "G2650CO12DBCatS6.0P{}Con2.fit".format(minPix)
		TB7_9= "G2650CO12DBCatS7.0P{}Con2.fit".format(minPix)


		TBFiles=[TB2_16,TB25_9,TB3_9, TB35_9, TB4_9, TB45_9,TB5_9, TB55_9, TB6_9 , TB65_9, TB7_9, TB75_9   ]
		TBList=[]

		Nlist=[]

		sigmas=[2,2.5,3,3.5,4,4.5, 5, 5.5, 6, 6.5,7,7.5]

		for eachTBF in TBFiles:
			#
			tb=Table.read(eachTBF)
			TBList.append( tb )
			goodT=tb

			goodT=goodT[ goodT["pixN"]>=16 ]
			goodT=goodT[ goodT["peak"]>=1.5 ]

			Nlist.append(len(goodT) )

		fig=plt.figure(figsize=(12,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })

		drawTB=Table.read( "Sigma1_P25FastDendro.fit" )


		axNumber=fig.add_subplot(1,2,1)

		axNumber.plot(sigmas,Nlist,'o-',color='blue')

		axNumber.set_ylabel(r"Total number of trunks")
		axNumber.set_xlabel(r"CO cutoff ($\sigma$)")
		axNumber.set_title("Plot of total trunk numbers with DBSCAN")

		plt.savefig( "dbscanNumber.png" ,  bbox_inches='tight',dpi=300)

	def getRealArea(self,TB):
		#print TB.colnames
		araeList=[]

		#print TB.colnames
		for eachR in TB:
			v = eachR["v_cen"]
			#dis= ( 0.033*v + 0.175)*1000 # pc
			dis= ( 0.033*v + 0.180)*1000 # pc

			if dis<0  or dis> 1500 :
				continue

			N=eachR["area_exact"]/0.25


			#print N,eachR["pixN"]
			length= dis*np.deg2rad(0.5/60) # square pc
			trueArea=length**2*N #eachR["pixN"]  #*10000
			#print N,  trueArea

			araeList.append(  trueArea )

		return np.asarray(araeList)

	def getMass(self,TB):

		massList=[]

		#print TB.colnames
		for eachR in TB:
			v = eachR["v_cen"]
			#dis= ( 0.033*v + 0.175)*1000 # pc
			dis= ( 0.033*v + 0.180)*1000 # pc

			if dis<0  or dis> 1500 :
				continue

			fluxSum= eachR["sum"]*0.2
			massSingle= self.calmassByXfactor(fluxSum,dis) #eachR["pixN"]  #*10000
			#print N,  trueArea

			massList.append(  massSingle )

		return np.asarray(massList)


	def calmassByXfactor(self, coInt, distance, xFactor=2.0e20):

		"""
		The unit of coInt must be K km/s,
		distance pc
		not calDis, calMass
		"""

		NH2 = coInt * xFactor # cm-2 #

		# distance=2200 # pc

		# parsecToMeter= 3.0857e16 #m
		# length1=np.radians(degSize1)*self.distance*self.parsecToMeter*100. #cm
		# length2=np.radians(degSize2)*self.distance*self.parsecToMeter*100.

		# should use single pix for 12CO

		length1 = np.radians(30. / 60. / 60.) * distance * self.parsecToMeter * 100.  # cm
		# length2=np.radians(degSize2)*self.distance*self.parsecToMeter*100.

		mu = 1.36

		Mh2 = 3.35e-27  # kg
		solarMass = 1.9891e30  # kg
		# s=np.pi*length1*length1
		s = length1 * length1

		coreSolar = s * NH2 * Mh2 * mu / solarMass
		return coreSolar


	def drawAreaDistribute(self,TBName,region="",algorithm='Dendrogram'):
		"""

		:return:
		"""

		TB=Table.read( TBName )

		TBLOcal=Table.read("DBSCAN35_9Sigma1_P1FastDBSCAN.fit")
		TBAll=vstack([TB,TBLOcal ])

		areaEdges=np.linspace(0,6,1000)
		areaCenter=self.getEdgeCenter( areaEdges )



		fig=plt.figure(figsize=(8,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })
		axArea=fig.add_subplot(1,1,1)

		##########
		goodT=TB

		if "pixN" in goodT.colnames:

			goodT=goodT[ goodT["pixN"]>=16 ]
			goodT=goodT[ goodT["peak"]>= self.rms*5. ]
		binN,binEdges=np.histogram(goodT["area_exact"]/3600., bins=areaEdges  )
		axArea.plot( areaCenter[binN>0],binN[binN>0], 'o-'  , markersize=1, lw=0.8,  alpha= 0.5, label= region)#   r"SCIMES,min3$\sigma$P16 12CO"  )
		print region
		self.getAlphaWithMCMC(  goodT["area_exact"]  , minArea= 50.*0.25/3600  , maxArea=None , physicalArea=False )
		print region,"Done"
		#a=np.linspace(1,3000,6000)

		#trueArea=1./a**2
		#self.getAlphaWithMCMC(  trueArea , minArea= 1e-7, maxArea=None , physicalArea=True )

		#print "Above???-2??"



		###############
 		goodT=TBLOcal

		if "pixN" in goodT.colnames:

			goodT=goodT[ goodT["pixN"]>=16 ]
			goodT=goodT[ goodT["peak"]>=1.5 ]
		binN,binEdges=np.histogram(goodT["area_exact"]/3600., bins=areaEdges  )
		axArea.plot( areaCenter[binN>0],binN[binN>0], 'o-'  , markersize=1, lw=0.8,  alpha= 0.5 ,label="Velocity range (0-30 km/s)12CO Raw"  )


		areaEdges=np.linspace(0,100,1000)
		areaCenter=self.getEdgeCenter( areaEdges )



		realArea=self.getRealArea(goodT)
		binN,binEdges=np.histogram( realArea  , bins=areaEdges  )

		axArea.plot( areaCenter[binN>0],binN[binN>0], 'o-'  , markersize=1, lw=0.8,  alpha= 0.5 ,label="Velocity range (0-30 km/s)12CO, distance Corrected"  )
		print min(realArea),"The minimum area?"
		self.getAlphaWithMCMC(  realArea  ,minArea= 0.42836824657505895 , maxArea=None,  physicalArea=True)

		areaEdges=np.linspace(0,6,1000)
		areaCenter=self.getEdgeCenter( areaEdges )


		############### Perseus
 		goodT=  Table.read("Local13DBSCAN3_9.fit")

		if "pixN" in goodT.colnames:

			goodT=goodT[ goodT["pixN"]>=16 ]
			goodT=goodT[ goodT["peak"]>=1.5 ]
		binN,binEdges=np.histogram(goodT["area_exact"]/3600., bins=areaEdges  )
		axArea.plot( areaCenter[binN>0],binN[binN>0], 'o-'  , markersize=1, lw=0.8,  alpha= 0.5 ,label=r"(26$^\circ$-50$^\circ$)13CO"  )


		############### Perseus
 		goodT=  Table.read("G210DBSCAN3_9.fit")

		if "pixN" in goodT.colnames:

			goodT=goodT[ goodT["pixN"]>=16 ]
			goodT=goodT[ goodT["peak"]>=1.5 ]
		binN,binEdges=np.histogram(goodT["area_exact"]/3600., bins=areaEdges  )
		axArea.plot( areaCenter[binN>0],binN[binN>0], 'o-'  , markersize=1, lw=0.8,  alpha= 0.5 ,label=r"G210(210$^\circ$-220$^\circ$)12CO"  )

		###############
 		goodT=  Table.read("DBSCAN3_9.fit")

		if "pixN" in goodT.colnames:

			goodT=goodT[ goodT["pixN"]>=16 ]
			goodT=goodT[ goodT["peak"]>=1.5 ]
		binN,binEdges=np.histogram(goodT["area_exact"]/3600., bins=areaEdges  )
		axArea.plot( areaCenter[binN>0],binN[binN>0], 'o-'  , markersize=1, lw=0.8,  alpha= 0.5 ,label="Velocity range (30-60 km/s)12CO"  )









		###############

		axArea.set_yscale('log')
		axArea.set_xscale('log')
		axArea.set_xlabel(r"Area (deg$^2$)")
		axArea.set_ylabel(r"Bin number of trunks ")



		axArea.legend()
		axArea.set_title("Plot of Area distribution with DBSCAN")

		plt.savefig( region+"dbscanArea.png" ,  bbox_inches='tight',dpi=300)
		plt.savefig( region+"dbscanArea.pdf" ,  bbox_inches='tight' )




	def getAlphaWithMCMC(self,areaArray,minArea=0.03,maxArea=1.,physicalArea=False,verbose=True,plotTest=False,saveMark="" ):
		"""
		areaArray should be in square armin**2
		:param areaArray:
		:param minArea:
		:param maxArea:
		:return:
		"""

		print "Fitting index with MCMC..."

		if not physicalArea:
			areaArray=areaArray/3600.

		if maxArea!=None:
			select=np.logical_and( areaArray>minArea, areaArray<maxArea)

		else:
			select= areaArray>minArea

		rawArea =   areaArray[ select ]

		if verbose:
			print "Run first chain for {} molecular clouds.".format( len( rawArea ) )
		part1=doG210.fitPowerLawWithMCMCcomponent1(rawArea, minV=minArea, maxV=maxArea)
		if verbose:
			print "Run second chain for {} molecular clouds.".format( len(rawArea) )

		part2=doG210.fitPowerLawWithMCMCcomponent1(rawArea, minV=minArea, maxV=maxArea)

		allSample=np.concatenate(  [ part1 , part2 ]    )



		#test plot
		if plotTest:
			fig = plt.figure(figsize=(12, 6))
			ax0 = fig.add_subplot(1, 1, 1)
			# fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
			rc('text', usetex=True)
			rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})

			ax0.scatter(part1,part2,s=10 )

			plt.savefig("mcmcSampleTest.pdf"  , bbox_inches='tight')
			aaaaaa

		meanAlpha= np.mean( allSample)
		stdAlpha=  np.std(allSample,ddof=1)
		if verbose:
			print "Alpha Mean: {:.2f}; std: {:.2f}".format( meanAlpha,  stdAlpha)

		return round(meanAlpha,2) , round(stdAlpha,2)

	def drawSumDistribute(self,TBName,region=""):
		"""
		:return:
		"""

		TB=Table.read( TBName )

		TBLOcal=Table.read("DBSCAN35_9Sigma1_P1FastDBSCAN.fit")
		TBAll=vstack([TB,TBLOcal ])



		goodT=TB





		fig=plt.figure(figsize=(12,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })
		axArea=fig.add_subplot(1,1,1)

		##########
		pixNCol =goodT["flux"]

		logPixN=np.log10( pixNCol  )

		print min( logPixN  ),max( logPixN  )

		areaEdges=np.linspace( min( logPixN  ),max( logPixN  ) ,100)
		areaCenter=self.getEdgeCenter( areaEdges )

		binN,binEdges=np.histogram( logPixN , bins=areaEdges  )


		drawBind=binN[binN>0]
		drawCenter= areaCenter[binN>0]

		axArea.plot(  drawCenter , np.log10(drawBind) , 'o-'  , markersize=1, lw=0.8,  alpha= 0.5 ,label="Flux"  )

		select=np.logical_and( drawCenter<6, drawCenter>4 )
		x=drawCenter[ select ]   #np.log(drawCenter)
		y=np.log10(drawBind  )[ select]

		#print np.polyfit(x,y,1)

		###########################



		if 0:
			tbVox=goodT[ goodT["pixN"]>16  ]
			pixNCol =tbVox["pixN"]

			logPixN=np.log10( pixNCol  )

			print min( logPixN  ),max( logPixN  )

			areaEdges=np.linspace( min( logPixN  ),max( logPixN  ) ,100)
			areaCenter=self.getEdgeCenter( areaEdges )

			binN,binEdges=np.histogram( logPixN , bins=areaEdges  )


			drawBind=binN[binN>0]
			drawCenter= areaCenter[binN>0]

			axArea.plot(  drawCenter , np.log10(drawBind) , 'o-'  , markersize=1, lw=0.8,  alpha= 0.5 ,label="Voxel"  )

			select=np.logical_and( drawCenter<4, drawCenter>1.5 )
			x=drawCenter[ select ]   #np.log(drawCenter)
			y=np.log10(drawBind  )[ select]

			print np.polyfit(x,y,1)

		###########################




		##########
		tbVox=goodT
		pixNCol =tbVox["area_exact"]

		logPixN=np.log10( pixNCol  )
		print "draw areas"
		print min( logPixN  ),max( logPixN  )

		areaEdges=np.linspace( min( logPixN  ),max( logPixN  ) ,15)

		print  areaEdges

		areaCenter=self.getEdgeCenter( areaEdges )

		binN,binEdges=np.histogram( logPixN , bins=areaEdges  )


		drawBind=binN[binN>0]
		drawCenter= areaCenter[binN>0]

		axArea.plot(  drawCenter , np.log10(drawBind) , 'o-'  , markersize=1, lw=0.8,  alpha= 0.5 ,label="area exact aa "  )

		select=np.logical_and( drawCenter<4, drawCenter>1.5  )


		x=drawCenter[ select ]   #np.log(drawCenter)
		y=np.log10(drawBind  )[ select]

		a= np.polyfit(x,y,1)
		p=np.poly1d(a)
		axArea.plot(  drawCenter ,  p(drawCenter) , 'o-'  , markersize=1, lw=0.8,  alpha= 0.5   )

		print a
		###########################








		axArea.set_xlabel(r"Voxel Number")
		axArea.set_ylabel(r"Bin number of trunks ")



		axArea.legend()
		axArea.set_title("Plot of Pixel distribution with DBSCAN")

		plt.savefig( region+"dbscanTotalPixel.png" ,  bbox_inches='tight',dpi=300)





	def drawPixNDistribute(self,TBName,region=""):
		"""

		:return:
		"""

		TB=Table.read( TBName )

		TBLOcal=Table.read("DBSCAN35_9Sigma1_P1FastDBSCAN.fit")
		TBAll=vstack([TB,TBLOcal ])



		goodT=TB
		goodT=goodT[ goodT["pixN"]>=16 ]

		pixNCol =goodT["pixN"]

		logPixN=np.log10( pixNCol  )

		print logPixN



		areaEdges=np.linspace( min( logPixN  ),max( logPixN  ) ,100)
		areaCenter=self.getEdgeCenter( areaEdges )



		fig=plt.figure(figsize=(12,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })
		axArea=fig.add_subplot(1,1,1)

		##########

		binN,binEdges=np.histogram( logPixN , bins=areaEdges  )


		drawBind=binN[binN>0]
		drawCenter= areaCenter[binN>0]

		axArea.plot(  drawCenter , np.log10(drawBind) , 'o-'  , markersize=1, lw=0.8,  alpha= 0.5 ,label=" ??"  )

		select=np.logical_and( drawCenter<4, drawCenter>1.5 )
		x=drawCenter[ select ]   #np.log(drawCenter)
		y=np.log(drawBind  )[ select]

		print np.polyfit(x,y,1)

		###############


		axArea.set_xlabel(r"Voxel Number")
		axArea.set_ylabel(r"Bin number of trunks ")



		axArea.legend()
		axArea.set_title("Plot of Pixel distribution with DBSCAN")

		plt.savefig( region+"dbscanTotalPixel.png" ,  bbox_inches='tight',dpi=300)





	def roughFit(self,centers,bins ):
		"""

		:return:
		"""
		y= bins[bins>0 ]  # areaCenter[binN>0]
		x=  centers[ bins>0  ]




		x1= x[x<= 0.1]  # areaCenter[binN>0]
		y1=  y[x<= 0.1 ]


		x2= x1[x1>=0.005 ]  # areaCenter[binN>0]
		y2=  y1[x1>= 0.005 ]

		x=np.log10(x2)
		y=np.log10(y2)


		print x
		print y

		print np.polyfit(x,y,1)


		return




	def drawTrueArea(self):

		goodTB=Table.read( "/home/qzyan/WORK/projects/maddalena/dendroDisPath/G2650/G2650goodDisTB.fit"  )
		dendroTB=Table.read( "/home/qzyan/WORK/myDownloads/testScimes/mosaicV1NewTB.fit" )

		areas=[]

		for eachG in goodTB:

			d=eachG["distance"]
			ID=int( eachG["sourceName"].split('oud')[1]  )
			dendroRow=  dendroTB[ID-1 ]

			area=dendroRow["area_exact"]

			area/0.25*(d*np.radians( 0.5/60. ) )**2
			#print area


			areas.append(area )
		#plot
		fig=plt.figure(figsize=(12,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })

		drawTB=Table.read( "Sigma1_P25FastDendro.fit" )

		bins=np.linspace( np.min(areas),np.max(areas),5  )
		areaCenter=self.getEdgeCenter( bins )

		ax=fig.add_subplot(1,2,1)

		binN,binEdges=np.histogram(areas, bins=bins  )

		ax.scatter(areaCenter,binN   )

		ax2=fig.add_subplot(1,2,2)


		ax2.scatter( goodTB["vlsr"], goodTB["distance"]  )


		plt.savefig( "exactArea.png" ,  bbox_inches='tight',dpi=300)




	def removeWrongEdges(self,TB):


		if TB==None:
			return None
		processTB=TB.copy()

		#remove cloudsThat touches the noise edge of the fits


		#part1= processTB[ np.logical_and( processTB["x_cen"]>=2815 ,processTB["y_cen"]>= 1003  )   ] #1003, 3.25

		#part2= processTB[ np.logical_and( processTB["x_cen"]<= 55 ,processTB["y_cen"]>= 1063  )   ] #1003, 3.25

		if "peak" in TB.colnames: #for db scan table

			part1= processTB[ np.logical_or( processTB["x_cen"]>26.25 ,processTB["y_cen"] < 3.25  )   ] #1003, 3.25
			#part1= processTB[ np.logical_or( processTB["x_cen"]>26.24166667 ,processTB["y_cen"] < 3.25833333 )   ] #1003, 3.25

			part2= part1[ np.logical_or( part1["x_cen"]<49.25 ,part1["y_cen"]<  3.75 )   ] #1003, 3.25
			#part2= part1[ np.logical_or( part1["x_cen"]<49.24166667 ,part1["y_cen"]<  3.75833333 )   ] #1003, 3.25

			return part2
		else: #dendrogram tb

			part1= processTB[ np.logical_or( processTB["x_cen"]< 2815 ,processTB["y_cen"] < 1003  )   ] #1003, 3.25

			part2= part1[ np.logical_or( part1["x_cen"]>  55 ,part1["y_cen"]< 1063  )   ] #1003, 3.25

			return part2



	def removeAllEdges(self,TBList):
		"""

		:param TBList:
		:return:
		"""
		newList=[]

		for eachTB in TBList:
			newList.append( self.removeWrongEdges(eachTB) )

			
		return newList

	def getNList(self,TBList):

		Nlist=[]

		for eachTB in TBList:
			Nlist.append( len(eachTB) )
		return Nlist



	def getTotalFluxList(self,TBList):

		fluxlist=[]
		omega =1  # np.deg2rad(0.5 / 60) * np.deg2rad(0.5 / 60) / 4. / np.log(2)

		for eachTB in TBList:

			if "sum" in eachTB.colnames:
				toalFlux=np.nansum( eachTB["sum"]  )*0.2*omega # K km/s

			else:
				toalFlux=np.nansum( eachTB["flux"]  )*0.2/self.getSumToFluxFactor()*omega # K km/s

			fluxlist.append(  toalFlux  )
		return fluxlist

	def fluxAlphaDistribution(self, onlyLocal=False):
		"""

		# draw alph distribution of flux, for each DBSCAN,

		#  onlyLocal,

		:return:
		"""

		algDendro = "Dendrogram"
		tb8Den, tb16Den, label8Den, label16Den, sigmaListDen = self.getTBList(algorithm=algDendro)
		tb8Den = self.removeAllEdges(tb8Den)
		tb16Den = self.removeAllEdges(tb16Den)

		if onlyLocal:
			tb8Den = self.removePerseusList( tb8Den )
			tb16Den = self.removePerseusList( tb16Den )



		algDB = "DBSCAN"
		tb8DB, tb16DB, label8DB, label16DB, sigmaListDB = self.getTBList(algorithm=algDB)
		tb8DB = self.removeAllEdges(tb8DB)
		tb16DB = self.removeAllEdges(tb16DB)

		if onlyLocal:
			tb8DB = self.removePerseusList(tb8DB)
			tb16DB = self.removePerseusList(tb16DB)




		algSCI="SCIMES"
		tb8SCI,tb16SCI,label8SCI,label16SCI,sigmaListSCI=self.getTBList(algorithm= algSCI )
		tb8SCI=self.removeAllEdges(tb8SCI)
		tb16SCI=self.removeAllEdges(tb16SCI)
		

		if onlyLocal:
			tb8SCI = self.removePerseusList(tb8SCI)
			tb16SCI = self.removePerseusList(tb16SCI)




		fig = plt.figure(figsize=(8, 6))
		rc('text', usetex=True)
		rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]



		#############   plot dendrogram
		axDendro = fig.add_subplot(1, 1, 1)

		alphaDendro, alphaDendroError = self.getFluxAlphaList(tb16Den, sigmaListDen)


		if onlyLocal: #
			ebDendro=axDendro.errorbar(sigmaListDen, alphaDendro, yerr=alphaDendroError, c='g', marker='D', capsize=3 , elinewidth=1.3, lw=1, label=algDendro+r" ($|b|>1^{\circ}$)")
			
		else:
			ebDendro=axDendro.errorbar(sigmaListDen, alphaDendro, yerr=alphaDendroError, c='g', marker='D', capsize=3 , elinewidth=1.3, lw=1, label=algDendro)


		axDendro.set_ylabel(r"$\alpha$ (flux)")
		axDendro.set_xlabel(r"CO cutoff ($\sigma$)")

		##############   plot DBSCAN

		alphaDB, alphaDBError =   self.getFluxAlphaList(tb16DB, sigmaListDB)
		
		if onlyLocal: #
			ebDBSCAN = axDendro.errorbar(sigmaListDB, alphaDB, yerr=alphaDBError, c='b', marker='^',linestyle=":" , capsize=3 , elinewidth=1.0, lw=1, label=algDB+r" ($|b|>1^{\circ}$)" ,  markerfacecolor='none' )
		else:
			ebDBSCAN = axDendro.errorbar(sigmaListDB, alphaDB, yerr=alphaDBError, c='b', marker='^',linestyle=":" , capsize=3 , elinewidth=1.0, lw=1, label=algDB ,  markerfacecolor='none' )

		ebDBSCAN[-1][0].set_linestyle(':')

		##########plot SCIMES

		alphaSCI, alphaDBError =   self.getFluxAlphaList(tb16SCI, sigmaListSCI)



		if onlyLocal: #
			ebSCISCAN = axDendro.errorbar(sigmaListSCI, alphaSCI, yerr=alphaDBError, c='purple', marker='s',linestyle="--", capsize=3 , elinewidth=1.0, lw=1, label=algSCI+r" ($|b|>1^{\circ}$)",  markerfacecolor='none' )

		else:
			ebSCISCAN = axDendro.errorbar(sigmaListSCI, alphaSCI, yerr=alphaDBError, c='purple', marker='s',linestyle="--", capsize=3 , elinewidth=1.0, lw=1, label=algSCI,  markerfacecolor='none' )



		ebSCISCAN[-1][0].set_linestyle('--')


		if 1: #plot average alpha
			allAlpha = alphaDB + alphaDendro
			allAlphaError = alphaDBError + alphaDendroError

			print "Average error, ", np.mean(allAlphaError)

			errorAlpha = np.mean(allAlphaError) ** 2 + np.std(allAlpha, ddof=1) ** 2
			errorAlpha = np.sqrt(errorAlpha)

			alphaMean = np.mean(allAlpha)

			print "The mean alpha of flux distribution is {:.2f}, error is {:.2f}".format(alphaMean, errorAlpha)
			labelAverage = r"{:.2f}$\pm${:.2f}".format(alphaMean, errorAlpha)
			axDendro.plot([min(sigmaListDB), max(sigmaListDB)], [alphaMean, alphaMean], '--', color='black', lw=1,label=labelAverage)

		plt.xticks( [2,3,4,5,6,7],["2 (1.0 K)", "3 (1.5 K)" , "4 (2.0 K)" , "5 (2.5 K)" , "6 (3.0 K)" , "7 (3.5 K)"    ] )


		axDendro.legend(loc=1)

		fig.tight_layout()



		if onlyLocal:

			plt.savefig("compareParaFluxAlphaOnlyLocal.pdf", bbox_inches='tight')
			plt.savefig("compareParaFluxAlphaOnlyLocal.png", bbox_inches='tight', dpi=300)

		else:
			plt.savefig("compareParaFluxAlpha.pdf", bbox_inches='tight')
			plt.savefig("compareParaFluxAlpha.png", bbox_inches='tight', dpi=300)

			


	def alphaDistribution(self,onlyLocal=False):
		"""

		# draw alph distribution, for each DBSCAN,

		:return:
		"""

		algDendro="Dendrogram"
		tb8Den,tb16Den,label8Den,label16Den,sigmaListDen=self.getTBList(algorithm=algDendro)
		tb8Den=self.removeAllEdges(tb8Den)
		tb16Den=self.removeAllEdges(tb16Den)

		if onlyLocal:
			tb8Den = self.removePerseusList(tb8Den)
			tb16Den = self.removePerseusList(tb16Den)



		algDB="DBSCAN"
		tb8DB,tb16DB,label8DB,label16DB,sigmaListDB=self.getTBList(algorithm=algDB)
		tb8DB=self.removeAllEdges(tb8DB)
		tb16DB=self.removeAllEdges(tb16DB)
		if onlyLocal:
			tb8DB = self.removePerseusList(tb8DB)
			tb16DB = self.removePerseusList(tb16DB)



		algSCI="SCIMES"
		tb8SCI,tb16SCI,label8SCI,label16SCI,sigmaListSCI=self.getTBList(algorithm= algSCI )
		tb8SCI=self.removeAllEdges(tb8SCI)
		tb16SCI=self.removeAllEdges(tb16SCI)

		if onlyLocal:
			tb8SCI = self.removePerseusList(tb8SCI)
			tb16SCI = self.removePerseusList(tb16SCI)


			
		fig=plt.figure(figsize=(8, 6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]


		#############   plot dendrogram
		axDendro=fig.add_subplot(1,1,1)

		alphaDendro , alphaDendroError = self.getAlphaList( tb16Den )

		if onlyLocal:
			axDendro.errorbar(sigmaListDen, alphaDendro, yerr= alphaDendroError , c='g', marker='D', capsize=3, elinewidth=1.3, lw=1, label= algDendro+r" ($|b|>1^{\circ}$)"  )
		else:
			axDendro.errorbar(sigmaListDen, alphaDendro, yerr= alphaDendroError , c='g', marker='D', capsize=3, elinewidth=1.3, lw=1, label= algDendro  )


		axDendro.set_ylabel(r"$\alpha$ (area)")
		axDendro.set_xlabel(r"CO cutoff ($\sigma$)")


		##############   plot DBSCAN

		#alphaDB,  alphaDBError = self.drawAlpha( axDB,tb8DB,tb16DB, label8DB,  label16DB ,sigmaListDB)

		alphaDB , alphaDBError = self.getAlphaList( tb16DB )

		if onlyLocal:
			ebDBSCAN=axDendro.errorbar(sigmaListDB, alphaDB, yerr= alphaDBError , c='b', marker='^',linestyle=":", capsize=3, elinewidth=1.0, lw=1,label= algDB+r" ($|b|>1^{\circ}$)" , markerfacecolor='none' )

		else:
			ebDBSCAN=axDendro.errorbar(sigmaListDB, alphaDB, yerr= alphaDBError , c='b', marker='^',linestyle=":", capsize=3, elinewidth=1.0, lw=1,label= algDB , markerfacecolor='none' )

		ebDBSCAN[-1][0].set_linestyle(':')

		##########plot SCIMES

		alphaSCI , alphaSCIError = self.getAlphaList( tb16SCI )


		if onlyLocal:
			ebSCISCAN=axDendro.errorbar(sigmaListSCI, alphaSCI, yerr= alphaSCIError , c='purple', marker='s', linestyle="--", capsize=3, elinewidth=1.0, lw=1,label= algSCI +r" ($|b|>1^{\circ}$)" ,  markerfacecolor='none' )
		else:
			ebSCISCAN=axDendro.errorbar(sigmaListSCI, alphaSCI, yerr= alphaSCIError , c='purple', marker='s', linestyle="--", capsize=3, elinewidth=1.0, lw=1,label= algSCI ,  markerfacecolor='none' )

		ebSCISCAN[-1][0].set_linestyle('--')
		plt.xticks( [2,3,4,5,6,7],["2 (1.0 K)", "3 (1.5 K)" , "4 (2.0 K)" , "5 (2.5 K)" , "6 (3.0 K)" , "7 (3.5 K)"    ] )

		#plot average
		allAlpha = alphaDB+alphaDendro
		allAlphaError = alphaDBError+alphaDendroError

		print "Average error, ",  np.mean(allAlphaError )

		errorAlpha= np.mean(allAlphaError )**2+  np.std(allAlpha,ddof=1)**2
		errorAlpha=np.sqrt( errorAlpha )

		alphaMean= np.mean( allAlpha)

		print "The mean alpha is {:.2f}, error is {:.2f}".format(alphaMean,errorAlpha )

		labelAverage=r"{:.2f}$\pm${:.2f}".format(alphaMean, errorAlpha  )
		axDendro.plot([min(sigmaListDB),max(sigmaListDB)],  [alphaMean, alphaMean],'--', color='black', lw=1,label=labelAverage )

		axDendro.legend(loc=1 )





		fig.tight_layout()

		if onlyLocal:
			plt.savefig( "compareParaAlphaOnlyLocal.pdf"  ,  bbox_inches='tight')
			plt.savefig( "compareParaAlphaOnlyLocal.png"  ,  bbox_inches='tight',dpi=300)

		else:
			plt.savefig( "compareParaAlpha.pdf"  ,  bbox_inches='tight')
			plt.savefig( "compareParaAlpha.png"  ,  bbox_inches='tight',dpi=300)






	def fluxDistribution(self):
		algDendro="Dendrogram"
		tb8Den,tb16Den,label8Den,label16Den,sigmaListDen=self.getTBList(algorithm=algDendro)
		tb8Den=self.removeAllEdges(tb8Den)
		tb16Den=self.removeAllEdges(tb16Den)


		algDB="DBSCAN"
		tb8DB,tb16DB,label8DB,label16DB,sigmaListDB=self.getTBList(algorithm=algDB)
		tb8DB=self.removeAllEdges(tb8DB)
		tb16DB=self.removeAllEdges(tb16DB)



		algSCI="SCIMES"
		tb8SCI,tb16SCI,label8SCI,label16SCI,sigmaListSCI=self.getTBList(algorithm= algSCI )
		tb8SCI=self.removeAllEdges(tb8SCI)
		tb16SCI=self.removeAllEdges(tb16SCI)







		compoleteFluxa= 9*(1500./250)**2*0.2*1.5 # 2 sigma


		#aaaaaa

		fig=plt.figure(figsize=(22,8))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 15,  'serif': ['Helvetica'] })
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]

		#plot dendrogram
		axDendro=fig.add_subplot(1,3,1)



		at = AnchoredText(algDendro, loc=3, frameon=False)
		axDendro.add_artist(at)

		self.drawFlux(axDendro,tb8Den,tb16Den, label8Den,label16Den, sigmaListDen)

		strXlabel=  r"Flux ($\rm K\ km\ s$$^{-1}$ $\Omega_{\rm A}$)"

		axDendro.set_xlabel( strXlabel )
		axDendro.set_ylabel(r"Number of trunks")

		plt.xticks( [2,3,4,5,6,7],["2 (1.0 K)", "3 (1.5 K)" , "4 (2.0 K)" , "5 (2.5 K)" , "6 (3.0 K)" , "7 (3.5 K)"    ] )

		axDendro.set_yscale('log')
		axDendro.set_xscale('log')
		#axDendro.plot( [compoleteFluxa,compoleteFluxa],[2,800],'--',color='black', lw=1  )



		axDendro.legend(loc=1, ncol=2 )
		#plot DBSCAN
		axDB=fig.add_subplot(1,3,2,sharex=axDendro,sharey=axDendro )

		self.drawFlux(axDB,tb8DB,tb16DB, label8DB,label16DB, sigmaListDB)


		at = AnchoredText(algDB, loc=3, frameon=False)
		axDB.add_artist(at)

		plt.xticks( [2,3,4,5,6,7],["2 (1.0 K)", "3 (1.5 K)" , "4 (2.0 K)" , "5 (2.5 K)" , "6 (3.0 K)" , "7 (3.5 K)"    ] )

		axDB.set_xlabel( strXlabel )
		#axDB.set_ylabel(r"Bin number of trunks ")

		axDB.set_yscale('log')
		axDB.set_xscale('log')

		axDB.legend(loc=1, ncol=2 )

		axDB.set_ylabel(r"Number of trunks")



		plt.xticks( [2,3,4,5,6,7],["2 (1.0 K)", "3 (1.5 K)" , "4 (2.0 K)" , "5 (2.5 K)" , "6 (3.0 K)" , "7 (3.5 K)"    ] )

		#plot SCIMES

		axSCI=fig.add_subplot(1,3,3,sharex=axDendro,sharey=axDendro )

		self.drawFlux(axSCI,tb8SCI,tb16SCI, label8SCI,label16SCI, sigmaListSCI )


		at = AnchoredText(algSCI, loc=3, frameon=False)
		axSCI.add_artist(at)

		axSCI.set_xlabel( strXlabel )
		#axSCI.set_ylabel(r"Bin number of trunks ")

		axSCI.set_yscale('log')
		axSCI.set_xscale('log')

		axSCI.legend(loc=1, ncol=2 )

		axSCI.set_ylabel(r"Number of clusters")



		fig.tight_layout()
		plt.savefig( "compareParaFlux.pdf"  ,  bbox_inches='tight')
		plt.savefig( "compareParaFlux.png"  ,  bbox_inches='tight',dpi=300)




	def areaDistribution(self):

		algDendro="Dendrogram"
		tb8Den,tb16Den,label8Den,label16Den,sigmaListDen=self.getTBList(algorithm=algDendro)
		tb8Den=self.removeAllEdges(tb8Den)
		tb16Den=self.removeAllEdges(tb16Den)


		algDB="DBSCAN"
		tb8DB,tb16DB,label8DB,label16DB,sigmaListDB=self.getTBList(algorithm=algDB)
		tb8DB=self.removeAllEdges(tb8DB)
		tb16DB=self.removeAllEdges(tb16DB)

		algSCI="SCIMES"
		tb8SCI,tb16SCI,label8SCI,label16SCI,sigmaListSCI=self.getTBList(algorithm= algSCI )
		tb8SCI=self.removeAllEdges(tb8SCI)
		tb16SCI=self.removeAllEdges(tb16SCI)



		#aaaaaa

		fig=plt.figure(figsize=(18,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 15,  'serif': ['Helvetica'] })


		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]



		#plot dendrogram
		axDendro=fig.add_subplot(1,3,1)

		#self.drawNumber(axDendro,tb8Den,tb16Den,sigmaListDen)
		at = AnchoredText(algDendro, loc=3, frameon=False)
		axDendro.add_artist(at)

		self.drawArea(axDendro,tb8Den,tb16Den, label8Den,label16Den, sigmaListDen)

		axDendro.set_xlabel(r"Angular area (deg$^2$)")
		axDendro.set_ylabel(r"Number of trunks")


		axDendro.set_yscale('log')
		axDendro.set_xscale('log')

		compoleteArea= 9*(1500./250)**2*0.25/3600. #0.0225


		axDendro.plot( [compoleteArea,compoleteArea],[2,2000],'--',color='black', lw=1  )

		axDendro.legend(loc=1, ncol=2 )
		#plot DBSCAN
		axDB=fig.add_subplot(1,3,2,sharex=axDendro,sharey=axDendro )
		#self.drawNumber(axDB,tb8DB,tb16DB,sigmaListDB)
		self.drawArea(axDB,tb8DB,tb16DB, label8DB,label16DB, sigmaListDB)


		at = AnchoredText(algDB, loc=3, frameon=False)
		axDB.add_artist(at)
		axDB.plot( [compoleteArea,compoleteArea],[2,2000],'--',color='black', lw=1  )
		axDB.set_xlabel(r"Angular area (deg$^2$)")
		#axDB.set_ylabel(r"Bin number of trunks ")
		axDB.set_ylabel(r"Number of trunks")

		axDB.set_yscale('log')
		axDB.set_xscale('log')

		axDB.legend(loc=1, ncol=2 )

		#plot SCIMES

		axSCI=fig.add_subplot(1,3,3,sharex=axDendro,sharey=axDendro )
		#self.drawNumber(axDB,tb8DB,tb16DB,sigmaListDB)
		self.drawArea(axSCI,tb8SCI,tb16SCI, label8SCI,label16SCI, sigmaListSCI )


		at = AnchoredText(algSCI, loc=3, frameon=False)
		axSCI.add_artist(at)
		axSCI.plot( [compoleteArea,compoleteArea],[2,2000],'--',color='black', lw=1  )
		axSCI.set_xlabel(r"Angular area (deg$^2$)")
		#axSCI.set_ylabel(r"Bin number of trunks ")
		axSCI.set_ylabel(r"Number of clusters")





		axSCI.set_yscale('log')
		axSCI.set_xscale('log')

		axSCI.legend(loc=1, ncol=2 )

		########




		fig.tight_layout()
		plt.savefig( "compareParaArea.pdf"  ,  bbox_inches='tight')
		plt.savefig( "compareParaArea.png"  ,  bbox_inches='tight',dpi=300)




	def totaFluxDistribution(self):
		"""
		Compare the change of molecular cloud numbers with
		:return:
		"""
		algDendro="Dendrogram"
		tb8Den,tb16Den,label8Den,label16Den,sigmaListDen=self.getTBList(algorithm=algDendro)
		tb8Den=self.removeAllEdges(tb8Den)
		tb16Den=self.removeAllEdges(tb16Den)


		algDB="DBSCAN"
		tb8DB,tb16DB,label8DB,label16DB,sigmaListDB=self.getTBList(algorithm=algDB)
		tb8DB=self.removeAllEdges(tb8DB)
		tb16DB=self.removeAllEdges(tb16DB)


		algSCI="SCIMES"
		tb8SCI,tb16SCI,label8SCI,label16SCI,sigmaListSCI=self.getTBList(algorithm= algSCI )
		tb8SCI=self.removeAllEdges(tb8SCI)
		tb16SCI=self.removeAllEdges(tb16SCI)



		fig=plt.figure(figsize=(8,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 10.5,  'serif': ['Helvetica'] })
		rc('text.latex',  preamble=r'\usepackage{upgreek}')

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]



		#plot dendrogram
		axDendro=fig.add_subplot(1,1,1)

		#self.drawTotalFlux(axDendro,tb8Den,tb16Den, label8Den,label16Den, sigmaListDen)

		#Nlist8Den=self.getTotalFluxList(tb8List)
		fluxList16Den=self.getTotalFluxList(tb8Den)
		axDendro.plot(sigmaListDen,fluxList16Den,'D-' , color="green",markersize=4, lw=1.0,label="Flux (dendrogram)",markerfacecolor='none' )


		axDendro.set_ylabel(r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_\mathrm{A}$)")
		axDendro.set_xlabel(r"CO cutoff ($\sigma$)")



		#plot DBSCAN
		fluxList16DB=self.getTotalFluxList( tb8DB )

		axDendro.plot(sigmaListDB,fluxList16DB,'^--' ,linestyle=':', color="blue",markersize=3, lw=1.0,label= "Flux (DBSCAN)", markerfacecolor='none' )


		#plot SCIMES
		fluxList16SCI=self.getTotalFluxList( tb8SCI )

		axDendro.plot(sigmaListSCI,fluxList16SCI,'s--' ,  color="purple",markersize=3, lw=1.0,label= "Flux (SCIMES)", markerfacecolor='none' )



		#drawTotal Flux

		#axDendro.set_yscale('log')
		#maskCO,_=myFITS.readFITS("G2650CO12MaskedCO.fits")
		#use dendro masked

		maskCO,_=myFITS.readFITS("G2650CO12DendroMaskedCO.fits")



		totalFluxList=[]

		omega = 1 #np.deg2rad(0.5 / 60) * np.deg2rad(0.5 / 60) / 4. / np.log(2)

		for eachS in sigmaListDen:
			maskCO[maskCO<eachS*self.rms]=0
			sumCO= np.sum( maskCO)*0.2* omega
			totalFluxList.append( sumCO )

		print totalFluxList[0], "maximum flux"

		totalFluxRatioDendro=np.asarray( fluxList16Den ) /np.asarray( totalFluxList  )
		totalFluxRatioDB=np.asarray( fluxList16DB ) /np.asarray( totalFluxList  )

		totalFluxRatioSCI=np.asarray( fluxList16SCI ) /np.asarray( totalFluxList  )

		tabBlue='tab:blue'
		axDendro.plot(sigmaListDB,totalFluxList,'o-' ,   color="black",markersize=3, lw=1.0,label= "Flux (total flux)", markerfacecolor='none' )

		z=np.polyfit(sigmaListDB, totalFluxList ,1 )

		p=np.poly1d(  z )
		print z
		#print p(0.0)/p(3), p(0.0)/p(2)
		print p(3)/p(0), p(2)/p(0) #to find the ratio to the infinite sensitivity

		#axDendro.plot(sigmaListDB, p( sigmaListDB )  ,'r-' ,   color="red",markersize=3, lw=1.0,label= "Total flux(cutoff mask)", markerfacecolor='none' )


		axRatio = axDendro.twinx()  # instantiate a second axes that shares the same x-axis

		axRatio.plot(sigmaListDen, totalFluxRatioDendro,'D-',color=tabBlue ,  markersize=3, lw=1.0, label= "Ratio to the total flux (dendrogram)" , markerfacecolor='none' )

		axRatio.plot(sigmaListDB, totalFluxRatioDB,'^--',linestyle=':',color= tabBlue ,  markersize=3, lw=1.0,label= "Ratio to  the total flux (DBSCAN)",  markerfacecolor='none' )

		axRatio.plot(sigmaListSCI, totalFluxRatioSCI,'s--', color= tabBlue ,  markersize=3, lw=1.0,label= "Ratio to  the total flux (SCIMES)",  markerfacecolor='none' )


		axRatio.set_ylabel('Ratio to  the total flux', color= tabBlue )

		#draw



		axDendro.legend(loc=6)
		plt.xticks( [2,3,4,5,6,7],["2 (1.0 K)", "3 (1.5 K)" , "4 (2.0 K)" , "5 (2.5 K)" , "6 (3.0 K)" , "7 (3.5 K)"    ] )

		leg=axRatio.legend(loc=7)
		for text in leg.get_texts():
			plt.setp(text, color=tabBlue)

		axRatio.set_ylim([0.4, 1.06])


		axRatio.yaxis.label.set_color( tabBlue )
		axRatio.spines["right"].set_edgecolor( tabBlue )
		axRatio.tick_params(axis='y', colors= tabBlue )


		fig.tight_layout()
		plt.savefig( "compareParaTotalFlux.pdf"  ,  bbox_inches='tight')
		plt.savefig( "compareParaTotalFlux.png"  ,  bbox_inches='tight',dpi=300)



	def printCatNumbers(self):
		"""
		print number of catalog
		:return:
		"""
		algDendro="Dendrogram"
		tb8Den,tb16Den,label8Den,label16Den,sigmaListDen=self.getTBList(algorithm=algDendro)
		tb8Den=self.removeAllEdges(tb8Den)
		tb16Den=self.removeAllEdges(tb16Den)



		algDB="DBSCAN"
		tb8DB,tb16DB,label8DB,label16DB,sigmaListDB=self.getTBList(algorithm=algDB)
		tb8DB=self.removeAllEdges(tb8DB)
		tb16DB=self.removeAllEdges(tb16DB)



		algSCI="SCIMES"
		tb8SCI,tb16SCI,label8SCI,label16SCI,sigmaListSCI=self.getTBList(algorithm= algSCI )
		tb8SCI=self.removeAllEdges(tb8SCI)
		tb16SCI=self.removeAllEdges(tb16SCI)

		######

		#for alg in [algDendro,algDB,algSCI]:


		############## DENDROGRAM
		p8Head = "{:^11s} & {:^2d} & ".format(  algDendro,   8 )
		print    p8Head,

		totalN=len(tb8Den)

		for i in range(totalN):
				if i<totalN-1:
					print len(tb8Den[i]), "&",
				else:
					print len(tb8Den[i]), "\\\\",

		print
		p16Head= "{:^11s} & {:^2d} & ".format(  algDendro,  16 )

		print p16Head,

		for i in range(totalN):
			if i < totalN - 1:
				print len(tb16Den[i]), "&",
			else:
				print len(tb16Den[i]), "\\\\",
		print


		############## DBSCAN
		p8Head = "{:^11s} & {:^2d} & ".format(  algDB,   8 )
		print    p8Head,

		totalN=len(tb8Den)

		for i in range(totalN):
				if i<totalN-1:
					print len(tb8DB[i]), "&",
				else:
					print len(tb8DB[i]), "\\\\",

		print
		p16Head= "{:^11s} & {:^2d} & ".format(  algDB,  16 )

		print p16Head,

		for i in range(totalN):
			if i < totalN - 1:
				print len(tb16DB[i]), "&",
			else:
				print len(tb16DB[i]), "\\\\",



		print


		############## SCIMES
		p8Head = "{:^11s} & {:^2d} & ".format(  algSCI,   8 )
		print    p8Head,

		totalN=len(tb8Den)

		for i in range(totalN):
				if i<totalN-1:
					print len(tb8SCI[i]), "&",
				else:
					print len(tb8SCI[i]), "\\\\",

		print
		p16Head= "{:^11s} & {:^2d} & ".format(  algSCI,  16 )

		print p16Head,

		for i in range(totalN):
			if i < totalN - 1:
				print len(tb16SCI[i]), "&",
			else:
				print len(tb16SCI[i])







	def numberDistribution(self):
		"""
		Compare the change of molecular cloud numbers with
		:return:
		"""
		algDendro="Dendrogram"
		tb8Den,tb16Den,label8Den,label16Den,sigmaListDen=self.getTBList(algorithm=algDendro)
		tb8Den=self.removeAllEdges(tb8Den)
		tb16Den=self.removeAllEdges(tb16Den)


		algDB="DBSCAN"
		tb8DB,tb16DB,label8DB,label16DB,sigmaListDB=self.getTBList(algorithm=algDB)
		tb8DB=self.removeAllEdges(tb8DB)
		tb16DB=self.removeAllEdges(tb16DB)



		algSCI="SCIMES"
		tb8SCI,tb16SCI,label8SCI,label16SCI,sigmaListSCI=self.getTBList(algorithm= algSCI )
		tb8SCI=self.removeAllEdges(tb8SCI)
		tb16SCI=self.removeAllEdges(tb16SCI)





		fig=plt.figure(figsize=(20,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 16,  'serif': ['Helvetica'] })
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]
		#plot dendrogram
		axDendro=fig.add_subplot(1,3,1)

		self.drawNumber(axDendro,tb8Den,tb16Den, label8Den,label16Den, sigmaListDen)
		at = AnchoredText(algDendro, loc=1, frameon=False)
		axDendro.add_artist(at)
		axDendro.set_ylabel(r"Total number of trunks")
		axDendro.set_xlabel(r"CO cutoff ($\sigma$)")
		axDendro.legend(loc=2)
		plt.xticks( [2,3,4,5,6,7],["2 (1.0 K)", "3 (1.5 K)" , "4 (2.0 K)" , "5 (2.5 K)" , "6 (3.0 K)" , "7 (3.5 K)"    ] )


		#plot DBSCAN
		axDB=fig.add_subplot(1,3,2,sharex=axDendro,sharey=axDendro)
		self.drawNumber(axDB,tb8DB,tb16DB, label8DB,  label16DB ,sigmaListDB)
		at = AnchoredText(algDB, loc=1, frameon=False)
		axDB.add_artist(at)

		#axDB.set_ylabel(r"Total number of clusters")

		axDB.set_xlabel(r"CO cutoff ($\sigma$)")
		axDB.set_ylabel(r"Total number of trunks")

		axDB.legend(loc=2)
		plt.xticks( [2,3,4,5,6,7],["2 (1.0 K)", "3 (1.5 K)" , "4 (2.0 K)" , "5 (2.5 K)" , "6 (3.0 K)" , "7 (3.5 K)"    ] )

		#plot SCIMES

		axSCI=fig.add_subplot(1,3,3,sharex=axDendro,sharey=axDendro)
		self.drawNumber(axSCI,tb8SCI,tb16SCI, label8SCI,  label16SCI ,sigmaListSCI)
		at = AnchoredText(algSCI, loc=1, frameon=False)
		axSCI.add_artist(at)
		plt.xticks( [2,3,4,5,6,7],["2 (1.0 K)", "3 (1.5 K)" , "4 (2.0 K)" , "5 (2.5 K)" , "6 (3.0 K)" , "7 (3.5 K)"    ] )

		#axSCI.set_ylabel(r"Total number of clusters")

		axSCI.set_xlabel(r"CO cutoff ($\sigma$)")
		axSCI.set_ylabel(r"Total number of clusters")

		axSCI.legend(loc=3)


		fig.tight_layout()
		plt.savefig( "compareParaNumber.pdf"  ,  bbox_inches='tight')
		plt.savefig( "compareParaNumber.png"  ,  bbox_inches='tight',dpi=300)


	def drawFlux(self,ax,tb8List,tb16List, label8,label16, sigmaListDen ):

		#areaEdges=np.linspace(0.25/3600.,150,10000)
		#areaCenter=self.getEdgeCenter( areaEdges )

		areaEdges=np.linspace(8,1e5,1000)
		areaCenter=self.getEdgeCenter( areaEdges )

		NUM_COLORS = 12
		clrs = sns.color_palette('husl', n_colors=NUM_COLORS)  # a list of RGB tuples

		for i in range( len(tb8List) ):

			eachTB8 = tb8List[i]
			eachTB16 = tb16List[i]

			if "sum" not in  eachTB8.colnames:
				#dendrogra
				sum8=eachTB8["flux"]/self.getSumToFluxFactor()*0.2 # K km/s
				sum16=eachTB16["flux"]/self.getSumToFluxFactor()*0.2 # K km/s

			else: #dbscan

				sum8 = eachTB8["sum"]*0.2 # K km/s
				sum16 = eachTB16["sum"]*0.2 # K km/s

			print i,"P8, minFlux:{:.2f}, P16, maxFlux:{:.2f}".format( np.min(sum8), np.max(sum16)  )
			print i, "P16, minFlux:{:.2f}, P16, maxFlux:{:.2f}".format( np.min(sum16), np.max(sum16)  )

			binN8,binEdges8 = np.histogram( sum8 , bins=areaEdges  )
			binN16,binEdges16 = np.histogram( sum16 , bins=areaEdges  )

			plot8 = ax.plot( areaCenter[binN8>0],binN8[binN8>0], 'o-'  , markersize=1, lw=0.8,label=label8[i] ,alpha= 0.5,color=  clrs[i])

			ax.plot( areaCenter[binN16>0],binN16[binN16>0], '^--'  , markersize=1, lw=0.8,label=label16[i] ,alpha= 0.5,color= clrs[i] )

			completeFlux=324*self.rms*0.2*sigmaListDen[i]*3    # K km/s, the last 2 is the two channels
			print "Complete, ",completeFlux,self.rms,sigmaListDen[i]
			ax.plot( [completeFlux,completeFlux], [2,3000]   ,'--', markersize=1, lw=0.8  ,alpha= 0.5,color= clrs[i] )




	def drawArea(self,ax,tb8List,tb16List, label8,label16, sigmaListDen ):


		areaEdges=np.linspace(0.25/3600.,150,10000)
		areaCenter=self.getEdgeCenter( areaEdges )


		NUM_COLORS = 12
		clrs = sns.color_palette('husl', n_colors=NUM_COLORS)  # a list of RGB tuples


		for i in range( len(tb8List) ):

			eachTB8 = tb8List[i]
			eachTB16 = tb16List[i]


			print  i, "8 pix: MinArea:{:.2f}, MaxArea:{:.2f}".format(   np.min(eachTB8["area_exact"]), np.max(eachTB8["area_exact"]/3600.) )

			binN8 , binEdges8 =np.histogram(eachTB8["area_exact"]/3600., bins=areaEdges  )
			binN16 , binEdges16 =np.histogram(eachTB16["area_exact"]/3600., bins=areaEdges  )

			print  i, "16 pix: MinArea:{:.2f}, MaxArea:{:.2f}".format(   np.min(eachTB16["area_exact"]), np.max(eachTB16["area_exact"]/3600.) )

			ax.plot( areaCenter[binN8>0],binN8[binN8>0], 'o-'  , markersize=1, lw=0.8,label=label8[i] ,alpha= 0.5,color=clrs[i] )
			ax.plot( areaCenter[binN16>0],binN16[binN16>0], '^--'   , markersize=1, lw=0.8,label=label16[i] ,alpha= 0.5,color=clrs[i] )


	def drawPhysicalAreaSingle(self, ax, tbDendro2, physicalEdges, physicalCenter,completeArea,  label=None ):

		# physicalEdges  physicalCenter
		#areaEdges=np.linspace(0,100,1000)
		#areaCenter=self.getEdgeCenter( areaEdges )

		realArea=self.getRealArea(tbDendro2)
		binN,binEdges=np.histogram( realArea , bins=physicalEdges )

		#calculate alpha

		meanA, stdA = self.getAlphaWithMCMC(realArea, minArea=completeArea, maxArea=None, physicalArea=True)

		if label!=None:
			label=label+r": $\alpha={:.2f}\pm{:.2f}$".format(meanA,stdA)

		stepa=ax.plot(physicalCenter[binN > 0], binN[binN > 0], 'o-', markersize=1, lw=0.8, label=label, alpha=0.5 )



		return stepa[-1].get_color()
	
	
	def drawNumber(self,ax,tb8List,tb16List,label8,label16,sigmaListDen ):
		Nlist8Den=self.getNList(tb8List)
		Nlist16Den=self.getNList(tb16List)



		ax.plot(sigmaListDen,Nlist8Den,'o-',label="min\_npix = 8",color="blue",lw=1 , markersize=3 )
		ax.plot(sigmaListDen,Nlist16Den,'o-',label="min\_npix = 16",color="green", lw=1, markersize=3 )



	def drawTotalFlux(self,ax,tb8List,tb16List,label8,label16,sigmaListDen ):
		#Nlist8Den=self.getTotalFluxList(tb8List)
		Nlist16Den=self.getTotalFluxList(tb16List)


		#ax.plot(sigmaListDen,Nlist8Den,'o-',label="min\_nPix = 8",color="blue",lw=0.5)
		ax.plot(sigmaListDen,Nlist16Den,'o-' , color="green", lw=0.5)




	def getAlphaList(self,tbList, minArea=0.01,colName="area_exact" ):
		# calculate alpha and  error for each alpha for each tb

		alphaList=[]
		errorList=[]

		for eachTB in tbList:

			area= eachTB[ colName ]

			meanA,stdA=self.getAlphaWithMCMC(  area ,  minArea= minArea ,  maxArea=None , physicalArea=False )

			alphaList.append(meanA)
			errorList.append( stdA)

		return  alphaList,  errorList


	def getFluxCol(self,TB):

		if "sum" in TB.colnames:
			return TB["sum"]*0.2 #K km/s


		else:

			return TB["flux"]*0.2/self.getSumToFluxFactor()



	def getFluxAlphaList(self,tbList, sigmaList  ):
		# calculate alpha and  error for each alpha for each tb

		alphaList=[]
		errorList=[]
		completeList=[]
		for i in range( len(sigmaList) ):
			eachTB = tbList[i]

			eachSigma = sigmaList[i]

			flux= self.getFluxCol(eachTB  )
			#minFlux=324*self.rms*0.2*eachSigma*3    # K km/s, the last 2 is the two channels
			minFlux=144*self.rms*0.2*eachSigma*3    # K km/s, the last 2 is the two channels

			meanA,stdA=self.getAlphaWithMCMC(  flux,  minArea= minFlux ,  maxArea=None , physicalArea=True )

			alphaList.append(meanA)
			errorList.append( stdA)
			completeList.append(minFlux)

		return  alphaList,  errorList,

	###########################################################
	def getMassAlphaList(self,tbList, sigmaList  ):
		# calculate alpha and  error for each alpha for each tb

		alphaList=[]
		errorList=[]
		completeList=[]
		for i in range( len(sigmaList) ):
			eachTB = tbList[i]

			eachSigma = sigmaList[i]

			massArray = self.getMass(eachTB  )
			#minFlux=324*self.rms*0.2*eachSigma*3    # K km/s, the last 2 is the two channels
			minFlux=144*self.rms*0.2*eachSigma*3    # K km/s, the last 2 is the two channels
			minMass=self.calmassByXfactor(minFlux,1500)

			print minMass,"The complete mass"

			meanA,stdA=self.getAlphaWithMCMC(  massArray,  minArea= minMass ,  maxArea=None , physicalArea=True )

			alphaList.append(meanA)
			errorList.append( stdA)
			completeList.append(minMass)

		return  alphaList,  errorList, completeList






	def drawAlpha(self,ax,tb8List,tb16List, label8, label16,sigmaListDen ): #


		#fitting alpha and draw

		alpha8List, alpha8ErrorList = self.getAlphaList(tb16List)
		#alpha16List, alpha16ErrorList = self.getAlphaList(tb16List)

		#ax.plot(sigmaListDen,alpha8List,'o-',label="MinPix = 8",color="blue", markersize= 3, lw=1)
		#ax.plot(sigmaListDen,alpha16List,'o-',label="MinPix = 16",color="green", lw=0.5, markersize= 2.5  ,  alpha=0.8 )

		ax.errorbar(sigmaListDen, alpha8List, yerr= alpha8ErrorList , c='b', marker='o', capsize=1.5, elinewidth=0.8, lw=1,label=r"min\_nPix = 16" )

		return alpha8List,alpha8ErrorList

	def areaAndNumberDistribution(self, algorithm="Dendrogram" ):
		"""
		#draw the area the
		:return:
		"""




		#first, get TBList

		tb8,tb16,label8,label16,sigmaList=self.getTBList(algorithm=algorithm)



		#need to constrain the minP and PeakN, PeakSigma=lower sigma cut + 3 sigma, as the minDelta,



		fig=plt.figure(figsize=(12,6))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })

		axNumber=fig.add_subplot(1,2,1)

		Nlist8=self.getNList(tb8)
		Nlist16=self.getNList(tb16)


		axNumber.plot(sigmaList,Nlist8,'o-',label="MinPix = 8",color="blue",lw=0.5)
		axNumber.plot(sigmaList,Nlist16,'o-',label="MinPix = 16",color="green", lw=0.5)





		#axArea.set_xlabel(r"Area (deg$^2$)")
		axNumber.set_ylabel(r"Total number of trunks")
		axNumber.set_xlabel(r"CO cutoff ($\sigma$)")

		axNumber.legend()

		at = AnchoredText(algorithm, loc=3, frameon=False)
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axNumber.add_artist(at)


		################ Area ##############


		axArea= fig.add_subplot(1,2,2)



		areaEdges=np.linspace(0,150,10000)
		areaCenter=self.getEdgeCenter( areaEdges )


		totalTB=tb8+tb16
		labelStr=label8+label16

		for i in range( len(totalTB) ):

			eachTB = totalTB[i]

			binN,binEdges=np.histogram(eachTB["area_exact"]/3600., bins=areaEdges  )


			axArea.plot( areaCenter[binN>0],binN[binN>0], 'o-'  , markersize=1, lw=0.8,label=labelStr[i] ,alpha= 0.5 )



		#set ticikes of Area


		axArea.set_yscale('log')
		axArea.set_xscale('log')

		axArea.set_xlim( [ 0.005,150 ] )


		if algorithm=="DBSCAN":
			axArea.set_ylim( [ 0.8,50000 ] )

		else:
			axArea.set_ylim( [ 0.8,10000 ] )





		axArea.set_xlabel(r"Area (deg$^2$)")
		axArea.set_ylabel(r"Bin number of trunks ")


		axArea.legend(ncol=2)


		#draw scimes

		scimesTB=Table.read("ClusterCat_3_16Ve20.fit")

		binN,binEdges=np.histogram(scimesTB["area_exact"]/3600., bins=areaEdges  )


		axArea.plot( areaCenter[binN>0],binN[binN>0], 'o-'  ,  color='red',  markersize=1, lw=0.8,label=r"Scimes,3.0$\sigma$, P16" ,alpha= 0.5 )

		at = AnchoredText("Red: SCIMES,3.0$\sigma$, P16", loc=4, frameon=False)
		#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
		axArea.add_artist(at)






		plt.savefig( "comparePara_{}.pdf".format(algorithm) ,  bbox_inches='tight')

		plt.savefig( "comparePara_{}.png".format(algorithm) ,  bbox_inches='tight',dpi=300)



	def getTBList(self, algorithm="DBSCAN"):
		"""
		return a list of table,
		:param minP:
		:param algorithm:
		:return:
		"""

		if algorithm=="DBSCAN":
			#ignore minP, only has 8


			TBList=[]
			TBList16=[]

			TBLabelsP8=[]
			TBLabelsP16=[]
			minPix=8

			#DbscanSigmaList= np.arange(2,6.5,0.5)
			DbscanSigmaList= np.arange(2,7.5,0.5)

			for sigmas in DbscanSigmaList:
				tbName= "G2650CO12DBCatS{:.1f}P{}Con2.fit".format(sigmas, minPix)
				ttt8=Table.read(tbName)
				ttt8=self.cleanDBTB(ttt8,pixN=8,minV=sigmas,minDelta=3)
				
				TBList.append(ttt8  )
				ttt16=ttt8[ttt8["pixN"]>=16]
				ttt16=self.cleanDBTB(ttt16,pixN=16,minV=sigmas,minDelta=3)


				TBList16.append(ttt16  )
				TBLabelsP8.append(  r"{:.1f}$\sigma$, P8".format( sigmas)   )
				TBLabelsP16.append( r"{:.1f}$\sigma$, P16".format( sigmas)   )

			


			return TBList,TBList16,TBLabelsP8,TBLabelsP16,DbscanSigmaList


		elif algorithm=="dendrogram" or algorithm=="Dendrogram" :

			TBListP8=[]
			TBListP16=[]

			TBLabelsP8=[]
			TBLabelsP16=[]

			#dendroSigmaList=[2,2.5 , 3, 3.5, 4,4.5,5, 5.5, 6]

			dendroSigmaList=[2,2.5 , 3, 3.5, 4,4.5,5, 5.5, 6,6.5,7 ]

			for sigmas in dendroSigmaList:
				tbName8= "minV{}minP{}_dendroCatTrunk.fit".format(sigmas, 8)
				tbName16= "minV{}minP{}_dendroCatTrunk.fit".format(sigmas, 16)

				TBListP8.append(Table.read(tbName8)  )
				TBListP16.append(Table.read(tbName16)  )


				TBLabelsP8.append(  r"{:.1f}$\sigma$, P8".format( sigmas)   )
				TBLabelsP16.append( r"{:.1f}$\sigma$, P16".format( sigmas)   )



			return TBListP8,TBListP16,TBLabelsP8,TBLabelsP16,dendroSigmaList

		elif algorithm=="SCIMES":

			TBListP8=[]
			TBListP16=[]

			TBLabelsP8=[]
			TBLabelsP16=[]

			#dendroSigmaList=[2,2.5 , 3, 3.5, 4,4.5,5, 5.5, 6,6.5,7]

			dendroSigmaList=[2,2.5 , 3, 3.5, 4,4.5,5, 5.5,6 ,6.5,7]
			path="./scimesG2650/"
			for sigmas in dendroSigmaList:
				tbName8= path+"ClusterCat_{}_{}Ve20.fit".format(sigmas, 8)
				tbName16= path+"ClusterCat_{}_{}Ve20.fit".format(sigmas, 16)

				TBListP8.append(Table.read(tbName8)  )
				TBListP16.append(Table.read(tbName16)  )


				TBLabelsP8.append(  r"{:.1f}$\sigma$, P8".format( sigmas)   )
				TBLabelsP16.append( r"{:.1f}$\sigma$, P16".format( sigmas)   )



			return TBListP8,TBListP16,TBLabelsP8,TBLabelsP16,dendroSigmaList



	def setMinVandPeak(self,cloudLabelFITS,COFITS, peakSigma=3,minP=8):
		"""
		:param cloudLabelFITS:
		:param peakSigma:
		:return:
		"""
		#reject cloud that peak values area less than peakSigma, and total peak number are less then minP



		dataCluster, head= myFITS.readFITS(cloudLabelFITS)

		dataCO,headCO= myFITS.readFITS(COFITS)

		clusterIndex1D= np.where( dataCluster>0)
		clusterValue1D=  dataCluster[clusterIndex1D ]

		Z0,Y0,X0 = clusterIndex1D

		allClouds,counts=np.unique( clusterValue1D , return_counts=True)


		widgets = ['Fast Dendro WithDBSCAN: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),  ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options

		pbar = ProgressBar(widgets=widgets, maxval=len(allClouds))
		pbar.start()


		for i in range( len(allClouds) ):


			pbar.update(i)

			cloudID =  allClouds[i]

			pixN=counts[i]
			cloudIndex=self.getIndices(Z0,Y0,X0,clusterValue1D,cloudID)

			if pixN<minP:#reject clouds

				dataCluster[cloudIndex]=0

				continue




			coValues=  dataCO[ cloudIndex ]

			if np.nanmax(coValues) <  peakSigma*self.rms: #reject
				cloudIndex=self.getIndices(Z0,Y0,X0,clusterValue1D,cloudID)
				dataCluster[cloudIndex]=0

				continue
		#relabelclouds
		pbar.finish()
		#dataCluster[dataCluster>0]=1
		s=generate_binary_structure(3,1)

		newDataCluster= dataCluster>0



		labeled_redo, num_features=label(newDataCluster, structure=s) #first label core, then expand, otherwise, the expanding would wrongly connected

		print "Total number of clouds? ",  num_features

		tbDendro=Table.read( "minV5minP8_dendroCatTrunk.fit" )
		print "The dendrogramN is ",len(tbDendro)

		#save the fits

		fits.writeto("relabelFastDendrominPeak{}_P{}.fits".format( peakSigma, minP ), labeled_redo, header=headCO,overwrite=True)


	def fastDendro(self,COFITS,minDelta=3,minV=3,minP=8):

		COData,COHead=myFITS.readFITS( COFITS)

		print np.max(COData)
		#first create dendrogram
		self.computeDBSCAN(COData,COHead, min_sigma=minV,min_pix=3,connectivity=1,region="fastDendroTest")

		dbFITS="fastDendroTestdbscanS{}P3Con1.fits".format(minV)



		self.setMinVandPeak(dbFITS,COFITS, peakSigma=minDelta+minV,minP=minP)

	def cleanLabelFITS(self,dataCluster,TB):
		"""

		only keep clusters in TB file

		:param data:
		:param TB:
		:return:
		"""

		noiseLabel=np.min( dataCluster[0] )
		clusterIndex1D= np.where( dataCluster>0 )
		clusterValue1D=  dataCluster[clusterIndex1D ]
		Z0,Y0,X0 = clusterIndex1D

		#

		widgets = ['Cleaning label fits:', Percentage(), ' ', Bar(marker='0',left='[',right=']'),  ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options

		pbar = ProgressBar(widgets=widgets, maxval=len(TB))
		pbar.start()

		indexRun=0

		returnCluster = np.zeros_like(dataCluster)+noiseLabel

		for eachDBRow in TB:
			indexRun=indexRun+1
			pbar.update(indexRun)
			cloudID=  eachDBRow["_idx"]


			cloudIndex = self.getIndices(Z0, Y0, X0, clusterValue1D, cloudID)
			returnCluster[cloudIndex] = cloudID # remove this cluster and do not record this cluster



		pbar.finish()

		return returnCluster


	def clearnDBAssign(self,DBLabelFITS,DBTableFile	,pixN=8,minDelta=3,minV=2 , MinPts=8 ,  minAreaPix=9, prefix="" ):

		minPeak=(minV+minDelta)*self.rms
		saveName=prefix+"DBCLEAN{}_{}Label.fits".format( minV, pixN )
		saveNameTB=prefix+"DBCLEAN{}_{}TB.fit".format( minV, pixN )

		DBTable=Table.read( DBTableFile )




		dataCluster,headCluster=myFITS.readFITS(DBLabelFITS )


		noiseLabel=np.min( dataCluster[0] )


		clusterIndex1D= np.where( dataCluster>0 )
		clusterValue1D=  dataCluster[clusterIndex1D ]

		Z0,Y0,X0 = clusterIndex1D
		#cloudIndex = self.getIndices(Z0, Y0, X0, clusterValue1D, newID)

		emptyTB= Table( DBTable[0] )
		emptyTB.remove_row(0)
		print "Cleaning DBSCAN table..."

		widgets = ['Recalculating cloud parameters: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),  ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options

		pbar = ProgressBar(widgets=widgets, maxval=len(DBTable))
		pbar.start()

		indexRun=0
		for eachDBRow in DBTable:
			indexRun=indexRun+1
			pbar.update(indexRun)
			cloudID=  eachDBRow["_idx"]

			pixNCloud=int(  eachDBRow["pixN"]  )
			peakCloud= eachDBRow["peak"]

			area =  eachDBRow["area_exact"] # by fefault, the area should be in arcmin2

			if peakCloud < minPeak or pixNCloud< pixN or  area< minAreaPix*self.pixelArea : #

				cloudIndex = self.getIndices(Z0, Y0, X0, clusterValue1D, cloudID)
				dataCluster[cloudIndex] = noiseLabel # remove this cluster and do not record this cluster

				continue



			emptyTB.add_row( eachDBRow  )
		pbar.finish()
		#save
		fits.writeto(saveName,dataCluster,header=headCluster,overwrite=True)

		emptyTB.write( saveNameTB,overwrite=True  )

		return saveName, saveNameTB


	def getCropDataAndHead(self,rawFITSFile,drawChannel,lRange,bRange):

		tmpFITS="checkCloudTmp.fits"

		tmpData, tmpHead = myFITS.readFITS(rawFITSFile)
		save2D=tmpData[drawChannel]

		fits.writeto(tmpFITS, save2D,header=tmpHead ,overwrite=True)


		cropTMP = "checkCloudTmpCrop.fits"

		doFITS.cropFITS2D(tmpFITS,cropTMP, Lrange=lRange, Brange=bRange , overWrite=True   )



		return myFITS.readFITS(cropTMP)


	def drawCloudMap(self,drawChannel=98,lRange=None,bRange=None ):
		"""
		#draw small clouds to check if the are real...

		one color for DBSCAN
		one color for dendrogram,

		draw 2sigma, because they would provide the smallest area of clouds,

		:return:
		"""

		xRange=[]
		yRange=[]

		#first

		#axCO.set_xlim( [1000, 1700] )
		#axCO.set_ylim( [350, 950 ] )


		COFITS="G2650Local30.fits"

		#data,head=myFITS.readFITS(COFITS)

		data,head=self.getCropDataAndHead(   COFITS,drawChannel=drawChannel,lRange=lRange,bRange= bRange )


		WCSCO=WCS(head)

		channelRawCO=data #[drawChannel]

		DBLabelFITS = "DBCLEAN2.0_8Label.fits"
		DBTableFile= "DBCLEAN2.0_8TB.fit"
		drawDBSCANtb=Table.read( DBTableFile )


		#relabelDB,newDBTable= self.clearnDBAssign( DBLabelFITS,DBTableFile	,pixN=16, minDelta=3, minV=2  )


		drawDBSCANtb=self.cleanDBTB( drawDBSCANtb, minDelta=3,minV=2,pixN=8)


		drawDBSCANtb=self.removeWrongEdges(drawDBSCANtb)



		drawDBSCANData,drawDBSCANHead=self.getCropDataAndHead(   DBLabelFITS,drawChannel=drawChannel,lRange=lRange,bRange= bRange )
		#drawDBSCANData,drawDBSCANHead = myFITS.readFITS(DBLabelFITS)
		WCSCrop = WCS( drawDBSCANHead )

		channelDBSCAN =drawDBSCANData  #drawDBSCANData[drawChannel]



		dbClouds=np.unique(channelDBSCAN)



		drawDENDROtb=Table.read("minV2minP8_dendroCatTrunk.fit")

		drawDENDROtb=self.removeWrongEdges(drawDENDROtb)


		#drawDENDROData,drawDENDROHead = myFITS.readFITS("minV2minP16_TrunkAsign.fits")
		#drawDENDROData=drawDENDROData-1
		#channelDENDRO =  drawDENDROData[drawChannel]


		drawDendrogram, drawDendrogramHead = self.getCropDataAndHead(   "minV2minP16_TrunkAsign.fits",drawChannel=drawChannel,lRange=lRange,bRange= bRange )


		channelDENDRO=drawDendrogram-1


		dendroClouds=np.unique(channelDENDRO)



		#scimes
		#drawSCIMESData, drawSCIMESHead = myFITS.readFITS("./scimesG2650/ClusterAsgn_2_8Ve20.fits")
		#channelSCIMES =  drawSCIMESData[drawChannel]

		drawSCIMES, drawSCIMESHead =self.getCropDataAndHead(   "./scimesG2650/ClusterAsgn_2_16Ve20.fits", drawChannel=drawChannel,lRange=lRange,bRange= bRange )
		channelSCIMES = drawSCIMES  #drawSCIMESData[drawChannel]

		maximumArea= 144 *0.25 #arcmin^2

		#
		#print drawDBSCANtb.colnames



		fig = plt.figure(1, figsize=(10, 8) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica']})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

		rc('text', usetex=True)

		axCO= pywcsgrid2.subplot(221, header=   WCSCrop  )
		axCO.imshow(channelRawCO,origin='lower',cmap="bone",vmin=0 ,vmax=3,interpolation='none')



		#draw Dendrogram.............
		#the trunk assign of Dendrogrom is wong, the cloud0 ara all missed, so we ignore them

		runIndex=0

		for eachDRC in dendroClouds:
			break
			if runIndex == 0:
				labelDendro = "Dendrogram"


			else:
				labelDendro = None



			eachDRC=int(eachDRC)
			if eachDRC==0:
				continue

			cRow=  drawDENDROtb[drawDENDROtb["_idx"]==eachDRC  ]

			area=cRow["area_exact"]

			if area>maximumArea:
				continue

			else:
				#draw

				if np.isnan(cRow["x_cen"] ):
					continue
				#print eachDRC

				axCO.scatter(cRow["x_cen"], cRow["y_cen"], s=13,facecolors='none',edgecolors='r', linewidths=0.3,  label= labelDendro  )
				runIndex=runIndex+1


		#draw DBSCAN.............
		runIndex=0
		for eachDBC in dbClouds:
			break
			if runIndex == 0:
				labelDB = "DBSCAN"


			else:
				labelDB = None

			if eachDBC==0:
				continue



			cRow=  drawDBSCANtb[drawDBSCANtb["_idx"]==eachDBC  ]



			if len(cRow)==0: #may be edge sources
				continue

			area=cRow["area_exact"]


			if area>maximumArea:
				continue

			else:
				#draw

				if np.isnan(cRow["x_cen"] ):
					continue


				axCO["gal"].scatter(cRow["x_cen"], cRow["y_cen"]  , s=10,facecolors='none',edgecolors='b',linewidths=0.3, label= labelDB )
				runIndex=runIndex+1

		axCO.legend(loc=3)

		axCO.set_ticklabel_type("absdeg", "absdeg")
		axCO.axis[:].major_ticks.set_color("w")

		#draw DBSCAN labels

		cmap = plt.cm.gist_rainbow

		cmap = plt.cm.jet

		cmap.set_bad('black', 1.)

		axDBSCAN= pywcsgrid2.subplot(222, header=   WCSCO  ,sharex=axCO,sharey=axCO )

		newLabels=self.showLabels(axDBSCAN, channelDBSCAN  )

		self.contourLabels(axDBSCAN, channelDBSCAN)


		at = AnchoredText("DBSCAN", loc=1, frameon=False,prop={"color":"w"} )
		axDBSCAN.add_artist(at)
		#draw Dendrogram labels
		axDendrogram= pywcsgrid2.subplot(223, header=   WCSCO ,sharex=axCO,sharey=axCO )


		#draw dnedrogram contours
		newLabels=self.showLabels(axDendrogram, channelDENDRO  )
		self.contourLabels(axDendrogram, channelDENDRO)

		#self.contourLabels(axDendrogram,newLabels )


		at = AnchoredText("Dendrogram", loc=1, frameon=False,prop={"color":"w"} )
		axDendrogram.add_artist(at)

		#draw SCIMES

		axSCIMES= pywcsgrid2.subplot(224, header=   WCSCO ,sharex=axCO,sharey=axCO )

		newLabels=self.showLabels(axSCIMES, channelSCIMES  )

		self.contourLabels(axSCIMES,channelSCIMES )


		at = AnchoredText("SCIMES", loc=1, frameon=False,prop={"color":"w"} )
		axSCIMES.add_artist(at)




		fig.tight_layout()
		plt.savefig("checkCloud.pdf", bbox_inches="tight")

		plt.savefig("checkCloud.png", bbox_inches="tight",dpi=600)

	def showLabels(self,ax,labelFITS2DArray):
		"""
		Use jet to show labels, and
		:param ax:
		:param labelFITS2DArray:
		:return:
		"""

		labelFITS2DArray=np.nan_to_num(labelFITS2DArray)

		minV=np.nanmin(labelFITS2DArray )
		print minV
		#labelFITS2DArray = labelFITS2DArray.astype(np.float)
		#labelFITS2DArray[labelFITS2DArray==minV]= np.NaN # only float is allowed to have NaN values

		cmap = plt.cm.jet
		cmap.set_bad('black', 1. )

		dbClouds = np.unique( labelFITS2DArray )
		newLabels = np.arange( len(dbClouds) )+1

		np.random.shuffle(newLabels)

		#reassign fits

		newChannelMap = labelFITS2DArray.copy()


		#
		index1D=np.where(labelFITS2DArray > minV)
		values1D=labelFITS2DArray[index1D]

		Y0, X0 = index1D



		for i in np.arange(len(dbClouds)):

			oldID = dbClouds[i]
			newID= newLabels[i]

			if oldID==minV:
				continue

			indicesL=self.getIndices2D(Y0,X0,values1D,oldID)

			newChannelMap[  indicesL  ] = newID

		newChannelMap = newChannelMap.astype(np.float)
		newChannelMap[   labelFITS2DArray == minV   ] = np.NaN


		ax.imshow(newChannelMap,origin='lower',cmap=cmap, vmin=minV,vmax=len(newLabels),  interpolation='none')


		return newChannelMap






	def contourLabels(self,ax ,  labelFITS2DArray ):
		"""

		:param ax:
		:param labelFITS2D:
		:return:
		"""

		dbClouds=np.unique( labelFITS2DArray )

		noiseV=np.min( dbClouds)

		for eachDendroC in dbClouds:
			#pass
			if eachDendroC == noiseV:
				continue

			cloudIndex= labelFITS2DArray == eachDendroC
			cloudIndex=cloudIndex.astype(int)


			ax.contour(cloudIndex, colors="white", linewidths=0.2, origin="lower", levels=[1] )







	def getLVFITSByDBMASK(self,DBlabel,CO12FITS,PVHeadTempFITS,algorithm="dendroMask"):
		dataDB, headDB = myFITS.readFITS( DBlabel )
		minMask=np.min(dataDB[0])

		self.maskEdges(dataDB,minMask)


		dataCO,headCO= myFITS.readFITS( CO12FITS )

		pvData,pvHead= myFITS.readFITS( PVHeadTempFITS )

		mask=dataDB>minMask
		mask=mask+0
		coMask=dataCO*mask
		Nz,Ny,Nx=dataCO.shape
		pvData=np.nansum(coMask, axis=1 )/Ny

		fits.writeto("G2650PV_{}.fits".format( algorithm ),pvData,header=pvHead, overwrite=True)

	def cleanAllDBfits(self):

		DbscanSigmaList = np.arange(2, 7.5, 0.5)


		for sigmas in DbscanSigmaList:

			for minPix in [8,16]:

				tbName = "G2650CO12DBCatS{:.1f}P{}Con2.fit".format(sigmas, 8)
				fitsName = "G2650CO12dbscanS{:.1f}P{}Con2.fits".format(sigmas, 8 )
				self.clearnDBAssign( fitsName,tbName, pixN=minPix, minV=sigmas, minDelta= 3 )



	def splitEdges(self,TB):
		"""

		:param TB:
		:return: good sources , and edge sources

		"""
 
		processTB=TB

		dbClusterTB=processTB
		if "peak" in dbClusterTB.colnames: #for db scan table

			select1=np.logical_and( processTB["x_cen"]<= 26.25 ,processTB["y_cen"] >= 3.25  )
			select2=np.logical_and( processTB["x_cen"]>=49.25 ,processTB["y_cen"]>=  3.75 )
			allSelect= np.logical_or( select1,select2 )

			badSource=dbClusterTB[allSelect]

			part1= processTB[ np.logical_or( processTB["x_cen"]>26.25 ,processTB["y_cen"] < 3.25  )   ] #1003, 3.25
			part2= part1[ np.logical_or( part1["x_cen"]<49.25 ,part1["y_cen"]<  3.75 )   ] #1003, 3.25


			return part2,badSource



		else: #dendrogram tb


			select1=np.logical_and( processTB["x_cen"]>= 2815 ,processTB["y_cen"] >= 1003  )
			select2=np.logical_and( processTB["x_cen"]<=  55 ,processTB["y_cen"]>= 1063  )
			allSelect= np.logical_or( select1,select2 )


			badSource=dbClusterTB[allSelect]
			part1= processTB[ np.logical_or( processTB["x_cen"]< 2815 ,processTB["y_cen"] < 1003  )   ] #1003, 3.25
			part2= part1[ np.logical_or( part1["x_cen"]>  55 ,part1["y_cen"]< 1063  )   ] #1003, 3.25



			return part2,badSource


		#used to remove edge sources, and




	def produceMask(self, COFITS, LabelFITS ,labelTB, region=""):

		#mask pixels that have been not ben labeled by other

		dataCluster, headCluster = myFITS.readFITS(LabelFITS)

		dataCO, headCO = myFITS.readFITS( COFITS )
		rmsData,rmsHead=myFITS.readFITS( "/home/qzyan/WORK/myDownloads/testScimes/RMS_G2650CO12.fits" )



		dbClusterTB=Table.read( labelTB )

		saveFITS=region+"MaskedCO.fits"
		minV= np.nanmin(dataCluster[0])
		print "The noise is maked with ",minV

		coGood=  dataCluster>minV


		clusterIndex1D= np.where( dataCluster>minV )
		clusterValue1D=  dataCluster[clusterIndex1D ]

		Z0,Y0,X0 = clusterIndex1D



		dataCO[ ~coGood ]=0#including nan


		#getedigeTBlist
		processTB=dbClusterTB
		
		goodSource,badSource=self.splitEdges( processTB )


		for eachBadC in badSource:

			badID=eachBadC["_idx"]
			cloudIndex = self.getIndices(Z0, Y0, X0, clusterValue1D, badID)
			
			dataCO[ cloudIndex ]=0




		fits.writeto(saveFITS,dataCO, header=headCO,   overwrite = True )

	def getDiffCO(self,maskedFITS,labelFITS, cutoff=3, ):


		dataCluster, headCluster = myFITS.readFITS(labelFITS)

		dataCO, headCO = myFITS.readFITS( maskedFITS )

		minV= np.nanmin(dataCluster)

		dataCO[dataCO<cutoff*self.rms]=0

		dataCO[dataCluster> minV]=0

		fits.writeto( "DiffTest.fits" , dataCO, header=headCluster,overwrite=True   )


	def getDiffLabel(self,label1FITS,label2FITS):
		"""
		#produce the region that has no overlapping between the two labeled fits

		:param label1FITS:
		:param label2FITS:
		:return:
		"""

		dataCluster1, headCluster1 = myFITS.readFITS(label1FITS)
		minV1=np.min( dataCluster1[0] )

		dataCluster1Copy= dataCluster1.copy()



		dataCluster2, headCluster2 = myFITS.readFITS(label2FITS)
		minV2=np.min( dataCluster2[0] )
		dataCluster2Copy= dataCluster2.copy()



		dataCluster1Copy[ dataCluster2 > minV2 ]=0

		dataCluster2Copy[ dataCluster1 > minV1 ]=0


		#save two copies,
		fits.writeto("Only_"+label1FITS,  dataCluster1Copy ,header= headCluster1 ,overwrite=True)

		fits.writeto("Only_"+label2FITS,  dataCluster2Copy ,header= headCluster2 ,overwrite=True)






	def extendAllScimes(self  ): #only extend,no others
		"""
		Extend the scimes Asgncube
		:return:
		"""
		#the SCIMES assign cubes stars from -1
		#scimesPath="./scimesG2650/"
		#minV=-1
		G2650MaskCO = "G2650CO12MaskedCO.fits"
		scimesPath="/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/"
		#dendroSigmaList = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7], negative values found in fits, no idea why
		dendroSigmaList = [  3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7]

		for eachSigma  in dendroSigmaList :
			for eachPix in [8, 16]:
				scimesAsgnFITS=scimesPath+"ClusterAsgn_{}_{}Ve20.fits".format( eachSigma ,eachPix )

				if os.path.isfile( scimesAsgnFITS):

					saveName="ClusterAsgn_{}_{}Ve20_extend.fits".format( eachSigma ,eachPix )
					self.myDilation(scimesAsgnFITS, None,startSigma=15, saveName=saveName, endSigma=eachSigma,maskCOFITS=G2650MaskCO ,savePath= scimesPath)


	def testFluxOfUM(self):

		UMfitsDB="./ursaMajor/UMCO12dbscanS2P16Con2.fits"
		dataDB,headDB=myFITS.readFITS(UMfitsDB)
		UMCO= "/home/qzyan/WORK/projects/NewUrsaMajorPaper/OriginalFITS/myCut12CO.fits"
		dataCO, headCO =myFITS.readFITS( UMCO )


		dbMask= dataDB>0
		dataCO[~dbMask]=0

		fluxList=[]
		sigmaList=np.arange(2,8,0.5)
		for sigmas in sigmaList:
			print sigmas
			dataCO[dataCO<sigmas*self.rms]=0

			totalFlux=np.sum(dataCO )
			fluxList.append( totalFlux*0.16)

		#plot
		fig=plt.figure(figsize=(10, 8))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })

		#############   plot dendrogram
		axUM=fig.add_subplot(1,1,1)

		axUM.plot(sigmaList, fluxList, 'o--' ,linestyle=':', color="blue", lw=1.5 )


		axUM.set_ylabel(r"Total flux ($\rm K\ km\ s$$^{-1}$ $\Omega_\mathrm{A}$)")
		axUM.set_xlabel(r"CO cutoff ($\sigma$)")

		plt.savefig( "ursaMajorFlux.pdf"  ,  bbox_inches='tight')
		plt.savefig( "ursaMajorFlux.png"  ,  bbox_inches='tight', dpi=300)



	def getAllSpectral(self  ):
		"""

		:param labelFITS:
		:param TB:
		:param savePath:
		:return:
		"""
		savePath="./spectralLineOfAll/"
		#labeledFITS= "G2650minV3minP16_TrunkAsignMask0.fits"
		labeledFITS = "G2650minV2minP8_TrunkAsignMask0.fits"

		#drawRegion="3sigma16pDendro"
		drawRegion="2sigma8pDendro"

		#first save TB to the path
		algDendro="Dendrogram"
		tb8Den,tb16Den,label8Den,label16Den,sigmaListDen=self.getTBList(algorithm=algDendro)
		
		tb8Den=self.removeAllEdges(tb8Den)
		tb16Den=self.removeAllEdges(tb16Den)
		
		TB= tb8Den[0]
		TB.write(savePath+"sourceTB.fit", overwrite=True  )

		#

		#getVdispersion(self, testFITS, testTB, regionName="", vCenPix=False, saveSpectral=False, savePath=None):

		vdisDendro3 = self.getVdispersion( labeledFITS , TB , regionName= drawRegion ,  vCenPix=True, savePath=savePath, saveSpectral=True )




	def drawVeDistribution(self, useGauss=True):
		"""
		:return:
		"""
		####
		#
		#first check if the trunk assign of dendrogram is right


		algDendro="Dendrogram"
		tb8Den,tb16Den,label8Den,label16Den,sigmaListDen=self.getTBList(algorithm=algDendro)
		tb8Den=self.removeAllEdges(tb8Den)
		tb16Den=self.removeAllEdges(tb16Den)


		algDB="DBSCAN"
		tb8DB,tb16DB,label8DB,label16DB,sigmaListDB=self.getTBList(algorithm=algDB)
		tb8DB=self.removeAllEdges(tb8DB)
		tb16DB=self.removeAllEdges(tb16DB)

		algSCI="SCIMES"
		tb8SCI,tb16SCI,label8SCI,label16SCI,sigmaListSCI=self.getTBList(algorithm= algSCI )
		tb8SCI=self.removeAllEdges(tb8SCI)
		tb16SCI=self.removeAllEdges(tb16SCI)

		velEdges = np.linspace(0,15,300)
		areaCenter = self.getEdgeCenter( velEdges )


		fig=plt.figure(figsize=(16,5))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]



		tbIndex1=0
		tbIndex2=4
		tbIndex3=8

		#Dendrogram....
		axDendro=fig.add_subplot(1,3,1)
		if useGauss:
			vdisDendro3=self.getVdispersion( "G2650minV2minP8_TrunkAsignMask0.fits" ,tb8Den[tbIndex1],regionName="2sigma8pDendro", vCenPix=True )
			print np.min(vdisDendro3), np.max( vdisDendro3 )
		else:
			vdisDendro3=None
		color2 = self.drawVelTBSingle(axDendro, tb8Den[tbIndex1] , velEdges ,areaCenter,label= label8Den[tbIndex1] ,inputVelDispersion=vdisDendro3 , useInputVel=useGauss )

		if useGauss:
			vdisDendro5 =self.getVdispersion( "G2650minV4minP8_TrunkAsignMask0.fits" ,tb8Den[tbIndex2],regionName="4sigma8pDendro", vCenPix=True )
			print np.min(vdisDendro5), np.max( vdisDendro5 )
			
		else:
			vdisDendro5=None

		color6 = self.drawVelTBSingle(axDendro,  tb8Den[tbIndex2] , velEdges ,areaCenter,label= label8Den[tbIndex2]  ,inputVelDispersion=vdisDendro5 , useInputVel=useGauss )

		if useGauss:
			vdisDendro7 =self.getVdispersion( "G2650minV6minP8_TrunkAsignMask0.fits" ,tb8Den[tbIndex3],regionName="6sigma8pDendro", vCenPix=True )
			print np.min(vdisDendro7), np.max( vdisDendro7 )

		else:
			vdisDendro7=None

			
		color10 = self.drawVelTBSingle(axDendro,  tb8Den[tbIndex3] , velEdges ,areaCenter,label= label8Den[tbIndex3]  ,inputVelDispersion=vdisDendro7 , useInputVel=useGauss )





		l=axDendro.legend(loc=1)

		colorsDendro=[color2, color6, color10 ]

		for text,color in zip(l.get_texts(), colorsDendro):
			text.set_color(color)





		at = AnchoredText("Dendrogram", loc=4, frameon=False ,pad=1)
		axDendro.add_artist(at)


		axDendro.set_ylabel(r"Number of trunks")

		axDendro.set_xlabel(r"Equivalent linewidth ($\rm km\ s$$^{-1}$)")

		#DBSCAN........
		axDBSCAN=fig.add_subplot(1,3,2,sharex=axDendro,sharey=axDendro)

		#getGaussFitting
		if useGauss:
			vdisDBSCAN3=self.getVdispersion( "G2650CO12dbscanS2.0P8Con2.fits" ,tb8DB[tbIndex1],regionName="2sigma8pDBSCAN" )
			print np.min(vdisDBSCAN3), np.max( vdisDBSCAN3 )
		else:
			vdisDBSCAN3=None
		color2=self.drawVelTBSingle(axDBSCAN, tb8DB[tbIndex1] , velEdges ,areaCenter,label= label8DB[tbIndex1] ,inputVelDispersion=vdisDBSCAN3 , useInputVel=useGauss )

		if useGauss:
			vdisDBSCAN5=self.getVdispersion( "G2650CO12dbscanS4.0P8Con2.fits" ,tb8DB[tbIndex2],regionName="4sigma8pDBSCAN" )
			print np.min(vdisDBSCAN5), np.max( vdisDBSCAN5 )
		else:
			vdisDBSCAN5=None
		color6=self.drawVelTBSingle(axDBSCAN,  tb8DB[tbIndex2] , velEdges ,areaCenter,label= label8DB[tbIndex2]   ,inputVelDispersion=vdisDBSCAN5 , useInputVel=useGauss )
		if useGauss:
			vdisDBSCAN7=self.getVdispersion( "G2650CO12dbscanS6.0P8Con2.fits" ,tb8DB[tbIndex3],regionName="6sigma8pDBSCAN" )
			print np.min(vdisDBSCAN7), np.max( vdisDBSCAN7 )
		else:
			vdisDBSCAN7=None

		color10=self.drawVelTBSingle(axDBSCAN,  tb8DB[tbIndex3] , velEdges ,areaCenter,label= label8DB[tbIndex3]  ,inputVelDispersion=vdisDBSCAN7 , useInputVel=useGauss )
		l=axDBSCAN.legend(loc=1)
 
		colorsDendro=[color2, color6, color10 ]

		for text,color in zip(l.get_texts(), colorsDendro):
			text.set_color(color)



		at = AnchoredText("DBSCAN", loc=4, frameon=False,pad=1)
		axDBSCAN.add_artist(at)
		axDBSCAN.set_xlabel(r"Equivalent linewidth ($\rm km\ s$$^{-1}$)")
		axDBSCAN.set_ylabel(r"Number of trunks")


		##SCIMES

		axSCIMES=fig.add_subplot(1,3,3,sharex=axDendro,sharey=axDendro)
		if useGauss:
			vdisSCIMES3=self.getVdispersion( "/home/qzyan/WORK/myDownloads/MWISPcloud/ClusterAsgn_2_8Ve20_mannual.fits" ,tb8SCI[tbIndex1],regionName="2sigma8pSCIMES", vCenPix=True )
			print np.min(vdisSCIMES3), np.max( vdisSCIMES3 )
		else:
			vdisSCIMES3=True

		color2=self.drawVelTBSingle(axSCIMES, tb8SCI[tbIndex1] , velEdges ,areaCenter,label= label8SCI[tbIndex1]  ,inputVelDispersion = vdisSCIMES3 , useInputVel=useGauss )

		if useGauss:
			vdisSCIMES5 = self.getVdispersion( "/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/ClusterAsgn_4_8Ve20.fits", tb8SCI[tbIndex2], regionName="4sigma8pSCIMES", vCenPix=True)
			print np.min(vdisSCIMES5), np.max( vdisSCIMES5 )
		else:
			vdisSCIMES5=True


		color6=self.drawVelTBSingle(axSCIMES,  tb8SCI[tbIndex2] , velEdges ,areaCenter,label= label8SCI[tbIndex2]    ,inputVelDispersion = vdisSCIMES5 , useInputVel=useGauss )
		
		if useGauss:
			vdisSCIMES7=self.getVdispersion( "/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/ClusterAsgn_6_8Ve20.fits" ,tb8SCI[tbIndex3],regionName="6sigma8pSCIMES", vCenPix=True )
			print np.min(vdisSCIMES7), np.max( vdisSCIMES7 )
		else:
			vdisSCIMES7=None
		color10=self.drawVelTBSingle(axSCIMES,  tb8SCI[tbIndex3] , velEdges ,areaCenter,label= label8SCI[tbIndex3]   ,inputVelDispersion = vdisSCIMES7 , useInputVel=useGauss )


		
		l=axSCIMES.legend(loc=1)
		colorsDendro=[color2, color6, color10 ]

		for text,color in zip(l.get_texts(), colorsDendro):
			text.set_color(color)

		at = AnchoredText("SCIMES", loc=4, frameon=False ,pad=1)
		axSCIMES.add_artist(at)


		axSCIMES.set_xlabel(r"Equivalent linewidth ($\rm km\ s$$^{-1}$)")

		axSCIMES.set_ylabel(r"Number of clusters")


		axSCIMES.set_xlim([0,10])


		#axDBSCAN.set_xlabel(r"Velocity dispersion $\sigma_V$ ($\rm km\ s$$^{-1}$)")

		fig.tight_layout()
		plt.savefig( "velDistribute.pdf"  ,  bbox_inches='tight')
		plt.savefig( "velDistribute.png"  ,  bbox_inches='tight', dpi=300)


	def drawVelTBSingle(self,ax,testTB,velEdges,areaCenter,label=None,inputVelDispersion=None,useInputVel=False):


		#testTB=self.removeWrongEdges(testTB)
		if not useInputVel :
			vDisperse = testTB["v_rms"]
		else:
			vDisperse =  inputVelDispersion

		maxV= np.max(vDisperse)
		meanV= np.mean(vDisperse)
		medianV= np.median(vDisperse)


		binN8 , binEdges8 =np.histogram(vDisperse, bins=velEdges  )

		peakV=areaCenter[ binN8.argmax() ]

		suffix="\nPeak: {:.2f}, Mean: {:.2f}, Median: {:.2f}".format( peakV, meanV, medianV )


		stepa=ax.step( areaCenter, binN8, lw=1.0 , label=label+suffix )





		return stepa[-1].get_color()

	def drawPeakTBSingle(self,ax,testTB,velEdges,areaCenter,label=None):

		#testTB=self.removeWrongEdges(testTB)
		vDisperse = testTB["peak"]

		peakV= np.max(vDisperse)
		meanV= np.mean(vDisperse)
		medianV= np.median(vDisperse)

		suffix="\nMean: {:.1f}, Median: {:.1f}".format(   meanV, medianV )

		print np.min(vDisperse),np.max(vDisperse)

		binN8 , binEdges8 =np.histogram(vDisperse, bins=velEdges  )
		stepa = ax.step( areaCenter, binN8, lw=1.0 , label=label+suffix )
		#ax.set_yscale('log')



		return stepa[-1].get_color()

	def drawPeakDistribution(self):
		"""
		:return:
		"""
		####

		# first check if the trunk assign of dendrogram is right

		algDendro="Dendrogram"
		tb8Den,tb16Den,label8Den,label16Den,sigmaListDen=self.getTBList(algorithm=algDendro)
		tb8Den=self.removeAllEdges(tb8Den)
		tb16Den=self.removeAllEdges(tb16Den)


		algDB="DBSCAN"
		tb8DB,tb16DB,label8DB,label16DB,sigmaListDB=self.getTBList(algorithm=algDB)
		tb8DB=self.removeAllEdges(tb8DB)
		tb16DB=self.removeAllEdges(tb16DB)

		algSCI="SCIMES"
		tb8SCI,tb16SCI,label8SCI,label16SCI,sigmaListSCI=self.getTBList(algorithm= algSCI )
		tb8SCI=self.removeAllEdges(tb8SCI)
		tb16SCI=self.removeAllEdges(tb16SCI)

		velEdges = np.linspace(0,40,200)
		areaCenter = self.getEdgeCenter( velEdges )


		fig=plt.figure(figsize=(16,5))
		rc('text', usetex=True )
		rc('font', **{'family': 'sans-serif',  'size'   : 13,  'serif': ['Helvetica'] })

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]



		#Dendrogram....
		axDendro=fig.add_subplot(1,3,1)
		#self.drawVelTBSingle(axDendro, tb16Den[2] , velEdges ,areaCenter,label= label16Den[2]  )
		#self.drawVelTBSingle(axDendro,  tb16Den[6] , velEdges ,areaCenter,label= label16Den[6]   )
		#self.drawVelTBSingle(axDendro,  tb16Den[10] , velEdges ,areaCenter,label= label16Den[10]   )

		#to be modified
		tbDendro2=Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV2minP8dendroMannualCat.fit")
		tbDendro6=Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV4minP8dendroMannualCat.fit")
		tbDendro10=Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV6minP8dendroMannualCat.fit")

		index1=0
		index2=4
		index3=8


		tbDendro2, tbDendro6,tbDendro10 =self.removeAllEdges( [ tbDendro2, tbDendro6,tbDendro10 ] )



		color2  = self.drawPeakTBSingle(axDendro,tbDendro2   , velEdges ,areaCenter,label= label8Den[index1]  )
		color6  = self.drawPeakTBSingle(axDendro,  tbDendro6 , velEdges ,areaCenter,label= label8Den[index2]   )
		color10 = self.drawPeakTBSingle(axDendro,  tbDendro10 , velEdges ,areaCenter,label= label8Den[index3]   )


		l=axDendro.legend(loc=1)
		colorsDendro=[color2, color6, color10 ]
		for text,color in zip(l.get_texts(), colorsDendro):
			text.set_color(color)



		at = AnchoredText("Dendrogram", loc=4, frameon=False,pad=1 )
		axDendro.add_artist(at)


		axDendro.set_ylabel(r"Number of trunks")

		
		axDendro.set_xlabel(r"Peak brightness temperature (K)")

		#DBSCAN........
		axDBSCAN=fig.add_subplot(1,3,2,sharex=axDendro,sharey=axDendro)
		color2= self.drawPeakTBSingle(axDBSCAN, tb8DB[index1] , velEdges ,areaCenter,label= label8DB[index1]  )
		color6= self.drawPeakTBSingle(axDBSCAN,  tb8DB[index2] , velEdges ,areaCenter,label= label8DB[index2]   )
		color10= self.drawPeakTBSingle(axDBSCAN,  tb8DB[index3] , velEdges ,areaCenter,label= label8DB[index3]   )



		l=axDBSCAN.legend(loc=1)
		colorsDendro=[color2, color6, color10 ]
		for text,color in zip(l.get_texts(), colorsDendro):
			text.set_color(color)




		at = AnchoredText("DBSCAN", loc=4, frameon=False, pad=1 )
		axDBSCAN.add_artist(at)
		#axDBSCAN.set_xlabel(r"Peak values (K)")
		axDBSCAN.set_xlabel(r"Peak brightness temperature (K)")
		axDBSCAN.set_ylabel(r"Number of trunks")

		##SCIMES

		axSCIMES=fig.add_subplot(1,3,3,sharex=axDendro,sharey=axDendro)

		tbSCIMES2=Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV2P8scimesMannual.fit")
		tbSCIMES6=Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV4P8scimesMannual.fit")
		tbSCIMES10=Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV6P8scimesMannual.fit")

		tbSCIMES2, tbSCIMES6,tbSCIMES10 =self.removeAllEdges( [ tbSCIMES2, tbSCIMES6,tbSCIMES10 ] )



		color2=self.drawPeakTBSingle(axSCIMES, tbSCIMES2 , velEdges ,areaCenter,label= label8SCI[index1]  )
		color6=self.drawPeakTBSingle(axSCIMES,  tbSCIMES6 , velEdges ,areaCenter,label= label8SCI[index2]   )
		color10=self.drawPeakTBSingle(axSCIMES,  tbSCIMES10 , velEdges ,areaCenter,label= label8SCI[index3]   )
		l=axSCIMES.legend(loc=1)


		colorsDendro=[color2, color6, color10 ]
		for text,color in zip(l.get_texts(), colorsDendro):
			text.set_color(color)




		at = AnchoredText("SCIMES", loc=4, frameon=False, pad=1 )
		axSCIMES.add_artist(at)

		axSCIMES.set_ylabel(r"Number of clusters")

		axSCIMES.set_xlabel(r"Peak brightness temperature (K)")

		axSCIMES.set_xlim(0,15)
		fig.tight_layout()
		plt.savefig("peakDistribute.pdf", bbox_inches='tight')
		plt.savefig("peakDistribute.png", bbox_inches='tight', dpi=300)




	def produceSCIMECat(self):
		"""
		To get the peak of clusters
		:return:
		"""
		sicmesPath="/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/"
		labelFITS = sicmesPath+"ClusterAsgn_3_16Ve20.fits"
		savename = sicmesPath+"ClusterAsgn_3_16Ve20ManualCat"
		doDBSCAN.getCatFromLabelArray(G2650CO12FITS, labelFITS, doDBSCAN.TBModel, minPix=16, rms=3,   saveMarker=savename)

		#######################
		sicmesPath="/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/"
		labelFITS = sicmesPath+"ClusterAsgn_5_16Ve20.fits"
		savename = sicmesPath+"ClusterAsgn_5_16Ve20ManualCat"
		doDBSCAN.getCatFromLabelArray(G2650CO12FITS, labelFITS, doDBSCAN.TBModel, minPix=16, rms=5,   saveMarker=savename)

		#######################
		sicmesPath="/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/"
		labelFITS = sicmesPath+"ClusterAsgn_7_16Ve20.fits"
		savename = sicmesPath+"ClusterAsgn_7_16Ve20ManualCat"
		doDBSCAN.getCatFromLabelArray(G2650CO12FITS, labelFITS, doDBSCAN.TBModel, minPix=16, rms=7,   saveMarker=savename)



	def produceDENDROCat(self):
		"""
		To get the peak of clusters
		:return:
		"""
		#sicmesPath="/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/"
		labelFITS =  "G2650minV7minP16_TrunkAsignMask0.fits"
		savename =  "Dendro_7_16Ve20ManualCat"
		doDBSCAN.getCatFromLabelArray(G2650CO12FITS, labelFITS, doDBSCAN.TBModel, minPix=16, rms=3,   saveMarker=savename)

		#######################
		#sicmesPath="/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/"
		labelFITS =  "G2650minV5minP16_TrunkAsignMask0.fits"
		savename =  "Dendro_5_16Ve20ManualCat"
		doDBSCAN.getCatFromLabelArray(G2650CO12FITS, labelFITS, doDBSCAN.TBModel, minPix=16, rms=5,   saveMarker=savename)

		#######################
		#sicmesPath="/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/"
		labelFITS =   "G2650minV3minP16_TrunkAsignMask0.fits"
		savename =  "Dendro_3_16Ve20ManualCat"
		doDBSCAN.getCatFromLabelArray(G2650CO12FITS, labelFITS, doDBSCAN.TBModel, minPix=16, rms=7,   saveMarker=savename)


	def removePerseusSingle(self,eachT,upTo=1  ):
		"""
		By default, upTo 0.5 degree
		:param mcTB:
		:return:
		"""
		#get upto pixel
		data,head=doFITS.readFITS("G2650Local30.fits")

		wcsCO= WCS(head )
		_,indexLow,_ = wcsCO.wcs_world2pix(40,-upTo, 0 ,0 )
		_,indexUp,_ = wcsCO.wcs_world2pix(40, upTo , 0 ,0 )



		if np.max(eachT["y_cen"]) > 10: #pixel

			select = np.logical_or(eachT["y_cen"] < indexLow, eachT["y_cen"]  > indexUp)
			locateTB = eachT[select]


		else:
			select = np.logical_or(eachT["y_cen"] < -upTo, eachT["y_cen"] > upTo )
			locateTB = eachT[select]

		return locateTB

	def removePerseusList(self,TBList):

		newList=[]

		for eachT in TBList:
			newList.append( self.removePerseusSingle(eachT)  )

		return newList




	def physicalAreaDistributionLocal(self): # ###
		"""
		#draw the physical Area Distribution distribution of molecular clouds
		:return:
		"""

		#9pixels at 1500 pc
		completeArea=0.428 #36824657505895 #pc^2 should equal to

		#use 3,5,7#because they all have observed area

		algDendro = "Dendrogram"
		tb8Den, tb16Den, label8Den, label16Den, sigmaListDen = self.getTBList(algorithm=algDendro)
		tb8Den = self.removeAllEdges(tb8Den)
		tb16Den = self.removeAllEdges(tb16Den)

		tb8Den = self.removePerseusList(tb8Den)
		tb16Den = self.removePerseusList(tb16Den)




		algDB = "DBSCAN"
		tb8DB, tb16DB, label8DB, label16DB, sigmaListDB = self.getTBList(algorithm=algDB)
		tb8DB = self.removeAllEdges(tb8DB)
		tb16DB = self.removeAllEdges(tb16DB)

		tb8DB = self.removePerseusList(tb8DB)
		tb16DB = self.removePerseusList(tb16DB)





		algSCI = "SCIMES"
		tb8SCI, tb16SCI, label8SCI, label16SCI, sigmaListSCI = self.getTBList(algorithm=algSCI)
		tb8SCI = self.removeAllEdges(tb8SCI)
		tb16SCI = self.removeAllEdges(tb16SCI)


		tb8SCI = self.removePerseusList(tb8SCI)
		tb16SCI = self.removePerseusList(tb16SCI)




		index1=0
		index2=4
		index3=8

		## physicalEdges  physicalCenter
		physicalEdges=np.linspace(0,100,1000) #square pc^2
		physicalCenter=self.getEdgeCenter( physicalEdges )


		velEdges = np.linspace(0, 40, 200)
		areaCenter = self.getEdgeCenter(velEdges)

		fig = plt.figure(figsize=(16, 5))
		rc('text', usetex=True)
		rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]



		# Dendrogram....
		axDendro = fig.add_subplot(1, 3, 1)
		# self.drawVelTBSingle(axDendro, tb16Den[2] , velEdges ,areaCenter,label= label16Den[2]  )
		# self.drawVelTBSingle(axDendro,  tb16Den[6] , velEdges ,areaCenter,label= label16Den[6]   )
		# self.drawVelTBSingle(axDendro,  tb16Den[10] , velEdges ,areaCenter,label= label16Den[10]   )

		# to be modified
		tbDendro2 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV2minP8dendroMannualCat.fit")
		tbDendro6 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV4minP8dendroMannualCat.fit")
		tbDendro10 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV6minP8dendroMannualCat.fit")

		tbDendro2, tbDendro6, tbDendro10 = self.removeAllEdges([tbDendro2, tbDendro6, tbDendro10])
		tbDendro2, tbDendro6, tbDendro10 = self.removePerseusList([tbDendro2, tbDendro6, tbDendro10])




		color2 = self.drawPhysicalAreaSingle(axDendro, tbDendro2, physicalEdges, physicalCenter, completeArea,label=label8Den[index1])
		color6 = self.drawPhysicalAreaSingle(axDendro, tbDendro6, physicalEdges, physicalCenter,completeArea, label=label8Den[index2])
		color10 = self.drawPhysicalAreaSingle(axDendro, tbDendro10, physicalEdges, physicalCenter, completeArea, label=label8Den[index3])



		l = axDendro.legend(loc=1)
		colorsDendro = [color2, color6, color10]
		for text, color in zip(l.get_texts(), colorsDendro):
			text.set_color(color)


		#draw complete line
		axDendro.plot( [completeArea,completeArea],[2,2000],'--',color='black', lw=1  )


		at = AnchoredText("Dendrogram ($|b|>1^{\circ}$)", loc=3, frameon=False, pad=0.2)
		axDendro.add_artist(at)


		xLabelStr= r"Phyiscal area ($\rm pc^{2}$)"
		axDendro.set_ylabel(r"Number of trunks")

		axDendro.set_xlabel( xLabelStr )


		axDendro.set_yscale('log')
		axDendro.set_xscale('log')


		# DBSCAN........
		axDBSCAN = fig.add_subplot(1, 3, 2, sharex=axDendro, sharey=axDendro)



		color2 = self.drawPhysicalAreaSingle(axDBSCAN, tb8DB[index1], physicalEdges, physicalCenter, completeArea,label=label8DB[index1])
		color6 = self.drawPhysicalAreaSingle(axDBSCAN, tb8DB[index2], physicalEdges, physicalCenter,completeArea, label=label8DB[index2])
		color10 = self.drawPhysicalAreaSingle(axDBSCAN, tb8DB[index3], physicalEdges, physicalCenter, completeArea, label=label8DB[index3])


		l = axDBSCAN.legend(loc=1)
		colorsDendro = [color2, color6, color10]
		for text, color in zip(l.get_texts(), colorsDendro):
			text.set_color(color)

		at = AnchoredText(r"DBSCAN ($|b|>1^{\circ}$)", loc=3, frameon=False, pad=0.2)
		axDBSCAN.add_artist(at)
		# axDBSCAN.set_xlabel(r"Peak values (K)")
		axDBSCAN.set_xlabel( xLabelStr )
		axDBSCAN.set_ylabel(r"Number of trunks")
		axDBSCAN.plot( [completeArea,completeArea],[2,2000],'--',color='black', lw=1  )


		##SCIMES

		axSCIMES = fig.add_subplot(1, 3, 3, sharex=axDendro, sharey=axDendro)

		tbSCIMES2 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV4P8scimesMannual.fit")
		tbSCIMES6 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV6P8scimesMannual.fit")
		tbSCIMES10 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV2P8scimesMannual.fit")

		tbSCIMES2, tbSCIMES6, tbSCIMES10 = self.removeAllEdges([tbSCIMES2, tbSCIMES6, tbSCIMES10])
		tbSCIMES2, tbSCIMES6, tbSCIMES10 = self.removePerseusList([tbSCIMES2, tbSCIMES6, tbSCIMES10])


		color2 = self.drawPhysicalAreaSingle(axSCIMES, tbSCIMES2, physicalEdges, physicalCenter, completeArea,label=label8SCI[index1])
		color6 = self.drawPhysicalAreaSingle(axSCIMES, tbSCIMES6 , physicalEdges, physicalCenter,completeArea, label=label8SCI[index2])
		color10 = self.drawPhysicalAreaSingle(axSCIMES,  tbSCIMES10 , physicalEdges, physicalCenter, completeArea, label=label8SCI[index3])



		l = axSCIMES.legend(loc=1)

		colorsDendro = [color2, color6, color10]
		for text, color in zip(l.get_texts(), colorsDendro):
			text.set_color(color)

		at = AnchoredText("SCIMES ($|b|>1^{\circ}$)", loc=3, frameon=False, pad=0.2)
		axSCIMES.add_artist(at)



		axSCIMES.set_ylabel(r"Number of clusters")
		axSCIMES.set_xlabel(xLabelStr)
		axSCIMES.plot( [completeArea,completeArea],[2,2000],'--',color='black', lw=1  )
		axSCIMES.set_yscale('log')
		axSCIMES.set_xscale('log')



		fig.tight_layout()
		plt.savefig("physicalAreaDistributeOnlyLocal.pdf", bbox_inches='tight')
		plt.savefig("physicalAreaDistributeOnlyLocal.png", bbox_inches='tight', dpi=300)




	def physicalAreaDistribution(self):
		"""
		#draw the physical Area Distribution distribution of molecular clouds
		:return:
		"""

		#9pixels at 1500 pc
		completeArea=0.428 # 9 pixels #pc^2 should equal to

		#use 3,5,7#because they all have observed area

		algDendro = "Dendrogram"
		tb8Den, tb16Den, label8Den, label16Den, sigmaListDen = self.getTBList(algorithm=algDendro)
		tb8Den = self.removeAllEdges(tb8Den)
		tb16Den = self.removeAllEdges(tb16Den)

		algDB = "DBSCAN"
		tb8DB, tb16DB, label8DB, label16DB, sigmaListDB = self.getTBList(algorithm=algDB)
		tb8DB = self.removeAllEdges(tb8DB)
		tb16DB = self.removeAllEdges(tb16DB)

		algSCI = "SCIMES"
		tb8SCI, tb16SCI, label8SCI, label16SCI, sigmaListSCI = self.getTBList(algorithm=algSCI)
		tb8SCI = self.removeAllEdges(tb8SCI)
		tb16SCI = self.removeAllEdges(tb16SCI)

		index1=0
		index2=4
		index3=8

		## physicalEdges  physicalCenter
		physicalEdges=np.linspace(0,100,1000) #square pc^2
		physicalCenter=self.getEdgeCenter( physicalEdges )


		velEdges = np.linspace(0, 40, 200)
		areaCenter = self.getEdgeCenter(velEdges)

		fig = plt.figure(figsize=(16, 5))
		rc('text', usetex=True)
		rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]


		# Dendrogram....
		axDendro = fig.add_subplot(1, 3, 1)
		# self.drawVelTBSingle(axDendro, tb16Den[2] , velEdges ,areaCenter,label= label16Den[2]  )
		# self.drawVelTBSingle(axDendro,  tb16Den[6] , velEdges ,areaCenter,label= label16Den[6]   )
		# self.drawVelTBSingle(axDendro,  tb16Den[10] , velEdges ,areaCenter,label= label16Den[10]   )

		# to be modified
		tbDendro2 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV2minP8dendroMannualCat.fit")
		tbDendro6 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV4minP8dendroMannualCat.fit")
		tbDendro10 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV6minP8dendroMannualCat.fit")

		tbDendro2, tbDendro6, tbDendro10 = self.removeAllEdges([tbDendro2, tbDendro6, tbDendro10])

		color2 = self.drawPhysicalAreaSingle(axDendro, tbDendro2, physicalEdges, physicalCenter, completeArea,label=label8Den[index1])
		color6 = self.drawPhysicalAreaSingle(axDendro, tbDendro6, physicalEdges, physicalCenter,completeArea, label=label8Den[index2])
		color10 = self.drawPhysicalAreaSingle(axDendro, tbDendro10, physicalEdges, physicalCenter, completeArea, label=label8Den[index3])



		l = axDendro.legend(loc=1)
		colorsDendro = [color2, color6, color10]
		for text, color in zip(l.get_texts(), colorsDendro):
			text.set_color(color)


		#draw complete line
		axDendro.plot( [completeArea,completeArea],[2,2000],'--',color='black', lw=1  )


		at = AnchoredText("Dendrogram", loc=3, frameon=False, pad=0.2)
		axDendro.add_artist(at)


		xLabelStr= r"Phyiscal area ($\rm pc^{2}$)"
		axDendro.set_ylabel(r"Number of trunks")

		axDendro.set_xlabel( xLabelStr )


		axDendro.set_yscale('log')
		axDendro.set_xscale('log')


		# DBSCAN........
		axDBSCAN = fig.add_subplot(1, 3, 2, sharex=axDendro, sharey=axDendro)



		color2 = self.drawPhysicalAreaSingle(axDBSCAN, tb8DB[index1], physicalEdges, physicalCenter, completeArea,label=label8DB[index1])
		color6 = self.drawPhysicalAreaSingle(axDBSCAN, tb8DB[index2], physicalEdges, physicalCenter,completeArea, label=label8DB[index2])
		color10 = self.drawPhysicalAreaSingle(axDBSCAN, tb8DB[index3], physicalEdges, physicalCenter, completeArea, label=label8DB[index3])


		l = axDBSCAN.legend(loc=1)
		colorsDendro = [color2, color6, color10]
		for text, color in zip(l.get_texts(), colorsDendro):
			text.set_color(color)

		at = AnchoredText("DBSCAN", loc=3, frameon=False, pad=0.2)
		axDBSCAN.add_artist(at)
		# axDBSCAN.set_xlabel(r"Peak values (K)")
		axDBSCAN.set_xlabel( xLabelStr )
		axDBSCAN.set_ylabel(r"Number of trunks")
		axDBSCAN.plot( [completeArea,completeArea],[2,2000],'--',color='black', lw=1  )


		##SCIMES

		axSCIMES = fig.add_subplot(1, 3, 3, sharex=axDendro, sharey=axDendro)

		#tbSCIMES2 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV4P8scimesMannual.fit")
		#tbSCIMES6 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV6P8scimesMannual.fit")
		#tbSCIMES10 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV2P8scimesMannual.fit")

		tbSCIMES2 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV2P8scimesMannual.fit")
		tbSCIMES6 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV4P8scimesMannual.fit")
		tbSCIMES10 = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/minV6P8scimesMannual.fit")


		tbSCIMES2, tbSCIMES6, tbSCIMES10 = self.removeAllEdges([tbSCIMES2, tbSCIMES6, tbSCIMES10])


		color2 = self.drawPhysicalAreaSingle(axSCIMES, tbSCIMES2, physicalEdges, physicalCenter, completeArea,label=label8SCI[index1])
		color6 = self.drawPhysicalAreaSingle(axSCIMES, tbSCIMES6 , physicalEdges, physicalCenter,completeArea, label=label8SCI[index2])
		color10 = self.drawPhysicalAreaSingle(axSCIMES,  tbSCIMES10 , physicalEdges, physicalCenter, completeArea, label=label8SCI[index3])



		l = axSCIMES.legend(loc=1)

		colorsDendro = [color2, color6, color10]
		for text, color in zip(l.get_texts(), colorsDendro):
			text.set_color(color)

		at = AnchoredText("SCIMES", loc=3, frameon=False, pad=0.2)
		axSCIMES.add_artist(at)



		axSCIMES.set_ylabel(r"Number of clusters")
		axSCIMES.set_xlabel(xLabelStr)
		axSCIMES.plot( [completeArea,completeArea],[2,2000],'--',color='black', lw=1  )
		axSCIMES.set_yscale('log')
		axSCIMES.set_xscale('log')



		fig.tight_layout()
		plt.savefig("physicalAreaDistribute.pdf", bbox_inches='tight')
		plt.savefig("physicalAreaDistribute.png", bbox_inches='tight', dpi=300)







	def momentAllWithMaskedCO(self,maskedCOFITS):


		data,head=myFITS.readFITS( maskedCOFITS )
		intData=np.sum(data,axis=0)*0.2 #K km/s
		fits.writeto( "dendroMaskedInt.fits", intData, header=head,overwrite =True )



	def drawOverallMoment(self, fitsFile="normalInt.fits" ):
		"""
		draw a simple map for over all moment of local molecular clouds
		:return:
		"""
		 # what is the unit of this integrated intensity?
		import matplotlib
		from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
		from mpl_toolkits.axes_grid1.axes_grid import AxesGrid
		from matplotlib.colors import LogNorm


		dataCO,headCO=myFITS.readFITS(fitsFile)

		wcsCO=WCS(headCO)



		fig = plt.figure(1 , figsize=(16, 7.8) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'], 'size':20})
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

		rc('text', usetex=True)
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]
		# grid helper
		grid_helper = pywcsgrid2.GridHelper(wcs=wcsCO)

		# AxesGrid to display tow images side-by-side
		fig = plt.figure(1, (6, 3.5))


		grid = ImageGrid(fig, (1, 1, 1), nrows_ncols=(1, 1),
						 cbar_mode="single", cbar_pad="0.3%",
						 cbar_location="right", cbar_size="1%",
						 axes_class=(pywcsgrid2.Axes, dict(header= wcsCO)))

		main_axes = grid[0]
		main_axes.locator_params(nbins=10)
		cb_axes = grid.cbar_axes[0]  # colorbar axes

		a=np.log(dataCO)

		a=a[a>0]
		print np.min(a),np.max(a), "minimum and maxinum after "


		#im = main_axes.imshow(np.sqrt(dataCO),origin='lower',cmap="jet",  vmin=np.sqrt(1.9*1.5), vmax= np.sqrt( 81) , interpolation='none')

		#im = main_axes.imshow( dataCO ,origin='lower',cmap="jet",  vmin=-10, vmax= 80 , interpolation='none')


		#im = main_axes.imshow(np.log(dataCO),origin='lower',cmap="jet",   interpolation='none')

		im = main_axes.imshow( np.log(dataCO) ,origin='lower',cmap="jet",  vmin=np.log(1.9*2.5), vmax= np.log( 110) , interpolation='none')

		#main_axes.set_facecolor('silver')
		main_axes.axis[:].major_ticks.set_color("purple")
		cb=cb_axes.colorbar(im)
		#cb_axes.axis["right"].toggle(ticklabels=True)
		cb_axes.set_ylabel(r"K km s$^{-1}$")
		cb_axes.set_xlabel("")

		#print dir(cb_axes.axes)

		#tickesArray=np.asarray( [0,0.1,0.5,1,2,3,4,5] )
		tickesArray=np.asarray( [ 5,10,20,40,80 ] )
		#tickesArray=tickesArray**2
		cb.ax.set_yticks(   np.log(tickesArray)  )
		cb.ax.set_yticklabels( map(str,tickesArray) )

		#cbar.ax.set_xticksset_xticklabels(['Low', 'Medium', 'High'])

		#add well-known star forming regions

		#W40 LBN 028.77+03.43
		fontAlpha=1
		fontSize=15
		fontColor="white"
		lw=1.2
		import matplotlib.patheffects as path_effects
		#W40
		#main_axes["gal"].text(28.77, 03.43, r"\textbf{W40}",  color='white', alpha=fontAlpha,fontsize=fontSize ,horizontalalignment='center',   verticalalignment='center' )
		#main_axes["gal"].text(28.77, 03.43, "W40",  color=fontColor, alpha=fontAlpha,fontsize=fontSize ,horizontalalignment='center',   verticalalignment='center' )
		text=main_axes["gal"].text(28.77, 03.43, "W40",  color=fontColor, alpha=fontAlpha,fontsize=fontSize ,horizontalalignment='center',   verticalalignment='center' )
		text.set_path_effects([path_effects.Stroke(linewidth=lw, foreground="black"), path_effects.Normal()])

		#W49  G043.2+0.01  W49 2013ApJ...775...79Z
		text=main_axes["gal"].text(43.2, 0.01, "W49",  color=fontColor, alpha=fontAlpha,fontsize=fontSize,horizontalalignment='center',   verticalalignment='center' )
		text.set_path_effects([path_effects.Stroke(linewidth=lw, foreground="black"), path_effects.Normal()])

		##G032.79+00.19(H)
		#main_axes["gal"].text(32.79,   0.19 , "1",  color=fontColor, alpha=fontAlpha,fontsize=fontSize,horizontalalignment='center',   verticalalignment='center' )

		###  049.1405 -00.6028
		#main_axes["gal"].text(49.14, -00.60 , "W51",  color=fontColor, alpha=fontAlpha,fontsize=fontSize ,horizontalalignment='center',   verticalalignment='center' )

		#Aquila Rift
		rA=30
		startL=38
		startB=0.0
		dL=startL-30.5

		text=main_axes["gal"].text(startL , startB ,  "Aquila Rift",  color=fontColor, alpha=fontAlpha,fontsize=fontSize , horizontalalignment='center', rotation = rA,  verticalalignment='center' )
		text.set_path_effects([path_effects.Stroke(linewidth=lw, foreground="black"), path_effects.Normal()])


		text=main_axes["gal"].text(startL-dL , startB+dL*np.tan(np.deg2rad(rA)) ,  "Aquila Rift",  color=fontColor, alpha=fontAlpha,fontsize=fontSize , horizontalalignment='center', rotation = rA,  verticalalignment='center' )
		text.set_path_effects([path_effects.Stroke(linewidth=lw, foreground="black"), path_effects.Normal()])

		dL2=startL-34.5
		text=main_axes["gal"].text(startL-dL2 , startB+dL2*np.tan(np.deg2rad(rA)) ,  "Aquila Rift",  color=fontColor, alpha=fontAlpha,fontsize=fontSize , horizontalalignment='center', rotation = rA,  verticalalignment='center' )
		text.set_path_effects([path_effects.Stroke(linewidth=lw, foreground="black"), path_effects.Normal()])

		#Serpens NE
		text=main_axes["gal"].text(31.5768569 , 3.0808449 ,  "Serpens NE",  color=fontColor, alpha=fontAlpha,fontsize=fontSize , horizontalalignment='center',    verticalalignment='center' )
		text.set_path_effects([path_effects.Stroke(linewidth=lw, foreground="black"), path_effects.Normal()])

		   

		#LDN 566
		#main_axes["gal"].text( 030.3739 , -01.0749,  r"\textbf{LDN 566}",  color='black', alpha=fontAlpha,fontsize=fontSize , ha='center',    va='center' )
		text=main_axes["gal"].text( 030.3739 , -01.0749,  "LDN 566",  color=fontColor, alpha=fontAlpha,fontsize=fontSize , ha='center',    va='center' )
		text.set_path_effects([path_effects.Stroke(linewidth=lw, foreground="black"), path_effects.Normal()])
		#LDN 617
		text=main_axes["gal"].text( 34.5724  ,-00.862,  "LDN 617",  color=fontColor, alpha=fontAlpha,fontsize=fontSize , horizontalalignment='center',    verticalalignment='center' )
		text.set_path_effects([path_effects.Stroke(linewidth=lw, foreground="black"), path_effects.Normal()])

		#LDN 673
		text=main_axes["gal"].text( 46.263 , -01.3303,  "LDN 673",  color=fontColor, alpha=fontAlpha,fontsize=fontSize , horizontalalignment='center',    verticalalignment='center' )
		text.set_path_effects([path_effects.Stroke(linewidth=lw, foreground="black"), path_effects.Normal()])

		##LDN 603
		#main_axes["gal"].text( 33.163 , +01.496 ,  "LDN 603",  color=fontColor, alpha=fontAlpha,fontsize=fontSize , horizontalalignment='center',    verticalalignment='center' )




		#####serpense NE cluster

		#2019AJ....157..200Z # 5 sources
		#G031.24-00.11(H)
		#main_axes["gal"].text(31.24, -0.11, "(2)",  color=fontColor, alpha=fontAlpha,fontsize=fontSize,horizontalalignment='center',   verticalalignment='center' )




		#G040.42+00.70(C)
		#main_axes["gal"].text(40.42,0.70, "(4)",  color=fontColor, alpha=fontAlpha,fontsize=fontSize ,horizontalalignment='center',   verticalalignment='center' )



		#G042.03+00.19(C)
		#main_axes["gal"].text(43.16, 0.01 , "(5)",  color=fontColor, alpha=fontAlpha,fontsize=fontSize ,horizontalalignment='center',   verticalalignment='center' )

		#G049.26+00.31(C)
		#main_axes["gal"].text(49.26, 0.31 , "(6)",  color=fontColor, alpha=fontAlpha,fontsize=fontSize ,horizontalalignment='center',   verticalalignment='center' )

 
		#####################################################
		#G035.19-00.74 2009ApJ...693..419Z
		#main_axes["gal"].text(35.19,-00.74, "(7)",  color=fontColor, alpha=fontAlpha,fontsize=fontSize ,horizontalalignment='center',   verticalalignment='center' )


		######################### two


		if 1:#Galactic longgitude
			xLocs = np.arange(26, 50,2)
			xLabels = map(int, xLocs)
			xLabels = map(str, xLabels)

			main_axes.set_ticklabel1_type("manual", locs=xLocs, labels=xLabels)

			main_axes.set_xlabel(r"Galactic Longitude ($^{\circ}$)")

		if 1:#Galactic latitudes
			xLocs = np.arange(-5,6,1)
			xLabels = map(int, xLocs)
			xLabels = map(str, xLabels)

			main_axes.set_ticklabel2_type("manual", locs=xLocs, labels=xLabels)

			main_axes.set_ylabel(r"Galactic Latitude ($^{\circ}$)")

		#at = AnchoredText(algDendro, loc=3, frameon=False)
		#axDendro.add_artist(at)
		#at = AnchoredText(r"$^{12}\mathrm{CO}~(J=1\rightarrow0)$", loc=4, frameon=True,  prop={"color": "black","size":13 },pad=0.1)
		#main_axes.add_artist(at )

		text=main_axes["gal"].text( 26 ,4.2,  r"$^{12}\mathrm{CO}~(J=1\rightarrow0)$",  color="black", alpha=fontAlpha,fontsize=8.3 , ha='center',    va='center',rotation=90 )
		#text.set_path_effects([path_effects.Stroke(linewidth=2, foreground="black"), path_effects.Normal()])





		fig.tight_layout()
		plt.savefig("localM0.pdf", bbox_inches='tight')
		plt.savefig("localM0.png", bbox_inches='tight', dpi=300)

	def drawLV(self,  vRange=[-6, 30]):
		"""
		drawPVDiagram of CO12
		:param LVFITS:s
		:return:
		"""

		from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
		import mpl_toolkits.axes_grid1.axes_grid as axes_grid

		#LVFITS= "/home/qzyan/WORK/myDownloads/MWISPcloud/G2650PV_DBMASK.fits"
		#use dendro mask, because we think dendro 2sigma, 8 pix, is the best way
		
		#LVFITS= "/home/qzyan/WORK/myDownloads/MWISPcloud/G2650PV_dendroMask.fits"
		LVFITS= "/home/qzyan/WORK/myDownloads/MWISPcloud/G2650PV_DBSCANS2P4CON1.fits"

		pvData, pvHead = doFITS.readFITS(LVFITS)

		wcs = WCS(pvHead)

		_, v0 = wcs.wcs_world2pix(40, vRange[0] * 1000, 0)

		_, v1 = wcs.wcs_world2pix(40, vRange[1] * 1000, 0)

		v0, v1 = map(round, [v0, v1])
		v0, v1 = map(int, [v0, v1])

		fig = plt.figure(figsize=(10, 1.35 ))
		#fig = plt.figure( )

		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'],"size":11.5,"weight":20})
		rc('text', usetex=True)
		#rc('mathtext', fontset='stixsans')
		#rc('text.latex', preamble=r'\usepackage{cmbright}')
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]
		#mpl.rcParams['text.usetex'] = True
		#mpl.rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
		#mpl.rcParams['font.family'] = 'scans-serif'
		#mpl.rcParams['font.scans-serif'] = 'cm'


		#rc('text.latex', preamble=r'\renewcommand{\familydefault}{\sfdefault}')
		#rc('text.latex', preamble=r'\renewcommand{\familydefault}{\sfdefault} \usepackage{helvet}')



		#mpl.rcParams['text.usetex'] = True
		#mpl.rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
		#mpl.rcParams['font.family'] = 'sans-serif'
		#mpl.rcParams['font.sans-serif'] = 'cm'


		# readfits

		# pvData=pvData*30./3600. #muliple degree

		grid = ImageGrid(fig, (1, 1, 1), nrows_ncols=(1, 1),
						 cbar_mode="each", cbar_pad="0.5%", cbar_size="1.8%",
						 cbar_location="right",
						 axes_class=(pywcsgrid2.Axes, dict(header=WCS(pvHead))))

		axPV = grid[0]

		aaaa=pvData[pvData>0]
		print np.min(aaaa),np.max(aaaa),"Maximum and minimum value"


		imCO12 = axPV.imshow(np.log(pvData), origin="lower", cmap="jet", aspect=1.5, vmin=np.log(0.01), vmax=np.log(2.2),
							 interpolation='none')  # draw color bar

		# axPV.set_facecolor('black')
		axPV.set_facecolor('silver')
		#axPV.set_xlabel(r"$l$($^{\circ}$)")

		cb_axes = grid.cbar_axes[0]

		#ticksP = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56]

		ticksP = [ 0.01,  0.1,  1]

		labels = map(str, ticksP)


		ticksPL = [0,0.1, 0.2, 0.4, 0.8, 1.6, 3.2]

		labelsL = map(str, ticksPL)

		#labels=[]
		#for eachV in ticksP:
			#labels.append( r"${:.1f}$".format(eachV)   )




		ticksP = np.log(ticksP)

		cb_axes.colorbar(imCO12, ticks=ticksP)
		cb_axes.set_yticklabels(labels)

		cb_axes.set_ylabel(r"K")

		if 1:  # velocity
			vRange = [0, 100]
			vInterval = 10  # km/s

			yLocs = np.arange(int(vRange[0]), int(vRange[1]) + vInterval, vInterval)

			yLabels = map(int, yLocs)
			yLabels = map(str, yLabels)

			axPV.set_ticklabel2_type("manual", locs=yLocs * 1000., labels=yLabels)
			axPV.set_ylabel(r"$V_{\rm LSR}$ (km s$^{-1}$)")

		if 1:#Galactic latitudes
			xLocs = np.arange(26, 50,2)
			xLabels = map(int, xLocs)
			xLabels = map(str, xLabels)

			axPV.set_ticklabel1_type("manual", locs=xLocs, labels=xLabels)
			axPV.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
			#axPV.set_xlabel(r"$l\left(^{\circ}\right)$")

		axPV.axis[:].major_ticks.set_color("w")
		axPV.set_ylim(v0, v1)
		fig.tight_layout()
		fig.savefig("G2650CO12LV.pdf", bbox_inches="tight" )

		fig.savefig("G2650CO12LV.png", bbox_inches="tight", dpi=600)


	def convertLBrangeToIndexRange(self,wcsCO,lRange,bRange):

		"""

		:param wcs:
		:param lRange:
		:param bRange:
		:return:
		"""



		startLIndex,startBIndex ,_  = wcsCO.wcs_world2pix( np.max( lRange ) , np.min(bRange) , 0, 0)

		endLIndex, endBIndex, _  = wcsCO.wcs_world2pix( np.min( lRange ) ,  np.max(bRange) , 0, 0)

		lIndexRange=[startLIndex,   endLIndex ]

		bIndexRange= [startBIndex,  endBIndex ]

		lIndexRange=map(round, lIndexRange)
		bIndexRange=map(round, bIndexRange)

		lIndexRange=map(int, lIndexRange)
		bIndexRange=map(int , bIndexRange)

		return lIndexRange, bIndexRange







	def drawCheckCloudsOneChannel(self,channelNumber=63):
		"""
		only take one channel from
		:param channelNumber:
		:return:
		"""

		# drawLrange,drawBrange=gaiaDis.box(38.9496172,0.1091115,3263.632 ,2991.628 ,4.9024796e-06)
		# region 2
		# drawLrange,drawBrange=gaiaDis.box(44.8346022,0.8519293,2361.857 ,2141.563 ,4.9024796e-06)
		#drawLrange, drawBrange = gaiaDis.box(42.8611667, 0.1834138, 3122.226, 2901.286, 4.9024796e-06)

		drawLrange, drawBrange = gaiaDis.box(38.0855621,-1.5002922,18000.000 , 14400.000 , 0)





		print drawLrange
		print drawBrange
		# vRange=[1.8,9.8] #km/s
		vRange = [-6, 6	]  # km/s

		# first crop fits

		rawCO ="G2650Local30.fits"

		# rawCO="G2650CO12MaskedCO.fits"

		# labelDendroFITS = "G2650minV3minP16_TrunkAsignMask0.fits"
		# labelDBSCANFITS = "DBCLEAN3.0_16Label.fits"
		# labelSCIMESFITS = "./scimesG2650/ClusterAsgn_3_16Ve20.fits"


		if 0:

			labelDendroFITS = "G2650minV2minP8_TrunkAsignMask0.fits"
			labelDBSCANFITS = "DBCLEAN2.0_8Label.fits"
			labelSCIMESFITS = "./scimesG2650/ClusterAsgn_2_8Ve20.fits"



			labelDendroData,labelDendroHead=myFITS.readFITS( labelDendroFITS  )
			labelDBSCANData,labelDBSCANHead=myFITS.readFITS( labelDBSCANFITS  )
			labelSCIMESData,labelSCIMESHead=myFITS.readFITS( labelSCIMESFITS  )


		if 1:
			labelCon1FITS = "/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/G2650LocaldbscanS2P4Con1_Clean.fits"
			labelCon2FITS = "/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/G2650LocaldbscanS2P8Con2_Clean.fits"
			labelCon3FITS = "/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/G2650LocaldbscanS2P11Con3_Clean.fits"

			labelDendroData, labelDendroHead = myFITS.readFITS(labelCon1FITS)
			labelDBSCANData, labelDBSCANHead = myFITS.readFITS(labelCon2FITS)
			labelSCIMESData, labelSCIMESHead = myFITS.readFITS(labelCon3FITS)

		labelDendro = labelDendroData[channelNumber]
		labelDBSCAN =  labelDBSCANData[channelNumber]
		labelSCIMES =  labelSCIMESData[channelNumber]

		WCSCrop =WCS(labelDendroHead)

		lIndexRange,bIndexRange = self.convertLBrangeToIndexRange(WCSCrop, drawLrange, drawBrange )

		#get velocity range

		_,_,channelV = WCSCrop.wcs_pix2world(0,0,channelNumber,0)
		channelV=channelV/1000.


		fig = plt.figure(1, figsize=(10, 7.5) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'], "size" :13 })
		# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

		rc('text', usetex=True)
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]


		axCO= pywcsgrid2.subplot(221, header=   WCSCrop  )


		dataCO ,headCO =myFITS.readFITS(rawCO)

		intData = dataCO[channelNumber]    #np.sum(dataCO ,axis=0 ) *0.2

		cmapCO = plt.cm.bone

		cmapCO.set_bad('black' )

		axCO.imshow( np.sqrt( intData ) ,origin='lower' ,cmap=cmapCO ,vmin=0 ,vmax=3, interpolation='none')

		at = AnchoredText(r"$^{{12}}\mathrm{{CO}}~(J=1\rightarrow0)$, $V_{{\rm LSR}}={}\ \rm km\ s^{{-1}}$".format(channelV) , loc=3, frameon=False,
						  prop={"color": "w", "size": 12})
		axCO.add_artist(at)

		axCO.set_ticklabel_type("absdeg", "absdeg")
		axCO.axis[:].major_ticks.set_color("w")
		axCO.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
		axCO.set_ylabel(r"Galactic Latitude ($^{\circ}$)")

		###dendrogram

		axDendrogram = pywcsgrid2.subplot(222, header=WCSCrop, sharex=axCO, sharey=axCO)
		# self.showLabels(axDendrogram,labelDendro  )

		self.showLabels(axDendrogram, labelDendro)
		# axDendrogram.imshow(labelDendro,origin='lower',cmap="jet", interpolation='none')

		at = AnchoredText("Connectivity 1", loc=3, frameon=False, prop={"color": "w", "size": 12})
		axDendrogram.add_artist(at)
		axDendrogram.set_ticklabel_type("absdeg", "absdeg")
		axDendrogram.axis[:].major_ticks.set_color("w")
		axDendrogram.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
		axDendrogram.set_ylabel(r"Galactic Latitude ($^{\circ}$)")

		###DBSCAN
		#labelDBSCAN, labelHead = self.getIntLabel(cropDBSCAN)

		axDBSCAN = pywcsgrid2.subplot(223, header=WCSCrop, sharex=axCO, sharey=axCO)
		self.showLabels(axDBSCAN, labelDBSCAN)
		# axDBSCAN.imshow(labelDBSCAN,origin='lower',cmap="jet", interpolation='none')

		at = AnchoredText("Connectivity 2", loc=3 , frameon=False, prop={"color": "w", "size": 12})
		axDBSCAN.add_artist(at)
		axDBSCAN.set_ticklabel_type("absdeg", "absdeg")
		axDBSCAN.axis[:].major_ticks.set_color("w")

		axDBSCAN.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
		axDBSCAN.set_ylabel(r"Galactic Latitude ($^{\circ}$)")



		##SCIMES
		#labelSCIMES, labelHead = self.getIntLabel(cropSCIMES)

		axSCIMES = pywcsgrid2.subplot(224, header=WCSCrop, sharex=axCO, sharey=axCO)
		self.showLabels(axSCIMES, labelSCIMES)

		# axSCIMES.imshow(labelSCIMES,origin='lower',cmap="jet", interpolation='none')

		at = AnchoredText("Connectivity 3", loc=3 , frameon=False, prop={"color": "w", "size": 12})
		axSCIMES.add_artist(at)
		axSCIMES.set_ticklabel_type("absdeg", "absdeg")
		axSCIMES.axis[:].major_ticks.set_color("w")


		axSCIMES.set_xlim(lIndexRange)
		axSCIMES.set_ylim(bIndexRange)

		axSCIMES.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
		axSCIMES.set_ylabel(r"Galactic Latitude ($^{\circ}$)")




		fig.tight_layout(pad=0.2)

		plt.savefig("checkCloudChannel{}.pdf".format(channelNumber ) , bbox_inches="tight")
		plt.savefig("checkCloudChannel{}.png".format( channelNumber) , bbox_inches="tight", dpi=600)




	def drawCheckClouds(self):
		"""
		compare the result of molecular clouds
		:return:
		"""
		#drawLrange,drawBrange=gaiaDis.box(38.9496172,0.1091115,3263.632 ,2991.628 ,4.9024796e-06)
		#region 2
		#drawLrange,drawBrange=gaiaDis.box(44.8346022,0.8519293,2361.857 ,2141.563 ,4.9024796e-06)
		drawLrange,drawBrange=gaiaDis.box(42.8611667,0.1834138,3122.226 ,2901.286 ,4.9024796e-06)

		print drawLrange
		print drawBrange
		#vRange=[1.8,9.8] #km/s
		vRange = [ -6,  6	 ]  # km/s

		#first crop fits

		rawCO="G2650Local30.fits"
		
		#rawCO="G2650CO12MaskedCO.fits"

		
		#labelDendroFITS = "G2650minV3minP16_TrunkAsignMask0.fits"
		#labelDBSCANFITS = "DBCLEAN3.0_16Label.fits"
		#labelSCIMESFITS = "./scimesG2650/ClusterAsgn_3_16Ve20.fits"

		labelDendroFITS = "G2650minV2minP8_TrunkAsignMask0.fits"
		labelDBSCANFITS = "DBCLEAN2.0_8Label.fits"
		labelSCIMESFITS = "./scimesG2650/ClusterAsgn_2_8Ve20.fits"





		cropRawcoFITS =  self.tmpPath + "cropRawco.fits"
		cropDendroFITS = self.tmpPath + "cropDendro.fits"
		cropSCIMES = self.tmpPath + "cropScimes.fits"
		cropDBSCAN=  self.tmpPath + "cropDbscan.fits"


		doFITS.cropFITS(rawCO,outFITS=cropRawcoFITS,Vrange=vRange,Brange=drawBrange,Lrange=drawLrange , overWrite=True )

		#cropFITS 3D
		doFITS.cropFITS(labelDendroFITS,outFITS=cropDendroFITS,Vrange=vRange,Brange=drawBrange,Lrange=drawLrange , overWrite=True )
		doFITS.cropFITS(labelDBSCANFITS,outFITS=cropDBSCAN,Vrange=vRange,Brange=drawBrange,Lrange=drawLrange , overWrite=True )
		

		doFITS.cropFITS(labelSCIMESFITS,outFITS=cropSCIMES,Vrange=vRange,Brange=drawBrange,Lrange=drawLrange , overWrite=True )


		#check uniqueness along each points
		if 1:
			self.checkUniqueness(cropDendroFITS )
			self.checkUniqueness(cropDBSCAN )
			self.checkUniqueness(cropSCIMES )

		labelDendro,labelHead=self.getIntLabel(cropDendroFITS )

		WCSCrop=WCS(labelHead)

		fig = plt.figure(1, figsize=(10, 9) )
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'], "size":13 })
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

		rc('text', usetex=True)

		axCO= pywcsgrid2.subplot(221, header=   WCSCrop  )


		dataCO,headCO=myFITS.readFITS(cropRawcoFITS)

		intData=np.sum(dataCO,axis=0)*0.2

		axCO.imshow(np.sqrt( intData),origin='lower',cmap="bone",vmin=0 ,vmax=4, interpolation='none')

		at = AnchoredText(r"$^{12}\mathrm{CO}~(J=1\rightarrow0)$", loc=1, frameon=False,  prop={"color": "w","size":13 })
		axCO.add_artist(at)

		axCO.set_ticklabel_type("absdeg", "absdeg")
		axCO.axis[:].major_ticks.set_color("w")


		###dendrogram

		axDendrogram= pywcsgrid2.subplot(222, header=   WCSCrop ,sharex=axCO,sharey=axCO )
		#self.showLabels(axDendrogram,labelDendro  )

		self.showLabels(axDendrogram,labelDendro  )
		#axDendrogram.imshow(labelDendro,origin='lower',cmap="jet", interpolation='none')


		at = AnchoredText("Dendrogram", loc=1, frameon=False,  prop={"color": "w","size":13 })
		axDendrogram.add_artist(at)
		axDendrogram.set_ticklabel_type("absdeg", "absdeg")
		axDendrogram.axis[:].major_ticks.set_color("w")

		###DBSCAN
		labelDBSCAN,labelHead=self.getIntLabel(cropDBSCAN )

		axDBSCAN = pywcsgrid2.subplot(223, header=WCSCrop, sharex=axCO, sharey=axCO)
		self.showLabels(axDBSCAN, labelDBSCAN)
		#axDBSCAN.imshow(labelDBSCAN,origin='lower',cmap="jet", interpolation='none')

		
		at = AnchoredText("DBSCAN", loc=1, frameon=False,  prop={"color": "w","size":13 })
		axDBSCAN.add_artist(at)
		axDBSCAN.set_ticklabel_type("absdeg", "absdeg")
		axDBSCAN.axis[:].major_ticks.set_color("w")


		##SCIMES
		labelSCIMES,labelHead=self.getIntLabel(cropSCIMES )

		axSCIMES = pywcsgrid2.subplot(224, header=WCSCrop, sharex=axCO, sharey=axCO)
		self.showLabels(axSCIMES, labelSCIMES)
		
		#axSCIMES.imshow(labelSCIMES,origin='lower',cmap="jet", interpolation='none')

		
		at = AnchoredText("SCIMES", loc=1, frameon=False, prop={"color": "w","size":13 })
		axSCIMES.add_artist(at)
		axSCIMES.set_ticklabel_type("absdeg", "absdeg")
		axSCIMES.axis[:].major_ticks.set_color("w")



		fig.tight_layout(pad=0.2)
		plt.savefig("checkCloud.pdf", bbox_inches="tight")

		plt.savefig("checkCloud.png", bbox_inches="tight",dpi=600)




	def getIntLabel(self,cropLabelFITS):
		labelDendroData, labelHead = doFITS.readFITS(cropLabelFITS)

		minVDendro = np.nanmin(labelDendroData)

		uniqueValues = np.unique(labelDendroData)
		# print 0 in uniqueValues

		labelDendroData[labelDendroData == minVDendro] = np.NaN

		meanLabel = np.nanmean(labelDendroData, axis=0)

		#meanLabel=np.nan_to_num(meanLabel)

		return meanLabel,labelHead



	def checkUniqueness(self,cropDendroFITS):

		labelDendroData,_ = doFITS.readFITS( cropDendroFITS )

		minVDendro = np.nanmin( labelDendroData )

		uniqueValues=np.unique(labelDendroData)
		#print 0 in uniqueValues
		 

		labelDendroData[ labelDendroData == minVDendro ] = np.NaN


		meanLabel=np.nanmean( labelDendroData,axis=0 )

		#meanLabel=np.nan_to_num(meanLabel)

		uniqueValuesV2=np.unique( meanLabel )
		uniqueValuesV2=uniqueValuesV2[uniqueValuesV2>minVDendro-1]


		print len( uniqueValues ), "==",  len( uniqueValuesV2 )+1








	#def


	def getVdispersion(self,testFITS,testTB,regionName="",vCenPix=False,saveSpectral=False,savePath=None,saveTBAs=None ):
		"""
		The weighted velooicty may be be seen as the veloicty dispersion
		:return:
		"""
		#testFITS="DBCLEAN3.0_16Label.fits"

		#testTB=Table.read("DBCLEAN3.0_16TB.fit")

		COFITS=  "G2650Local30.fits"
		dataLabel,headLabel = myFITS.readFITS(  testFITS )

		dataCO,headCO = myFITS.readFITS(  COFITS   )

		cubeCO=SpectralCube.read(COFITS )

		vAxis=cubeCO.spectral_axis
		#print vAxis.real.to_value()
		wcsCO=WCS(headCO)

		vAxis= vAxis.value/1000. #convert to rms

		noiseV=np.nanmin( dataLabel[0] )
		index1D=np.where(dataLabel > noiseV)
		values1D=dataLabel[index1D]

		Z0,Y0,X0 =  index1D

		dataZero=np.zeros_like(dataLabel)

		fitter = modeling.fitting.LevMarLSQFitter()

		lineWdith=[]


		widgets = ['Geting line width: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),  ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
		pbar = ProgressBar(widgets=widgets, maxval=len(testTB))
		pbar.start()

		try:
			testTB["lineWidth"]=testTB["v_rms"]
		except:
			pass

		i=0
		for eachR in testTB:


			testID=  int(   eachR["_idx"] )



			testIndices=self.getIndices(Z0,Y0,X0,values1D,testID)
			singleZ0,singleY0,singleX0=testIndices

			dataZero[testIndices]= dataCO[testIndices]

			#cropThe cloudRange
			minY=np.min(singleY0  )
			maxY=np.max(singleY0  )
			###########
			minX=np.min(singleX0  )
			maxX=np.max(singleX0  )

			###########
			minZ=np.min( singleZ0  )
			maxZ=np.max( singleZ0  )

			#########

			cloudCropSpectra=dataZero[:,minY:maxY+1,minX:maxX+1]

			cloudCropCube=dataZero[minZ:maxZ+1,minY:maxY+1,minX:maxX+1]


			averageSpectraCrop= np.nansum( cloudCropSpectra,axis=(1,2) )


			intCloud =  np.nansum( cloudCropCube,axis=0 )

			#count the number spectra

			totalSpectral=len( intCloud[intCloud>0] )

			meanSpectral = averageSpectraCrop/1./totalSpectral

			if saveSpectral:
				savefileName=savePath+"{}_{}Spectral".format( regionName , testID)

				np.save( savefileName, [  vAxis , meanSpectral  ]  )
			#model = modeling.models.Gaussian1D(  amplitude=np.max(vAxis), mean=vCen  ,stddev= testRow["v_rms"].data[0]  )  # depending on the data you need to give some initial values
			#fitted_model = fitter(model, vAxis, meanSpectral)

			#diff=  fitted_model.stddev.value - testRow["v_rms"].data[0]

			#print  diff, testID,  fitted_model.stddev.value

			spectraPeak= np.max( meanSpectral  )

			area=(vAxis[1]-vAxis[0])*np.sum( meanSpectral )

			eqLineWidth= area/spectraPeak
			dataZero[testIndices]=  0

			eachR["lineWidth"] = eqLineWidth


			lineWdith.append( eqLineWidth  )

			i = i+1

			pbar.update(i)
			#plot
			#if abs(diff) >0.5 :
			if  0 :
				figVVV = plt.figure(figsize=(12, 6))
				ax = figVVV.add_subplot(1, 1, 1)
				ax.plot(vAxis,meanSpectral)
				ax.plot(vAxis, fitted_model(vAxis),'r--',label="Diff: {:.3f} km/s".format(diff) )

				ax.legend()

				# fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
				rc('text', usetex=True)
				rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})
				plt.savefig("./spectralTestFigures/testSpectral{}_{}.png".format(regionName,testID) , bbox_inches='tight')
				#print dir(ax)
				#ax.cla()
				#
		pbar.finish()

		if saveTBAs!=None:
			testTB.write(saveTBAs,overwrite=True)

		return np.asarray( lineWdith )

	def getExamineSpectral(self,savePath="./spectralLineOfAll/" ):

		#savePath="./spectralLineOfAll/" #this path saves a tb file, and alll spectral for each clouds

		allFiles=glob.glob(savePath+"*.npy")

		figVVV = plt.figure(figsize=(10, 6))
		ax = figVVV.add_subplot(1, 1, 1)
		rc('text', usetex=True)
		rc('font', **{'family': 'sans-serif', 'size': 14, 'serif': ['Helvetica']})



		#should add spectral here



		for eachF in allFiles:
			pass
			saveFigureName=eachF[:-4]+".png"

			#get cloud ID
			currentID= int(eachF.split("_")[1].split("Spectral.npy")[0])


			velocity, tempearature = np.load( eachF )


			ax.step(velocity, tempearature, color='b', label="cloud ID: {}".format(currentID))
			#ax.plot(velocity, fitted_model(vAxis), 'r--', label="cloud ID: ".format(currentID))

			ax.legend()

			ax.set_xlabel(r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)")
			ax.set_ylabel(r"Brightness temperature ($\rm K$)")

			plt.savefig( saveFigureName   , bbox_inches='tight')
			# print dir(ax)
			ax.cla()



			#print saveFigureName

			#break

	def getVelTemByID(self,cloudID,savePath):
		files=glob.glob(savePath+"*_{}*.npy".format(cloudID))

		return np.load(files[0] )



	def drawSpectraExample(self):

		savePath="./spectralLineOfAll/" #this path saves a tb file, and alll spectral for each clouds
		#cloudTB=Table.read( "/home/qzyan/WORK/myDownloads/MWISPcloud/spectralLineOfAll/sourceTB.fit" )

		cloudTB = Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/G2650LocaldbscanS2P4Con1_CleanWithLW.fit")






		data,header=myFITS.readFITS("G2650Local30.fits")
		wcsCO=WCS(header)

		#draw nine cases
		#18225, 4components
		#16417, one strong, one weak, 2 components, 7864?5135
 		#12580, two equally strong componennts ,close, 6305
		#11132, more than 4 components,9530

		#16171, three componennts,7703 4288
		#11810, very large velocity range

		# 22547, small, but multiple components

		#15497 ,#very small dispersion  gaussion, 7111
		# 11505,#large velocity dispersion? 8009?

		# 9218, line wing, like outflow 1876, 227, very long wind
		# 5401, multicomponents, with small veloicty dispersion
		#8195, single line winge,6898,5595
		#3753, no clear componnets shown

		###########draw cases
		#1. 18202, 1 componentns
		#2. 7864, 2 componnents, seperatted
		#3, 16171, 3 componennts,
		#4, 18225, 4 componnents

		#5, 11132, morethan 4 componennts

		#6 3753 #  unregular

		#7, 11505, large velocity dispersion
		#8, 9218, line wine
		#9, 5595, single wind

		dv= 18
		fig = plt.figure(figsize=(18, 15 ))

		rc('text', usetex=True)
		rc('font', **{'family': 'sans-serif', 'size': 22, 'serif': ['Helvetica']})
		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]



		yLabelStr=  r"Brightness temperature ($\rm K$)"
		xLabelStr= r"$V_{\rm LSR}$ ($\rm km\ s^{-1}$)" 

		######draw Spectral 1
		ax0 = fig.add_subplot(3, 3, 1)
		self.drawOneCloud(ax0,133683,dv, cloudTB,wcsCO,first=True) #cloud ID: 18202, 1 componentns
		ax0.set_ylabel( yLabelStr )


		######draw Spectral 2
		ax1 = fig.add_subplot(3, 3, 2 ,sharey=ax0)
		self.drawOneCloud(ax1,74029,dv, cloudTB,wcsCO) #cloud ID: 18202, 1 componentns
		ax1.set_ylim(0, 1.8)

		######draw Spectral 3
		ax1 = fig.add_subplot(3, 3, 3 ,sharey=ax0)
		self.drawOneCloud(ax1,107557 ,dv, cloudTB,wcsCO) #cloud ID: 18202, 1 componentns



		######draw Spectral 4
		ax1 = fig.add_subplot(3, 3, 4 ,sharey=ax0)
		self.drawOneCloud(ax1,125004,dv, cloudTB,wcsCO) #cloud ID: 18202, 1 componentns
		ax1.set_ylabel( yLabelStr )
		######draw Spectral 5
		ax1 = fig.add_subplot(3, 3, 5 ,sharey=ax0)
		self.drawOneCloud(ax1,87629,dv, cloudTB,wcsCO) #cloud ID: 18202, 1 componentns


		####################

		######draw Spectral 6
		ax1 = fig.add_subplot(3, 3, 6 ,sharey=ax0)
		self.drawOneCloud(ax1,64578,dv, cloudTB,wcsCO) #cloud ID: 18202, 1 componentns

		######draw Spectral 7
		ax1 = fig.add_subplot(3, 3, 7 ,sharey=ax0)
		self.drawOneCloud(ax1,128307,dv, cloudTB,wcsCO) #cloud ID: 18202, 1 componentns
		ax1.set_ylabel( yLabelStr )
		ax1.set_xlabel( xLabelStr )
		######draw Spectral 8
		ax1 = fig.add_subplot(3, 3, 8 ,sharey=ax0)
		self.drawOneCloud(ax1,113981,dv, cloudTB,wcsCO) #cloud ID: 18202, 1 componentns
		ax1.set_xlabel( xLabelStr )
		######draw Spectral 9
		ax1 = fig.add_subplot(3, 3, 9 ,sharey=ax0)
		self.drawOneCloud(ax1,53488,dv, cloudTB,wcsCO) #cloud ID: 18202, 1 componentns
		ax1.set_xlabel( xLabelStr )



		fig.tight_layout()

		plt.savefig("spectraExamples.pdf", bbox_inches='tight')
		plt.savefig("spectraExamples.png", bbox_inches='tight')

	def drawOneCloud(self,ax1,cloud1ID,dv,cloudTB,wcsCO,first=False):
		#cloud1ID=18202 #case 1,

		#savePath="./spectralLineOfAll/" #this path saves a tb file, and alll spectral for each clouds
		savePath="./spectraSave/" #this path saves a tb file, and alll spectral for each clouds



		vCenRow="v_cen"
		xCenRow="x_cen"
		yCenRow="y_cen"
		
		cloud1Row= cloudTB[cloudTB["_idx"] == cloud1ID ][0]
		vCen= cloud1Row[ vCenRow ]


		xCen=   cloud1Row[ xCenRow ]
		yCen=  cloud1Row[ yCenRow ]

		#xCen,yCen,vCen= wcsCO.wcs_pix2world( xCen,yCen,vCen,0 )
		#vCen=vCen/1000.
		V1,T1=self.getVelTemByID(cloud1ID,savePath)

		ax1.step(V1, T1, color='b' )
		if first:
			at = AnchoredText(r"$\left(l,\ b\right)=\left({:.2f}^\circ, \ {:.2f}^\circ\right)$".format( xCen,yCen), loc=1, frameon=False)
		else:
			at = AnchoredText(r"$\left({:.2f}^\circ, \ {:.2f}^\circ\right)$".format( xCen,yCen), loc=1, frameon=False)

		ax1.add_artist(at)



		ax1.set_xlim(vCen-dv/2.,vCen+dv/2.)


	def maskEdges(self,data,minV):
		"""
		mask the data, tha is one the edge
		:param data:
		:return:
		"""
		data[:, 1063:, 0: 55]=minV #left up corner

		data[:, 1003: , 2815:]=minV


	def getFalseRate(self,labelFITS,noiseThreshold=0.5  ):



		"""
		fits bad pixels in the labelled fits, to get the false rate
		the noise threshold is 0.5 K, 1 simga
		:param label1FITS:
		:return:
		"""
		dataCluster,headCluster=myFITS.readFITS(labelFITS)

		noiseMask=np.nanmin( dataCluster[0] )

		self.maskEdges( dataCluster ,  noiseMask )

		dataCO, headCO=myFITS.readFITS(self.rawCOFITS)

		dataCOFront=dataCO.copy()

		#find pixels, that in spectra are isolated  channels, which could be seen as noise

		signalMask= dataCluster>noiseMask

		dataCOFront[0]=False
		dataCOFront[-1]=False #remove the edge

		signalMask[0] = False
		signalMask[-1] = False



		#fits, get those pixels, that the pixel in front of it is less than the noise threshold
		coFrontToCurrent= np.roll( dataCOFront,1,axis=0  )
		coBackToCurrent=np.roll( dataCOFront,-1,axis=0  )



		frontGoodPixels= coFrontToCurrent>=noiseThreshold
		coFrontToCurrent[frontGoodPixels]  = 0
		coFrontToCurrent[~frontGoodPixels] = 1


		backGoodPixels = coBackToCurrent>=noiseThreshold
		coBackToCurrent[ backGoodPixels ] = 0 #to collect, how many pixels are less tha the threshold, if both front and back are less than the threshold, it looks like this is a noise pixels
		coBackToCurrent[ ~backGoodPixels] = 1

		coFrontToCurrent[ ~signalMask ] = False
		coBackToCurrent[ ~signalMask ]  = False



		badCube=np.logical_and(coBackToCurrent,coFrontToCurrent)

		#badCube[badCube<2]=0

		#badCube[badCube==2]=1

		badVoxelNumber= np.sum( badCube ) #number of bad pixels
		totalVoxelNumber=np.sum( signalMask)

		print "Bad voxels number {}, total voxelsNumber {}, ratio: {:f}".format(badVoxelNumber, totalVoxelNumber,badVoxelNumber/1./totalVoxelNumber  )
		#






	def getCorrectRate(self,labelFITS,signalThreshold=1.5):



		"""
		fits bad pixels in the labelled fits, to get the  correct
		# at least one point adjacent has to be larger than 3 sigma, which is 1.5 K
		#
		#
		the noise threshold is 0.5 K, 1 simga
		:param label1FITS:
		:return:
		"""
		dataCluster,headCluster=myFITS.readFITS(labelFITS)

		noiseMask=np.nanmin( dataCluster[0] )

		self.maskEdges( dataCluster ,  noiseMask )

		dataCO, headCO=myFITS.readFITS(self.rawCOFITS)
		dataCOFront=dataCO.copy()

		#find pixels, that in spectra are isolated  channels, which could be seen as noise

		signalMask= dataCluster>noiseMask

		dataCOFront[0]=False
		dataCOFront[-1]=False #remove the edge

		signalMask[0] = False
		signalMask[-1] = False



		#fits, get those pixels, that the pixel in front of it is less than the noise threshold
		coFrontToCurrent= np.roll( dataCOFront,1,axis=0  )
		coBackToCurrent=np.roll( dataCOFront,-1,axis=0  )



		frontGoodPixels= coFrontToCurrent >= signalThreshold
		coFrontToCurrent[frontGoodPixels]  = 1
		coFrontToCurrent[~frontGoodPixels] = 0


		backGoodPixels = coBackToCurrent >= signalThreshold
		coBackToCurrent[ backGoodPixels ] = 1 #to collect, how many pixels are less tha the threshold, if both front and back are less than the threshold, it looks like this is a noise pixels
		coBackToCurrent[ ~backGoodPixels] = 0



		coFrontToCurrent[ ~signalMask ] = False
		coBackToCurrent[ ~signalMask ]  = False



		goodCube=np.logical_or(coBackToCurrent,coFrontToCurrent)

		#badCube[badCube<2]=0

		#badCube[badCube==2]=1

		goodVoxelNumber= np.sum( goodCube ) #number of bad pixels
		totalVoxelNumber=np.sum( signalMask)

		print "Good voxels number {}, total voxelsNumber {}, ratio: {:f}".format(goodVoxelNumber, totalVoxelNumber,goodVoxelNumber/1./totalVoxelNumber  )
		#



	def produceSingleCubesForEachCloud(self, labelsFITS , cloudTBFile ):
		"""

		#output all data cubes for each cloud

		:return:
		"""

		#################

		savePath="./cloudSubCubes/"

		cloudTB=Table.read(cloudTBFile)

		cloudTB=self.removeWrongEdges(cloudTB)
		print len( cloudTB),"Number of molecular clouds"

		dataCluster,headCluster=myFITS.readFITS( labelsFITS )
		dataCO,headCO=myFITS.readFITS(self.rawCOFITS)
		#print cloudTB

		minV = np.nanmin(dataCluster[0])
		wcsCloud = WCS( headCluster )
		clusterIndex1D = np.where( dataCluster>minV )
		clusterValue1D =  dataCluster[clusterIndex1D ]
		Z0,Y0,X0=clusterIndex1D

		fitsZero=np.zeros_like(dataCluster)
		#print cloudTB.colnames
		for eachC in cloudTB:


			cloudID=eachC["_idx"]
			saveName= "cloud{}cube.fits".format( cloudID )



			cloudIndex = self.getIndices(Z0, Y0, X0, clusterValue1D, cloudID)
			fitsZero[cloudIndex]= dataCO[cloudIndex]

			cloudZ0,cloudY0,cloudX0 = cloudIndex

			minZ = np.min( cloudZ0 )
			maxZ = np.max( cloudZ0 )

			minY = np.min( cloudY0 )
			maxY = np.max( cloudY0 )

			minX = np.min( cloudX0 )
			maxX = np.max( cloudX0 )

			cropWCS=wcsCloud[minZ:maxZ+1, minY:maxY+1, minX:maxX+1    ]

			cropData=fitsZero[minZ:maxZ+1, minY:maxY+1, minX:maxX+1    ]


			fits.writeto(savePath+ saveName , cropData, header=cropWCS.to_header(),overwrite=True)

			fitsZero[cloudIndex]=0



	def getPublicCatalog(self):
		"""

		:param TBFileName:
		:return:
		"""



		#tb=Table.read(TBFileName)
		#better use

		rawCatalog="minV2minP8dendroMannualCat_LineWidth.fit"
		#rawCatalog="minV2minP8_dendroCatTrunk.fit"

		
		saveFileName="cloudCatalogDendrgram.fit"

		tbRaw = Table.read( rawCatalog )
		
		print len(tbRaw)
		tbRaw=self.removeWrongEdges(tbRaw)

		print len(tbRaw)

		publishCat = Table()
		tbRaw["name"] = tbRaw["_idx"].astype(str)
		#add cloudName to the row
		#print tbRaw
		for eachR in tbRaw:

			l = np.float(eachR["x_cen"])
			b = np.float(eachR["y_cen"])
			v = np.float(eachR["v_cen"])    # to km/s
			#print l,b,v


			cloudName=self.getCloudNameByLB(l,b)

			eachR["name"]=cloudName


		print tbRaw.colnames
		#addName
		publishCat.add_column( tbRaw["name"] )

		#add l, b, v, values
		colL = tbRaw["x_cen"]
		colL.name="l"
		colL.unit= u.deg
		publishCat.add_column( colL )


		colb = tbRaw["y_cen"]
		colb.name="b"
		colb.unit= u.deg
		publishCat.add_column( colb )

		colv = tbRaw["v_cen"]
		colv.name="Vlsr"
		colv.unit= u.km/u.s
		publishCat.add_column( colv )

		#######################################

		#lbv, dispersion

		colLrms = tbRaw["l_rms"]
		colLrms.name="l_sigma"
		colLrms.unit= u.deg
		publishCat.add_column( colLrms )

		colBrms = tbRaw["b_rms"]
		colBrms.name = "b_sigma"
		colBrms.unit = u.deg
		publishCat.add_column(colBrms)

		colVrms = tbRaw["lineWidth"]
		#colBrms.name = "b_sigma"
		colVrms.unit = u.km/u.s
		publishCat.add_column(colVrms)


		colVoxN=tbRaw["pixN"].astype(int)
		colVoxN.name="voxelN"
		colVoxN.unit=None
		publishCat.add_column(colVoxN)

		colPeak = tbRaw["peak"]
		#colBrms.name = "b_sigma"
		colPeak.unit = u.K
		publishCat.add_column(colPeak)

		##############
		colArea = tbRaw["area_exact"]
		#colBrms.name = "b_sigma"
		#colArea.unit = u.K
		publishCat.add_column(colArea)

		#########################################

		colFlux = tbRaw["sum"]*0.2*0.25
		colFlux.name = "flux"
		colFlux.unit = u.K*u.km/u.s*u.arcmin**2 #*u.def_unit("Omega A")
		publishCat.add_column(colFlux)






		publishCat.write( saveFileName, overwrite=True  )

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


	def contourDisClouds(self,fitsFile="normalInt.fits"):


		# what is the unit of this integrated intensity?
		import matplotlib
		from mpl_toolkits.axes_grid1.axes_grid import ImageGrid
		from mpl_toolkits.axes_grid1.axes_grid import AxesGrid
		from matplotlib.colors import LogNorm

		dataCO, headCO = myFITS.readFITS(fitsFile)

		wcsCO = WCS(headCO)

		fig = plt.figure(1, figsize=(16, 7.8))
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'], 'size': 20})
		# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

		rc('text', usetex=True)

		# grid helper
		grid_helper = pywcsgrid2.GridHelper(wcs=wcsCO)

		# AxesGrid to display tow images side-by-side
		fig = plt.figure(1, (6, 3.5))

		grid = ImageGrid(fig, (1, 1, 1), nrows_ncols=(1, 1),
						 cbar_mode="single", cbar_pad="0.3%",
						 cbar_location="right", cbar_size="1%",
						 axes_class=(pywcsgrid2.Axes, dict(header=wcsCO)))

		main_axes = grid[0]
		main_axes.locator_params(nbins=10)
		cb_axes = grid.cbar_axes[0]  # colorbar axes

		a = np.log(dataCO)

		a = a[a > 0]
		print np.min(a), np.max(a), "minimum and maxinum after "

		# im = main_axes.imshow(np.sqrt(dataCO),origin='lower',cmap="jet",  vmin=np.sqrt(1.9*1.5), vmax= np.sqrt( 81) , interpolation='none')

		# im = main_axes.imshow( dataCO ,origin='lower',cmap="jet",  vmin=-10, vmax= 80 , interpolation='none')

		# im = main_axes.imshow(np.log(dataCO),origin='lower',cmap="jet",   interpolation='none')

		im = main_axes.imshow(np.log(dataCO), origin='lower', cmap="bone", vmin=np.log(1.9  ), vmax=np.log(110),
							  interpolation='none')

		main_axes.set_facecolor('black')
		main_axes.axis[:].major_ticks.set_color("white")


		# print dir(cb_axes.axes)

		# tickesArray=np.asarray( [0,0.1,0.5,1,2,3,4,5] )
		tickesArray = np.asarray([5, 10, 20, 40, 80])
		# tickesArray=tickesArray**2
		#cb.ax.set_yticks(np.log(tickesArray))
		#cb.ax.set_yticklabels(map(str, tickesArray))

		#contour good Clouds
		import matplotlib as mpl
		goodCloudTB="/home/qzyan/WORK/projects/maddalena/dendroDisPath/G2650/G2650goodDisTB.fit"
		goodCloudTB=Table.read( goodCloudTB )
		import matplotlib.cm as cm
		norm = mpl.colors.Normalize(vmin=0.2, vmax=1.5)
		cmap = cm.jet
		m = cm.ScalarMappable(norm=norm, cmap=cmap)


		for eachC in  goodCloudTB:
			#pass
			cloudName= eachC["sourceName"]

			cloudID=int( cloudName.split("oud")[1])
			maskFITS=glob.glob("/home/qzyan/WORK/myDownloads/testScimes/G2650Formal/cloudInt/Cloud{}_mask.fits".format(cloudID))[0]

			dataMask,headMask=myFITS.readFITS( maskFITS )
			distance=eachC["distance"]/1000.
			cloudColor = m.to_rgba( distance)

			c=main_axes.contour(dataMask, colors=[cloudColor], linewidths=0.6, origin="lower", levels=[1])

			#trueCloud Name

			main_axes["gal"].text(eachC["l"],  eachC["b"],  eachC["Note"],color=cloudColor,fontsize=8 ,horizontalalignment='center',   verticalalignment='center' )

		cb = mpl.colorbar.ColorbarBase(cb_axes,norm=norm,cmap=cmap)

		#cb_axes = fig.colorbar(c, ax=main_axes, cmap=cmap, norm=norm, pad=0.01)
		# cb_axes.axis["right"].toggle(ticklabels=True)
		cb_axes.set_ylabel("Distance (kpc)")
		#cb_axes.set_xlabel("")


		fig.tight_layout()
		plt.savefig("contourDisClouds.pdf", bbox_inches='tight')
		plt.savefig("contourDisClouds.png", bbox_inches='tight', dpi=300)


	def testAngularAndLineWith(self):
		"""
		
		:return:
		"""
		algDendro="Dendrogram"
		tb8Den,tb16Den,label8Den,label16Den,sigmaListDen=self.getTBList(algorithm=algDendro)
		tb8Den=self.removeAllEdges(tb8Den)
		tb16Den=self.removeAllEdges(tb16Den)


		algDB="DBSCAN"
		tb8DB,tb16DB,label8DB,label16DB,sigmaListDB=self.getTBList(algorithm=algDB)
		tb8DB=self.removeAllEdges(tb8DB)
		tb16DB=self.removeAllEdges(tb16DB)

		algSCI="SCIMES"
		tb8SCI,tb16SCI,label8SCI,label16SCI,sigmaListSCI=self.getTBList(algorithm= algSCI )
		tb8SCI=self.removeAllEdges(tb8SCI)
		tb16SCI=self.removeAllEdges(tb16SCI)

		fig = plt.figure(figsize=(12, 6))
		ax1 = fig.add_subplot(1, 2, 1)
		# fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
		rc('text', usetex=True)
		rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})

		ax1.scatter( tb8Den[0]["area_exact"],  tb8Den[0]["v_rms"], s=9,color='blue',label=r"2$\sigma$, P8" )

		ax1.scatter( tb8Den[-1]["area_exact"],  tb8Den[-1]["v_rms"], s=9,color='red',label=r"7$\sigma$, P8" )

		ax1.legend()



		ax1.set_xlim(0,1000)
		plt.savefig("testFarComponnents.png" , bbox_inches='tight')

	def testLargestFluxRatio(self):
		"""
		print the largest ratio of flux

		
		:return:
		"""
		algDendro="Dendrogram"
		tb8Den,tb16Den,label8Den,label16Den,sigmaListDen=self.getTBList(algorithm=algDendro)
		tb8Den=self.removeAllEdges(tb8Den)
		tb16Den=self.removeAllEdges(tb16Den)


		algDB="DBSCAN"
		tb8DB,tb16DB,label8DB,label16DB,sigmaListDB=self.getTBList(algorithm=algDB)
		tb8DB=self.removeAllEdges(tb8DB)
		tb16DB=self.removeAllEdges(tb16DB)

		algSCI="SCIMES"
		tb8SCI,tb16SCI,label8SCI,label16SCI,sigmaListSCI=self.getTBList(algorithm= algSCI )
		tb8SCI=self.removeAllEdges(tb8SCI)
		tb16SCI=self.removeAllEdges(tb16SCI)


		#ratio of clouds with 0.5 degree of the Galactic latitudes, which is used to estimate the conamination of Perseus arm
		#0.5 672
		#-0.5, 552

		#-1, 492
		#1, 732

		cutB=1.
		cutIndex1=492
		cutIndex2=732

		print algDendro
		for eachT in tb8Den:
			print "Ratio of max flux",np.max( eachT["flux"])/np.sum(   eachT["flux"] )

			if np.max( eachT["y_cen"] )>10:

				select= np.logical_and( eachT["y_cen"]>=cutIndex1, eachT["y_cen"]<= cutIndex2  )
				perSeusTB=eachT[select]
				print "Ratio of Perseus Clouds", len(perSeusTB)/1./len(eachT)


			else:
				select= np.logical_and( eachT["y_cen"]>=-cutB, eachT["y_cen"]<= cutB  )
				perSeusTB=eachT[select]
				print "Ratio of Perseus Clouds", len(perSeusTB)/1./len(eachT)
			#print eachT["y_cen"].unit



		for eachT in tb8DB:
			print np.max( eachT["sum"])/np.sum(   eachT["sum"] )


			if np.max( eachT["y_cen"] )>10:

				select= np.logical_and( eachT["y_cen"]>=cutIndex1, eachT["y_cen"]<= cutIndex2  )
				perSeusTB=eachT[select]
				print "Ratio of Perseus Clouds", len(perSeusTB)/1./len(eachT)


			else:
				select= np.logical_and( eachT["y_cen"]>=- cutB , eachT["y_cen"]<= cutB  )
				perSeusTB=eachT[select]
				print "Ratio of Perseus Clouds", len(perSeusTB)/1./len(eachT)
			#print eachT["y_cen"].unit
		print algSCI

		for eachT in tb8SCI:
			print np.max( eachT["flux"])/np.sum(   eachT["flux"] )

			if np.max( eachT["y_cen"] )>10:

				select= np.logical_and( eachT["y_cen"]>=cutIndex1 , eachT["y_cen"]<= cutIndex2  )
				perSeusTB=eachT[select]
				print "Ratio of Perseus Clouds", len(perSeusTB)/1./len(eachT)


			else:
				select= np.logical_and( eachT["y_cen"]>=- cutB , eachT["y_cen"]<=cutB  )
				perSeusTB=eachT[select]
				print "Ratio of Perseus Clouds", len(perSeusTB)/1./len(eachT)
			#print eachT["y_cen"].unit

	def getTablesByPath(self,CORawFITS,processPath):

		"""

		:param processPath:
		:return:
		"""

		labelFITS=glob.glob(processPath+"*.fits")

		for eachLabel in labelFITS:

			#fileName= os.path.split(eachLabel )[-1][:-5]

			self.getCatFromLabelArray(CORawFITS, eachLabel,  self.TBModel, saveMarker=eachLabel[0:-5]  )

	def getDBSCANTBList(self,conType,PixList,processPath):

		TBList= []
		TBCount= []

		for eachPix in PixList:

			searchStr =  processPath+"*P{}Con{}.fit".format(eachPix,conType)
			tbFile=glob.glob( searchStr  )[0]
			tb=  Table.read(tbFile)
			TBList.append( tb )

			TBCount.append( len(tb) )

		#print len(Table.read(tbFile)   )


		return TBList,TBCount

	def drawTestDBSCAN(self, processPath ):

		"""

		:return:
		"""
		TableFiles=glob.glob(processPath+"*.fits")

		connect2PixList= [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
		connect1PixList=   [3,4,5,6,7 ]
		connect3PixList= [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]

		TBCon1,countCon1 = self.getDBSCANTBList(  1, connect1PixList,  processPath)
		TBCon2,countCon2 = self.getDBSCANTBList(  2, connect2PixList,  processPath)
		TBCon3 ,countCon3= self.getDBSCANTBList(  3, connect3PixList,  processPath)


		#draw

		fig = plt.figure(figsize=(10, 6))
		ax = fig.add_subplot(1, 1, 1)
		# fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
		rc('text', usetex=True)
		rc('font', **{'family': 'sans-serif', 'size': 15, 'serif': ['Helvetica']})

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]

		scatterSize=14

		ax.plot( connect1PixList, countCon1,'b--' ,lw=1.5)
		ax.scatter( connect1PixList, countCon1, color='b',s= scatterSize )

		
		ax.plot( connect2PixList, countCon2,'g--',lw=1.5)
		ax.scatter( connect2PixList, countCon2, color='g',s= scatterSize )



		ax.plot( connect3PixList, countCon3,'r--',lw=1.5)
		ax.scatter( connect3PixList, countCon3, color='r',s= scatterSize )
		print connect3PixList
		print countCon3
		ax.set_xlim(2,27)
		#ax.set_ylim(-1,12)

		ax.set_ylabel("Fasely detected numbers")

		ax.set_xlabel("minPts")
		
		plt.savefig("dbscanTestParameter.png" , bbox_inches='tight')

	def countTBNlist(self,TBList):

		nList=[]

		for eachTB in TBList:
			nList.append( len(eachTB) )

		return nList

	def cleanUMMCTB(self,inputTB,rmsData,minValue=2.0, minPix=3, minDelta=1):
		"""
		clearn DBSCAN table according to minPix and minDelta
		:param TB:
		:param minPix:
		:param minDelta:
		:return:
		"""



		if type(inputTB) is list:

			newTBlist = []

			for eachTB in inputTB:
				TB = eachTB.copy()

				# need to convert peakvalue into sigma

				peakBIndex = TB["peakB"]  # the unit of peakB is
				peakLIndex = TB["peakL"]  # the unit of peakL is

				peakBIndex = map(int, peakBIndex)
				peakLIndex = map(int, peakLIndex)

				cooridnates = tuple([peakBIndex, peakLIndex])

				peakSigma = TB["peak"] / rmsData[cooridnates]

				minPeakSigma = minValue + minDelta

				select1 = TB[peakSigma >= minPeakSigma]  # based on minDelta
				select2 = select1[select1["pixN"] >= minPix]

				newTBlist.append( select2 )
			return newTBlist
		else:

			TB=inputTB.copy()

			#need to convert peakvalue into sigma

			peakBIndex=TB["peakB"] #the unit of peakB is
			peakLIndex=TB["peakL"] #the unit of peakL is

			peakBIndex= map(int, peakBIndex)
			peakLIndex= map(int, peakLIndex)

			cooridnates= tuple( [ peakBIndex ,peakLIndex  ] )

			peakSigma= TB["peak"] /rmsData[cooridnates]

			minPeakSigma=minValue+minDelta

			select1 = TB[ peakSigma>= minPeakSigma] #based on minDelta

			select2= select1[ select1["pixN"]>=minPix]
			return select2

	def relabelDBWithRMSFITS(self,FITSLabel,FITSTB, rmsFITS,  minPix=8,minDelta=3,minValue=2):
		"""

		:param FITSLabel:
		:param FITSTB:
		:param minPix:
		:param minDelta:
		:param minValue:
		:return:
		"""
		#first, get the sub TB
		rmsData,rmsHead = myFITS.readFITS(rmsFITS   )
		rawTB=Table.read(  FITSTB )

		goodCloudList= self.cleanUMMCTB( rawTB , rmsData,  minValue=minValue, minPix=minPix,minDelta=minDelta ) #self,inputTB,rmsData,minValue=2.0, minPix=3, minDelta=1):


		dataCluster,headCluster=myFITS.readFITS(FITSLabel)

		clusterIndex1D= np.where( dataCluster>0 )
		clusterValue1D=  dataCluster[clusterIndex1D ]

		Z0,Y0,X0 = clusterIndex1D
		


		newLabel=np.zeros_like( dataCluster )

		print "Cleaning DBSCAN UMMC table..."

		widgets = ['Recalculating cloud parameters: ', Percentage(), ' ', Bar(marker='0',left='[',right=']'),  ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options

		pbar = ProgressBar(widgets=widgets, maxval=len(goodCloudList))
		pbar.start()

		indexRun=0
		for eachDBRow in goodCloudList:
			indexRun=indexRun+1
			pbar.update(indexRun)

			cloudID=  eachDBRow["_idx"]

			cloudIndex = self.getIndices(Z0, Y0, X0, clusterValue1D, cloudID)
			newLabel[cloudIndex] =  cloudID

		pbar.finish()

		fits.writeto( "UMMCFormaClouds.fits",newLabel,header=headCluster,overwrite=True)



	def testMinDelta(self,processPath):
 		"""
 		test the effect of minDelta, i.e., the peak value
 		:param processPath:
 		:return:
 		"""

		CO12RMSFITS = "/home/qzyan/WORK/projects/NewUrsaMajorPaper/rmsGood_UMMC_CO12_CropEmpyt.fits"
		rmsData, rmsHead = myFITS.readFITS(CO12RMSFITS)

		connect2PixList = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
		connect1PixList = [3, 4, 5, 6, 7]
		connect3PixList = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27]

		TBCon1, countCon1 = self.getDBSCANTBList(1, connect1PixList, processPath)
		TBCon2, countCon2 = self.getDBSCANTBList(2, connect2PixList, processPath)
		TBCon3, countCon3 = self.getDBSCANTBList(3, connect3PixList, processPath)

		# draw
		fig = plt.figure(figsize=(10, 6))
		ax = fig.add_subplot(1, 1, 1)
		# fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
		rc('text', usetex=True)
		rc('font', **{'family': 'sans-serif', 'size': 15, 'serif': ['Helvetica']})

		mpl.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',  # helvetica font
			r'\usepackage{sansmath}',  # math-font matching  helvetica
			r'\sansmath'  # actually tell tex to use it!
			r'\usepackage{siunitx}',  # micro symbols
			r'\sisetup{detect-all}',  # force siunitx to use the fonts
		]

		scatterSize = 14

		testDeltaList=  [  2,3]

		minPixTest=15

		if 1:#test connectivity

			for eachDelta in  testDeltaList :

				tempTB=self.cleanUMMCTB( TBCon1,rmsData,minValue=2.0, minPix=minPixTest, minDelta=eachDelta)
				tempTBCount= self.countTBNlist( tempTB  )
				print tempTBCount
				ax.plot(connect1PixList, tempTBCount,'--', marker=11, lw=1.5)
				#ax.scatter(connect1PixList, tempTBCount, color='b', s=scatterSize)

			for eachDelta in testDeltaList:
				tempTB = self.cleanUMMCTB(TBCon2, rmsData, minValue=2.0, minPix=minPixTest, minDelta=eachDelta)
				tempTBCount = self.countTBNlist(tempTB)
				print tempTBCount
				ax.plot(connect2PixList, tempTBCount, '-', marker=3, lw=1.5)
			# ax.scatter(connect1PixList, tempTBCount, color='b', s=scatterSize)

			for eachDelta in testDeltaList:
				tempTB = self.cleanUMMCTB(TBCon3, rmsData, minValue=2.0, minPix=minPixTest, minDelta=eachDelta)
				tempTBCount = self.countTBNlist(tempTB)
				print tempTBCount
				ax.plot(connect3PixList, tempTBCount, '-', marker=5, lw=1.5)
			# ax.scatter(connect1PixList, tempTBCount, color='b', s=scatterSize)


		if 0:
			ax.plot(connect1PixList, countCon1, 'b--', lw=1.5)
			ax.scatter(connect1PixList, countCon1, color='b', s=scatterSize)

			ax.plot(connect2PixList, countCon2, 'g--', lw=1.5)
			ax.scatter(connect2PixList, countCon2, color='g', s=scatterSize)

			ax.plot(connect3PixList, countCon3, 'r--', lw=1.5)
			ax.scatter(connect3PixList, countCon3, color='r', s=scatterSize)


		ax.set_xlim(2, 27)
		# ax.set_ylim(-1,12)

		ax.set_ylabel("Fasely detected numbers")

		ax.set_xlabel("minPts")

		plt.savefig("dbscanTestMinDelta.png", bbox_inches='tight')


	def DBSCANAllInOne(self, rawCOFITS, saveTag, rmsFITS=None, minDelta = 3, minValue=2, minPix=8, MinPts=8, sigma=0.5 ,    minAreaPix=9  ,  connectivity=2,outPath="./" ,  ):
		"""

		:param rawCOFITS:
		:param saveTag:
		:param rmsFITS:
		:param minDelta:
		:param minValue:
		:param minPix:
		:param sigma:
		:return:

		"""

		#used to save
		tmpPath="./dbscanTmp/"

		dataCO, headCO= myFITS.readFITS( rawCOFITS )

		#step 1 compute DBSCAN
		print "Step 1: computing DBSCAN......"

		dbscanLabelFITS= self.computeDBSCAN(dataCO, headCO,  min_sigma= minValue , min_pix=MinPts , connectivity=connectivity, region= tmpPath+saveTag , rmsFITS= None, inputRMS= sigma  )

		print "Step 2: computing DBSCAN table......"
		rawDBSCANTBFile = self.getCatFromLabelArray(rawCOFITS,  dbscanLabelFITS , self.TBModel, minPix=minPix, rms=minValue , saveMarker=tmpPath+saveTag )

		print "Step 2: clean DBSCAN clusters..."

		self.clearnDBAssign( dbscanLabelFITS, rawDBSCANTBFile, pixN=minPix, minDelta=minDelta, minV=minValue , minAreaPix=minAreaPix,   prefix=outPath+saveTag+"_" )

	def cleanAndGetCatalog(self,rawCOFITS, labelFITS,saveTag, minDelta=3, rms=0.5,minPix=8,minValue=2):
		"""
		# labelFITS, is the fits created by HDBSCAN or DBSCAN or any other programs

		:param labelFITS:
		:param saveTag:
		:param minDelta:
		:param minPix, only used for saving
		:return:
		"""

		rawTBFile = self.getCatFromLabelArray(rawCOFITS, labelFITS, self.TBModel, minPix=minPix,  rms=rms, saveMarker= saveTag)
		self.clearnDBAssign( labelFITS, rawTBFile, pixN=minPix, minDelta=minDelta, minV=minValue ,prefix=  saveTag+"_" )



def ZZZZZZ(self):
		pass
		"""
		draw a simple map for over all moment of local molecular clouds
		:return:
		"""


#############################

if 0:

	doDBSCAN = myDBSCAN()

	#



if 0: # UMMC DBSCAN
	doDBSCAN = myDBSCAN()
	#UMMCCO12Data,UMMCCO12Head=myFITS.readFITS("/home/qzyan/WORK/projects/NewUrsaMajorPaper/OriginalFITS/UMMC_CO12_CropEmpty.fits")
	UMMCCO12Data,UMMCCO12Head=myFITS.readFITS("/home/qzyan/WORK/projects/NewUrsaMajorPaper/OriginalFITS/UMMC_CO12_Crop.fits")

	#CO12RMSFITS = "/home/qzyan/WORK/projects/NewUrsaMajorPaper/rmsGood_UMMC_CO12_CropEmpyt.fits"
	UMMC12COrms=0.16
	doDBSCAN.rms = UMMC12COrms

	savePathUMMC = "./UMMCTest/"
	#for minPixs in [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]:
		#doDBSCAN.computeDBSCAN(UMMCCO12Data, UMMCCO12Head, min_sigma=3, min_pix=minPixs, connectivity=2, region=savePathUMMC+"UMMC_DBSCAN_" ,inputRMS= UMMC12COrms )

	#for minPixs in [3,4,5,6,7 ]:

		#doDBSCAN.computeDBSCAN(UMMCCO12Data, UMMCCO12Head, min_sigma=3, min_pix=minPixs, connectivity=1, region=savePathUMMC+"UMMC_DBSCAN_" ,inputRMS= UMMC12COrms )

	#for minPixs in [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]:
		#doDBSCAN.computeDBSCAN(UMMCCO12Data, UMMCCO12Head, min_sigma=3, min_pix=minPixs, connectivity=3, region=savePathUMMC+"UMMC_DBSCAN_",inputRMS= UMMC12COrms )


	sys.exit()




if 0:
	G2650CO12FITS="/home/qzyan/WORK/myDownloads/testFellwalker/WMSIPDBSCAN/G2650Local30.fits"
	#doDBSCAN.getLVFITSByDBMASK( "G2650CO12dbscanS2.0P8Con2.fits", G2650CO12FITS, "/home/qzyan/WORK/myDownloads/testScimes/G2650PV.fits"  )
	#doDBSCAN.getLVFITSByDBMASK( "G2650minV2minP8_TrunkAsignMask0.fits", G2650CO12FITS, "/home/qzyan/WORK/myDownloads/testScimes/G2650PV.fits"  )
	doDBSCAN = myDBSCAN()
	doDBSCAN.getLVFITSByDBMASK( "/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/G2650LocaldbscanS2P4Con1_Clean.fits", G2650CO12FITS, "/home/qzyan/WORK/myDownloads/testScimes/G2650PV.fits",algorithm="DBSCANS2P4CON1"  )
	
	sys.exit()





# G2650Figures
if 0:
	doDBSCAN = myDBSCAN()
	# doDBSCAN.drawOverallMoment()
	doDBSCAN.drawLV()

	# doDBSCAN.numberDistribution()
	# doDBSCAN.drawCheckCloudsOneChannel()
	#doDBSCAN.drawSpectraExample()

	# doDBSCAN.drawVeDistribution(useGauss=True)
	# doDBSCAN.drawPeakDistribution()

	# doDBSCAN.areaDistribution()

	#doDBSCAN.alphaDistribution()
	# doDBSCAN.alphaDistribution( onlyLocal=True )

	# doDBSCAN.fluxAlphaDistribution()
	# doDBSCAN.fluxAlphaDistribution(onlyLocal=True)

	# doDBSCAN.fluxDistribution()
	# doDBSCAN.totaFluxDistribution()
	#doDBSCAN.physicalAreaDistribution()

	#doDBSCAN.physicalAreaDistributionLocal()

	sys.exit()

	pass



if 0:#get catalog for HDBSCAN
	saveTag = "./testHdbscan/ExtendTestCatalog"
	doDBSCAN = myDBSCAN()

	rawCOFITS = "HDBSCANTestLargeSubExtendCh1.fits"
	dbscanLabelFITS="/home/qzyan/WORK/myDownloads/MWISPcloud/testHdbscan/ExtendTesthdbscanProduceLocalS3P16.fits"
	rawDBSCANTBFile = doDBSCAN.getCatFromLabelArray(rawCOFITS, dbscanLabelFITS, doDBSCAN.TBModel, minPix=16, rms=0.5, saveMarker= saveTag, )

	sys.exit()





if 0:
	G2650CO12FITS="/home/qzyan/WORK/myDownloads/testFellwalker/WMSIPDBSCAN/G2650Local30.fits"
	doDBSCAN = myDBSCAN()
	G2650MaskCO = "/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/G2650LocaldbscanS2P4Con1_CO_Masked.fits" #expand edge of labels
	dendroLabel = "ClusterAsgn_3_1000Ve20.fits"
	extendedLabel = "G2650DisCloudVe20_extend.fits"
	doDBSCAN.myDilation(dendroLabel, G2650CO12FITS, startSigma=10, endSigma=2, saveName="G2650DisCloudVe20", maskCOFITS=G2650MaskCO)
	sys.exit()
	
if 0:#get catalog for HDBSCAN
	saveTag = "./testHdbscan/NorMalTestCatalog"
	doDBSCAN = myDBSCAN()

	rawCOFITS = "HDBSCANTestLargeSub.fits"
	dbscanLabelFITS="/home/qzyan/WORK/myDownloads/MWISPcloud/testHdbscan/NormalTesthdbscanProduceLocalS3P16.fits"
	rawDBSCANTBFile = doDBSCAN.getCatFromLabelArray(rawCOFITS, dbscanLabelFITS, doDBSCAN.TBModel, minPix=16, rms=0.5, saveMarker= saveTag, )

	sys.exit()





if 0:

	doDBSCAN=myDBSCAN()

	G2650CO12FITS="/home/qzyan/WORK/myDownloads/testFellwalker/WMSIPDBSCAN/G2650Local30.fits"
	DBMaskFITS= "/home/qzyan/WORK/myDownloads/testFellwalker/G2650DB_1_25.fits"
	TaurusCO12FITS="/home/qzyan/WORK/dataDisk/Taurus/t12_new.fits"
	PerCO12="/home/qzyan/WORK/dataDisk/MWISP/G2650/merge/G2650Per3060.fits"

	localCO13="/home/qzyan/WORK/dataDisk/MWISP/G2650/merge/G2650Local30CO13.fits"

	G210CO12="/home/qzyan/WORK/myDownloads/newMadda/data/G210CO12sm.fits"
	G210CO13="/home/qzyan/WORK/myDownloads/newMadda/data/G210CO13sm.fits"

	ursaMajor=""
	G2650MaskCO = "G2650CO12MaskedCO.fits"

	G2650MaskCODendro = "G2650CO12DendroMaskedCO.fits"








if 0:  #
	pass
	saveTag = "./sparseTest/HDBTESTCH_sparse"
	rawCO12fits = "sparseHDBSCAN.fits"  # "testHDBSCANCH_3.fits"
	labelFITS = "./sparseTest/HDBTESTCH_sparse.fits"

	doDBSCAN.cleanAndGetCatalog(rawCO12fits, labelFITS, saveTag)
	sys.exit()

if 0:  # #DBSCAN, find clouds at a spare region

	rawCO12fits = "sparseHDBSCAN.fits"  #"testHDBSCANCH_3.fits"

	#co12FITSInRMS = "/home/qzyan/WORK/projects/NewUrsaMajorPaper/UMMCCO12InRmsUnit.fits"
	saveTag = "DBSCANCH_sparse"
	doDBSCAN.DBSCANAllInOne(rawCO12fits, saveTag, rmsFITS=None, minDelta=3, minValue=2, minPix=8, MinPts=8, sigma=0.5, connectivity=2, outPath="./sparseTest/")

	sys.exit()







if 0:# #DBSCAN, find clouds

	
	rawCO12fits =  "/home/qzyan/WORK/projects/NewUrsaMajorPaper/OriginalFITS/UMMC_CO12_Crop.fits"

	co12FITSInRMS = "/home/qzyan/WORK/projects/NewUrsaMajorPaper/UMMCCO12InRmsUnit.fits"
	saveTag="DBSCANInRMSUNIT"


	doDBSCAN.DBSCANAllInOne(co12FITSInRMS, saveTag, rmsFITS=None, minDelta=3, minValue=2, minPix=8, MinPts=4, sigma=1 , connectivity=1, outPath="./UMMCFormal/"  )

	sys.exit()


if 0:#[30, 60] km/s DBSCAN, 12CO, sagitarius

	#rawCOFITS="/home/qzyan/WORK/myDownloads/MWISPcloud/G2650V3060/G2650V3060Sub.fits"
	#saveTag="G2650TestSagitarius"

	rawCOFITS="/home/qzyan/WORK/myDownloads/MWISPcloud/G2650V3060/G2650V3060CO13.fits"
	saveTag="G2650V3060CO13"

	doDBSCAN.DBSCANAllInOne(rawCOFITS, saveTag, rmsFITS=None, minDelta=3, minValue=2, minPix=8, MinPts=8, sigma=0.3 , connectivity=2, outPath="./G2650V3060/"  )

	sys.exit()


if 0:#test DBSCAN

	rawCOFITS="/home/qzyan/WORK/myDownloads/MWISPcloud/G2650V3060/G2650V3060Sub.fits"
	saveTag="G2650V3060MimicDendro"

	doDBSCAN.DBSCANAllInOne(rawCOFITS, saveTag, rmsFITS=None, minDelta=3.00, minValue=2, minPix=8, MinPts=3, sigma=0.5 , connectivity=1   )

	sys.exit()



if 0:#[30, 60] km/s DBSCAN

	#rawCOFITS="/home/qzyan/WORK/myDownloads/MWISPcloud/G2650V3060/G2650V3060Sub.fits"
	#saveTag="G2650TestSagitarius"

	rawCOFITS="/home/qzyan/WORK/myDownloads/MWISPcloud/G2650V3060/G2650Per3060.fits"
	saveTag="G2650V3060"

	doDBSCAN.DBSCANAllInOne(rawCOFITS, saveTag, rmsFITS=None, minDelta=3, minValue=2, minPix=8, MinPts=8, sigma=0.5 , connectivity=2, outPath="./G2650V3060/"  )

	sys.exit()

if 0: #testDBSCAN with empty cubes, which
	rawCOfits= "/home/qzyan/WORK/projects/NewUrsaMajorPaper/OriginalFITS/UMMC_CO12_Crop.fits"
	UMMCCO12DataTrue,UMMCCO12HeadTrue=myFITS.readFITS( rawCOfits )
	CO12RMSFITS = "/home/qzyan/WORK/projects/NewUrsaMajorPaper/rmsGood_UMMC_CO12_CropEmpyt.fits"
	savePathUMMC = "./UMMCFormal/"
	#doDBSCAN.computeDBSCAN(UMMCCO12DataTrue, UMMCCO12HeadTrue, min_sigma=2, min_pix=8, connectivity=2,  region=savePathUMMC + "UMMCSignal_", rmsFITS=CO12RMSFITS)

	#doDBSCAN.getCatFromLabelArray( rawCOfits ,  "/home/qzyan/WORK/myDownloads/MWISPcloud/UMMCFormal/UMMCSignal_dbscanS2P8Con2.fits",   doDBSCAN.TBModel, saveMarker= "UMMCCO12RawTB")


	labelFITS=  "/home/qzyan/WORK/myDownloads/MWISPcloud/UMMCFormal/UMMCSignal_dbscanS2P8Con2.fits"
	rawTBFile= "UMMCCO12RawTB.fit" 
	doDBSCAN.relabelDBWithRMSFITS( labelFITS ,   rawTBFile ,  CO12RMSFITS  )


	sys.exit()

if 0:# test minDelta

	processPath= "./UMMCTest/"
	doDBSCAN.testMinDelta( processPath )
	sys.exit()




if 0:
	CO12RMSFITS = "/home/qzyan/WORK/projects/NewUrsaMajorPaper/rmsGood_UMMC_CO12_CropEmpyt.fits"
	rmsData,rmsHead=myFITS.readFITS( CO12RMSFITS )
	testTB=Table.read("./UMMCTest/UMMCNoise_dbscanS2P8Con2.fit")

	doDBSCAN.cleanUMMCTB(testTB,rmsData,minDelta=3,minValue=2,minPix=3)

	sys.exit()




if 0: #draw TB test results
	processPath= "./UMMCTest/"

	doDBSCAN.drawTestDBSCAN( processPath )

	sys.exit()

if 0:#Calculate all tables in a path folder
	CORawFITS="/home/qzyan/WORK/projects/NewUrsaMajorPaper/OriginalFITS/UMMC_CO12_CropEmpty.fits"

	processPath= "./UMMCTest/"
	doDBSCAN.getTablesByPath(CORawFITS,processPath)

	sys.exit()



if 0: #testDBSCAN with empty cubes, which

	UMMCCO12Data,UMMCCO12Head=myFITS.readFITS("/home/qzyan/WORK/projects/NewUrsaMajorPaper/OriginalFITS/UMMC_CO12_CropEmpty.fits")

	CO12RMSFITS = "/home/qzyan/WORK/projects/NewUrsaMajorPaper/rmsGood_UMMC_CO12_CropEmpyt.fits"
	savePathUMMC = "./UMMCTest/"
	for minPixs in [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]:
		doDBSCAN.computeDBSCAN(UMMCCO12Data, UMMCCO12Head, min_sigma=2, min_pix=minPixs, connectivity=2, region=savePathUMMC+"UMMCNoise_",rmsFITS= CO12RMSFITS)

	for minPixs in [3,4,5,6,7 ]:
		doDBSCAN.computeDBSCAN(UMMCCO12Data, UMMCCO12Head, min_sigma=2, min_pix=minPixs, connectivity=1, region=savePathUMMC+"UMMCNoise_" ,rmsFITS= CO12RMSFITS)

	for minPixs in [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]:
		doDBSCAN.computeDBSCAN(UMMCCO12Data, UMMCCO12Head, min_sigma=2, min_pix=minPixs, connectivity=3, region=savePathUMMC+"UMMCNoise_" ,rmsFITS= CO12RMSFITS)



	sys.exit()



if 0: # DBSCAN test completeness

	UMMCCO12Data,UMMCCO12Head=myFITS.readFITS("/home/qzyan/WORK/projects/NewUrsaMajorPaper/OriginalFITS/UMMC_CO12_CropEmpty.fits")

	CO12RMSFITS = "/home/qzyan/WORK/projects/NewUrsaMajorPaper/rmsGood_UMMC_CO12_CropEmpyt.fits"


	savePathUMMC = "./UMMCTest/"
	for minPixs in [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]:
		doDBSCAN.computeDBSCAN(UMMCCO12Data, UMMCCO12Head, min_sigma=2, min_pix=minPixs, connectivity=2, region=savePathUMMC+"UMMC_completeness_" )

	#for minPixs in [3,4,5,6,7 ]:
		#doDBSCAN.computeDBSCAN(UMMCCO12Data, UMMCCO12Head, min_sigma=2, min_pix=minPixs, connectivity=1, region=savePathUMMC+"UMMC_" )

	#for minPixs in [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]:
		#doDBSCAN.computeDBSCAN(UMMCCO12Data, UMMCCO12Head, min_sigma=2, min_pix=minPixs, connectivity=3, region=savePathUMMC+"UMMC_" )


	sys.exit()




if 0:
	TBName="/home/qzyan/WORK/myDownloads/Q3/Q3LocalDendroCO12minV2minP8_dendroCatTrunk.fit"
	doDBSCAN.drawAreaDistribute(TBName, region="Q3")

	sys.exit()

if 0:
	doDBSCAN.testLargestFluxRatio()
	sys.exit()

if 0:
	doDBSCAN.testAngularAndLineWith()

	#doDBSCAN.contourDisClouds()
	sys.exit()






if 0:

	#label1FITS="minV2minP8_TrunkAsign.fits"
	#label2FITS="minV2minP16_TrunkAsign.fits"

	#label1FITS="minV2minP8_TrunkAsign.fits"
	#label2FITS="minV3minP8_TrunkAsign.fits"
	# compare,
	#label1FITS="minV2minP8_TrunkAsign.fits"
	#label2FITS="DBCLEAN2.0_8Label.fits"


	#doDBSCAN.getDiffLabel(label1FITS, label2FITS)


	#doDBSCAN.getFalseRate( "Only_minV2minP8_TrunkAsign.fits", noiseThreshold=0.5  )
	#doDBSCAN.getFalseRate( "Only_DBCLEAN2.0_8Label.fits" , noiseThreshold=0.5  )

	#doDBSCAN.getFalseRate("minV2.5minP8_TrunkAsign.fits", noiseThreshold = 0.5)
	#doDBSCAN.getFalseRate("minV2.5minP16_TrunkAsign.fits", noiseThreshold = 0.5)

	#doDBSCAN.getFalseRate( "DBCLEAN2.0_8Label.fits", noiseThreshold = 0.5  )


	#doDBSCAN.getCorrectRate("minV2minP8_TrunkAsign.fits", signalThreshold = 1.5)
	#doDBSCAN.getCorrectRate( "DBCLEAN2.0_8Label.fits", signalThreshold = 1.5 )

	#doDBSCAN.getCorrectRate("testMinDelta2DBCLEAN2_8Label.fits", signalThreshold = 1.5)

	#doDBSCAN.getCorrectRate("DBCLEAN2.0_8Label.fits", signalThreshold = 1.5)



	#doDBSCAN.getFalseRate("testMinDelta2DBCLEAN2_8Label.fits", noiseThreshold = 0.5)


	#doDBSCAN.getFalseRate("testMinDelta2DBCLEAN2_8Label.fits", noiseThreshold=0.5)

	doDBSCAN.getFalseRate("G2650minV2minP8_TrunkAsignMask0.fits", noiseThreshold=0.5)
	doDBSCAN.getFalseRate( "DBCLEAN2.0_8Label.fits", noiseThreshold=0.5  )

	sys.exit()




if 0:  # compare distributions

	doDBSCAN.printCatNumbers()




	#doDBSCAN.fluxDistribution()




	sys.exit()



if 0:

	doDBSCAN.getPublicCatalog()
	sys.exit()


if 0:

	#doDBSCAN.getAllSpectral()

	#doDBSCAN.getExamineSpectral()

	sys.exit()


if 0:


	doDBSCAN.physicalAreaDistribution()
	sys.exit()




if 0: #test get velocity dispersion

	inputTB=Table.read( "minV2minP8dendroMannualCat.fit" )
	doDBSCAN.getVdispersion( "G2650minV2minP8_TrunkAsignMask0.fits", inputTB, saveTBAs="minV2minP8dendroMannualCat_LineWidth.fit")

	inputTB=Table.read( "minV4minP8dendroMannualCat.fit" )
	doDBSCAN.getVdispersion( "G2650minV4minP8_TrunkAsignMask0.fits", inputTB, saveTBAs="minV4minP8dendroMannualCat_LineWidth.fit")

	inputTB=Table.read( "minV6minP8dendroMannualCat.fit" )
	doDBSCAN.getVdispersion( "G2650minV6minP8_TrunkAsignMask0.fits", inputTB, saveTBAs="minV6minP8dendroMannualCat_LineWidth.fit")



	inputTB=Table.read( "minV2P8scimesMannual.fit" )
	doDBSCAN.getVdispersion( "ClusterAsgn_2_8Ve20_mannual.fits", inputTB, saveTBAs="minV2minP8ScimesMannualCat_LineWidth.fit")

	inputTB=Table.read( "minV4P8scimesMannual.fit" )
	doDBSCAN.getVdispersion( "./scimesG2650/ClusterAsgn_4_8Ve20.fits", inputTB, saveTBAs="minV4minP8ScimesMannualCat_LineWidth.fit")

	inputTB=Table.read( "minV6P8scimesMannual.fit" )
	doDBSCAN.getVdispersion( "./scimesG2650/ClusterAsgn_6_8Ve20.fits", inputTB, saveTBAs="minV6minP8ScimesMannualCat_LineWidth.fit")




	sys.exit()




if 0: #Generate manucate for peak distribution, because SCIMES and dendro has not peak and intensity

	#dendrogram
	#doDBSCAN.getCatFromLabelArray(doDBSCAN.rawCOFITS, "G2650minV2minP8_TrunkAsignMask0.fits", doDBSCAN.TBModel, minPix=8, rms=2, saveMarker="minV2minP8dendroMannualCat")
	#doDBSCAN.getCatFromLabelArray(doDBSCAN.rawCOFITS, "G2650minV4minP8_TrunkAsignMask0.fits", doDBSCAN.TBModel, minPix=8, rms=4, saveMarker="minV4minP8dendroMannualCat")
	#doDBSCAN.getCatFromLabelArray(doDBSCAN.rawCOFITS, "G2650minV6minP8_TrunkAsignMask0.fits", doDBSCAN.TBModel, minPix=8, rms=6, saveMarker="minV6minP8dendroMannualCat")

	#scimes
	#doDBSCAN.getCatFromLabelArray(doDBSCAN.rawCOFITS, "./scimesG2650/ClusterAsgn_2_8Ve20.fits", doDBSCAN.TBModel, minPix=8, rms=2, saveMarker="scimesMannual")
	#doDBSCAN.getCatFromLabelArray(doDBSCAN.rawCOFITS, "./scimesG2650/ClusterAsgn_4_8Ve20.fits", doDBSCAN.TBModel, minPix=8, rms=4, saveMarker="scimesMannual")
	#doDBSCAN.getCatFromLabelArray(doDBSCAN.rawCOFITS, "./scimesG2650/ClusterAsgn_6_8Ve20.fits", doDBSCAN.TBModel, minPix=8, rms=6, saveMarker="scimesMannual")

	doDBSCAN.getCatFromLabelArray(doDBSCAN.rawCOFITS, "ClusterAsgn_2_8Ve20_mannual.fits", doDBSCAN.TBModel, minPix=8, rms=2, saveMarker="minV2P8scimesMannual")
	#doDBSCAN.getCatFromLabelArray(doDBSCAN.rawCOFITS, "./scimesG2650/ClusterAsgn_4_8Ve20.fits", doDBSCAN.TBModel, minPix=8, rms=4, saveMarker="minV4P8scimesMannual")
	#doDBSCAN.getCatFromLabelArray(doDBSCAN.rawCOFITS, "./scimesG2650/ClusterAsgn_6_8Ve20.fits", doDBSCAN.TBModel, minPix=8, rms=6, saveMarker="minV6P8scimesMannual")


	sys.exit()







if 0:
	algSCI = "SCIMES"
	tb8SCI, tb16SCI, label8SCI, label16SCI, sigmaListSCI = doDBSCAN.getTBList(algorithm=algSCI)
	tb8SCI = doDBSCAN.removeAllEdges(tb8SCI)
	tb16SCI = doDBSCAN.removeAllEdges(tb16SCI)

	doDBSCAN.getVdispersion("/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/ClusterAsgn_2_8Ve20.fits", tb8SCI[0], regionName="2sigma8pSCIMES", vCenPix=True)


	sys.exit()





if 0:

	doDBSCAN.produceSingleCubesForEachCloud( "G2650minV2minP8_TrunkAsignMask0.fits" ,  "minV2minP8_dendroCatTrunk.fit" )






if 0:

	#doDBSCAN.drawCheckClouds()
	sys.exit()








#veloicty distance, relation
# 13.46359868  4.24787753

if 0: #produce mask with dendrogram, because dendrogram gives the lowest false rate
		#maskCO="G2650CO12MaskedCO.fits"
		#doDBSCAN.computeDBSCAN(COData, COHead, min_sigma=2, min_pix=8, connectivity=2, region="G2650CO12Mask")

		doDBSCAN.produceMask( G2650CO12FITS, "G2650minV2minP8_TrunkAsignMask0.fits", "minV2minP8_dendroCatTrunk.fit",  region="G2650CO12Dendro")
		#doDBSCAN.getDiffCO("G2650CO12MaskedCO.fits","minV4minP16_TrunkAsign.fits",cutoff=4)


		sys.exit()





if 0: #produce new mask
		maskCO="G2650CO12MaskedCO.fits"
		#doDBSCAN.computeDBSCAN(COData, COHead, min_sigma=2, min_pix=8, connectivity=2, region="G2650CO12Mask")

		doDBSCAN.produceMask( G2650CO12FITS, "G2650CO12dbscanS2.0P8Con2.fits", "G2650CO12DBCatS2.0P8Con2.fit",  region="G2650CO12")
		#doDBSCAN.getDiffCO("G2650CO12MaskedCO.fits","minV4minP16_TrunkAsign.fits",cutoff=4)


		sys.exit()



if 0: #
	tbName = "G2650CO12DBCatS{:.1f}P{}Con2.fit".format(2, 8)
	fitsName = "G2650CO12dbscanS{:.1f}P{}Con2.fits".format(2, 8)
	doDBSCAN.clearnDBAssign(fitsName, tbName, pixN=8, minV=2, minDelta=3 , prefix="testMinDelta2" )
	sys.exit()






if 0:
	#doDBSCAN.getVdispersion()
	sys.exit()





if 0:
	doDBSCAN.produceDENDROCat()
	sys.exit()


##
if 0:#Scimes pipe line

	doDBSCAN.produceSCIMECat()
	#get catalog from scimes fits

	sys.exit()



if 0:
	#drawLrange,drawBrange=gaiaDis.box(38.4031840,0.3732480,22768.128 ,18363.802 ,0)
	drawLrange,drawBrange=gaiaDis.box(28.9977093,-0.0541066,22680.000 ,18662.400 ,0)

	doDBSCAN.drawCloudMap( lRange=drawLrange,bRange=drawBrange,  drawChannel= 91 )

	#doDBSCAN
	#doDBSCAN.drawAreaDistribute("ClusterCat_3_16Ve20.fit", region="scimes")
	#doDBSCAN.drawAreaDistribute("taurusDB3_8.fit" , region="scimes" )

	sys.exit()




if 0:#use DBSCAN to produe dendrotrunks
	COData, COHead = myFITS.readFITS(G2650CO12FITS)
	#DBLabelFITS="G2650CO12DendroByDBSCANdbscanS7P3Con1.fits"
	#DBTableFile="G2650CO12DendroByDBSCANCatS7P3Con1.fit"

	#doDBSCAN.clearnDBAssign(  DBLabelFITS, DBTableFile, pixN=16, minDelta=3, minV=7,   prefix="mimicDendro" )

	#DbscanSigmaList = np.arange(6)
	#for sigmas in [7]:
		#print "Calculating dengram with dBSCAN sigma:", sigmas
		#doDBSCAN.computeDBSCAN(COData, COHead, min_sigma=sigmas, min_pix= 3 , connectivity=1, region="G2650CO12DendroByDBSCAN" , mimicDendro= False )
	#sys.exit()
	if 1: #step2, calculate all catalog
		DbscanSigmaList = np.arange(2, 7.5, 0.5)
		for sigmas in [7]:
			labelFITS="G2650CO12DendroByDBSCANdbscanS7P3Con1.fits"
			savename="G2650CO12DendroByDBSCANCatS{}P{}Con1".format(sigmas,3)

			doDBSCAN.getCatFromLabelArray( G2650CO12FITS,  labelFITS  ,  doDBSCAN.TBModel, minPix=3, rms=sigmas, saveMarker= savename )

	sys.exit()

if 0: #small test
	pass

	TBListP8,TBListP16,TBLabelsP8,TBLabelsP16,dendroSigmaList=doDBSCAN.getTBList("SCIMES")

	print TBListP8
	sys.exit()







if 0: #high Galacticlatitudes, ursa major
	doDBSCAN.rms=0.16
	coFITS12UM="/home/qzyan/WORK/projects/NewUrsaMajorPaper/OriginalFITS/myCut12CO.fits"
	COData,COHead=myFITS.readFITS( coFITS12UM)
	#doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=2, min_pix=16, connectivity=2, region="UMCO12",savePath="./ursaMajor/")

	doDBSCAN.testFluxOfUM(  )

	#doDBSCAN.getCatFromLabelArray(coFITS12UM,"UMCO12dbscanS2P8Con2.fits",doDBSCAN.TBModel,saveMarker="UMCO12_2_8")
	#doDBSCAN.drawAreaDistribute("UMCO12_2_8.fit" , region="Taurus" )

	sys.exit()






if 0: #produce new mask
		maskCO="G2650CO12MaskedCO.fits"
		#doDBSCAN.computeDBSCAN(COData, COHead, min_sigma=2, min_pix=8, connectivity=2, region="G2650CO12Mask")

		doDBSCAN.produceMask( G2650CO12FITS, "G2650CO12dbscanS2.0P8Con2.fits", "G2650CO12DBCatS2.0P8Con2.fit",  region="G2650CO12")
		#doDBSCAN.getDiffCO("G2650CO12MaskedCO.fits","minV4minP16_TrunkAsign.fits",cutoff=4)


		sys.exit()







if 0: #DBSCAN PipeLine

	if 0: #produce all DBSCAN cases step1
		COData, COHead = myFITS.readFITS(G2650CO12FITS)

		DbscanSigmaList = np.arange(2, 7.5, 0.5)
		for sigmas in DbscanSigmaList:
			print "Calculating ",sigmas
			doDBSCAN.computeDBSCAN(  COData,COHead, min_sigma=sigmas, min_pix=8, connectivity=2, region="G2650CO12")

	if 0: #step2, calculate all catalog
		DbscanSigmaList = np.arange(2, 7.5, 0.5)
		for sigmas in DbscanSigmaList:
			labelFITS="G2650CO12dbscanS{}P8Con2.fits".format( sigmas )
			savename="G2650CO12DBCatS{}P{}Con2".format(sigmas,8)

			doDBSCAN.getCatFromLabelArray( G2650CO12FITS,  labelFITS  ,  doDBSCAN.TBModel, minPix=8, rms=sigmas, saveMarker= savename )




	if 0:

		doDBSCAN.cleanAllDBfits()
	sys.exit()


if 0:#distance pipeline

	#step1, extend 3,sigma,1000 pixels, to 2 sigma
	if 0: # use masked fits, a cutoff to
		dendroLabel= "ClusterAsgn_3_1000Ve20.fits"
		extendedLabel="G2650DisCloudVe20_extend.fits"

		doDBSCAN.myDilation( dendroLabel, G2650CO12FITS, startSigma=10,endSigma=2, saveName="G2650DisCloudVe20",maskCOFITS= G2650MaskCO )



	if 1: # get catalog
		dendro1000ExtendFITS = "G2650DisCloudVe20_extend.fits"
		savename="G2650CloudForDisCatDendro1000"
		doDBSCAN.getCatFromLabelArray(G2650CO12FITS, dendro1000ExtendFITS, doDBSCAN.TBModel, minPix=1000, rms=3, saveMarker=savename)


	if 0:
		pass
		#then produce int figures, usemyScimes.py








if 0: #Taurus
	COData,COHead=myFITS.readFITS( TaurusCO12FITS)
	doDBSCAN.rms=0.3
	#doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=3,min_pix=8,connectivity=2,region="Taurus")

	#doDBSCAN.getCatFromLabelArray(TaurusCO12FITS,"TaurusdbscanS3P8Con2.fits",doDBSCAN.TBModel,saveMarker="taurusDB3_8")
	#doDBSCAN.drawAreaDistribute("taurusDB3_8.fit" , region="Taurus" )

	


	#doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=1,min_pix=25,connectivity=3)

	#doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=4,min_pix=9,connectivity=2)
	#doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=5,min_pix=9,connectivity=2)
	#doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=6,min_pix=9,connectivity=2)

	sys.exit()


if 0:


	doDBSCAN.drawDBSCANNumber()
	doDBSCAN.drawDBSCANArea()



if 0: #dilation SCIMES

	#scimesFITS= "/home/qzyan/WORK/myDownloads/MWISPcloud/ClusterAsgn_ComplicateVe.fits"
	#rawFITS="/home/qzyan/WORK/myDownloads/testScimes/complicatedTest.fits"

	scimesFITS= "ClusterAsgn_3_16Ve20.fits"
	rawFITS= G2650CO12FITS  #"/home/qzyan/WORK/myDownloads/testScimes/complicatedTest.fits"

	doDBSCAN.myDilation( scimesFITS , rawFITS, saveName="G2650SCIMES_3_16Ve20", startSigma=15 )

	sys.exit()




if 0: #test Fast Dendrogram


	doDBSCAN.fastDendro("testDendro.fits",minV=5,minP=8)


	sys.exit()

if 0: #test Fast Dendrogram
	COData,COHead=myFITS.readFITS( G2650CO12FITS)


	#doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=5,min_pix=3,connectivity=1,region="G2650CO12DBDendro")

	#should peak sigma be larger?
	doDBSCAN.setMinVandPeak( "G2650CO12DBDendrodbscanS5P3Con1.fits" ,G2650CO12FITS,minP=8,peakSigma=8 )

	sys.exit()






if 0: # get catalog from extended fits


	doDBSCAN.getCatFromLabelArray(G2650CO12FITS,"G2650DisCloudVe20_extend.fits",doDBSCAN.TBModel,saveMarker="G2650CloudForDisCat")
	sys.exit()







if 0:
	#doDBSCAN.getCatFromLabelArray(G2650CO12FITS,"G2650CO12dbscanS2P16Con2.fits",doDBSCAN.TBModel, saveMarker="G2650CO12DBCatS2P16Con2" )
	for i in np.arange(2 ,8,0.5):
		savename="G2650CO12DBCatS{}P{}Con2".format(i,8)
		doDBSCAN.getCatFromLabelArray(G2650CO12FITS,"G2650CO12dbscanS{}P8Con2.fits".format(i),doDBSCAN.TBModel,saveMarker=savename)










if 0: #
	G214COFITS="G214CO12.fits"
	COData,COHead=myFITS.readFITS( G214COFITS)

	doDBSCAN.computeDBSCAN(COData, COHead,region="G214")
	doDBSCAN.slowDBSCAN(COData, COHead,region="G214")




	sys.exit()






if 0:# DBSCAN for G210
	region="G210CO13"

	processFITS=G210CO13

	doDBSCAN.rms=0.25


	if 1:#find clouds
		COData,COHead=myFITS.readFITS( processFITS)

		for sigmas in [3,4,5]:
			saveName=doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=sigmas, min_pix=8,connectivity=2, region=region)

			doDBSCAN.getCatFromLabelArray(processFITS, saveName , doDBSCAN.TBModel  ,  rms=1,minPix=1 , saveMarker=region+"DBSCAN{}_8".format(sigmas)   )

	sys.exit()




if 0:# DBSCAN for CO13
	region="Local13"
	doDBSCAN.rms=0.25


	if 0:#find clouds
		COData,COHead=myFITS.readFITS( localCO13)

		for sigmas in [3,4,5]:
			doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=sigmas, min_pix=9,connectivity=2, region=region)

		sys.exit()

	else:#calcatelog

		doDBSCAN.getCatFromLabelArray(localCO13,  "Local13dbscanS3P9Con2.fits" , doDBSCAN.TBModel  ,  rms=1,minPix=1 , saveMarker=region+"DBSCAN3_9"  )
		doDBSCAN.getCatFromLabelArray(localCO13,  "Local13dbscanS4P9Con2.fits" , doDBSCAN.TBModel  ,  rms=1,minPix=1 , saveMarker=region+"DBSCAN4_9"  )
		doDBSCAN.getCatFromLabelArray(localCO13,  "Local13dbscanS5P9Con2.fits" , doDBSCAN.TBModel  ,  rms=1,minPix=1 , saveMarker=region+"DBSCAN5_9"  )



if 0:
	#draw perseus
	#doDBSCAN.drawAreaDistribute("DBSCAN3_9.fit"  )

	doDBSCAN.drawAreaDistribute("ClusterCat_3_16Ve20.fit" , region="scimes" )


	#doDBSCAN.drawAreaDistribute("minV3minP16_dendroCatTrunk.fit" , region="Perseus" )

	#doDBSCAN.drawSumDistribute("DBSCAN3_9Sigma1_P1FastDBSCAN.fit"  )


	#doDBSCAN.drawSumDistribute("DBSCAN3_9.fit"  )



	#doDBSCAN.drawSumDistribute("minV3minP16_dendroCatTrunk.fit"  )
	#doDBSCAN.drawDBSCANArea()

	sys.exit()



if 0:# DBSCAN for perseus
	region="PerG2650"
	PerCO12="/home/qzyan/WORK/dataDisk/MWISP/G2650/merge/G2650Per3060.fits"

	if 0:#find clouds


		COData,COHead=myFITS.readFITS( PerCO12)
		doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=4,min_pix=9,connectivity=2, region=region)
		doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=5,min_pix=9,connectivity=2, region=region)

		sys.exit()
	else:#calcatelog

		doDBSCAN.getCatFromLabelArray(PerCO12,  "PerG2650dbscanS3P9Con2.fits" , doDBSCAN.TBModel  ,  rms=1,minPix=1 , saveMarker="DBSCAN3_9"  )
		doDBSCAN.getCatFromLabelArray(PerCO12,  "PerG2650dbscanS4P9Con2.fits" , doDBSCAN.TBModel  ,  rms=1,minPix=1 , saveMarker="DBSCAN4_9"  )
		doDBSCAN.getCatFromLabelArray(PerCO12,  "PerG2650dbscanS5P9Con2.fits" , doDBSCAN.TBModel  ,  rms=1,minPix=1 , saveMarker="DBSCAN5_9"  )

if 0:#Example
	doDBSCAN.rms=0.5
	COData,COHead=myFITS.readFITS( CO12FITS)
	doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=i,min_pix=9,connectivity=2)




if 0:
	ModelTB="minV3minP16_dendroCatTrunk.fit"

	#doDBSCAN.getCatFromLabelArray(CO12FITS, "dbscanS1P25Con3.fits",  ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN1_25"  )
	#doDBSCAN.getCatFromLabelArray(CO12FITS, "dbscanS2P16Con2.fits",  ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN2_16"  )
	#doDBSCAN.getCatFromLabelArray(CO12FITS,  "dbscanS3P9Con2.fits" ,   ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN3_9"  )
	#doDBSCAN.getCatFromLabelArray(CO12FITS,  "dbscanS4P9Con2.fits" ,   ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN4_9"  )
	#doDBSCAN.getCatFromLabelArray(CO12FITS,  "dbscanS5P9Con2.fits" ,   ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN5_9"  )
	#doDBSCAN.getCatFromLabelArray(CO12FITS,  "dbscanS6P9Con2.fits" ,   ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN6_9"  )

	doDBSCAN.getCatFromLabelArray(CO12FITS,  "dbscanS2.5P9Con2.fits" ,   ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN25_9"  )
	doDBSCAN.getCatFromLabelArray(CO12FITS,  "dbscanS3.5P9Con2.fits" ,   ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN35_9"  )
	doDBSCAN.getCatFromLabelArray(CO12FITS,  "dbscanS4.5P9Con2.fits" ,   ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN45_9"  )
	doDBSCAN.getCatFromLabelArray(CO12FITS,  "dbscanS5.5P9Con2.fits" ,   ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN55_9"  )
	doDBSCAN.getCatFromLabelArray(CO12FITS,  "dbscanS6.5P9Con2.fits" ,   ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN65_9"  )
	doDBSCAN.getCatFromLabelArray(CO12FITS,  "dbscanS7.5P9Con2.fits" ,   ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN75_9"  )
	doDBSCAN.getCatFromLabelArray(CO12FITS,  "dbscanS7P9Con2.fits" ,   ModelTB,  rms=1,minPix=1 , saveMarker="DBSCAN7_9"  )

	import sys
	sys.exit()

if 0:
	COData,COHead=myFITS.readFITS( CO12FITS)

	#doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=1,min_pix=25,connectivity=3)
	doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=2,min_pix=9,connectivity=2)
	doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=3,min_pix=9,connectivity=2)
	doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=4,min_pix=9,connectivity=2)
	doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=5,min_pix=9,connectivity=2)
	doDBSCAN.computeDBSCAN(COData,COHead, min_sigma=6,min_pix=9,connectivity=2)







if 0:

	ModelTB="minV3minP16_dendroCatTrunk.fit"

	doDBSCAN.getCatFromLabelArray(CO12FITS,DBMaskFITS,  ModelTB,  rms=1,minPix=25 )