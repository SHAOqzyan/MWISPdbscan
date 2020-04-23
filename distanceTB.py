import numpy as np

from astropy.table import Table, Column
import matplotlib.pyplot as plt
from myMCMC import MCMC
from matplotlib import rc
from astropy import wcs
from astropy.io import fits

from matplotlib import gridspec

from myPYTHON import *

#from mpl_toolkits.axes_grid1.anchored_artists import AnchoredText
from matplotlib.offsetbox import AnchoredText

import os



from astropy.table import   Table
import numpy as np

class disTB: #do not use mysql
	
	
	
	fitsPathDefault="/home/qzyan/WORK/projects/gaiaDistance/data/"
	
	saveFITSNameDefault=fitsPathDefault+"disCollection.fit"
	

	
	sourceName="sourceName"
 
	l='l' #float DEFAULT 999.0,
	b='b' # float DEFAULT 999.0,
	
	boxCenterL='cloudBoxCenterL' # float DEFAULT 999.0,
	boxCenterB='cloudBoxCenterB' # float DEFAULT 999.0,
	boxSizeL='cloudBoxSizeL' # float DEFAULT 999.0,
	boxSizeB ='cloudBoxSizeB' #float DEFAULT 999.0,
	
	paraErrorCut="paraErrorCut"
 
	noiseLevel='noiseLevel' # float DEFAULT 0
	signalLevel='signalLevel' # float DEFAULT 0

	cloudVlsr='vlsr'
 
	Note="Note"
	cutDistanceLower="cutDistanceLower"
	cutDistanceUpper="cutDistanceUpper"
	
	fitsFile="fitsFile"
 
	 
	AG1="AG1"
	AG1Std="AG1Std" #new  #hpd

	AG1HPDLow="AG1HPDLow" #new #hpd
	AG1HPDUp="AG1HPDUp" #new #hpd
 
	AG2="AG2"
	AG2Std="AG2Std" #new  #hpd

	
	AG2HPDLow="AG2HPDLow" #new #hpd
	AG2HPDUp="AG2HPDUp" #new #hpd
 
 
 
	distance='distance' # float DEFAULT 0,
	disStd="disStd" #new std
	disHPDLow="disHPDLow" #new hpd
	disHPDUp="disHPDUp" #new





	sigma1="sigma1" #new
	sigma1Std="sigma1Std" #new

	Sigma1HPDLow="Sigma1HPDLow" #new
	Sigma1HPDUp="Sigma1HPDUp" #new


	sigma2="sigma2" #new
	sigma2Std="sigma2Std" #new

	Sigma2HPDLow="Sigma2HPDLow" #new
	Sigma2HPDUp="Sigma2HPDUp" #new

 


	onStarN="onStarN"
	offStarN="offStarN"

	COsum="COsum" # sum of CO fits in the fits FILE
	
	colnames=[ sourceName, l,b,boxCenterL,boxCenterB,boxSizeL,boxSizeB ,paraErrorCut,noiseLevel,signalLevel,cloudVlsr,cutDistanceLower,cutDistanceUpper,AG1,AG1Std,AG1HPDLow,\
			AG1HPDUp, AG2,AG2Std,AG2HPDLow,AG2HPDUp, distance,disStd,disHPDLow,disHPDUp,sigma1,sigma1Std,Sigma1HPDLow,Sigma1HPDUp,sigma2,\
			sigma2Std,Sigma2HPDLow,Sigma2HPDUp,onStarN,offStarN,COsum,fitsFile,Note]
 
	dataTypes=[ 'S20', float,float,float,float,float,float ,float,float,float,float,float,float,float,float,\
			float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,\
			float,float,int,int,float,'S200','S200']
	
	parsecToMeter= 3.0857e16 #m 
	
	
	
	def __init__(self, TBName=None):
		
		

		

		if TBName!=None:
			
			#print  TBName ,"????????"

			self.fitsPath= os.path.split(TBName)[0]+"/"  # fitsPathDefault="/home/qzyan/WORK/projects/gaiaDistance/data/"
		
			self.saveFITSName =  TBName    # saveFITSNameDefault=fitsPath+"disCollection.fit"
		else:
			
			self.fitsPath= self.fitsPathDefault #="/home/qzyan/WORK/projects/gaiaDistance/data/"
		
			self.saveFITSName= self.saveFITSNameDefault #=fitsPath+"disCollection.fit"
			
			
		if TBName!=None and  not os.path.isfile(TBName):
			#get emptyp te 
			
			empTB=self.getEmptyTB()
			empTB.write( TBName)

	def calMassByXfactor(self,coInt,distance,xFactor=2.0e20):
		
		"""
		The unit of coInt must be K km/s, 
		distance pc
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


	def calMassByRow(self,dataRow,xFactor=2.0e20):
		
		"""
		The unit of coInt must be K km/s, 
		distance pc
		"""
		
		distance=dataRow[self.distance]
		coInt=dataRow[self.COsum]


		#print "coInt",coInt,"???????????????"


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


	def addNewColumn(self,newColName,dtype=float,afterCol=None):
	
		"""
		add a new column with 0 values
		
		if afterCol is none, after the last column
		
		"""
		#first get tabe
		oldTB=self.getTB()
		
		#create a new column
		oldColNames=oldTB.colnames
		N=len( oldTB)
		
		newCol=Column( np.zeros(N),dtype=dtype,name= newColName)	
		
		

		insertIndex=len(oldColNames )
		
		try:
			insertIndex= oldColNames.index(afterCol )
			insertIndex=insertIndex+1
		except:
			pass
		
 
		
		oldTB.add_column(newCol,insertIndex) #
		
		#print oldTB.colnames
		
		
		self.writeTB(oldTB,overwrite=True)
		



	def getEmptyTB(self):
		
		"""
		build a TB with columns
		"""
		
		return Table(names=self.colnames,dtype=self.dataTypes)
		
		
	def getDefaultValues(self):

 

		Dvalues=[]
		
		for eachDtype in self.dataTypes:
			
			Dvalues.append(  0 )
 
		
		return Dvalues

		
	def getEmptyRow(self):
		
		"""
		provide an empyt row,
		"""
		emptyTB=self.getEmptyTB()
		rowData=self.getDefaultValues()
		
		emptyTB.add_row(rowData)
		emptyTB[-1][self.cutDistanceLower]=1.
		emptyTB[-1][self.paraErrorCut]=0.2

		
		return emptyTB[-1]

 
	def weighted_avg_and_std(self,values, weights):
	    """
	    Return the weighted average and standard deviation.
	
	    values, weights -- Numpy ndarrays with the same shape.
	    """
	    average = np.average(values, weights=weights)
	    # Fast and numerically precise:
	    variance = np.average((values-average)**2, weights=weights)
	    
	    
	    std=np.sqrt(variance)*np.sqrt( len(values)/ (len(values)-1 ) )
 
	    
	    return (average,std )

	
	##def updateRow(self,newRow):
		##"""
		#"""

	def addRowToTB(self,newRow):
		"""
		Add a row to the fits collection
		"""
		
		#readExistingFITS 
		
		if os.path.exists(self.saveFITSName):
		
		
			preTB=Table.read(self.saveFITSName)
			os.remove(self.saveFITSName)
		else:
			preTB=self.getEmptyTB()
			
		
		checkRow=preTB[preTB[self.sourceName]==newRow[self.sourceName]]  

		if len(checkRow)==0:
			preTB.add_row(newRow)
			print 'Source {} saved!'.format(newRow[self.sourceName])

		if len(checkRow)==1:
			
			for eachR in preTB:
				if eachR[self.sourceName] == newRow[self.sourceName] :
					for colS in self.colnames:
						eachR[colS]=newRow[colS]
 
					print 'Source {} updated!'.format(newRow[self.sourceName])
					break
		preTB.write(self.saveFITSName)
 
		print "Distances are saved as "+self.saveFITSName


	def getTB(self):
		return Table.read(self.saveFITSName) 

	def writeTB(self,newTB,overwrite=False):
		
		if newTB ==None:
			return 
		
		if  os.path.isfile( self.saveFITSName ) and   overwrite :
			os.remove( self.saveFITSName) 
			
		
		
		
		if overwrite:
			newTB.write( self.saveFITSName )

			return 
			
 
			
		if  os.path.isfile( self.saveFITSName ) and not overwrite :
			return 
			
		if  os.path.isfile( self.saveFITSName ) and   overwrite :
 			
			newTB.write( self.saveFITSName )
		
		
		
		return  

		
		
		
	def getRowByName(self,sName):
		"""
		Add a row to the fits collection
		"""
		
		#readExistingFITS 

		if str.strip(sName)=="":
			return None
		if not os.path.exists(self.saveFITSName):
			return None
			
			
		preTB=Table.read(self.saveFITSName)

		
		checkRow=preTB[preTB[self.sourceName]==sName]  

		if len(checkRow)==0:
			#print sName,"does not exist. Returing None..."
			return None
		if len(checkRow)==1:
			#print "Source", sName," found"
 
			return    checkRow[0]
			
		return None
	def updateRow(self,dataRow):
		
		"""
		This function is used to update the distance result 
		"""
		oldTB=self.getTB()

		if dataRow[self.sourceName] not in oldTB[self.sourceName]:
			try:
				os.remove(self.saveFITSName)
			except:
				pass
			oldTB.add_row(dataRow)
			oldTB.write(self.saveFITSName )
			return

		for i,eachC in enumerate(oldTB):
			
			if eachC[self.sourceName] ==  dataRow[self.sourceName]:
				#update
				#eachC=dataRow
				
				#updating fits 
				
				try:
					os.remove(self.saveFITSName)
				except:
					pass
				oldTB[i]=dataRow
				oldTB.write(self.saveFITSName )
					
				#print dataRow
				print eachC[self.sourceName], "updated!!!"
				break
				
		return 

	def getMask(self,sName):
		"""
		
		#produce mask file according to lb range, this should not be used to astrodendro sources, because astrodendro provides mask fits naturely
		
		"""

		sourceRow=self.getRowByName( sName )

		lRange,bRange=self.getBoxRange( sourceRow)
		#print lRange
		fitsFile=sourceRow[self.fitsFile]
		
		hdu=  myFITS.downTo2D(fitsFile)
		
		data=hdu.data
		head=hdu.header
		wcsCloud=wcs.WCS(head)
		
		leftXIndex,lowerYIndex = wcsCloud.wcs_world2pix(lRange[1],bRange[0],0   )
		rightXIndex,upYIndex = wcsCloud.wcs_world2pix(lRange[0],bRange[1],0   )

		leftXIndex,lowerYIndex , rightXIndex,upYIndex  = map(round, [leftXIndex,lowerYIndex , rightXIndex,upYIndex  ] )

		leftXIndex,lowerYIndex , rightXIndex,upYIndex  = map(int, [leftXIndex,lowerYIndex , rightXIndex,upYIndex  ] )

		#print leftXIndex,lowerYIndex , rightXIndex,upYIndex
		
		maskData=data*0+1
		
		maskData[ :, :leftXIndex ]=0
		
		maskData[:, rightXIndex:   ]=0
		maskData[  : lowerYIndex, :]=0
		maskData[   upYIndex: ,: ]=0

		return maskData,data,head
 		
	def getBoxRange(self,Row):
		
		"""
		return the cloud box range of the row
		"""

 

		minL=Row[self.boxCenterL  ] - Row[ self.boxSizeL]/2.
		maxL=Row[self.boxCenterL  ] + Row[ self.boxSizeL]/2.

		minB=Row[self.boxCenterB  ] - Row[ self.boxSizeB]/2.
		maxB=Row[self.boxCenterB  ] + Row[ self.boxSizeB]/2.
		
		return [minL,maxL], [minB,maxB]
		

	def ZZZ(self):
		pass
		
 