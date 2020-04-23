
from astropy.io import fits
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

import sys
from astropy.table import Table
from starlink import kappa
from starlink import convert

from starlink import cupid

import os

#this file is used to examine the molecular cores found with starlink



sys.path.append("/home/qzyan/WORK/projects/maddalena/")

#from madda import myG210

#doMadd=myG210()


class myCore():


		corePath="./maddalenaCore/"

		CO13FileFITS="/home/qzyan/WORK/projects/maddalena/data/mosaic_L.fits" #the 13CO FITS
		CO13FileSDF="/home/qzyan/WORK/projects/maddalena/data/mosaic_L.sdf" #the 13CO FITS

		#CO13FileSDF_madd="/home/qzyan/WORK/projects/maddalena/data/madd_crop.sdf" #the 13CO FITS
		CO13FileSDF_madd="/home/qzyan/WORK/projects/maddalena/G216_13.sdf" #the 13CO FITS
		#CO12FileSDF_madd="/home/qzyan/WORK/projects/maddalena/G216_CO12.sdf" #the 13CO FITS

		CO13FileSDF_W345="/home/qzyan/WORK/projects/maddalena/data/cropW345CO13.sdf" #the 13CO FITS



		CO12FileSDF="/home/qzyan/WORK/projects/maddalena/data/mosaic_U.sdf" #the 13CO FITS


		G216SDF="G216.sdf"

		RMS12=0.5 #kelvin
		
		RMS13=0.25 #kelvin


		fellwakerCoresSDF= corePath+"fellwakerCores.sdf"


		fellWalker="fellwalker"
		GaussClumps="GaussClumps"


		test2DFITS= corePath+"CO12RGB_R.fits"
		test2DSDF= corePath+"CO12RGB_R.sdf"


		#default values for Gaussclump
		
		gaussPara={}
		
		#gaussPara["MAXBAD"]=0.05
		gaussPara["MAXNF"]= 100
		gaussPara["MAXSKIP"]= 10
		gaussPara["MODELLIM"]= 0.5
		#gaussPara["NPAD"]= 10
		gaussPara["S0"]= 1.0
		gaussPara["SA"]= 1.0
		gaussPara["SC"]= 1.0
		gaussPara["SB"]= 0.1
		#gaussPara["THRESH"]= 2.0 #minimum peak
		gaussPara["WMIN"]= 0.05 #minimum peak
		gaussPara["WWIDTH"]= 2.0 #minimum peak



		
		

		def __init__(self):
			pass
		
		


		def randomPara(self,startN,endN):
			
			for i in np.arange(startN,endN+1):
				
				print "Radomly testing the {}th case".format(i)
				self.writeConfigRandom(i)

				saveName="./gaussClumpTest/"+"RANDOM_W345_{}".format(i)
				doCore.testCoreMethod(doCore.CO13FileSDF_W345,doCore.GaussClumps, self.RMS13, i, saveName=saveName, useConfig=True)

				saveName="./gaussClumpTest/"+"RANDOM_G216_{}".format(i)
				doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13,i, saveName=saveName, useConfig=True)
				
				
				os.system("cp GaussClumps_config_{}.conf ".format(i) + "./gaussClumpTest/RANDOM_{}.conf".format(i) )


 

		def testPara(self,paraName,factor=1.0):
		
		
			defaultValue=self.gaussPara[paraName]
 
			testValue=defaultValue*factor

			
 			self.writeConfig(paraName,testValue) #the minimum value to be considered #default 0.5
 			
 			
 			valueStr=str(testValue).replace(".","_")
 			
 			
 			saveName="./gaussClumpTest/"+"TEST_{}{}".format(paraName,valueStr)
 			
 			#print saveName
			self.testCoreMethod(self.CO13FileSDF_madd,self.GaussClumps, self.RMS13, 1, saveName=saveName, useConfig=True)


 			#self.writeConfig(paraName,0.5) #the minimum value to be considered #default 0.5
			#doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13, saveName="TEST_MODELLIM0_5", useConfig=True)
 


		def drawSensitivePara(self):
 
			
			factors=np.arange(0.5,1.6,0.1)

			clumpsNumberDict={}
			
			


 			for k,v in self.gaussPara.iteritems():
				
				clumpNList=[]
				
				
				for factor in factors:
					factor=round(factor,2)
					value=v*factor
					factorStr=str(value).replace(".","_")
					
					TableName="./gaussClumpTest/TEST_{}{}GaussClumpscoreCat_1.fit".format(k,factorStr  )
					
					N= len(Table.read(TableName))
					clumpNList.append(N)
				clumpsNumberDict[k]=clumpNList
 
			#draw all parameters
			fig, axs = plt.subplots(nrows=3, ncols=3, constrained_layout=True)
			
			indexStart=0
			
 			for k,v in clumpsNumberDict.iteritems():
				#
				
				ax=axs.flat[indexStart]
				
				ax.plot(v,'b.')
				
				indexStart=indexStart+1
				ax.set_ylabel(k)
			plt.savefig("sensitive.eps",bbox_inches='tight')
			plt.savefig("sensitive.png",bbox_inches='tight',dpi=300)

		def findDeepSensitivePara(self):
			"""
			"""
			i=0
			
			
			factors=np.arange(0.5,1.6,0.1)

 			for k,v in self.gaussPara.iteritems():
				
				
				if i%2==1:
					i=i+1
					continue
 				
				print i
				for factor in factors:
					factor=round(factor,2)
					
 
					self.testPara(k,factor=factor)
				i=i+1
					 
					#self.testPara(k,factor=factor)
				
			

		def findSensitivePara(self):
			#the default one 			
			#self.testPara("MAXBAD",factor=1.)


 			for k,v in self.gaussPara.iteritems():
				
				self.testPara(k,factor=1.5)
				self.testPara(k,factor=0.5)





		def convertToSDF(self,FITSFile,overWrite=False):
			"""
			
			A method of convering FITS files to SDF files
			"""
		
			SDFFile=FITSFile.replace('.fits','.sdf') #assume they are in the same path?
		
			if os.path.isfile(SDFFile) and not overWrite:
				print "File exist...doing nothing"
				return 

			if os.path.isfile(SDFFile) and   overWrite:
				print "Overwriting....."
 


			convert.fits2ndf(FITSFile,SDFFile)
			



		def convertToFITS(self,SDFFile,overWrite=True):
			"""
			
			A method of convering FITS files to SDF files
			"""
		
			FITSFile=SDFFile.replace('.sdf','.fits' ) #assume they are in the same path?
		
			if os.path.isfile(FITSFile) and not overWrite:
				print "File exist...doing nothing"
				return 

			if os.path.isfile(FITSFile) and overWrite:
				print "Overwriting....."
				os.system('rm '+FITSFile)
 


			convert.fits2ndf(SDFFile,FITSFile )
			







		def getCore13(self,method="fellwalker",outcat="fellwakerCores.fit"  ):
			"""
			"""
			
			
			cores=cupid.findclumps( self.CO13FileSDF, self.fellwakerCoresSDF, method= method ,rms=self.RMS13,outcat=outcat,wcspar=True)
 		


 




 		
 		def rerunCoreMethod(self, source,i):
 			
 			if "345" in source:
 				SDFFile=self.CO13FileSDF_W345
 			if "216" in source:
 				SDFFile=self.CO13FileSDF_madd
 			
 			rms=0.25
 			
			saveName="./gaussClumpTest/"+"RANDOM_{}_{}".format(source,i)

 			
			method="GaussClumps"
			configFile=method+"_config_{}.conf".format(i) 
 
			CoreSDF=saveName+method+"core_{}.sdf".format(1)
			CoreFITS=saveName+method+"core_{}.fits".format(1)
			
			CoreCat=saveName+method+"coreCat_{}.fit".format(1)
			
			print "Using config file: {}".format(configFile)
			cores=cupid.findclumps( SDFFile, CoreSDF, method= method ,rms=rms,outcat=CoreCat,wcspar=True,deconv=True,config="^"+configFile)
 

			coreTB=Table.read(CoreCat)
			
			print "{} cores found with {}_1".format(len(coreTB),method, i   )
			

			#cenvert to coreFITS
			try :
				os.remove(CoreFITS)
			except:
				pass
		
			
			convert.ndf2fits(CoreSDF,CoreFITS)
			
			#read the core file 
			

			#
			return coreTB
 			
		
 			
 		
 
		def testCoreMethod(self, SDFFile, method,rms,configIndex,  useConfig=True,saveName=""):
			"""
			"""
					
			print "Searching cores...with file: {},using method: {}, RMS: {}".format(SDFFile,method,rms)

			configFile=method+"_config_{}.conf".format(configIndex) 
			
			

			
			CoreSDF=saveName+method+"core_{}.sdf".format(1)
			CoreFITS=saveName+method+"core_{}.fits".format(1)
			
			CoreCat=saveName+method+"coreCat_{}.fit".format(1)
			
			print configFile
			print CoreSDF
			print CoreFITS
			print CoreCat
			
			
			if os.path.isfile(configFile) and useConfig:
				print "Using config file: {}".format(configFile)
				cores=cupid.findclumps( SDFFile, CoreSDF, method= method ,rms=rms,outcat=CoreCat,wcspar=True,deconv=True,config="^"+configFile)
			else:
				print "No config file provided."
				cores=cupid.findclumps( SDFFile, CoreSDF, method= method ,rms=rms,outcat=CoreCat,wcspar=True,deconv=True,)



			coreTB=Table.read(CoreCat)
			
			print "{} cores found with {} ".format(len(coreTB),method  )
			

			#cenvert to coreFITS
			try :
				os.remove(CoreFITS)
			except:
				pass
		
			
			convert.ndf2fits(CoreSDF,CoreFITS)
			
			#read the core file 
			

			#
			return coreTB
			#cores=cupid.findclumps( self.CO13FileFITS, self.fellwakerCoresSDF, method= method ,rms=0.3,outcat=outcat,wcspar=True)
 		
 		
 		def writeConfigRandom(self,i):
 			
 			"""
 			write the configue value to configue file
 			
 			
 			"""
 		
			#first echo 

			#if parameterName=="MODELLIM":
				#os.system('echo "GAUSSCLUMPS.{}={}" > GaussClumps_config_1.conf'.format(parameterName,value))
				#return 
				
			
			#os.system('echo "GAUSSCLUMPS.MODELLIM=1.5" > GaussClumps_config_1.conf') #this is the only default value
			os.system('echo "GAUSSCLUMPS.{}={}" > GaussClumps_config_{}.conf'.format("NPAD",10,i)) #the true function is to remove previsou written parameters

			np.random.seed()

 			for k,v in self.gaussPara.iteritems():
 
				randomValue=np.random.uniform(0.5*v, 1.5*v )
				
				if k=="MAXNF":
					
					randomValue=np.arange(int(round(0.5*v)),  int(round(1.5*v))  )
					randomValue= np.random.choice(randomValue)
				

				os.system('echo "GAUSSCLUMPS.{}={}" >> GaussClumps_config_{}.conf'.format(k,randomValue,i))

 		
 		
 		def writeConfig(self,parameterName,value):
 			
 			"""
 			write the configue value to configue file
 			
 			
 			"""
 		
			#first echo 

			#if parameterName=="MODELLIM":
				#os.system('echo "GAUSSCLUMPS.{}={}" > GaussClumps_config_1.conf'.format(parameterName,value))
				#return 
				
			
			#os.system('echo "GAUSSCLUMPS.MODELLIM=1.5" > GaussClumps_config_1.conf') #this is the only default value
 
			os.system('echo "GAUSSCLUMPS.{}={}" > GaussClumps_config_1.conf'.format(parameterName,value))

			

		def getRMSs(self,fitsName):
			
			"""
			find the RMS of FITS file, the method is fitting a trucated Gaussian using MCMC, using maximum likelihood
			
			only negative values are concerned
			
			"""

		


 		
 		def testGaussianCores(self):
 			
 			"""
 			findcores with GaussianClumps, and change one parameters at a time to see the result of cores 
 			
 			#examine the peak residual of each core, which should given us some clues about hte Gaussian fitting
 			
 			"""
 			# examine the peak residual of 
 			
 			#test1
 			
 			if 0: # use 0.2 
				testValues=np.arange(0.5,4.2,0.1)
				
				
				print testValues
				for value in testValues:
					
					value=float(  "{:.1f}".format(value) )
					
					valueStr= "TEST_MODELLIM{:.1f}".format(value).replace(".","_")
	 			
		 			self.writeConfig("MODELLIM",value) #the minimum value to be considered
					doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13,saveName=valueStr,useConfig=True)
					

 			#self.writeConfig("MODELLIM",1.8) #the minimum value to be considered
			#doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13,saveName="TEST_MODELLIM1.8",useConfig=True)
			
 			
	
 			self.writeConfig("MODELLIM",0.5) #the minimum value to be considered #default 0.5
			doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13, saveName="TEST_MODELLIM0_5", useConfig=True)
 
 			
 			if 0:
 				
	 			self.writeConfig("MODELLIM",2.0) #the minimum value to be considered
				doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13,saveName="TEST_MODELLIM2_0",useConfig=True)
				
	 			self.writeConfig("MODELLIM",1.5) #the minimum value to be considered #default 0.5
				doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13, saveName="TEST_MODELLIM1_5", useConfig=True)
	 
	 			self.writeConfig("MODELLIM",1.0) #the minimum value to be considered
				doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13, saveName="TEST_MODELLIM1_0",useConfig=True)
	 		

	 		
	 		
	 			self.writeConfig("THRESH",2.5) #the minimum value to be considered #default 2.0
				doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13,saveName="TEST_THRESH2_5",useConfig=True)
	 		
	 			self.writeConfig("THRESH",1.5) #the minimum value to be considered
				doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13,saveName="TEST_THRESH1_5",useConfig=True)
	 		
	 			self.writeConfig("THRESH",4) #the minimum value to be considered
				doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13,saveName="TEST_THRESH4_0",useConfig=True)
	 		
 			#self.writeConfig("THRESH",8) #the minimum value to be considered
			#doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13,saveName="TEST_THRESH8_0",useConfig=True)
	
	 			self.writeConfig("MODELLIM",0.5) #the minimum value to be considered #default 0.5
				doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13, saveName="TEST_MODELLIM0_5", useConfig=True)
	 
	 			self.writeConfig("MODELLIM",2.5) #the minimum value to be considered #default 0.5
				doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13, saveName="TEST_MODELLIM2_5", useConfig=True)
	 
	 			self.writeConfig("MODELLIM",3.0) #the minimum value to be considered #default 0.5
				doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13, saveName="TEST_MODELLIM3_0", useConfig=True)
	 
	 			self.writeConfig("MODELLIM",3.5) #the minimum value to be considered #default 0.5
				doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, self.RMS13, saveName="TEST_MODELLIM3_5", useConfig=True)
 			#self.writeConfig("MODELLIM",3.5) #the minimum value to be considered #default 0.5
			#doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.fellWalker, self.RMS13, saveName="TEST_FellwalerG216", useConfig=True)
 
			#doCore.testCoreMethod(doCore.CO12FileSDF_madd,doCore.fellWalker, self.RMS12, saveName="FellwalerG216CO12", useConfig=True)


		def findCoreForAllG200220(self,saveName='allCoresG200220'):
			"""
			
			this function is used to find cores with the whold 13CO
			using default values
			"""

			SDFFile='./data/G210200merge13.sdf'
			method=self.GaussClumps

			rms=self.RMS13
 
			print "Searching cores...with file: {},using method: {}, RMS: {}".format(SDFFile,method,rms)

			#configFile=method+"_config_{}.conf".format(configIndex) 
			
 
			
			CoreSDF=saveName+method+"core.sdf" 
			CoreFITS=saveName+method+"core.fits" 
			
			CoreCat=saveName+method+"coreCat.fit" 
 
			print "No config file is using...."
			cores=cupid.findclumps( SDFFile, CoreSDF, method= method ,rms=rms,outcat=CoreCat,wcspar=True,deconv=True )
			
 
			coreTB=Table.read(CoreCat)
			
			print "{} cores found with {} ".format(len(coreTB),method  )
			

			#cenvert to coreFITS
			try :
				os.remove(CoreFITS)
			except:
				pass
 
			convert.ndf2fits(CoreSDF,CoreFITS)
			
			#read the core file 
			

			#
			return coreTB
			#cores=cupid.findclumps( self.CO13FileFITS, self.fellwakerCoresSDF, method= method ,rms=0.3,outcat=outcat,wcspar=True)
 		



		def findCoreForAll(self,saveName='allCores'):
			"""
			
			this function is used to find cores with the whold 13CO
			using default values
			"""

			SDFFile='./data/mosaic_L.fits'
			method=self.GaussClumps

			rms=self.RMS13
 
			print "Searching cores...with file: {},using method: {}, RMS: {}".format(SDFFile,method,rms)

			#configFile=method+"_config_{}.conf".format(configIndex) 
			
 
			
			CoreSDF=saveName+method+"core.sdf" 
			CoreFITS=saveName+method+"core.fits" 
			
			CoreCat=saveName+method+"coreCat.fit" 
 
			print "No config file is using...."
			cores=cupid.findclumps( SDFFile, CoreSDF, method= method ,rms=rms,outcat=CoreCat,wcspar=True,deconv=True )
			
 
			coreTB=Table.read(CoreCat)
			
			print "{} cores found with {} ".format(len(coreTB),method  )
			

			#cenvert to coreFITS
			try :
				os.remove(CoreFITS)
			except:
				pass
 
			convert.ndf2fits(CoreSDF,CoreFITS)
			
			#read the core file 
			

			#
			return coreTB
			#cores=cupid.findclumps( self.CO13FileFITS, self.fellwakerCoresSDF, method= method ,rms=0.3,outcat=outcat,wcspar=True)


		def getTBFITSNameListByPara(self, parameter):
			factors=np.arange(0.5,1.6,0.1)

			#clumpsNumberDict={}
			
			TBNameList=[]
			FitsNameList=[]
 			for k,v in self.gaussPara.iteritems():
				
				clumpNList=[]
				
				
				if k!=parameter:
					continue
				
				for factor in factors:
					
					factor=round(factor,2)
					value=v*factor
					factorStr=str(value).replace(".","_")
					
					TableName="./gaussClumpTest/TEST_{}{}GaussClumpscoreCat_1.fit".format(k,factorStr  )
					FITSName="./gaussClumpTest/TEST_{}{}GaussClumpscore_1.fits".format(k,factorStr  )
					
					
					TBNameList.append( TableName )
					FitsNameList.append( FITSName )
 
			return TBNameList,	FitsNameList		
			


		def getTestTBByPara(self,parameter):
			"""
			
			get TB list by Paraname
			
			"""
			
			TBNameList,	FitsNameList		=self.getTBFITSNameListByPara(parameter)
 
			
			TBList=[]


 			for TableName in TBNameList:
 
					TBList.append(  Table.read(TableName)  )
 
 
			return TBList


		def calG214Age(self):
			
			"""
			"""
			
			blueFITS="cropblue.fits"
			
			#blueSDF=blueFITS.replace(".fits",".sdf")
			
			#doconvert
			
			#blueCoreSDF=blueFITS.replace(".fits","core.sdf")

			#blueCore=blueFITS.replace(".fits","Cat.fit")
			#find the center with GaussClump
			
			#cores=cupid.findclumps( blueSDF,blueCoreSDF, method= self.GaussClumps,rms=0.6,outcat=blueCore,wcspar=True)
			
			peakBlueL=214.4833333
			peakBlueB=-1.8166667
			
			S287ClumpL=0
			S287ClumpB=0
			
			




		def findTaurusCores(self):

			"""
			"""
			TaurusCO13SDF="./data/t13_new.sdf"
			
			TaurusCO13CoreSDF="./data/t13_newCore.sdf"

			TaurusCO13CoreCat="./data/t13_newCoreCat.fit"
			rms=0.13
			
			
			cores=cupid.findclumps( TaurusCO13SDF, TaurusCO13CoreSDF, method= self.GaussClumps,rms=rms,outcat=TaurusCO13CoreCat,wcspar=True,deconv=True  )
 

		def zzz():
			pass




doCore=myCore()

if 0:
	doCore.findTaurusCores()


if 0:
	doCore=myCore()

	doCore.calG214Age()

if 0:
	TBList=doCore.getTestTBByPara( "S0" )
	print TBList[0]
	print TBList[1]

if 0:
	doCore.drawSensitivePara()


if 0:
	doCore=myCore()

	doCore.findDeepSensitivePara()


if 0:
	
	doCore=myCore()
	
	
	#number 210, is the default value 
	
	doCore.rerunCoreMethod("W345",210)
	doCore.rerunCoreMethod("G216",210)


if 0:
	"""
	"""
	
	doCore=myCore()
	
	
	#doCore.randomPara(0,100) #G216 and W345 should be calculated simultaneously
	#doCore.randomPara(152,152) #G216 and W345 should be calculated simultaneously
	doCore.randomPara(153,153) #G216 and W345 should be calculated simultaneously

	


			
if 0:
	
	"""
	"""
	doCore=myCore()
	
	doCore.findCoreForAllG200220()



	
	#doCore.randomPara(11,40)
	#doCore.randomPara(6,10)


	#doCore.randomPara(101,150) #G216 and W345 should be calculated simultaneously

	#doCore.randomPara(151,200) #G216 and W345 should be calculated simultaneously

	#doCore.randomPara(71,100)

	
	#doCore.findSensitivePara()






if 0:
	doCore=myCore()
	doCore.findCoreForAll()
	#


#doCore=myCore()





if 0:
	"""
	"""
	doCore=myCore()

	doCore.testGaussianCores()
	

if 0:
	"""
	"""
	doCore=myCore()
	doCore.testCoreMethod(doCore.CO13FileSDF,doCore.GaussClumps, 0.25 ,useConfig=True,saveName="G210ALL")

#doCore.testCoreMethod(doCore.G216SDF,doCore.GaussClumps, 0.3,useConfig=True,saveName="G216")


#use 12CO fellwaker to split the fits into clumps
#doCore.testCoreMethod(doCore.CO12FileSDF,doCore.fellWalker, 0.5,useConfig=True)




#do not use 12CO
#doCore.testCoreMethod(doCore.test2DSDF,doCore.GaussClumps,  1.16,useConfig=False)


#doCore.testCoreMethod(doCore.CO13FileSDF,doCore.GaussClumps, 0.3,useConfig=True)
#doCore.testCoreMethod(doCore.CO13FileSDF_madd,doCore.GaussClumps, 0.3,useConfig=False)
#doCore.testCoreMethod(doCore.CO13FileSDF ,doCore.GaussClumps, 0.3,useConfig=True)


#do not use fellwaker
#doCore.testCoreMethod(doCore.test2DSDF,doCore.fellWalker, 1.16)
#doCore.testCoreMethod(doCore.CO13FileSDF,doCore.fellWalker, 0.3)

#doMadd.drawCMF("GaussClumpscoreCat_1.fit")
#doMadd.drawCMF("fellwalkercoreCat_1.fit")


#doMadd.drawCMF("fellwalker1.fit")

if 0:

	convert.fits2ndf("CO12RGB_R.fits","CO12RGB_R.sdf")
	
	statsinfo = kappa.stats('CO12RGB_R.sdf')
	
	
	print statsinfo