from astropy.table import   Table
import numpy as np

from gaiaTB import  GAIATB
from progressbar import * 
import pymc3 as pm
from myPYTHON import *
from astropy.table import Column
from scipy.stats import norm,expon
from scipy  import special

def getDisAndErrorMCMCTrucatedGau6p(dataDis, dataAG,disError, AGError,mu1=0,distance0=-1,maxDistanceSearch=1000,distanceError=200,sampleN=1000,burn_in=500):
	
	"""
	do MCMC in this single function
	"""
	#first
	print "Using truncated Gaussion"
 
	#print "Calculating distances with MCMC...total smaple number is {}....".format(sampleN)
	ns=np.arange(int(np.min(dataDis)),int(np.max(dataDis)),1)

	minDis=int(np.min(dataDis))
	
	minDis=max(minDis,distance0)
	
	maxDis=int(np.max(dataDis))

		
	mu10=np.mean(dataAG[0:50])
	mu20=np.mean(dataAG[-50])

	imu2sigmaSigmaTest=np.std(dataAG[-200:],ddof=1)

	c0=[  mu20, 500 ]

	Nsample=len(dataAG)

	pns=ns #np.arange(100,1000,1)
 
	mu2, disCloud  =c0
	#mu1=0.1
	

	disCloud= 50 # minDis+(maxDis-minDis)*1.5/3
	disSigma =disCloud*0.1 # 50 #pc #disCloud*0.1 #minDis+(maxDis-minDis)*2./3

	imu1sigma,imu2sigma=[2,2]
  
	aaa=0.5*np.log(2*np.pi)
	
	sqrt2pi=np.sqrt(2*np.pi)
	sqrt2=np.sqrt(2)
	 
	logsqrt2pi=0.5*np.log(2*np.pi)
  

	agmin= 0 #np.min(dataAG)
	agmax= 3.609 #np.max(dataAG)
	nominator=np.sqrt(2) #*sigma

	AGErrorSquare=AGError*AGError
	sqrtPi2=np.sqrt(2*np.pi)*0.5

	#a=  special.erf((agmax-goodMu)*nominatorInv*igoodSigma)+ special.erf((goodMu-agmin)*nominatorInv*igoodSigma)  

	#= np.exp(-0.5*np.sum((dataAG-goodMu)**2*igoodSigma**2))/a**Nsample*igoodSigma**Nsample #true probability
	#return -np.sum( conUdis-0.5*(dataAG-goodMu)**2/errorVar - np.log(a)-0.5*np.log(errorVar)    )
	
	
	trueError=np.sqrt(disError**2+disSigma**2)
	
	weight= norm.cdf(disCloud,dataDis,trueError)
	goodMu=mu2+(mu1-mu2)*weight #true norm
	igoodSigma=imu2sigma+(imu1sigma-imu2sigma)*weight #true norm
 
	#a= special.erf((agmax-goodMu)/nominator*igoodSigma)+ special.erf((goodMu-agmin)/nominator*igoodSigma) 

	errorVar= 1./igoodSigma**2+AGErrorSquare

	convolveSTD=  1./np.sqrt(errorVar)	
	a=  special.erf((agmax-goodMu)/nominator*convolveSTD)+ special.erf((goodMu-agmin)/nominator*convolveSTD)  
	p0= np.sum(-0.5*(dataAG-goodMu)**2/errorVar - np.log(a)-0.5*np.log( errorVar)  )
	
	 
	
	modeGibbs=0

	disList=[]
	mu1List=[]
	mu2List=[]
	imu1sigmaList=[]
	imu2sigmaList=[]
	disSigmaList=[]

	
	widgets = ['MCMCSmapleDistance: ', Percentage(), ' ', Bar(marker='>',left='|',right='|'),
	           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
	
	pbar = ProgressBar(widgets=widgets, maxval=sampleN+burn_in+1)
	pbar.start() 
 
	disCloud0=0
	sumAGError=np.sum(-aaa-np.log(AGError))

	print "Blind search for a reasonalbe distance.....",


	

	for i in range(200000):
		#c1=[np.random.normal(0.5,0.45), np.random.normal(0.8,0.45),np.random.choice(pns),np.random.normal(20,10)  ]
		if modeGibbs==0:
			#disCloud=np.random.normal(disList[-1],disList[-1]*0.1 )
			disCloud=np.random.uniform(minDis,maxDis) # maxDistanceSearch) #np.random.normal(distance0,distanceError) # #choice(pns)
			#disSigma= abs( np.random.normal(disCloud*0.1,disCloud*0.05)  #np.random.exponential(20)
			#disSigma=   np.random.exponential(disCloud*0.1) #np.random.exponential(20)

		if modeGibbs==1:

			mu2=abs( np.random.normal(mu20,0.45) )
			
			 
			#mu2=max([0.3,mu2])
		if modeGibbs==2:
			#mu1=abs( np.random.normal(0,0.45) )
			disSigma=abs(np.random.normal(disCloud*0.1,disCloud*0.05) )
			#disSigma=np.random.exponential(50)

			#mu1= np.random.exponential(0.1)   #abs( np.random.normal(0,0.45) )
			#mu1=np.random.normal(mu1List[-1],mu1List[-1]*0.2) 
 
		if modeGibbs==3:
			imu1sigma=abs(np.random.normal(2,2))
			#imu1simga=np.random.normal(0.45,0.3)
			
			#c1=[np.random.normal(0.5,0.45), c0[1],c0[2] ]
		if modeGibbs==4:
			#imu2sigma=np.random.normal(0.45,0.3)
			imu2sigma= abs( np.random.normal(2,2)  )
			
		if modeGibbs==4:
			mu1=abs( np.random.normal(mu10,0.45) )

			
			
			
		trueError=np.sqrt(disError**2+disSigma**2)

		weight= norm.cdf(disCloud,dataDis,trueError)
		goodMu=mu2+(mu1-mu2)*weight #true norm
		igoodSigma=imu2sigma+(imu1sigma-imu2sigma)*weight #true norm
	
		#p1=np.sum(-logsqrt2pi - np.log(goodSigma) -0.5*(goodMu-dataAG)**2/goodSigma**2 )
		errorVar= 1./igoodSigma**2+AGErrorSquare

		convolveSTD=  1./np.sqrt(errorVar)
		
		a=  special.erf((agmax-goodMu)/nominator*convolveSTD)+ special.erf((goodMu-agmin)/nominator*convolveSTD)  

		#a=  special.erf((agmax-goodMu)/nominator*igoodSigma)+ special.erf((goodMu-agmin)/nominator*igoodSigma)  
		#p1= np.exp(-0.5*np.sum((dataAG-goodMu)**2*igoodSigma**2))/a**Nsample*igoodSigma**Nsample #true probability
		#p1= np.sum(-0.5*(dataAG-goodMu)**2*igoodSigma**2 - np.log(a)+ np.log(igoodSigma))
 

	#= np.exp(-0.5*np.sum((dataAG-goodMu)**2*igoodSigma**2))/a**Nsample*igoodSigma**Nsample #true probability
		p1= np.sum(-0.5*(dataAG-goodMu)**2/errorVar - np.log(a)-0.5*np.log( errorVar)  )
 
		randomR=np.random.uniform(0,1)
		
  		if p1>=p0 or p1-p0>np.log(randomR):
			p0=p1; 
			#disCloud0=disCloud
			if modeGibbs==0:
				#print disCloud,disSigma , mu2,imu1sigma,imu2sigma,"->" 
				
				disList.append(disCloud )
				mu1List.append(mu1 )
				disSigmaList.append(disSigma )
				mu2List.append(mu2 )
				imu1sigmaList.append(imu1sigma )
				imu2sigmaList.append(imu2sigma )

				#print disCloud,disSigma,mu2,imu1sigma,imu2sigma
 				if len(disList)>50:
					break
	 		modeGibbs=(modeGibbs+1)%5

				
	print "over, now sampling! "	
	
	cutIndex=len(disList)/3
	disList[-1]=np.mean(disList[-cutIndex:])
	disList=disList[-2:]
	mu1List=mu1List[-2:]
	mu2List[-1]=np.mean(mu2List[-cutIndex:])
	mu2List=mu2List[-2:] 
	mu1List=mu1List[-2:] 

	disSigmaList[-1]= np.mean(disSigmaList[-cutIndex:])
	disSigmaList=disSigmaList[-2:]

	imu1sigmaList[-1]=np.mean(imu1sigmaList[-cutIndex:])
	imu1sigmaList=imu1sigmaList[-2:]
		
	imu2sigmaList[-1]=np.mean(imu2sigmaList[-cutIndex:])
	imu2sigmaList=imu2sigmaList[-2:] 	
		

	#print  disList[-1],disSigmaList[-1],mu2List[-1],imu1sigmaList[-1],imu2sigmaList[-1],"is this correcet?"
	
	
	disErrorSQUARE=disError**2	
	for i in range(1000000):

 		#c1=[np.random.normal(0.5,0.45), np.random.normal(0.8,0.45),np.random.choice(pns),np.random.normal(20,10)  ]
		if modeGibbs==0:
			disCloud=np.random.normal(disList[-1],disList[-1]*0.1 )
			
			#disSigma=   np.random.exponential(disCloud*0.1) #np.random.exponential(20)

			
			#disCloud=np.random.uniform(minDis,maxDis) # maxDistanceSearch) #np.random.normal(distance0,distanceError) # #choice(pns)
		if modeGibbs==1:
			
			#mu1= np.random.normal(mu1List[-1],mu1List[-1]*0.1)
			disSigma= np.random.normal(disSigmaList[-1],disSigmaList[-1]*0.1)

			#mu1= np.random.exponential(0.1)   #abs( np.random.normal(0,0.45) )

			#mu1= np.random.normal(mu1List[-1],mu1List[-1]*0.1)

			#mu1= abs(np.random.normal(0,0.2)) 
		if modeGibbs==2:
			mu2=np.random.normal(mu2List[-1],mu2List[-1]*0.1)
			#mu2=  np.random.normal(mu20,0.45)  

			
  
		if modeGibbs==3:
			imu1sigma=abs(np.random.normal(imu1sigmaList[-1],imu1sigmaList[-1]*0.1))
			#imu1simga=np.random.normal(0.45,0.3)
			
			#c1=[np.random.normal(0.5,0.45), c0[1],c0[2] ]
		if modeGibbs==4:
			imu2sigma=abs(np.random.normal(imu2sigmaList[-1],imu2sigmaList[-1]*0.1))
			#imu2sigma=np.random.normal(0.45,0.3)

			#c1=[np.random.normal(0.5,0.45), c0[1],c0[2] ]


		if modeGibbs==5:
			mu1= np.random.normal(mu1List[-1],mu1List[-1]*0.1)

			#mu1=abs(np.random.normal(0.3,0.45))
			#imu2sigma=np.random.normal(0.45,0.3)

			#c1=[np.random.normal(0.5,0.45), c0[1],c0[2] ]
			

		trueError=np.sqrt(disErrorSQUARE+disSigma**2)
		weight= norm.cdf(disCloud,dataDis,trueError)
		goodMu=mu2+(mu1-mu2)*weight #true norm
		igoodSigma=imu2sigma+(imu1sigma-imu2sigma)*weight #true norm
	
		#p1=np.sum(-logsqrt2pi - np.log(goodSigma) -0.5*(goodMu-dataAG)**2/goodSigma**2 )
 
		#a=  special.erf((agmax-goodMu)/nominator*igoodSigma)+ special.erf((goodMu-agmin)/nominator*igoodSigma)   
 
		errorVar= 1./igoodSigma**2+AGErrorSquare

		convolveSTD=  1./np.sqrt(errorVar)
			
		a=  special.erf((agmax-goodMu)/nominator*convolveSTD)+ special.erf((goodMu-agmin)/nominator*convolveSTD)  



	 	p1= np.sum(-0.5*(dataAG-goodMu)**2/errorVar - np.log(a)-0.5*np.log(errorVar)  )
 
  

		randomR=np.random.uniform(0,1)
		
  
 		if p1>=p0 or p1-p0>np.log(randomR):
			p0=p1; 
	 		modeGibbs=(modeGibbs+1)%6
			#disCloud0=disCloud
			if modeGibbs==0:
				#print disCloud,mu1,mu2,imu1sigma,imu2sigma,"----->" 
				disList.append(disCloud )
				mu1List.append(mu1 )
				#print mu1
				
				
				mu2List.append(mu2 )
				imu1sigmaList.append(imu1sigma )
				imu2sigmaList.append(imu2sigma )
				disSigmaList.append(disSigma )

				
				pbar.update(len(disList)) #this adds a little symbol at each iteration
 				if len(disList)>sampleN+burn_in:
					break

	pbar.finish()


	disArray=np.array(disList[burn_in+1:]) 
 
	
	#for mu1 t	
	#mu1ArrayT=np.array(mu1List[burn_in:])
	#print "The modeled mu1 is ",np.mean(mu1ArrayT)

	disSigmaArray=np.array(disSigmaList[burn_in+1:])
	mu1Array=np.array(mu1List[burn_in+1:])
		
	mu2Array=np.array(mu2List[burn_in+1:])
		
	mu1SigmaArray=1./np.array( imu1sigmaList[burn_in+1:])
	mu2SigmaArray=1./np.array( imu2sigmaList[burn_in+1:])


	#print "Sampling number,",len(mu2SigmaArray)
 
	return  [disArray,disSigmaArray,mu1Array,mu1SigmaArray,mu2Array,mu2SigmaArray] #disEqual,errorEqual,round()
 

#
 


#this  module is used to calculate the distance with Gaia ag extinction

class GAIADIS:
	#name="GSF_Gaia2_all"
	
	
	coint="coint"
	
	gaiaDo=GAIATB()
	
	def calDisWith5p(self,oncloudStars):
		"""
		calculate the distance with oncloud stars, no figure is provided
		
		5 parameers are used 
		mu1 is assgned to be zero
		"""
		pass
	

	def getDisAndAGFromTB(self,oncloudStars):
 
		dataDis= oncloudStars["parallax"]*0
		disError=dataDis*0
		dataAG=oncloudStars["a_g_val"]
		
		#if 1 :#:#calculate parallax with mcmc
		for i in range(len(dataDis)):
			para=oncloudStars[i]["parallax"]
			paraError=oncloudStars[i]["parallax_err"] #+0.1 #   #
			
			dA=1./np.random.normal(para,paraError,10000)*1000
			
			dataDis[i]=np.mean(dA)
			disError[i]=np.std(dA,ddof=1)
		#no bining
		AGError= oncloudStars["agError"] 
		
		return dataDis, dataAG,disError, AGError

	#def calDisWith6p(self,oncloudStars):
		#"""
		#calculate the distance with oncloud stars, no figure is provided
		#
		#5 parameers are used 
		#mu1 is assgned to be zero
		#"""

		
		
		
		#return getDisAndErrorMCMCTrucatedGau6p(dataDis, dataAG,disError, AGError)


	def getDisAndHPD(self,disSample):
		
		

		disSample=np.array(disSample)
		
		meanDis=np.mean(disSample)
		
		#disHPD=pm.stats.hpd(disSample,alpha=0.1)
		disHPD=pm.stats.hpd(disSample,alpha=0.05)

		
		print disHPD,"????????95%HPD????????"
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
		 
		 
		if meanMu>10:
			return np.round(meanMu,1), np.round(lowerMu,1),np.round(upperMu,1) 

		return np.round(meanMu,3), np.round(lowerMu,3),np.round(upperMu,3) 


	def assignOnsourceGaia(self,GaiaTB,bgFITS):
 
		GaiaTB= GaiaTB.copy()
		#print "processing {} gaia sourxes".format(len(GaiaTB))
		#cloudMark="cloudMark"
		
		
		bgData,bgHead=myFITS.readFITS(bgFITS)
		
 
		
		if len(bgData.shape)==3:
			
			bgData=bgData[0]
			
			bgHead["NAXIS"]=2
			del bgHead["CRPIX3"] 
			del bgHead["CDELT3"] 
			
			del bgHead["CRVAL3"] 
			del bgHead["CTYPE3"] 
 
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

		maxY,maxX=bgData.shape


		for i in range(len(Xs)):
			
			if Ys[i] >=maxY-1 or  Xs[i]>=maxX-1:
				GaiaTB[i][self.coint]= -100 #self.getCointByLB(bgData, WCS(bgHead), eachStar["l"],eachStar["b"])

			else:
			
				GaiaTB[i][self.coint]= bgData[Ys[i]][Xs[i]]  #self.getCointByLB(bgData, WCS(bgHead), eachStar["l"],eachStar["b"])



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


	def ZZZ(self):
		pass

	#def testDis
