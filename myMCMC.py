# collection of MCMC functions, such as linear regression and gaussian, non-gaussian fitting


import numpy as np

from scipy  import special
from scipy.stats import norm
logRoot2pi= np.log( 2*np.pi )
from progressbar import *





def logProbGauss(parameters,meanList,stdList):
	
	
	#n=len(meanList)
	
	mu,sigma=parameters
	
	
	#sigmaCal=sigma**2+stdList**2
	
	sigmaCal=sigma**2   +stdList**2


	
	#return np.sum( -0.5*np.log(2*np.pi*sigma**2)-0.5*(meanList-mu)**2/sigma**2-0.5*np.log(2*np.pi*stdList**2)-0.5*(meanList-mu)**2/stdList**2  )

	
	
	return np.sum( -0.5*np.log(2*np.pi*sigmaCal)-0.5*(meanList-mu)**2/sigmaCal )
	
	
	



def logProbGaussTruncate(parameters, dataArray,vMax):
	
	"""

	By default, vMin= minus infinite
	
	"""
	
	mu,isigma=parameters # isigma, is the inverse of the sigma
	
	part1=np.log(isigma)-0.5*logRoot2pi-0.5*(dataArray-mu)**2*isigma**2
	
	
	
	denominator=np.log( 0.5*(  special.erf( (vMax-mu)*isigma/np.sqrt(2) )  +1 ) )
	
	
	
	return np.sum(part1) - denominator*len(dataArray) 
	
def getNextValues( runSamples, RMSs):
 
	currentValues=[ ]

	for j in range(len(runSamples)):
		currentValues.append(runSamples[j][-1])

	for i in range(len(currentValues) ):
		currentValues[i]=currentValues[i] +np.random.normal(0,RMSs[i])
	
	return currentValues


class MCMC():

	
	def __init__(self):
		pass
	
	def getMeanAndStd(self,meanList,stdList,burn_in=100,thinning=15,sampleN=500):
		"""
		 
		"""
		pass
		meanList=np.array(meanList)
		
		stdList=np.array(stdList)
		
		


		mu0=np.mean(meanList)
		
		sigma0=np.std(meanList,ddof=1)
		
		
		
		
 		
		RMSs=[0.5,  0.1]

		if sigma0>0:

			mu=   np.random.normal(mu0,sigma0 )
			sigma= np.random.exponential(  sigma0 )

		else:
			
			mu=   np.random.normal(mu0,0.1 )
			sigma= np.random.exponential(  0.1 )


		#isigma= np.random.exponential( 1./sigma0 )
		
		parameters=[mu,sigma]

		p0=logProbGauss(parameters,  meanList, stdList)
			
		runSamples= [[mu],[sigma] ]
		
		

		muList=[mu]
		sigmaList=[ sigma]
		
		
		
		widgets = ['MCMCSmapleDistance: ', Percentage(), ' ', Bar(marker='>',left='|',right='|'),
		           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
		
		pbar = ProgressBar(widgets=widgets, maxval=sampleN+burn_in+1)
		pbar.start() 
		
		
		
		for i in range(30000000):
			
			
			
			nextValue= getNextValues(runSamples,RMSs)
			mu,sigma=nextValue
			
 
			
			if sigma< 0 :
				continue
			
			p1=logProbGauss(nextValue,   meanList, stdList)
			
			if np.isinf( p1):
				continue
		
		
			randomR=np.random.uniform(0,1)
			
			
			if p1>=p0 or p1-p0>np.log(randomR):
			#print disCloud,  imu1sigma,mu2,imu2sigma,"--->",processID
		
				p0=p1;

				#print mu,  1./isigma,  p0
										
				runSamples[0][0]=mu
				runSamples[1][0]=sigma
			else:
				continue

			if i%thinning==0:
 
 
				muList.append(runSamples[0][-1] )
				sigmaList.append( runSamples[1][-1] )
				pbar.update(len(muList)) #this adds a little symbol at each iteration
				
				
			if len(muList)>burn_in+sampleN:
			
				break
 			
 		
		pbar.finish()
		muList=np.array(muList[burn_in:]) 
		
		#for mu1 t	
		sigmaList=np.array(sigmaList[burn_in:-1])
		
		#print  np.mean(muList),   np.std(muList,ddof=1), np.mean( sigmaList ),  np.std( sigmaList,ddof=1), 
		
		return muList,sigmaList
		
		
	
	def getRMSWithNegative(self,dataArray, vMax=0.0,burn_in=50,thinning=10,sampleN=100):
		"""
		fitting the mean and rms of the truncated data
		"""
 
		#maxMu=
		mu0=np.mean(dataArray)
		
		sigma0=np.std(dataArray,ddof=1)
		
		
		minIsigma=1./sigma0/10.
		
		
		
		RMSs=[1.,  sigma0/2.]
		
		mu=   np.random.normal(mu0,sigma0 )
		
		isigma= np.random.exponential( 1./sigma0 )
 
		parameters=[mu,isigma]
		

		
		p0=logProbGaussTruncate(parameters,  dataArray, vMax)
			
		runSamples= [[mu],[isigma] ]



		muList=[mu]
		sigmaList=[1./isigma]
		
		
		
		widgets = ['MCMCSmapleDistance: ', Percentage(), ' ', Bar(marker='>',left='|',right='|'),
		           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
		
		pbar = ProgressBar(widgets=widgets, maxval=sampleN+burn_in+1)
		pbar.start() 
		
		
		
		for i in range(3000000):
			nextValue= getNextValues(runSamples,RMSs)
			mu,isigma=nextValue
			
			
			if isigma< 0 :
				continue
			
			p1=logProbGaussTruncate(nextValue,  dataArray, vMax)
			
			if np.isinf( p1):
				continue
		
		
			randomR=np.random.uniform(0,1)
			
			
			if p1>=p0 or p1-p0>np.log(randomR):
			#print disCloud,  imu1sigma,mu2,imu2sigma,"--->",processID
		
				p0=p1;

				#print mu,  1./isigma,  p0
										
				runSamples[0][0]=mu
				runSamples[1][0]=isigma
			else:
				continue

			if i%thinning==0:
 
 
				muList.append(runSamples[0][-1] )
				sigmaList.append(1./runSamples[1][-1] )
				pbar.update(len(muList)) #this adds a little symbol at each iteration
				
				
			if len(muList)>burn_in+sampleN:
			
				break
 			
 		
		pbar.finish()
		muList=np.array(muList[burn_in:]) 
		
		#for mu1 t	
		sigmaList=np.array(sigmaList[burn_in:-1])
		
		print np.mean(muList), np.mean( sigmaList )
		
		
		
	
		
		
		
if 0:
			
	testMCMC=MCMC()
	
	
	np.random.seed()
	a=np.random.normal(0,3,3000)
	vMax=0.
	a=a[a<vMax] 
	
	testMCMC.getRMSWithNegative(a,vMax=vMax)