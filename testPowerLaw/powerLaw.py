import numpy as np
import random as random
import matplotlib.pyplot as plt
import time
from progressbar import * 
# set the seed



def random_PowerLaw(N, alpha, M_min, M_max):
    """
    Draw random samples from a power-law defined over M_min, M_max.
        dN/dM = Z x M ** -alpha
    INPUTS
    ------
    N: int
        number of samples.
    alpha: float
        power-law index.
    M_min: float
        lower bound of mass interval.
    M_max: float
        upper bound of mass interval.
    OUTPUTS
    -------
    masses: ndarray[float, ndim=1]
        list of N masses drawn from the power-law
    """
    beta = 1. - alpha
    x = np.random.uniform(0., 1., N)
    if beta == 0:
        y = M_min * np.exp( np.log(M_max / M_min) * x )
    else:
        y = ((M_max ** beta - M_min ** beta) * x + M_min ** beta) ** (1. / beta)
    return y







class myPowerLaw(): 
 
	def __init__(self,alpha, M_min, M_max): 
		self.alpha=alpha
		self.M_min=M_min
		self.M_max=M_max
		self.beta=1.-self.alpha

 


	def testPLPDF(self):
		a=np.linspace(self.M_min,self.M_max,1000)
		pdfs=[]
		for eacha in a:
			pdfs.append(self.powerLawPDF(eacha)  )

		plt.figure(1)
		plt.plot(a,pdfs)
		plt.show()

	def powerLawPDF(self,x):
	
		if x<self.M_min or x>self.M_max:
			return 0.
			
		if self.alpha==1.:
			Z=1./(np.log(self.M_max/self.M_min))
			return 		Z*x**(-self.alpha)	
			
		 
	
		Z=self.M_max**self.beta-self.M_min**self.beta
		Z=self.beta/Z
		
		return Z*x**(-self.alpha)
		
		
	def getSample(self,N,burn_in=50,thinning=50):
		"""
		return an arrry of samples
		
		"""

		#sample this PDF with MCMC
		
		np.random.seed()
		X0=np.random.uniform(self.M_min,self.M_max) 
		
		
		P0=self.powerLawPDF(X0)
		
		
		#interval= self.M_max/10.
		interval=10
		
		
		xList=[X0]
		
		
		for i in range(10000000):
 
			xCand= xList[-1]+ np.random.normal(0,interval)
					
			P1= self.powerLawPDF(xCand)

			if P1==0.:

				xList.append(xList[-1])
				continue
 
			u=np.random.uniform(0,1)
			
			if P1>P0 or P1/P0>u:
				
				xList.append(xCand)
				P0=P1
 
			else:
				xList.append(xList[-1])
			
			if len(xList)/thinning+burn_in>N:
				break
		returnList=xList[::thinning]
		return returnList[burn_in:]
	
	
	
 	def getNextValues(self,runSamples, indexPara,RMSs):
 
		currentValues=[ ]
		
		for a in range(len(runSamples)):
			currentValues.append(runSamples[a][-1])
	
		currentValues[indexPara]=currentValues[indexPara] +np.random.normal(0,RMSs[indexPara])
		
		return currentValues

	#fitting function
	def calPowerLawProbLog3p(self,theta, dataArraym ,minV,maxV  ):
		
		
		alpha1,alpha2,Mt =theta #Mt is the turnover mass of molecular cores
		
		
		#should normallize them together
		
		part1=dataArraym[dataArraym<Mt]
		part2=dataArraym[dataArraym>=Mt]
		
		beta1=1-alpha1
		beta2=1-alpha2


 

 		normalFactor1=beta1/(  Mt**beta1-minV**beta1 )
 		normalFactor2=beta2/(  maxV**beta2-Mt**beta2 )
 
		return    len(part1)*np.log(normalFactor1)-alpha1*np.sum(np.log( part1 )) +len(part2)*np.log(normalFactor2)-alpha2*np.sum(np.log( part2 ))





	def calPowerLawProbLog5p(self,theta, dataArraym ,minV,maxV  ):
		
		
		alpha1,alpha2,alpha3,Mt1,Mt2 =theta #Mt is the turnover mass of molecular cores
		
		
		#should normallize them together
		
		part1=dataArraym[dataArraym<Mt1]
		
		part2=dataArraym[dataArraym>=Mt1]
		part2=part2[part2<Mt2]

		
		part3=dataArraym[dataArraym>=Mt2]
		
		
		#print len(part1),len(part2),len(part3),len(dataArraym),Mt1,Mt2
		
		
		beta1=1-alpha1
		beta2=1-alpha2
		beta3=1-alpha3


 		normalFactor1=beta1/(  Mt1**beta1-minV**beta1 )
 		normalFactor2=beta2/(  Mt2**beta2-Mt1**beta2 )
 		normalFactor3=beta3/(  maxV**beta3-Mt2**beta3 )

		return    len(part1)*np.log(normalFactor1)-alpha1*np.sum(np.log( part1 )) +len(part3)*np.log(normalFactor3)-alpha3*np.sum(np.log( part3 ))+len(part2)*np.log(normalFactor2)-alpha2*np.sum(np.log( part2 ))


	def fitPowerLawWithMCMC5p(self,dataArray,sampleN=200,burn_in=50,thinning=10):
		"""
		fit a power law with MCMC
		"""
		
		dim=5
		minV= min(dataArray)
		maxV= max(dataArray)
		
		#mass,massN=self.madeCount(dataArray)
		
		print "Moddling 3 parameters."
		
		print "The minimum and maximum masses are (in solar mass)",minV,maxV
		
		mass=np.sort(dataArray)
		
		
		maxV=mass[-10]
		
		
		np.random.seed()

		#alpha1=np.random.exponential(1)
		#alpha2=np.random.exponential(1.5)
		#alpha3=np.random.exponential(2)

		alpha1=np.random.uniform(0,3)
		alpha2=np.random.uniform(alpha1,3)
		alpha3=np.random.uniform(alpha2,3)

		


		Mt1=np.random.uniform(minV,maxV)
		Mt2=np.random.uniform(Mt1,maxV)

		
		#Mt=7000
		
		
		theta=[ alpha1,alpha2,alpha3,Mt1,Mt2] #

		p0= self.calPowerLawProbLog5p(theta,mass,minV,maxV)

 
		
		alpha1List=[ alpha1]
		alpha2List=[ alpha2]
		alpha3List=[ alpha3]

		Mt1List=[ Mt1]
		Mt2List=[ Mt2]

		runSamples=[ alpha1List,alpha2List,alpha3List, Mt1List, Mt2List]
		RMSs=[0.5,0.5, 0.5, maxV/20 , maxV/20]
		
		acceptN=0
		widgets = ['MCMC Smaple CMF: ', Percentage(), ' ', Bar(marker='>',left='|',right='|'),
		           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
		
		pbar = ProgressBar(widgets=widgets, maxval=thinning*(sampleN+burn_in))
		pbar.start() 
		
		
		for i in range(10000):
			
			
			
			paraJ=0
			while paraJ<dim:

			#for paraJ in range(3):
				valueCand=self.getNextValues(runSamples,paraJ,RMSs)
 
			#theta=[   ]
				alpha1,alpha2,alpha3,Mt1,Mt2=valueCand
 
				if  alpha1<0 or alpha2<0 or alpha3<0 or Mt1<minV or Mt2>maxV or Mt1>Mt2 :
					#paraJ=paraJ-1
					continue

				if alpha1>alpha2 or alpha2>alpha3:
					continue

				p1= self.calPowerLawProbLog5p(valueCand,mass,minV,maxV)

				randomR=np.random.uniform(0,1)
 
				if p1>p0 or p1-p0>np.log(randomR):
					#print disCloud,  imu1sigma,mu2,imu2sigma,"--->",processID
					acceptN=acceptN+1
					p0=p1;
					runSamples[paraJ].append( valueCand[paraJ] )
				
				else:
					runSamples[paraJ].append( runSamples[paraJ][-1] )
				paraJ=paraJ+1

			if len(runSamples[0])>thinning*(sampleN+burn_in):
				break
			pbar.update(len(runSamples[0]))
			

			if i %20 ==0:
				print runSamples[0][-1],runSamples[1][-1],runSamples[2][-1],runSamples[3][-1],runSamples[4][-1] 

		pbar.finish()


		alpha1Array=runSamples[0][::thinning] 
		alpha1Array=alpha1Array[burn_in:]
		


		alpha2Array=runSamples[1][::thinning] 
		alpha2Array=alpha2Array[burn_in:]
		
		mtArray=runSamples[2][::thinning] 
		mtArray=mtArray[burn_in:]
 
		print  "alpha1, alpha2, Mt:",np.mean(alpha1Array),np.mean(alpha2Array),np.mean(mtArray)
		#sample theta
		return [np.mean(alpha1Array),np.mean(alpha2Array),np.mean(mtArray)]
	 










	#fitting function
	def calPowerLawProbLog1p(self,theta, dataArraym ,minV,maxV  ):
		
	 
		
		alpha1  =theta[0] #Mt is the turnover mass of molecular cores

		beta1=1. - alpha1

 

 		normalFactor=beta1/(  maxV**beta1-minV**beta1 )
  
		return    len(dataArraym)*np.log(normalFactor)-alpha1*np.sum(np.log( dataArraym)  )



  
	def fitPowerLawWithMCMC1p(self,dataArray,sampleN=5000,burn_in=50,thinning=10):
		"""
		fit a power law with MCMC
		"""
		mass=np.sort(dataArray)

		dim=1
		minV= min(dataArray)
		maxV= max(dataArray)
		np.random.seed()
		alpha1=np.random.exponential(2)
		theta=[ alpha1 ]

		p0= self.calPowerLawProbLog1p(theta,mass,minV,maxV)
		alpha1List=[ alpha1]
		runSamples=[ alpha1List  ]


		RMSs=[0.5,0.5,5]
		
		acceptN=0
		widgets = ['MCMCSmapleDistance: ', Percentage(), ' ', Bar(marker='>',left='|',right='|'),
		           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
		
		pbar = ProgressBar(widgets=widgets, maxval=thinning*(sampleN+burn_in))
		pbar.start() 
		
		
		for i in range(10000):

			
			paraJ=0
			while paraJ<dim:

			#for paraJ in range(3):
				valueCand=self.getNextValues(runSamples,paraJ,RMSs)
 
			#theta=[   ]
				alpha1 =valueCand[0]
 
				if  alpha1<0  :
					#paraJ=paraJ-1
					continue
 
				p1= self.calPowerLawProbLog1p(valueCand,mass,minV,maxV)

				randomR=np.random.uniform(0,1)
 
				if p1>=p0 or p1-p0>np.log(randomR):
					#print disCloud,  imu1sigma,mu2,imu2sigma,"--->",processID
					acceptN=acceptN+1
					p0=p1;
					runSamples[paraJ].append( valueCand[paraJ] )
				
				else:
					runSamples[paraJ].append( runSamples[paraJ][-1] )
				paraJ=paraJ+1
				
			if len(runSamples[0])>thinning*(sampleN+burn_in):
				break
			pbar.update(len(runSamples[0]))
			
 
		pbar.finish()


		alpha1Array=runSamples[0][::thinning] 
		alpha1Array=alpha1Array[burn_in:]
		
 
 
		print  "alpha1 :",np.mean(alpha1Array) 
		#sample theta
		return  np.mean(alpha1Array) 
 



	def fitPowerLawWithMCMC3p(self,dataArray,sampleN=200,burn_in=50,thinning=10):
		"""
		fit a power law with MCMC
		"""
		
		dim=3
		minV= min(dataArray)
		maxV= max(dataArray)
		
		#mass,massN=self.madeCount(dataArray)
		
		print "Moddling 3 parameters."
		
		print "The minimum and maximum masses are (in solar mass)",minV,maxV
		
		mass=np.sort(dataArray)
		
		
		maxV=mass[-10]
		
		
		np.random.seed()

		alpha1=np.random.exponential(1)
		alpha2=np.random.exponential(2)
		Mt=np.random.uniform(minV,maxV)
		Mt=7000
		
		
		theta=[ alpha1,alpha2,Mt]

		p0= self.calPowerLawProbLog3p(theta,mass,minV,maxV)

 
		
		alpha1List=[ alpha1]
		alpha2List=[ alpha2]
 
		MtList=[ Mt]
		
		runSamples=[ alpha1List,alpha2List, MtList ]
		RMSs=[0.5,0.5,  maxV/20 ]
		
		acceptN=0
		widgets = ['MCMC Smaple CMF: ', Percentage(), ' ', Bar(marker='>',left='|',right='|'),
		           ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
		
		pbar = ProgressBar(widgets=widgets, maxval=thinning*(sampleN+burn_in))
		pbar.start() 
		
		
		for i in range(10000):
			
			
			
			paraJ=0
			while paraJ<dim:

			#for paraJ in range(3):
				valueCand=self.getNextValues(runSamples,paraJ,RMSs)
 
			#theta=[   ]
				alpha1,alpha2,Mt=valueCand
 
				if  alpha1<0 or alpha2<0 or Mt<minV or Mt>maxV :
					#paraJ=paraJ-1
					continue
 
				p1= self.calPowerLawProbLog3p(valueCand,mass,minV,maxV)

				randomR=np.random.uniform(0,1)
 
				if p1>p0 or p1-p0>np.log(randomR):
					#print disCloud,  imu1sigma,mu2,imu2sigma,"--->",processID
					acceptN=acceptN+1
					p0=p1;
					runSamples[paraJ].append( valueCand[paraJ] )
				
				else:
					runSamples[paraJ].append( runSamples[paraJ][-1] )
				paraJ=paraJ+1

			if len(runSamples[0])>thinning*(sampleN+burn_in):
				break
			pbar.update(len(runSamples[0]))
			

			if i %20 ==0:
				print runSamples[0][-1],runSamples[1][-1],runSamples[2][-1] 

		pbar.finish()


		alpha1Array=runSamples[0][::thinning] 
		alpha1Array=alpha1Array[burn_in:]
		


		alpha2Array=runSamples[1][::thinning] 
		alpha2Array=alpha2Array[burn_in:]
		
		mtArray=runSamples[2][::thinning] 
		mtArray=mtArray[burn_in:]
 
		print  "alpha1, alpha2, Mt:",np.mean(alpha1Array),np.mean(alpha2Array),np.mean(mtArray)
		#sample theta
		return [np.mean(alpha1Array),np.mean(alpha2Array),np.mean(mtArray)]
	def zzz(self):
		pass




if 0:
	
	
	
	breakP=2


			
	myPL1=myPowerLaw(2.3,breakP,30)
 
	sample1=myPL1.getSample(1000,thinning=20)
	
	
	myPL2=myPowerLaw(1.3,0.5,breakP)
	sample2=myPL2.getSample(1000,thinning=20)
	
	myPL3=myPowerLaw(3.3,30,100)
 
	sample3=myPL3.getSample(1000,thinning=20)
	
	mixedData=np.concatenate( (sample2,sample1,sample3) )
	mixedData2=np.concatenate( (sample2,sample1) )

	
	#myPL1.fitPowerLawWithMCMC5p(mixedData)
	
	
	myPL2.fitPowerLawWithMCMC3p(mixedData2)
	if 0:
		plt.figure(1)
		
		
		plt.hist(  mixedData,bins=100)
		
		#plt.hist( sample2   ,bins=50)
		
		
		plt.show()
		
