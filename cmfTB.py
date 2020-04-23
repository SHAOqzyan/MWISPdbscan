import numpy as np

from astropy.table import Table
import matplotlib.pyplot as plt
from myMCMC import MCMC
from matplotlib import rc

from matplotlib import gridspec


from parameterTB import parameter

from mpl_toolkits.axes_grid.anchored_artists import AnchoredText


class myCMF():
	
	"""
	
	define the table of calculated CMF foe each case
	"""
	
	
	#colnames
	
	
	case="case"
	residual="residual"
	
	peakMin="peakMin"
	Nclump="Nclump"
	totalMass="totalMass"
	
	alpha1="alpha1"
	alpha1Std="alpha1Std"
	
	alpha2="alpha2"
	alpha2Std="alpha2Std"

	Mturn="Mturn"
	MturnStd="MturnStd"

	Mvirial="Mvirial"	
	MvirialRatio="MvirialRatio"

	note="note"
	
	sourceName="source"
	
	source=None
	
	
	names=( case,residual, peakMin, Nclump, totalMass, alpha1,alpha1Std, alpha2, alpha2Std, Mturn, MturnStd, Mvirial,MvirialRatio, note,sourceName )
	dtypes=('i8','f8','f8','i8','f8','f8','f8','f8','f8','f8','f8', 'f8', 'f8',  'S100','S10')
	
	


	#W345CMFFile="W345CMF.fit"
	#G216CMFFile="G216CMF.fit"


	W345CMFFile="W345CMFnoUpper.fit"
	G216CMFFile="G216CMFnoUpper.fit"



	def __init__(self,source="G216"):
		self.source=source
		
		self.doPara=parameter()
		
		self.paraTB=self.doPara.getAllParameter()
		
	def getEmptyTB(self):
		
		"""
		build a TB with columns
		"""
		
		return Table(names=self.names,dtype=self.dtypes)
		
		
	def getDefaultValues(self):
	
		return [0,0.,0.,0.,0,0.,0.,0.,0.,0.,0.,  0.,0.,"",self.source  ]
		
	def getEmptyRow(self):
		
		"""
		provide an empyt row,
		"""
		emptyTB=self.getEmptyTB()
		rowData=self.getDefaultValues()
		
		emptyTB.add_row(rowData)

		return emptyTB[-1]


	def printCMFTB(self):
		
		"""
		
		print CMF fitting result
		
		"""

	def checkParaBiase(self,source='G216',checkResult="alpha2"):
		
		fitTBG216=Table.read(self.G216CMFFile)
 
		fitTBW345=Table.read(self.W345CMFFile)
		
		from astropy.table import  hstack
		if source=="G216":
			bigTB=hstack([fitTBG216,self.paraTB])
		else:
			bigTB=hstack([fitTBW345,self.paraTB])

		
		#examine alpha2
		
		fig, axs = plt.subplots( figsize=(16,12) )  
		AX = gridspec.GridSpec(5,2)
		AX.update(wspace = 0.2, hspace = 0.4)
		indexDraw=0

		for i, drawP in enumerate(self.doPara.names):
			
			if drawP=="NPAD":
				#indexDraw=indexDraw+1
				continue
			
			
				
				
			
			ax  = plt.subplot(AX[indexDraw/2,indexDraw%2])
	
			ax.scatter(bigTB[drawP],  bigTB[checkResult] ,s=5 )
		
			#at = AnchoredText(drawP, loc=1, frameon=False,prop={"color":"black"} )
			#at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			#ax.add_artist(at)

			ax.set_ylabel(checkResult)
			ax.set_xlabel(drawP)

		
		
			indexDraw=indexDraw+1

	
		
		saveName="bias{}{}.eps".format(source,checkResult)
		
		plt.savefig(saveName,bbox_inches='tight')

		#plt.show()



	


	def fittingGaussian(self,source="G216"):
		
		"""
		sfsfddf  ssdfsdfsdf
		"""
 
		fitTBG216=Table.read(self.G216CMFFile)

		print fitTBG216


		#print fitTBG216
 
		fitTBW345=Table.read(self.W345CMFFile)
		
		print fitTBW345

		print "case number",len(fitTBG216)
			
			
		print "mean MvirialRatio G216",np.mean(fitTBG216["MvirialRatio"] )
		print  "mean MvirialRatio W345",np.mean(fitTBW345["MvirialRatio"] )
	 
		print  "total mass G216",np.mean(fitTBG216["totalMass"] )
		print  "total mass W345",np.mean(fitTBW345["totalMass"] )

		
		print  "total mass G216",np.mean(fitTBG216["Mvirial"] )
		print  "total mass W345",np.mean(fitTBW345["Mvirial"] )

			
			
		NclumpG216=fitTBG216[self.Nclump]
		NclumpW345=fitTBW345[self.Nclump]

		alpha2ListG216=fitTBG216[self.alpha2]
		alpha2ErrorListG216=fitTBG216[self.alpha2Std]

		alpha2ListW345=fitTBW345[self.alpha2]
		alpha2ErrorListW345=fitTBW345[self.alpha2Std]


			
		alpha1ListG216=fitTBG216[self.alpha1]
		alpha1ErrorListG216=fitTBG216[self.alpha1Std]

		alpha1ListW345=fitTBW345[self.alpha1]
		alpha1ErrorListW345=fitTBW345[self.alpha1Std]


		MturnListG216=fitTBG216[self.Mturn]
		MturnErrorListG216=fitTBG216[self.MturnStd]

		MturnListW345=fitTBW345[self.Mturn]
		MturnErrorListW345=fitTBW345[self.MturnStd]


		totalMassG216= fitTBG216[self.totalMass]
		totalMassW345= fitTBW345[self.totalMass]

		gBoundG216= fitTBG216[self.MvirialRatio]
		gBoundW345= fitTBW345[self.MvirialRatio]




		weightsG216=1./alpha2ErrorListG216**2
		weightsG216=weightsG216/np.sum(weightsG216)

		weightsW345=1./alpha2ErrorListW345**2
		weightsW345=weightsW345/np.sum(weightsW345)




		mean,Part1std=self.weighted_avg_and_std(alpha2ListG216, weightsG216)
		Part2std=np.average( alpha2ErrorListG216 , weights= weightsG216 )
		G216Alpha2mean=mean
		G216Alpha2error=np.sqrt(Part1std**2+ Part2std**2)
		print "G216, mean and std:",G216Alpha2mean,G216Alpha2error




		#Part2std=np.average( alpha2ErrorListG216**2, weights= weightsG216 )

		#print "G216, mean and std:",mean,np.sqrt(Part1std**2+ Part2std)


		#MeanalphaG216=mean
		#MeanalphaG216Error=np.sqrt(Part1std**2+ Part2std)

		#print "weighted mean and std (G216):",mean,Part1std
		
		
		#print np.mean( alpha2ListG216 )
		mean,Part1std=self.weighted_avg_and_std(alpha2ListW345, weightsW345)
		
		Part2std=np.average( alpha2ErrorListW345 , weights= weightsW345 )
		W345Alpha2mean=mean
		W345Alpha2error= np.sqrt(Part1std**2+ Part2std**2)
		#print "G216, mean and std:",G216Alpha2mean,G216Alpha2error

		
		print "W345, mean and std:",W345Alpha2mean,W345Alpha2error
		
		#Part2std=np.average( alpha2ErrorListW345**2, weights= weightsW345 )

		#print "W345, mean and std:",mean,np.sqrt(Part1std**2+ Part2std)

		#print len(weightsG216)

		#print "Mean, std",np.mean(alpha2List ),np.std(alpha2List,ddof=1 )
 




		


		#draw all four parameters
		
		
		#print fitTBG216

		#a= ( alpha2ListG216-np.mean(alpha2ListG216) ) *  ( alpha2ListW345-np.mean(alpha2ListW345) ) 
		
		#print np.mean(a)
		
		#domcmc.getMeanAndStd([1,4], [0.5,0.5] )
		#domcmc=MCMC()
		
		#domcmc.getMeanAndStd(np.arange(400)*0.+2, np.arange(400)*0.+2 )

		
		if 0:
			print "MCMC G216"
			muList,sigmaList=domcmc.getMeanAndStd(alpha2ListG216.data, alpha2ErrorListG216.data )
			
			print np.mean(muList),np.std(muList,ddof=1)
			
			print np.mean(sigmaList)*np.sqrt( len(alpha2ListG216) ),np.std(sigmaList,ddof=1)


			print "MCMC W345"
			muList,sigmaList=domcmc.getMeanAndStd(alpha2ListW345.data, alpha2ErrorListW345.data )
			
			print np.mean(muList),np.std(muList,ddof=1)
			print np.mean(sigmaList)*np.sqrt( len(alpha2ListW345) ),np.std(sigmaList,ddof=1)
			
			
			
			
		if  0:
			muList,sigmaList=domcmc.getMeanAndStd(np.arange(400)*0.+2, np.arange(400)*0.+1 )
	
			print np.mean(muList),np.std(muList,ddof=1)
			
			print np.mean(sigmaList)*20 ,np.std(sigmaList,ddof=1)



		
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


		
		rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'] })
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
 
		rc('text', usetex=True)
		
		
		fontLabel=12
		
		#plt.rcParams.update({'font.size': 22})
		#plt.rc('font', size=20)
		######################
		#axNclump= #plt.subplot2grid((2,2), (0, 0) ,rowspan=1,colspan=1)
		axNclump.scatter(NclumpG216,NclumpW345,color='g',s=5)
		#axNclump.set_aspect("equal")

		axNclump.set_xlabel("Clump number of G216" ,fontsize=fontLabel)
		axNclump.set_ylabel("Clump number of W345" ,fontsize=fontLabel)

		
		######################
		#axAlpha1=plt.subplot2grid((2,2), (0, 1) ,rowspan=1,colspan=1)
		axAlpha1.errorbar(alpha1ListG216,alpha1ListW345,xerr=alpha1ErrorListG216, capsize=1.0, yerr= alpha1ErrorListW345,marker='o', ms=2.5, ls = 'none',color='blue',lw=0.4)
		#axAlpha1.set_aspect("equal")
		axAlpha1.set_xlabel(r"$\alpha_1$ of G216" ,fontsize=fontLabel)
		axAlpha1.set_ylabel(r"$\alpha_1$ of W345" ,fontsize=fontLabel)




		#axAlpha2=plt.subplot2grid((2,2), (1, 0) ,rowspan=1,colspan=1)
		axAlpha2.errorbar(alpha2ListG216,alpha2ListW345,xerr=alpha2ErrorListG216, capsize=1.0, yerr= alpha2ErrorListW345,marker='o', ms=2.5, ls = 'none',color='blue',lw=0.4,zorder=1)
	#	axAlpha2.set_aspect("equal")
		axAlpha2.errorbar([G216Alpha2mean] ,[W345Alpha2mean],xerr=[G216Alpha2error], capsize=1.0, yerr= [W345Alpha2error],marker='o', ms=2.5, ls = 'none',color='red',lw=0.4)

	
		axAlpha2.set_xlabel(r"$\alpha_2$ of G216" ,fontsize=fontLabel)
		axAlpha2.set_ylabel(r"$\alpha_2$ of W345" ,fontsize=fontLabel)



		#axMturn=plt.subplot2grid((2,2), (1, 1) ,rowspan=1,colspan=1)
		axMturn.errorbar(MturnListG216,MturnListW345,xerr=MturnErrorListG216, capsize=1.0, yerr= MturnErrorListW345,marker='o', ms=2.5, ls = 'none',color='blue',lw=0.4)
		#axMturn.set_aspect("equal")
		axMturn.set_xlabel(r"$M_{\rm turn}$ of G216 ($M_\odot$)" ,fontsize=fontLabel)
		axMturn.set_ylabel(r"$M_{\rm turn}$ of W345 ($M_\odot$)" ,fontsize=fontLabel)
		

		#totalMassG216= fitTBG216[self.totalMass]
		#totalMassW345= fitTBW345[self.totalMass]

		#gBoundG216= fitTBG216[self.MvirialRatio]
		#gBoundW345= fitTBW345[self.MvirialRatio]
		
		#axTotal.scatter( )
		
		axTotal.scatter(totalMassG216/1e4,totalMassW345/1e4,color='g',s=5)
		
		#axNclump.set_aspect("equal")
		
		axTotal.set_xticks([1.0,1.5,2.0])
		axTotal.set_yticks([4,6,8,10,12])

		
		#locs, labels = plt.xticks()            # Get locations and labels
		#plt.xticks(ticks, [labels], **kwargs)  # Set locations and labels

		#xticks(np.arange(0, 1, step=0.2))


		axTotal.set_xlabel(r" Total $M_{\rm LTE}$ of G216 ($10^4\ M_\odot$)" ,fontsize=fontLabel)
		axTotal.set_ylabel(r" Total $M_{\rm LTE}$ of W345 ($10^4\ M_\odot$)" ,fontsize=fontLabel)


		axRatio.scatter(gBoundG216 ,gBoundW345  ,color='g',s=5,label='Ratio of mass in gravitationally bound')
		axRatio.set_xticks([0.1,0.2,0.3,0.4])
		axRatio.set_yticks([0.6,0.7,0.8,0.9,1.0])
		
		axRatio.set_xlabel(r"Ratio of $M_{\rm LTE}$ with $\alpha_{\rm vir}<2$ (G216)",fontsize=fontLabel)
		axRatio.set_ylabel(r"Ratio of $M_{\rm LTE}$ with $\alpha_{\rm vir}<2$ (W345)",fontsize=fontLabel)
		#axRatio.legend( loc=4)
		#plt.tight_layout()
		plt.savefig("G216W345Compare.pdf",bbox_inches='tight')

		
		
		return 
		
		f=plt.figure(figsize=(8,6))
		ax=f.add_subplot(111)

		#ax.errorbar(alpha2ListG216,alpha2ListW345,xerr=alpha2ErrorListG216, capsize=1.0, yerr= alpha2ErrorListW345,marker='o', ms=4, ls = 'none',color='blue',lw=0.8)

		#ax.errorbar(X,Y,yerr=[lowerY,upperY],c='b',marker='o',capsize=1.5,elinewidth=0.8,lw=1,label="{} clumps in total".format(len(coreTB)))

		#ax.errorbar(alpha1ListG216,alpha1ListW345,xerr=alpha1ErrorListG216, capsize=1.0, yerr= alpha1ErrorListW345,marker='o', ms=4, ls = 'none',color='blue',lw=0.8)
		#ax.errorbar(MturnListG216,MturnListW345,xerr=MturnErrorListG216, capsize=1.0, yerr= MturnErrorListW345,marker='o', ms=2.5, ls = 'none',color='blue',lw=0.5)


		#print np.corrcoef(alpha1ListG216, alpha1ListW345)
		#print np.corrcoef(alpha2ListG216, alpha2ListW345)
		#print np.corrcoef(MturnListG216, MturnListW345)
		print np.corrcoef(NclumpG216, NclumpW345)


		###

		ax.plot(NclumpG216,NclumpW345,'b.')
		
		#ax.set_xlim(1.5,3)
		#ax.set_ylim(0,200)

	#	ax.set_xlim(1.5,3)
		#ax.set_ylim(0.7,2.2)

		#ax.plot( alpha2ListG216, alpha2ListW345,'b.'  )
		
		plt.savefig("G216W345.eps",bbox_inches='tight')
		
		
		
		
		

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

	def ZZZ(self):
		pass
		
		
if 0: 
	
	docmf=myCMF("G216")
		
	for eachP in docmf.names:
		docmf.checkParaBiase(source="G216",checkResult= eachP)
		docmf.checkParaBiase(source="W345",checkResult= eachP)



if  0: 
	docmf=myCMF("G216")

	docmf.fittingGaussian()



	#docmf2=myCMF("W345")
	#newRow= docmf2.fittingGaussian(  source="W345") 


	#print newRow[docmf.sourceName]