import numpy as np

from astropy.table import Table
import matplotlib.pyplot as plt
from myMCMC import MCMC
from matplotlib import rc

from matplotlib import gridspec


from parameterTB import parameter

from mpl_toolkits.axes_grid.anchored_artists import AnchoredText



from rawTB import rawTable

  
import numpy as np

from astropy.table import Table
import matplotlib.pyplot as plt
from myMCMC import MCMC
from matplotlib import rc

from matplotlib import gridspec

 

from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

import os

 

class myHII(rawTable): #do not use mysql
	#fitsPath="/home/qzyan/WORK/projects/maddalena/HIICat/"

	#saveFITSName=fitsPath+


	#fitsPath="/home/qzyan/WORK/projects/maddalena/HIICat/"
	
	#saveFITSName=fitsPath+"HIIG210.fit"
	#sourceName="sourceName" #WISE or IRAS 

	l="l"
	
	b="b"
	
	radius='radius'
	
	
	CO13L="CO13L"
	CO13B="CO13B"
	CO13Radius="CO13Radius"  #the region of 13CO used to get the spectrum of 13CO to fit a gaussian
 
	
	
	velocity="velocity" # gauss fitting with 13CO
	
	vStd="vStd" # gauss fiting with 12CO

	CO13Peak="CO13Peak"
	CO13VFit="CO13VFit"

	arm="arm"
	
	membership='membership' #which ground does this source belongs
	catalog="catalog"

	
	note="note"
	
	colnames=[  rawTable.sourceName,l,b,  CO13L,CO13B, CO13Radius,  radius,velocity,CO13VFit,vStd,CO13Peak, arm,membership,catalog,note]
	
	dataTypes=[ 'S30',  float,float,float,float,float,float,float,float,float,float,float, float,'S20', 'S200']



	def __init__(self):
		rawTable.saveFITSName="/home/qzyan/WORK/projects/maddalena/HIICat/"+"HIIG210.fit"
 
		self.saveLatex= "/home/qzyan/WORK/projects/maddalena/HIICat/HIILatexG210.txt"

		pass

 
 
	def printHIILatex(self):
		
		"""
		
		print HII region Latex
		"""
		
		#get allHII
		
		printTB=self.getAllTB()
		
		N=len(printTB)
		
		#rowN=N/2+N%2
		rowN=N
		 
		 
		GoodHII=printTB[ printTB[self.vStd]>0  ]
		badHII=printTB[ printTB[self.vStd]<=0  ]

		GoodHII.sort(self.sourceName)
		
		GoodHII.reverse()
		
		badHII.sort(self.sourceName)
		
		badHII.reverse()
		
		from astropy.table import vstack
		
		
		
		printTB=vstack([ GoodHII,badHII])
		
		
		LocalN=0
		PerseusN=0
		OuterN=0
		NoCO13=0
		#
 		
		f=open(self.saveLatex, 'w')
 
		for eachHII in printTB:
		
		
			sourceName= eachHII[myHII.sourceName]  
			#print sourceName
			
			
			if not sourceName.startswith("G"):
				sourceName="IRAS"+sourceName
		
			nameStr= "{:>20}".format(sourceName)
		
			positionStr="{:>8.3f} & {:>8.3f} ".format(eachHII[myHII.l] ,eachHII[myHII.b]  )
	
			wisetRadiusStr="{:>4.1f}".format(eachHII[myHII.radius]*60 ) #arcmin
			if not sourceName.startswith("G"):

				wisetRadiusStr="{:>4}".format("--" ) #arcmin


			CO13RegionPosition= "{:>8.3f} & {:>8.3f} ".format(eachHII[myHII.CO13L] ,eachHII[myHII.CO13B]  )
			CO13RegionRadius=  "{:>4.1f}".format(eachHII[myHII.CO13Radius]*60 ) #arcmin


 
 
			gaussian= " {:>8.2f} & {:>8.2f} & {:>8.2f}  ".format(  eachHII[self.CO13VFit], eachHII[self.vStd], eachHII[self.CO13Peak]   ) 

			if eachHII[self.vStd]<=0:
				gaussian= " {:>8} & {:>8} & {:>8}  " .format( '--', '--', '--' )
				CO13RegionPosition= "{:>8} & {:>8} ".format(  '--', '--'  )
				CO13RegionRadius=  "{:>4}".format( '--'  ) #arcmin

			ArmStr="???????????"
 
			if   eachHII[self.CO13VFit] >38:
				ArmStr="Outer"

				OuterN=OuterN+1
			if  14<eachHII[self.CO13VFit] <=38:
				ArmStr="Perseus"
				PerseusN=PerseusN+1
			if  0< eachHII[self.CO13VFit] <=14:
				ArmStr="Local"
				LocalN=LocalN+1
				
				print nameStr,eachHII[self.CO13VFit]
				
			if eachHII[self.vStd]<=0:
				ArmStr="--"
				NoCO13=NoCO13+1
				
				
			ArmStr="{:>10}".format(ArmStr)

			classOfHII= eachHII[self.catalog]

			if classOfHII=="0":
				classOfHII="--"

			latexLine="{} & {}  & {}  & {} &{}  & {} & {} & {:>3} \\\\ ".format(nameStr, positionStr,wisetRadiusStr, CO13RegionPosition,CO13RegionRadius ,gaussian, ArmStr,classOfHII) 
	
			f.write(latexLine   + os.linesep)

		f.close()	

 


		print "Local Arm HII:",LocalN
		
		print "Perseus Arm HII:",PerseusN
		
		print "Outer Arm HII:",OuterN
		print "no Arm HII:",NoCO13
		print "In total", LocalN+PerseusN+OuterN+NoCO13,len(printTB)