import numpy as np

from astropy.table import Table
import matplotlib.pyplot as plt
from myMCMC import MCMC
from matplotlib import rc

from matplotlib import gridspec

 

from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

import os


from astropy.table import   Table
import numpy as np

class rawTable: #do not use mysql
	
	
 
	
	saveFITSName=""
	
	sourceName="sourceName" #sourceName has to be the id
	colnames=[  ]
 
	dataTypes=[  ]
 
	def __init__(self):
		pass
		
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
 
		return emptyTB[-1]
 

	def addRowToTB(self,newRow):
		"""
		Add a row to the fits collection
		"""
		
		#readExistingFITS 
		
		if os.path.exists(self.saveFITSName):
		
		
			preTB=Table.read(self.saveFITSName)
			os.remove(self.saveFITSName)
			
			nameList=list(preTB[ self.sourceName]  )
			
			if newRow[self.sourceName]  not in nameList:
				preTB.add_row(newRow)
				preTB.write(self.saveFITSName)
				print 'Source {} added!'.format(newRow[self.sourceName])
	
				
			else:
				indexOfSource=nameList.index( newRow[self.sourceName])
 
				for colS in self.colnames:
					preTB[indexOfSource][colS]=newRow[colS]
					
				preTB.write(self.saveFITSName)
				print 'Source {} updated!'.format(newRow[self.sourceName])
 
 
			
			
			
		else: #no sources
			preTB=self.getEmptyTB()
			preTB.add_row(newRow)
			preTB.write(self.saveFITSName)
			print 'Source {} added!'.format(newRow[self.sourceName])

	 
		print "Sources saved as: ",self.saveFITSName 
		
	def getRowByName(self,sName):
		"""
		Add a row to the fits collection
		"""
		
		#readExistingFITS 

		if str.strip(sName)=="":
			return None
		if not os.path.exists(self.saveFITSName):
			return self.getEmptyRow()
			
			
		preTB=Table.read(self.saveFITSName)

		
		checkRow=preTB[preTB[self.sourceName]==sName]  

		if len(checkRow)==0:
			print sName,"does not exist. Returing a new one..."
			return self.getEmptyRow()
		if len(checkRow)==1:
			print "Source", sName," found"
 
			return    checkRow[0]
			

  
	def getAllTB(self):
  	
		return Table.read(self.saveFITSName)


	def removeIndex(self,TB):
		
		TestTB=Table() 
		
		for eachCol in TB.colnames:
			
			
			TB[eachCol].mask=False #do not mask
			
			aaa=Column( TB[eachCol]  ,name=eachCol )	
			
 
			TestTB.add_column(aaa )
		return TestTB



	def ZZZ(self):
		pass
		
 