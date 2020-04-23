from astropy.table import Table
import numpy as np

from astropy.table import Column

class myTB:

 

	def __init__(self ):
		pass
		
		
		
	@staticmethod
	def filterByRange(TB,colName,vRange):
		"""
		
		filter the table by the name and vRange
		
		"""
		

		
		if vRange[0]==None:
			
			doTB=TB.copy()
			
			doTB.add_index(colName)
			
			returnTB=doTB.loc[colName,  :vRange[1]]
	
 
			TestTB=Table() 
			
			for eachCol in returnTB.colnames:
				aaa=Column( returnTB[eachCol],name=eachCol )	
 
				TestTB.add_column(aaa )
			#gaiaOffCloudStars.remove_column(self.coint)
			return TestTB
			
 
		
		if vRange[1]==None:
			
			doTB=TB.copy()
			
			doTB.add_index(colName)
			
			returnTB=doTB.loc[colName,  vRange[0]:]
	
	
 
			TestTB=Table() 
			
			for eachCol in returnTB.colnames:
				aaa=Column( returnTB[eachCol],name=eachCol )	
 
				TestTB.add_column(aaa )
			#gaiaOffCloudStars.remove_column(self.coint)
			return TestTB
			
 
			
			
		#add inDex
		doTB=TB.copy()
		
		doTB.add_index(colName)
		
		returnTB=doTB.loc[colName,  min(vRange):max(vRange)]

		TestTB=Table() 
		
		for eachCol in returnTB.colnames:
			aaa=Column( returnTB[eachCol],name=eachCol )	

			TestTB.add_column(aaa )
		#gaiaOffCloudStars.remove_column(self.coint)
		return TestTB
	