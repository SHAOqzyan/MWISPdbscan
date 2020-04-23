import numpy as np

from astropy.table import Table
 




class parameter():
	
	"""
	
	define the table of calculated CMF foe each case
	"""
	
	
	#colnames
	
  
	ID="ID"
	NPAD="NPAD"
	MODELLIM="MODELLIM"
	
	S0="S0"
	MAXSKIP="MAXSKIP"
	MAXNF="MAXNF"
	
	WMIN="WMIN"
	SC="SC"
	
	SB="SB"
	SA="SA"

 
	WWIDTH="WWIDTH"
 
 
	names=( ID,NPAD,MODELLIM, S0, MAXSKIP, MAXNF, WMIN,SC, SB, SA, WWIDTH  )
	dtypes=('i8','i8','f8','f8','f8','i8','f8','f8','f8','f8','f8'  )
	
 
 

	def __init__(self,source="G216"):
		self.source=source
		
		
	def getEmptyTB(self):
		
		"""
		build a TB with columns
		"""
		
		return Table(names=self.names,dtype=self.dtypes)
		
	def getDefaultValues(self):
	
		return [0,0,0.,0.,0.,0,0.,0.,0.,0.,0.  ]
		
	def getEmptyRow(self):
		
		"""
		provide an empyt row,
		"""
		emptyTB=self.getEmptyTB()
		rowData=self.getDefaultValues()
		
		emptyTB.add_row(rowData)

		return emptyTB[-1]
		
		
		
		

	def getOnRow(self,i):
		"""
		get the ith paramether and return the row
		"""


		#construct the file:
		
		paraRow=self.getEmptyRow()
		
		figureFile="./gaussClumpTest/RANDOM_{}.conf".format(i)
		
		#read this file
		
		file = open(figureFile, 'r') 
		for line in file: 
			name,value= str.strip(line).split("=")
			name=name.split('.')[1]
			
			paraRow[name]= float(value)



		paraRow[self.ID]=i
		
 
		
		return paraRow
		

		
	def getAllParameter(self,startI=0,endI=150):
	
		returnTB=self.getEmptyTB()
	
	
		for i in np.arange(startI,endI ):
		
			rowI=self.getOnRow(i)
		
			returnTB.add_row(rowI)
		
		
		#print returnTB
		
		return returnTB
		#print figureFile




