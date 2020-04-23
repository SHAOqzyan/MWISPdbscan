
#some table of gaia data base

import MySQLdb
from astropy.table import   Table
import numpy as np


class GAIATB:
	name="GSF_Gaia2_all" #table name
	

	
	GAIA_source_id ="source_id"


	GAIA_parallax ="parallax"

	GAIA_parallax_err ="parallax_err" 
	
	GAIA_l ="l"
	GAIA_b ="b"

	GAIA_distance ="distance"
	GAIA_distanceError ="distance_err"


	GAIA_a_g_val ="a_g_val"

	GAIA_a_g_percentile_lower ="a_g_percentile_lower"
	GAIA_a_g_percentile_upper ="a_g_percentile_upper"

	agError ="agError" # relative ag error


	relative_error ="relative_error" ##relative parallax error

 
	colnames=[GAIA_source_id, GAIA_parallax,GAIA_parallax_err, GAIA_l,  GAIA_b, GAIA_a_g_val,  GAIA_a_g_percentile_lower, GAIA_a_g_percentile_upper,relative_error,agError   ]
 

	dataTypes=[ float,float,float,float, float, float,float,float,float,float ] #all float
	
	GAIA_distance ="distance"
	GAIA_distanceError ="distance_err"
	

	def getDB(self):
		#read dbInfo.txt
 		
		#fileA = open('dbInfo.txt', 'r')
		#a = fileA.readlines()
		#db = MySQLdb.connect(str.strip(a[0]),str.strip(a[1]),str.strip(a[2]),str.strip(a[3]) )
		db = MySQLdb.connect("127.0.0.1","root","shao1234","gaia")

		#db = MySQLdb.connect(str.strip(a[0]),str.strip(a[1]),str.strip(a[2]),str.strip(a[3]) )
		return db



	def getByLBRange(self,Lrange,Brange,lowerPara=0.2,paraError=0.2,upperPara=None,calDis=False):
		
		"""
		dedicated to find other regions that are overlaping with the current source

		# lowerPara is the samallest parallax, corresponding to the farthest distance
		
		# upperPara, is largest parallax, corresponding to the nearest distance
		
		
		
		By default, paraError less than 20%, and disances less than 5 kpc
		
		"""
		 
 		
		
		sL= min(Lrange)
		
		eL= max(Lrange)
#
		sB=  min(Brange)
		eB=  max(Brange)


		db = self.getDB() #MySQLdb.connect("localhost","root","shao1234","gaia" )
		# prepare a cursor object using cursor() method
		cursor = db.cursor()
		
		
		
 		
		if upperPara is None:

			sqlCommand="select * from {} where  l > {} and l < {} and b > {} and b < {} and parallax>{} and parallax_err<parallax*{};".format(self.name,sL,eL,sB,eB,lowerPara,paraError)
		
		else:
			sqlCommand="select * from {} where  l > {} and l < {} and b > {} and b < {} and parallax>{} and parallax<{}  and parallax_err<parallax*{};".format(self.name,sL,eL,sB,eB,lowerPara,upperPara,paraError)

 
 
		# execute SQL query using execute() method.
		cursor.execute(sqlCommand)
		# Fetch a single row using fetchone() method.
 
		data = np.array(cursor.fetchall() )
		db.close()
 
		if len(data)==0  :
			#db.commit()
 
 
			print "No Gaia stars found."
			return None
 
		#t=  self.converSqlToTB(data)
		#construct a TB with data
		#db.commit()
		
		queryTB=Table(rows= data , names=self.colnames,dtype=('f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))

		if calDis:
			return self.addDisToTB(queryTB )
		
		return queryTB
		
		

	def addDisToTB(self,gaiaTB ):
		
		"""
		
		calculate distance, and add them to gaiaTB,
		
		usually only do once
		"""
		#copy file
		
		newTB=gaiaTB.copy()
 
		disCol= gaiaTB["parallax"]*0
		disCol.name= self.GAIA_distance
		
 
		disErrCol= gaiaTB["parallax"]*0
		disErrCol.name= self.GAIA_distanceError


		parallaxErrCol="parallax_err"

		if "parallax_error" in newTB.colnames: # for new gaia TB
			parallaxErrCol="parallax_error"


		
		newTB.add_columns( [disCol,  disErrCol  ])
		
		
		for eachRow in newTB: 
			para=eachRow["parallax"]
			paraError= eachRow[parallaxErrCol] #+0.1 #   #
			
			dA=1./np.random.normal(para,paraError,20000)*1000
			
			eachRow[self.GAIA_distance]=  round(np.mean(dA), 2) 
			
			eachRow[self.GAIA_distanceError]=   round(np.std(dA,ddof=1), 2) 
		

 		return newTB




	def converSqlToTB(self,sqlTB):
		
		"""
		"""
		
		pass
		
	def filterByLB(self,TB):
		
		"""
		"""
		
		pass