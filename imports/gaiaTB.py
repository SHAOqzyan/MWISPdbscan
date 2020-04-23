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

	agError ="agError"

	relative_error ="relative_error"

 
	colnames=[GAIA_source_id, GAIA_parallax,GAIA_parallax_err, GAIA_l,  GAIA_b, GAIA_a_g_val,  GAIA_a_g_percentile_lower, GAIA_a_g_percentile_upper,relative_error,agError   ]
 

	dataTypes=[ float,float,float,float, float, float,float,float,float,float ] #all float
	
	




	def searchByRadius(self,L,B,radius):
		"""
		return the radius serach with respect to L and B
		"""
		
		Lrange=[L-radius/2., L+radius/2.]
		Brange=[L-radius/2., L+radius/2.]
		
		sL= min(Lrange)
		eL= max(Lrange)
#
		sB=  min(Brange)
		eB=  max(Brange)
		
		db = MySQLdb.connect("localhost","root","shao1234","gaia" )
		# prepare a cursor object using cursor() method
		cursor = db.cursor()
		
		 
		sqlCommand="select * from {} where  l > {} and l < {} and b > {} and b < {};".format(self.name,sL,eL,sB,eB )

 

		# execute SQL query using execute() method.
		cursor.execute(sqlCommand)
		# Fetch a single row using fetchone() method.
 
		data = np.array(cursor.fetchall() )
		db.close()

		if len(data)==0  :
 
			print "No Gaia stars found."
			return None
 
		return Table(rows= data , names=self.colnames,dtype=('f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
		
 

	def getByLBRange(self,Lrange,Brange,lowerPara=0.2,paraError=0.2,upperPara=None):
		
		"""
		dedicated to find other regions that are overlaping with the current source
		
		
		By default, paraError less than 10%, and disances less than 5 kpc
		
		"""
		 
 		
		
		sL= min(Lrange)
		
		eL= max(Lrange)
#
		sB=  min(Brange)
		eB=  max(Brange)
		
  
		db = MySQLdb.connect("localhost","root","shao1234","gaia" )
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
		
		
		return Table(rows= data , names=self.colnames,dtype=('f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
		
 
  

	def converSqlToTB(self,sqlTB):
		
		"""
		"""
		
		pass