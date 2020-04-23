
#some table of gaia data base

import MySQLdb
from astropy.table import   Table
import numpy as np

#import psycopg2

from gaiaTB import GAIATB
 
doAG= GAIATB()

class GAIAAVFlag:
	name='public."MWISPGaiaNew"' #table name 
	
	
	#MWISPGaiaNew is the table with FLag
	
	
	GAIA_source_id ="source_id"


	GAIA_parallax ="parallax"

	GAIA_parallax_err ="parallax_error" 
	
	GAIA_l ="l"
	GAIA_b ="b"

	phot_g_mean_mag = "phot_g_mean_mag"

	GAIA_a_g_val ="a_g_val"

	GAIA_a_g_percentile_lower ="a_g_percentile_lower"
	GAIA_a_g_percentile_upper ="a_g_percentile_upper"

	agError ="agError" # relative ag error

	relative_error ="relative_error" ##relative parallax error

	av16="av16"
	av50="av50"
	av84="av84"

	dist16="dist16"
	dist50="dist50"
	dist84="dist84"

 

	sh_photoflag ="sh_photoflag"
	sh_parallaxflag  ="sh_parallaxflag"
	sh_gaiaflag  ="sh_gaiaflag" 
	sh_outflag  ="sh_outflag" 

	teff50 ="teff50"
	abp50  ="abp50"
	arp50  ="arp50" 
	
	




	colnames=[GAIA_source_id, GAIA_parallax,GAIA_parallax_err, GAIA_l,  GAIA_b, phot_g_mean_mag, GAIA_a_g_val,  GAIA_a_g_percentile_lower,\
	
	 GAIA_a_g_percentile_upper,    av16    , av50   ,   av84   , dist16   ,  dist50   , dist84, sh_photoflag, sh_parallaxflag,  sh_gaiaflag, sh_outflag, teff50, abp50, arp50]
 

	dataTypes= [float, float,float, float,  float, float, float, float,\
	 float,    float    , float   ,   float   , float   ,  float   , float, str, str, str, str, float,float,float  ]
	 
 

	def getDB(self):
		#read dbInfo.txt
 		
		#fileA = open('dbInfo.txt', 'r')
		#a = fileA.readlines()
		#db = MySQLdb.connect(str.strip(a[0]),str.strip(a[1]),str.strip(a[2]),str.strip(a[3]) )
		#db = MySQLdb.connect("127.0.0.1","root","shao1234","gaia")
		
 
		connection = psycopg2.connect(user = "postgres", password = "100425", host = "127.0.0.1", port = "5432", database = "gaia")


		#db = MySQLdb.connect(str.strip(a[0]),str.strip(a[1]),str.strip(a[2]),str.strip(a[3]) )
		return connection



	def getByLBRange(self,Lrange,Brange,lowerPara=0.2,paraError=0.2,upperPara=None,mimicAG=True):
		
		"""
		dedicated to find other regions that are overlaping with the current source

		# lowerPara is the samallest parallax, corresponding to the farthest distance
		
		# upperPara, is largest parallax, corresponding to the nearest distance
		
		
		
		By default, paraError less than 20%, and disances less than 5 kpc
		
		"""
		
		
		
		
		upperDis=1./lowerPara 
		
		lowerDis=0;
		
		if upperPara !=None:
			lowerDis= 1./upperPara 

 		
		
		sL= min(Lrange)
		
		eL= max(Lrange)
#
		sB=  min(Brange)
		eB=  max(Brange)


		db = self.getDB() #MySQLdb.connect("localhost","root","shao1234","gaia" )
		# prepare a cursor object using cursor() method
		cursor = db.cursor()
		
		
		
 
		#sqlCommand="select * from {} where  l > {} and l < {} and b > {} and b < {} and dist50>{} and dist50<{}  and parallax_error<parallax*{} and sh_outflag='00000' and sh_gaiaflag='000' ;".format(self.name,sL,eL,sB,eB,lowerDis,upperDis,paraError)

		#sqlCommand="select * from {} where  l > {} and l < {} and b > {} and b < {} and dist50>{} and dist50<{}  and parallax_error<parallax*{} and sh_outflag like '0%000' and sh_gaiaflag='000'   ;".format(self.name,sL,eL,sB,eB,lowerDis,upperDis,paraError)
		#sqlCommand="select * from {} where  l > {} and l < {} and b > {} and b < {} and dist50>{} and dist50<{}  and parallax_error<parallax*{} and  sh_outflag like '0%000' and sh_outflag != '01000'  and sh_gaiaflag='000'   ;".format(self.name,sL,eL,sB,eB,lowerDis,upperDis,paraError)
		#only good values
		sqlCommand="select * from {} where  l > {} and l < {} and b > {} and b < {} and dist50>{} and dist50<{}  and parallax_error<parallax*{} and    sh_outflag = '00000'  and sh_gaiaflag='000'   ;".format(self.name,sL,eL,sB,eB,lowerDis,upperDis,paraError)

		#sqlCommand="select * from {} where  l > {} and l < {} and b > {} and b < {} and dist50>{} and dist50<{}  and parallax_error<parallax*{} and (sh_outflag  = '00000' or sh_outflag  = '02000' )    and sh_gaiaflag='000'  and ( av84-av50)>0.02 ;".format(self.name,sL,eL,sB,eB,lowerDis,upperDis,paraError)

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
 
		
		returnTB=Table(rows= data , names=self.colnames,dtype= self.dataTypes )
		
		if mimicAG:
			return self.addAGCol(returnTB)
		
		return returnTB
		
	

	def addAGCol(self,TB):
		"""
		produce severl columns that mimic AG 
		"""
		newTB=TB.copy()
		
		disCol= ( TB[self.dist50])*1000 
		disCol.name= doAG.GAIA_distance
		
 
		disErrCol= ( TB[ self.dist84] -  TB[ self.dist16] )/2.*1000.
		disErrCol.name= doAG.GAIA_distanceError


		
		agErrorCol= ( TB[ self.av84] -  TB[ self.av16] )/2./TB[ self.av50]
		agErrorCol.name= doAG.agError
		
		
		
		
		agCol=  TB[ self.av50][:]
		#agCol.name= 
		newTB["gaiaAg"]=newTB["a_g_val"].copy() #keep the old ag value
		newTB[ doAG.GAIA_a_g_val] = agCol


		
		newTB.add_columns( [disCol,  disErrCol ,agErrorCol ])

		return newTB
		
		



	def converSqlToTB(self,sqlTB):
		
		"""
		"""
		
		pass
		
	def filterByLB(self,TB):
		
		"""
		"""
		
		pass