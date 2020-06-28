from myPYTHON import *
from mwispDBSCAN import MWISPDBSCAN

doMWdbscan= MWISPDBSCAN()



if 1: #an example
    rawCOFITS = "Q1Sub.fits"

    doMWdbscan.rawCOFITS =  rawCOFITS
    doMWdbscan.rmsFITS = None #you can provde rms fits if you hae one
    doMWdbscan.averageRMS = 0.5 # if you do not have an rm fits file, use an average rms

    doMWdbscan.setDBSCANParameters( cutoff_sigma=2,minPts=4,connectivity= 1 )
    doMWdbscan.processPath = './'

    doMWdbscan.computeDBSCAN()
    doMWdbscan.getCatFromLabelArray(doClean=True) #by cleaning, we remove noise clusters
    doMWdbscan.produceCleanFITS()


