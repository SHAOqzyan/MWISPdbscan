from myPYTHON import *
from mwispDBSCAN import MWISPDBSCAN

doMWdbscan= MWISPDBSCAN()


rawCOFITS="Q1Sub.fits"

if 1:
    doMWdbscan.rawCOFITS =  rawCOFITS
    doMWdbscan.rmsFITS = None #you can provde rms fits if you hae one
    doMWdbscan.averageRMS = 0.5 # if you do not have an rm fits file, use an average rms

    doMWdbscan.setDBSCANParameters( cutoff_sigma=2,minPts=4,connectivity= 1 )
    doMWdbscan.processPath = './'

    doMWdbscan.computeDBSCAN()
    doMWdbscan.getCatFromLabelArray(doClean=True)
    doMWdbscan.produceCleanFITS()


