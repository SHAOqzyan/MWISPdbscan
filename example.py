from mpPYTHON import *
from mwispDBSCAN import MWISPDBSCAN

doMWdbscan= MWISPDBSCAN()


rawCOFITS="Q1Sub.fits"

if 1:
doMWdbscan.rawCOFITS =  coFITS
    doMWdbscan.rmsFITS = rawCOFITS

    doMWdbscan.setDBSCANParameters( cutoff_sigma=2,minPts=4,connectivity= 1 )
    doMWdbscan.processPath = './'

    doMWdbscan.computeDBSCAN()
    doMWdbscan.getCatFromLabelArray(doClean=True)
    doMWdbscan.produceCleanFITS()


