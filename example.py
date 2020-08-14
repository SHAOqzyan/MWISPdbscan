from myPYTHON import *
from mwispDBSCAN import MWISPDBSCAN

doMWdbscan= MWISPDBSCAN()


def pipeLine(rawCOFITS,rmsFITS=None,averageRMS=0.5,processPath="./"):

    #rawCOFITS = "../fcrao_rep_m.fits"
    rawCOFITS =  rawCOFITS #"Q1Sub.fits"

    doMWdbscan.rawCOFITS =  rawCOFITS
    doMWdbscan.rmsFITS = rmsFITS #you can provde rms fits if you hae one
    doMWdbscan.averageRMS =averageRMS # if you do not have an rm fits file, use an average rms

    doMWdbscan.setDBSCANParameters( cutoff_sigma=2,minPts=4,connectivity= 1 )
    doMWdbscan.processPath =  processPath

    doMWdbscan.setCatalogSelectionCriteria( minVox=16,minChannel=3,hasBeam=1,minPeakSigma=5)

    #if rmsFITS is not None:
    #dataRMs,headRMs=myFITS.readFITS(rmsFITS)
    doMWdbscan.computeDBSCAN()
    doMWdbscan.labelFITSName= doMWdbscan.getLabelFITSName()
    doMWdbscan.getCatFromLabelArray(doClean=True) #by cleaning, we remove noise clusters
    doMWdbscan.produceCleanFITS()


    #doMWdbscan.produceMask(doMWdbscan.rawCOFITS, doMWdbscan.cleanFITSName ) # labelFITSName may be none if the computeDBSCAN() dot no process, you need to

    #doMWdbscan.produceIndividualClouds( doMWdbscan.rawCOFITS, doMWdbscan.cleanFITSName ,doMWdbscan.cleanCatName  )

if 1:
    rawCOFITS =  "Q1Sub.fits" #"Q1Sub.fits"

    doMWdbscan.rawCOFITS =  rawCOFITS

    doMWdbscan.getEquivalentLinewidth( "Q1SubdbscanS2P4Con1_Clean.fits", "Q1SubdbscanS2P4Con1_Clean.fit")
if 0: #an example

    pipeLine("Q1Sub.fits", averageRMS=0.5)



if 0:
    #doMWdbscan.produceIndividualClouds( "Q1Sub.fits", "Q1SubdbscanS2P4Con1_Clean.fits", "Q1SubdbscanS2P4Con1_Clean.fit"  )

    doMWdbscan.produceCloudIntFITS( "Q1Sub.fits", "Q1SubdbscanS2P4Con1_Clean.fits", "Q1SubdbscanS2P4Con1_Clean.fit",outputPath="./cloudIntPath", foreground=True  )


if 0: #an example for only produceIndividualCLouds

    rawCOFITS =  "Q1Sub.fits" #"Q1Sub.fits"

    doMWdbscan.rawCOFITS =  rawCOFITS
    doMWdbscan.rmsFITS = None #you can provde rms fits if you hae one
    doMWdbscan.averageRMS =0.5 # if you do not have an rm fits file, use an average rms

    doMWdbscan.setDBSCANParameters( cutoff_sigma=2,minPts=4,connectivity= 1 )
    doMWdbscan.processPath =  "./"

    #cloudSubCubes
    doMWdbscan.produceIndividualClouds( doMWdbscan.rawCOFITS, doMWdbscan.getLabelFITSName(doClean=True) ,doMWdbscan.getSaveCatName(doClean=True)  )
