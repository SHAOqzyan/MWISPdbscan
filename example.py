
from onlyDBSCAN import allDBSCAN


doAllDBSCAN=allDBSCAN()


if 1: # Q1subfits test
    ############################################parameter setings
    #rawFITS="/home/qzyan/WORK/dataDisk/MWISP/G120/100_150_U.fits" #input raw coFITS
    #saveTag = "Q2test"  # your project code

    if 1:#simple test

        rawFITS="Q1Sub.fits" #input raw coFITS
        saveTag="Q1SubTest" #your project code


    if 0:#Q2Per
        rawFITS="/sci/yansun/100_150/100_150_U_per.fits" #input raw coFITS
        saveTag="Q2Per" #your project code



    cutoff=2 #in units of sigma
    minPts=4
    contype=1
    fitsRMS=0.5
    ########################

    doAllDBSCAN.rmsCO12= fitsRMS # K, set the rms



    doAllDBSCAN.setfastDBSCANrms(fitsRMS)
    doAllDBSCAN.pureDBSCAN(rawFITS, cutoff , MinPts=minPts, saveTag= saveTag , connectivity= contype , inputRMS=None, redo=True, keepFITSFile=True)


    if 1:
        outputFITS= doAllDBSCAN.fitsPath+saveTag+"dbscanS{}P{}Con{}.fits".format(cutoff,minPts,contype)
        outputTable= doAllDBSCAN.fitsPath+saveTag+"dbscanS{}P{}Con{}.fit".format(cutoff,minPts,contype)

    doAllDBSCAN.getMaskByLabel(  outputFITS  , outputTable, onlyClean=True )

    TB=doAllDBSCAN.getMaskByLabel( outputFITS ,  outputTable , onlySelect=True ) #only Select=True, meas we only want to select the catlaog
    TB.write(doAllDBSCAN.fitsPath+saveTag+"dbscanS{}P{}Con{}_clean.fit".format(cutoff,minPts,contype), overwrite=True)

