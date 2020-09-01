
#This script is used reconstruct the script of DBSCAN.
#To make is more usefule to use. and can be imported septeratrely

import numpy as np
from astropy import units as u
from spectral_cube import SpectralCube
from astropy.table import Table,vstack
#need to download myPYTHON.py
from myPYTHON import *
import os
from progressbar import *

import sys
from skimage.morphology import erosion, dilation
from scipy.ndimage import label, generate_binary_structure,binary_erosion,binary_dilation
from sklearn.cluster import DBSCAN
import requests

doFITS=myFITS()

#print doFITS.weighted_avg_and_std

class  MWISPDBSCAN(object):


    processPath="./"#  by default save the product in the current path

    rawCOFITS=None
    rmsFITS=None
    labelFITSName=None
    catFITSName=None

    modelTBFile=None

    rmsData = None #must has the same shape with Ny,and Nz of the data file

    averageRMS = None

    cleanCatName = None
    cleanFITSName = None
    ################# DBSCAN parameters

    cutoff_sigma = 2
    minPts = 4
    connectivity = 1


    ###################### seting of catalog selection

    minVox=16
    minChannel=3
    hasBeam=1
    minPeakSigma=5
    ###



    def __init__(self):
        self.getModelTB()


    def setDBSCANParameters(self,cutoff_sigma=2,minPts=4,connectivity=1):

        self.cutoff_sigma=cutoff_sigma
        self.minPts=minPts
        self.connectivity=connectivity

    def setCatalogSelectionCriteria(self,minVox=16,minChannel=3,hasBeam=1,minPeakSigma=5):
        self.minVox=minVox
        self.minChannel=minChannel
        self.hasBeam=hasBeam
        self.minPeakSigma=minPeakSigma



    def sumEdgeByCon1(self, extendMask):  # 7 in total
        raw = extendMask[1:-1, 1:-1, 1:-1]

        leftShiftZ = extendMask[0:-2, 1:-1, 1:-1]
        rightShiftZ = extendMask[2:, 1:-1, 1:-1]

        leftShiftY = extendMask[1:-1, 0: -2, 1:-1]
        rightShiftY = extendMask[1:-1, 2:, 1:-1]

        leftShiftX = extendMask[1:-1, 1:-1, 0:-2]
        rightShiftX = extendMask[1:-1, 1:-1, 2:]

        sumAll = raw + leftShiftZ + rightShiftZ + leftShiftY + rightShiftY + leftShiftX + rightShiftX

        return sumAll

    def sumEdgeByCon2(self, extendMask):  # 27 in total
        sumAll = extendMask[1:-1, 1:-1, 1:-1] * 0
        Nz, Ny, Nx = sumAll.shape
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:

                    if np.sqrt(abs(i) + abs(j) + abs(k)) > 1.5:
                        continue

                    sumAll = sumAll + extendMask[1 + i:Nz + 1 + i, j + 1:Ny + 1 + j, k + 1: Nx + 1 + k]

        return sumAll

    def sumEdgeByCon3(self, extendMask):  # 27 in total
        raw = extendMask[1:-1, 1:-1, 1:-1]
        Nz, Ny, Nx = raw.shape
        sumAll = raw * 0
        for i in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for k in [-1, 0, 1]:
                    sumAll = sumAll + extendMask[1 + i:Nz + 1 + i, j + 1:Ny + 1 + j, k + 1: Nx + 1 + k]

        return sumAll


    ############

    def setOutputPath(self,targetPath):

        if os.path.isdir(targetPath):
            self.processPath=targetPath
        else:
            print "The target output path does not exist, please check your path: ", targetPath


    def getModelTB(self):
        """
        provide a TB file that can be used a model of catalog, new columns could be added to this table
        :return:
        """

        #targetFileName="https://raw.githubusercontent.com/SHAOqzyan/MWISPdbscan/master/minV3minP16_dendroCatTrunk.fit"

        targetFileName="minV3minP16_dendroCatTrunk.fit"

        if os.path.isfile(targetFileName):
            self.modelTBFile=targetFileName
        else:
            print "Downloading a modelTB file from the internet..."
            os.system("wget https://raw.githubusercontent.com/SHAOqzyan/MWISPdbscan/master/minV3minP16_dendroCatTrunk.fit")

        if os.path.isfile(targetFileName):
            self.modelTBFile=targetFileName

        else:
            print "Cannot download modelTB onthe internet, you need to manually download a table file named minV3minP16_dendroCatTrunk.fit (https://raw.githubusercontent.com/SHAOqzyan/MWISPdbscan/master/minV3minP16_dendroCatTrunk.fit)"




    def getDBSCANTag(self ):

        """
        return
        :param cutoff: #the sigma to which the overall datacube will be cut
        :param minpts:
        :param conType: #thre are
        :return:
        """

        saveName = "dbscanS{}P{}Con{}".format( self.cutoff_sigma,self.minPts,self.connectivity)

        return saveName



    def getUniformRMSData(self,averageRMS,Ny,Nx):



        return averageRMS+np.zeros(  (Ny,Nx) )

    def computeDBSCAN(  self  ):
        """
        There are two steps of computing DBSCAN, compute the DBSCAN label, 2 extract catlog, select post catalog and produce clean fits, would be doe seperatrely
        :param COdata:
        :param COHead:
        :param min_sigma: cutoff
        :param min_pix:
        :param connectivity:
        :param region:
        :param getMask:
        :param savePath:
        :param mimicDendro:
        :param rmsFITS:
        :param inputRMS:
        :return:
        """
        # pass

        if self.rawCOFITS is None:

            print "rawCOFITS need to be provided, stoping..."
            return

        COdata, COHead= doFITS.readFITS( self.rawCOFITS )


        if len(COdata.shape)==4:
            COdata=COdata[0]

        Nz, Ny, Nx = COdata.shape


        if self.averageRMS is None and self.rmsFITS is None:
            print "You need either use averaged rms, or use rmsFITS, which can be accurate to each spectrum. Please check you input"


        if self.rmsFITS == None: #use average rms, te rmsFITS has a higher priority
            rmsData= self.getUniformRMSData(self.averageRMS,Ny,Nx)

        else:
            print "Using rms fits...", self.rmsFITS
            rmsData, rmsHead = myFITS.readFITS(self.rmsFITS)
        self.rmsData=rmsData #record the rms Data

        rmsCOData = COdata / rmsData
        goodValues = rmsCOData >= self.cutoff_sigma
        extendMask = np.zeros([Nz + 2, Ny + 2, Nx + 2], dtype=int)


        extendMask[1:-1, 1:-1, 1:-1] = goodValues  # [COdata>=minValue]=1

        s = generate_binary_structure(3, self.connectivity )

        if self.connectivity == 1:
            coreArray = self.sumEdgeByCon1(extendMask)

        if self.connectivity == 2:
            coreArray = self.sumEdgeByCon2(extendMask)

        if self.connectivity == 3:
            coreArray = self.sumEdgeByCon3(extendMask)

        coreArray = coreArray >= self.minPts
        coreArray[~goodValues] = False  # nan could be, #remove falsely, there is a possibility that, a bad value may have lots of pixels around and clould be
        #coreArray = coreArray + 0

        labeled_core, num_features = label(coreArray,    structure=s)  # first label core, then expand, otherwise, the expanding would wrongly connected

        selectExpand = np.logical_and(labeled_core == 0, goodValues)
        # expand labeled_core
        # coreLabelCopy=labeled_core.copy()

        expandTry = dilation(labeled_core,   s)  # first try to expand, then only keep those region that are not occupied previously
        # it is possible  that a molecular cloud may have less pixN than 8, because of the closeness of two
        # only expanded, one time
        labeled_core[selectExpand] = expandTry[selectExpand]

        labeled_array = labeled_core



        print num_features, "features found!"

        saveLabelFITSName=self.getLabelFITSName()
        fits.writeto(  saveLabelFITSName , labeled_array, header=COHead, overwrite=True)
        self.labelFITSName=saveLabelFITSName

        return saveLabelFITSName



    def getLabelFITSName(self, doClean=False):
        """
        used to get the DBSCAN label fits name, when you do not want to rerun the DBSCAN process
        :return:
        """
        baseName  =  os.path.basename( self.rawCOFITS )
        saveTag = self.getDBSCANTag( )


        if doClean:

            saveLabelFITSName= os.path.join(self.processPath,baseName[0:-5]+saveTag+"_Clean.fits")

        else:
            saveLabelFITSName= os.path.join(self.processPath,baseName[0:-5]+saveTag+".fits")


        return saveLabelFITSName


    def getIndices(self, Z0, Y0, X0, values1D, choseID):

        cloudIndices = np.where(values1D == choseID)

        cX0 = X0[cloudIndices]
        cY0 = Y0[cloudIndices]
        cZ0 = Z0[cloudIndices]

        return tuple([cZ0, cY0, cX0])

    def getSaveCatName(self,doClean=False):


        processLabelFITSName=self.labelFITSName

        if processLabelFITSName is None:
            processLabelFITSName=self.getLabelFITSName(doClean=False)

        if not doClean: # _Clean.fit
            return processLabelFITSName[0:-1]

        else:
            return processLabelFITSName[0:-5]+"_Clean.fit"

    def rmsmap(self, outPUT=None, overwrite=True):
        """
        3d rms
        :param outPUT:
        :param overwrite:
        :return:
        """

        if self.rawCOFITS is None:
            print "rawCOFITS need to be provided, stoping..."
            return

        COdata, COHead = doFITS.readFITS(self.rawCOFITS)

        if len(COdata.shape) == 4:
            COdata = COdata[0]

        if outPUT is None:
            writeName = "rmsmap.fits"

        else:
            writeName = outPUT

        fileExist = os.path.isfile(writeName)

        if overwrite and fileExist:
            os.remove(writeName)

        Nz, Ny, Nx = COdata.shape


        rmsData = np.zeros_like(COdata, dtype=np.float32)

        for i in range(Nz):
            channelI = COdata[i]
            negativeValues = channelI[channelI < 0]

            sigma = np.std(negativeValues) / np.sqrt(1 - 2. / np.pi)

            print sigma
            rmsData[i, :, :] = sigma

        fits.writeto(writeName, rmsData, header=COHead)
        return fits.open(writeName)[0]

    def getCatFromLabelArray(self,doClean=True ):
        """
        Extract catalog from label fits, the minPix and rms is only used for saving
        #need a TB Model, if TB model does not exist, better downlooad from online
        This catalog has already been cleaned
        four criteria minPix"
        if you want to change the catlog selection criteria, please see setCatalogSelectionCariteria
        :param labelArray:
        :param head:
        :return:
        """

        if self.labelFITSName is None:

            print "The DBLABLE FITS is needed  "
            return

        dataCluster, headCluster = myFITS.readFITS(self.labelFITSName)


        #calculate the area of one pixel

        sizeL=  headCluster["CDELT1"]*60
        sizeB= headCluster["CDELT2"]*60

        sizeB=abs(sizeB) #in arcmin
        sizeL=abs(sizeL) #in arcmin

        pixelArea= sizeB*sizeL


        Nz, Ny, Nx = dataCluster.shape

        if self.rawCOFITS is None:

            print "A COFITS is needed  "
            return

        dataCO, headCO = myFITS.readFITS(self.rawCOFITS)


        if len(dataCO.shape)==4:
            dataCO=dataCO[0]


        if self.averageRMS is not None and self.rawCOFITS is not None:

            self.getUniformRMSData(self.averageRMS,Ny,Nx)


        if self.rmsData is None:

            print "A rmsData is needed to select peakSgima"
            return



        saveCatName=self.getSaveCatName(doClean=doClean)


        clusterTBOld = Table.read(self.modelTBFile)



        minV = np.nanmin(dataCluster[0]) #noise mask value



        wcsCloud = WCS(headCluster,naxis=3)

        clusterIndex1D = np.where(dataCluster > minV)
        clusterValue1D = dataCluster[clusterIndex1D]

        Z0, Y0, X0 = clusterIndex1D

        newTB = Table(clusterTBOld[0])
        newTB["sum"] = newTB["flux"]


        newTB["l_rms"] = newTB["v_rms"]*0
        newTB["b_rms"] = newTB["v_rms"]*0

        newTB["pixN"] = newTB["v_rms"]*0
        newTB["peak"] = newTB["v_rms"]*0

        newTB["peakL"] = newTB["v_rms"]*0
        newTB["peakB"] = newTB["v_rms"]*0
        newTB["peakV"] = newTB["v_rms"]*0
        newTB["area_accurate"] = newTB["area_exact"]*0  # the column of area_accurate need to cos(b) facor

        # newTB["Nchannel"]=newTB["v_rms"] #number of channels, whis is used to sellect

        newTB["allChannel"] = newTB["v_rms"]*0  # number channel involved
        newTB["has22"] = newTB["v_rms"]*0  # number channel involved


        #clear TB
        newTB["sum"].unit="K"
        del newTB["area_ellipse"]
        del newTB["flux"]
        del newTB["major_sigma"]
        del newTB["minor_sigma"]
        del newTB["radius"]
        del newTB["position_angle"]





        newTB["v_cen"].unit="km/s"
        newTB["x_cen"].unit="deg"
        newTB["y_cen"].unit="deg"
        newTB["l_rms"].unit="deg"
        newTB["b_rms"].unit="deg"
        newTB["peakB"].unit=""
        newTB["peakL"].unit=""
        newTB["peakV"].unit=""
        newTB["peak"].unit="K"

        newTB["allChannel"].unit=""
        newTB["has22"].unit=""






        zeroProjection = np.zeros((Ny, Nx))  # one zero channel, used to get the projection area and
        zeroProjectionExtend = np.zeros((Ny + 1, Nx + 1))

        idCol = "_idx"

        # count all clusters

        # ids,count=np.unique(dataCluster,return_counts=True )
        ids, count = np.unique(clusterValue1D, return_counts=True)
        GoodIDs = ids
        GoodCount = count
        print "Total number of trunks? ", len(ids)
        # print "Total number of Good Trunks? ",len(GoodIDs)

        # dataCO,headCO=doFITS.readFITS( CO12FITS )
        widgets = ['Recalculating cloud parameters: ', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(),
                   ' ', FileTransferSpeed()]  # see docs for other options

        catTB = newTB.copy()
        catTB.remove_row(0)

        # zeroP
        # remove any cloud with voxels less than minPix
        selectVox = GoodCount >= self.minVox  # #criteria 1

        if  doClean:
            GoodCount = GoodCount[selectVox]
            GoodIDs = GoodIDs[selectVox]
            print len(GoodIDs), "cluster has more than {} voxels.".format(self.minVox)

        pbar = ProgressBar(widgets=widgets, maxval=len(GoodIDs))
        pbar.start()

        ####

        for i in range(len(GoodIDs)):

            # i would be the newID
            newID = GoodIDs[i]

            pixN = GoodCount[i]

            newRow = newTB[0]

            newRow[idCol] = newID

            cloudIndex = self.getIndices(Z0, Y0, X0, clusterValue1D, newID)
            coValues = dataCO[cloudIndex]
            cloudV = cloudIndex[0]
            cloudB = cloudIndex[1]
            cloudL = cloudIndex[2]
            # step1, check peak
            # peak should be done first

            peakIndex = coValues.argmax()
            peakV = cloudV[peakIndex]
            peakB = cloudB[peakIndex]
            peakL = cloudL[peakIndex]
            peak=coValues[peakIndex]
            # #criteria 2

            peakSimga=0

            if len(self.rmsData.shape)==2:
                peakSimga =  self.minPeakSigma * self.rmsData[peakB,peakL]
            if len(self.rmsData.shape)==3:
                peakSimga =  self.minPeakSigma * self.rmsData[peakV,peakB,peakL]


            if peak <  peakSimga and doClean:  #  accurate to lines
                continue #do not consider the minimum peaks



            # #criteria 3, check consecutive channels

            ############### ###
            diffVs = np.unique(cloudV)

            if len(diffVs) < self.minChannel and doClean:  # reject all cloud that has channels less than 3 channels
                continue



            # get the exact peak position, which would be used to

            projectIndex = tuple([cloudB, cloudL])

            zeroProjection[projectIndex] = 1
            zeroProjectionExtend[0:-1, 0:-1] = zeroProjection
            sum22 = zeroProjectionExtend[0:-1, 0:-1] + zeroProjectionExtend[0:-1, 1:] + zeroProjectionExtend[1:,
                                                                                        0:-1] + zeroProjectionExtend[1:, 1:]

            # if any pixel>4:
            #

            if 4 in sum22:
                newRow["has22"] = 1

            else:  # Serious bug
                newRow["has22"] = 0



            if self.hasBeam and newRow["has22"]==0 and doClean: # #criteria 4
                zeroProjection[projectIndex] = 0
                zeroProjectionExtend[0:-1, 0:-1] = zeroProjection
                continue
            # newRow["has22"] = 0

            # calculate the accurate, this is only used for data at high Galactic latitude, and useless form the MWISP survey
            # to increase the speed, this part is commentted

            # indexB2D,indexL2D=np.where(zeroProjection==1 )
            # _,BS2D, LS2D = wcsCloud.wcs_pix2world(indexL2D,indexB2D,  0, 0)
            # area_accurate=np.sum( np.cos( np.deg2rad(BS2D) )    )*0.25
            # newRow["area_accurate"]= area_accurate

            sumCO = np.sum(coValues,dtype=np.float64) #float32 is not enough for large molecular clouds

            Vcen, Vrms = doFITS.weighted_avg_and_std(cloudV, coValues)
            Bcen, Brms = doFITS.weighted_avg_and_std(cloudB, coValues)
            Lcen, Lrms = doFITS.weighted_avg_and_std(cloudL, coValues)

            # calculate the exact area

            # LBcore = zip(cloudB, cloudL)
            # pixelsN= {}.fromkeys(LBcore).keys() #len( set(LBcore) )
            # area_exact=len(pixelsN)*0.25 #arc mins square
            area_exact = np.sum(zeroProjection) * pixelArea

            newRow["area_exact"] = area_exact

            # dataClusterNew[cloudIndex] =newID

            # save values
            newRow["pixN"] = pixN
            newRow["peak"] = peak

            newRow["peakV"] = peakV
            newRow["peakB"] = peakB
            newRow["peakL"] = peakL


            # newRow["peak2"]= peak2

            newRow["sum"] = sumCO

            newRow["x_cen"], newRow["y_cen"], newRow["v_cen"] = wcsCloud.wcs_pix2world(Lcen, Bcen, Vcen, 0)
            newRow["v_cen"] = newRow["v_cen"] / 1000.
            dv = headCluster["CDELT3"] / 1000.  # km/s

            dl = abs(headCluster["CDELT1"])  # deg

            newRow["v_rms"] = Vrms * dv
            newRow["l_rms"] = Lrms * dl
            newRow["b_rms"] = Brms * dl

            # _, Nchan=np.unique( cloudV, return_counts=True)

            # newRow["Nchannel"] =    np.max(P3Value)# if there is a three consecutive spectra in the cloud
            newRow["allChannel"] = len(diffVs)

            catTB.add_row(newRow)

            zeroProjection[projectIndex] = 0
            zeroProjectionExtend[0:-1, 0:-1] = zeroProjection



            pbar.update(i)

        pbar.finish()
        # save the clouds

        # fits.writeto(self.regionName+"NewCloud.fits", dataClusterNew,header=headCluster,overwrite=True   )
        catTB.write(saveCatName, overwrite=True)
        if doClean:
            self.cleanCatName=saveCatName
            #self.catFITSName=saveCatName

        else:
            self.catFITSName=saveCatName
        return saveCatName

    def getEmptyTB(self,tb):
        newTB=Table( tb[0] )

        newTB.remove_row(0)

        return newTB


    def cleanTB(self):

        """
        select clouds that satifies the selection criteria
        :return:
        """



        if self.catFITSName is None:
            print "no catalog file provided"

            return

        if self.rmsData is None:
            print "You need to provide the rms FITS data to select peak accordingly"


        saveCleanTBName=self.catFITSName[0:-4]+"_Clean.fit"

        rawTB=Table.read(self.catFITSName )

        #first voxel
        filterTB = rawTB[rawTB["pixN"] >= self.minVox ]


        filterTB = filterTB[filterTB["has22"] ==self.hasBeam ]

        #select by minCkannel
        filterTB = filterTB[filterTB["allChannel"] >= self.minChannel  ]

        #second peak

        peakB = map(int,filterTB["peakB"] )
        peakL = map(int,filterTB["peakL"] )
        peakCoordinate=tuple( [peakB, peakL] )


        rmsArray=self.rmsData[peakCoordinate]
        filterTB = filterTB[filterTB["peak"] >= self.minPeakSigma * rmsArray ]

        filterTB.write(saveCleanTBName, overwrite=True)
        self.cleanCatName=saveCleanTBName

        return  filterTB




    def produceCloudIntFITS(self,rawCOFITS,LabelCOFITS,tbFile, outputPath="./cloudIntPath",  useTB=False, pureInt=False, foreground=False):
        """
        the minimum Area of 0.015 square deg = 54 square arcmin, is the largest cloud we coud perform distance examination

        :param rawCOFITS:
        :param LabelCOFITS:
        :param tb:
        :param minimumArea:
        :return:
        """
        if not useTB:
            TB=Table.read(tbFile)
        else:
            TB=tbFile
        #TB=TB[TB["area_exact"]>=minimumArea] #only select large angular size molecular clouds

        TB.sort("area_exact")
        TB.reverse()

        dataCO,headCO=myFITS.readFITS(rawCOFITS)
        coSpec0, vaxis0 = doFITS.getSpectraByIndex( dataCO,headCO,0,0 )

        #Nz,Ny,Nx=dataCO.shape




        dataCluster,headCluster=myFITS.readFITS( LabelCOFITS  )

        ####################

        clusterIndex1D = np.where(dataCluster > 0)
        clusterValue1D = dataCluster[clusterIndex1D]
        Z0, Y0, X0 = clusterIndex1D

        velsolution = vaxis0[1]- vaxis0[0]

        ################
        projection0=np.zeros_like(dataCO[0])



        widgets = ['Integrating cloud fits:', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',   FileTransferSpeed()]  # see docs for other options

        pbar = ProgressBar(widgets=widgets, maxval=len(TB))
        pbar.start()
        indexRun = 0


        zeroCube = np.zeros_like(dataCluster, dtype=np.float32 )

        for eachDBRow in TB:
            indexRun = indexRun + 1
            pbar.update(indexRun)
            cloudID = eachDBRow["_idx"]

            cloudIndex = self.getIndices(Z0, Y0, X0, clusterValue1D, cloudID)
            iz,iy,ix=cloudIndex

            projectIndices=tuple( [iy,ix] )
            ############ determine v range
            centerV=eachDBRow["v_cen"]
            dv=eachDBRow["v_rms"]
            vRange= [centerV-dv*3, centerV+dv*3]
            startV0,endVindex= self.getVindexRange(vaxis0, vRange )

            cropCOCube=dataCO[startV0:endVindex+1]


            saveMaskFITSname= os.path.join( outputPath , "Cloud{}_mask.fits".format( cloudID  )  )
            savePureIntFITSname= os.path.join( outputPath , "Cloud{}_PureInt.fits".format( cloudID  ) )
            #produce pureint

            ###int fits
            sumCO=np.sum( cropCOCube, axis=0,dtype=float )*velsolution
            saveIntFITSname= os.path.join( outputPath , "Cloud{}_int.fits".format( cloudID  ) )
            fits.writeto(saveIntFITSname, sumCO, header=headCluster, overwrite=True)

            ##mask fits
            projection0[projectIndices] =1
            fits.writeto(saveMaskFITSname, projection0, header=headCluster, overwrite=True)

            ##foreground fits

            if foreground: #the way of generating foreground fits is different, the following code is used to produce foreground fits for Q2
                foreCOCube = dataCO[endVindex + 1:]

                sumForeground =   np.sum( foreCOCube , axis=0,dtype=float )*velsolution

                saveForeGroundFITSname= os.path.join( outputPath , "Cloud{}_fore.fits".format( cloudID  ) )
                fits.writeto(saveForeGroundFITSname, sumForeground, header=headCluster, overwrite=True)


            if pureInt: #usually do no use this
                startV0=min(iz)
                endV0= max(iz)

                zeroCube[cloudIndex] = dataCO[cloudIndex]
                cropCOCubePure=zeroCube[startV0:endV0+1]


                sumCOPure = np.sum(cropCOCubePure, axis=0, dtype=float) * velsolution
                fits.writeto(savePureIntFITSname, sumCOPure, header=headCluster, overwrite=True)


            projection0[projectIndices] =0
            zeroCube[cloudIndex] =  0

            #zeroCluster[cloudIndex] = cloudID  # remove this cluster and do not record this cluster

        pbar.finish()
        #fits.writeto(saveCleanFITSName,returnCluster,header=headCluster,overwrite=True)




    def getVindexRange(self,vaxis,vRange):
        """

        :param vaxis:
        :param vRange:
        :return:
        """

        minValue = min( vRange )
        maxValue = max( vRange )

        indexV0=doFITS.find_nearestIndex(vaxis,minValue)
        indexV1=doFITS.find_nearestIndex(vaxis,maxValue)

        return [indexV0, indexV1]



    def produceCleanFITS(self):
        """
        remove noise cluster according to the label fits and catFITS name
        :return:
        """
        if self.labelFITSName is None :

            print "You need to provde the label FITS  "

            return

        if self.cleanCatName is not None:

            targetTBFile= self.cleanCatName #table has alread been cleaned
        elif self.catFITSName is not None:
            #catalog has been calculated but not cleaned, need to clean table first, then do clean
            self.cleanTB()
            targetTBFile= self.cleanCatName #table has alread been cleaned

        elif self.catFITSName is None:
            "You need to provde a raw or a cleaned fits first to produce a clean fits"
            return

        saveCleanFITSName=targetTBFile+"s"

        TB=Table.read(targetTBFile)
        dataCluster,headCluster=doFITS.readFITS(self.labelFITSName)
        noiseLabel = np.min(dataCluster[0])
        clusterIndex1D = np.where(dataCluster > noiseLabel )
        clusterValue1D = dataCluster[clusterIndex1D]
        Z0, Y0, X0 = clusterIndex1D

        #

        widgets = ['Cleaning label fits:', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
                   FileTransferSpeed()]  # see docs for other options

        pbar = ProgressBar(widgets=widgets, maxval=len(TB))
        pbar.start()

        indexRun = 0

        returnCluster = np.zeros_like(dataCluster) + noiseLabel

        for eachDBRow in TB:
            indexRun = indexRun + 1
            pbar.update(indexRun)
            cloudID = eachDBRow["_idx"]

            cloudIndex = self.getIndices(Z0, Y0, X0, clusterValue1D, cloudID)
            returnCluster[cloudIndex] = cloudID  # remove this cluster and do not record this cluster

        pbar.finish()
        fits.writeto(saveCleanFITSName,returnCluster,header=headCluster,overwrite=True)

        self.cleanFITSName=saveCleanFITSName
        return saveCleanFITSName


    def produceMask(self,COFITS,labelFITS,outFITS=None):
        """
        produce a masked COFITS, with labelsFITS
        :param COFITS:
        :param labelFITS:
        :return:
        """


        dataCO,headCO =doFITS.readFITS(COFITS)

        dataLabel,headLabel=doFITS.readFITS(labelFITS)


        noiseLabel=np.min(dataLabel[0])


        if outFITS is None:
            outFITS="mask_"+COFITS


        #maskCO
        dataCO[dataLabel==noiseLabel] = np.nan

        fits.writeto(outFITS,dataCO,header=headCO,overwrite=True)

    def produceIndividualClouds(self, rawCOFITS, labelsFITS, cloudTBFile, savePath="./cloudSubCubes/" ,noiseMask=0 ):
        """

        #output all data cubes for each cloud

        :return:
        """

        #################

        # savePath=""

        if os.path.isdir(savePath):
            pass
        else:
            os.makedirs(savePath)


        cloudTB = Table.read(cloudTBFile)

        # cloudTB=self.removeWrongEdges(cloudTB)
        print len(cloudTB), " molecular clouds in total."

        dataCluster, headCluster = myFITS.readFITS(labelsFITS)
        dataCO, headCO = myFITS.readFITS( rawCOFITS )
        # print cloudTB

        minV = np.nanmin(dataCluster[0])
        wcsCloud = WCS(headCluster,naxis=3)
        clusterIndex1D = np.where(dataCluster > minV)
        clusterValue1D = dataCluster[clusterIndex1D]
        Z0, Y0, X0 = clusterIndex1D

        fitsZero = np.zeros_like(dataCluster,dtype=np.float32)
        fitsZero=fitsZero+noiseMask
         # print cloudTB.colnames
        for eachC in cloudTB:
            cloudID = eachC["_idx"]
            saveName = "cloud{}cube.fits".format(cloudID)

            cloudIndex = self.getIndices(Z0, Y0, X0, clusterValue1D, cloudID)
            fitsZero[cloudIndex] = dataCO[cloudIndex]

            cloudZ0, cloudY0, cloudX0 = cloudIndex

            minZ = np.min(cloudZ0)
            maxZ = np.max(cloudZ0)

            minY = np.min(cloudY0)
            maxY = np.max(cloudY0)

            minX = np.min(cloudX0)
            maxX = np.max(cloudX0)

            cropWCS = wcsCloud[minZ:maxZ + 1, minY:maxY + 1, minX:maxX + 1]

            cropData = fitsZero[minZ:maxZ + 1, minY:maxY + 1, minX:maxX + 1]

            fits.writeto(savePath + saveName, cropData, header=cropWCS.to_header(), overwrite=True)

            fitsZero[:] =  noiseMask
        print "Cloud fits writing done!"

    def getCloudIDByRow(self, eachC):


        if "_idx" in eachC.colnames:
            return eachC["_idx"]
        else:
            cloudName = eachC["sourceName"]
            id = cloudName.split("Cloud")[-1]
            return int(id)

    def getEquivalentLinewidth(self,labelFITSName, inputTBFile ,saveSpectral=False):
        """
        add a colnames of linewidth to the inputTB
        :param inputTB:
        :return:
        """

        #first read label and CO data
        dataCluster,headCluster=myFITS.readFITS( labelFITSName)
        dataCO,headCO= myFITS.readFITS( self.rawCOFITS   )


        filterTB= Table.read(inputTBFile)  # remove unrelated sources  #self.selectTB(rawTB)

        cubeCO = SpectralCube.read( labelFITSName )
        vAxis = cubeCO.spectral_axis


        vAxis = vAxis.value / 1000.  # convert to rms

        noiseV = np.nanmin(dataCluster[0])
        index1D = np.where(dataCluster > noiseV)
        values1D = dataCluster[index1D]

        Z0, Y0, X0 = index1D

        dataZero = np.zeros_like(dataCluster)

        widgets = ['Geting line equivalent width: ', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
                   FileTransferSpeed()]  # see docs for other options
        pbar = ProgressBar(widgets=widgets, maxval=len(filterTB))
        pbar.start()

        try:
            filterTB["lineWidth"] = filterTB["v_rms"]
        except:
            filterTB["lineWidth"] = filterTB["vlsr"]



        i = 0
        for eachR in filterTB:

            testID =  self.getCloudIDByRow( eachR )


            testIndices = self.getIndices(Z0, Y0, X0, values1D, testID)
            singleZ0, singleY0, singleX0 = testIndices

            dataZero[testIndices] = dataCO[testIndices]

            # cropThe cloudRange
            minY = np.min(singleY0)
            maxY = np.max(singleY0)
            ###########
            minX = np.min(singleX0)
            maxX = np.max(singleX0)

            ###########
            minZ = np.min(singleZ0)
            maxZ = np.max(singleZ0)

            #########

            cloudCropSpectra = dataZero[:, minY:maxY + 1, minX:maxX + 1]

            cloudCropCube = dataZero[minZ:maxZ + 1, minY:maxY + 1, minX:maxX + 1]

            averageSpectraCrop = np.nansum(cloudCropSpectra, axis=(1, 2))

            intCloud = np.nansum(cloudCropCube, axis=0)

            # count the number spectra

            totalSpectral = len(intCloud[intCloud > 0])

            meanSpectral = averageSpectraCrop / 1. / totalSpectral

            if saveSpectral:
                savefileName = saveSpectralPath + "{}_{}Spectral".format(regionName, testID)

                np.save(savefileName, [vAxis, meanSpectral])


            spectraPeak = np.max(meanSpectral)

            area = (vAxis[1] - vAxis[0]) * np.sum(meanSpectral)

            eqLineWidth = area / spectraPeak
            dataZero[testIndices] = 0

            eachR["lineWidth"] = eqLineWidth

            #lineWdith.append(eqLineWidth)

            i = i + 1

            pbar.update(i)

        pbar.finish()

        saveTBAs=  inputTBFile[0:-5]+"_LW.fit"
        filterTB.write(saveTBAs, overwrite=True)

    def selectTB(self, TB, rmsFITS=None,meanRMS=None,minVox=16,minChannel=3,hasBeam=1,minPeakSigma=5):

        """

        :param rmsFITS:
        :param meanRMS:
        :param minVox:
        :param minChannel:
        :param hasBeam:
        :param minPeakSigma:
        :return:
        """

        if rmsFITS is None and meanRMS is None:
            return


        if rmsFITS is not None:
            if rmsFITS is not None:
                rmsData, rmsHead = doFITS.readFITS(rmsFITS)

            else:
                rmsData =   meanrms



    def ZZZ(self):
        pass