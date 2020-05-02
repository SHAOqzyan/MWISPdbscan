
import numpy as np

from astropy.io import fits
#this script is used to exame all possible parameter spaces
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

from fastDBSCAN import myDBSCAN
import sys
from myPYTHON import *

import glob
from astropy.table import Table,vstack
import seaborn as sns
from progressbar import *

import matplotlib as mpl
import matplotlib.pyplot as plt
import os
from spectral_cube import SpectralCube
from  myGAIA import GAIADIS

doFITS=myFITS()

gaiaDis=GAIADIS()

doDBSCAN=myDBSCAN()



class allDBSCAN:

    #

    rootPath="./"

    emptyFITS = "G2650LocalEmpty.fits"
    rmsCO12=0.5 # K

    #tmpPath= rootPath+"tmpOUT/"  #"./DBSCANTest/"
    tmpPath="./tmpFiles/"
    testAllPath =tmpPath  # rootPath+"testAll/"

    con1PtsAll= np.arange(3,8,1)
    con2PtsAll= np.arange(3,20,1)
    con3PtsAll=  np.arange(3,28,1)

    G2650Local = "G2650Local30.fits"

    #con1G2650= np.arange( 3, 8, 1)
    #con2G2650= np.arange( 5, 20,1) #start from 5
    #con3G2650=  np.arange( 6, 28, 1) #start from 6

    #########
    con1G2650= np.arange( 4 , 8, 1)
    con2G2650= np.arange( 8 , 20, 1 ) #start from 5
    con3G2650= np.arange( 11 , 28, 1) #start from 6
    #########





    conTypeG2650={1:con1G2650,2:con2G2650,3:con3G2650}

    cutoffList= [2, 2.5,3.0,3.5,4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
    cutoffInKelvin=np.asarray( cutoffList )*rmsCO12

    #cutoffList=   np.arange( 2, 7.5, 0.5)

    saveTag=None



    #color for minPts

    NUM_COLORS = 25
    clrs = sns.color_palette('husl', n_colors=NUM_COLORS)  # a list of RGB tuples


    # norm = matplotlib.colors.BoundaryNorm(np.arange(0,30,0.1), cmap.N)
    normV = mpl.colors.Normalize(vmin=3, vmax=27)
    mMinPts = plt.cm.ScalarMappable(norm=normV, cmap=plt.cm.jet)

    con1Str = "Connectivity 1"
    con2Str = "Connectivity 2"
    con3Str = "Connectivity 3"


    #fitsPath="/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/"

    fitsPath="./DBSCANresults/"



    saveTag = "G2650Local"


    selectFormalCode = "selectFormal"
    selectVoxAndPeakCode = "voxAndPeak"
    selectAreaAndChannelCode = "selectAreaAndChannel"

    dendroSigmaList = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7]

    SCIMESStr="SCIMES"

    averageFluxAlpha=  1.79 #+./0.03
    averageAreaAlpha= 1.97

    def __init__(self):
        pass

        if not os.path.isdir(self.fitsPath):
            os.system("mkdir "+self.fitsPath)

        if not os.path.isdir(self.tmpPath):
            os.system("mkdir "+self.tmpPath)


    def getMiddleMinPts(self, minPtList ):
        """

        :param minPtList:
        :return:
        """

        minV=  np.min(minPtList)
        maxV= np.max(minPtList)

        midV= (minV+maxV)/1./2

        return int(round( midV) )


    def getColorByMinPts(self,minPts):

        return self.clrs[minPts-3]

    def setfastDBSCANrms(self,inputRMS):

        doDBSCAN.rms=inputRMS

    def  DBSCANTest(self, rawCOFITS, saveTag ,connecttivity=2,minValue=2,minPix=8,MinPts=8,minDelta=0,minAreaPix=0):
        """

        :param rawCOFITS:
        :param saveTag:
        :return:
        """

        # co12FITSInRMS = "/home/qzyan/WORK/projects/NewUrsaMajorPaper/UMMCCO12InRmsUnit.fits"
        doDBSCAN.DBSCANAllInOne(rawCOFITS, saveTag, rmsFITS=None, minDelta=minDelta, minValue=minValue, minPix=minPix, minAreaPix=minAreaPix,  MinPts=MinPts, sigma=self.rmsCO12,  connectivity=connecttivity, outPath= self.outPath )


    def pureDBSCAN(self, FITSfile, minValue, MinPts=3, saveTag="",  connectivity=2,inputRMS=None ,redo=True ,keepFITSFile=True, onlyGetFITS=False ):
        """
        #
        :param FITSfile:
        :param minValue:
        :param MinPts:
        :param saveTag:
        :param connectivity:
        :param inputRMS:
        :return:
        """

        dataCO,headCO = myFITS.readFITS( FITSfile )


        if not redo:#check if target fit fit file exisit
            searchStr=self.testAllPath+saveTag+"*S{}P{}Con{}.fit".format(minValue,MinPts,connectivity)
            isDone=glob.glob(searchStr)
            if len(isDone) ==1:
                print "*S{}P{}Con{}.fit".format(minValue,MinPts,connectivity),"done, skipping..."
                return
        # step 1 compute DBSCAN
        print "Step 1: computing DBSCAN......"

        dbscanLabelFITS = doDBSCAN.computeDBSCAN(dataCO, headCO, min_sigma=minValue, min_pix=MinPts,   connectivity=connectivity, region=self.fitsPath + saveTag, rmsFITS=None,     inputRMS=inputRMS)

        if onlyGetFITS:
            return dbscanLabelFITS


        saveMarker= dbscanLabelFITS[:-5]

        print "Step 2: computing DBSCAN table......"

        rawDBSCANTBFile = doDBSCAN.getCatFromLabelArray(FITSfile, dbscanLabelFITS, doDBSCAN.TBModel,  saveMarker= saveMarker )


        if not keepFITSFile:
            os.remove(dbscanLabelFITS)









    def selectTBByCode(self,cutOff,rawTB,selectionCode):
        """

        :param rawTB:
        :param selectCode:
        :return:
        """

        if selectionCode == self.selectFormalCode:  #
            filterTB = self.selectTBFormal(rawTB, cutOff=cutOff )

        if selectionCode == self.selectVoxAndPeakCode:  #
            filterTB = self.selectTBByPeak(rawTB, cutOff=cutOff )

        if selectionCode == self.selectAreaAndChannelCode:  #
            filterTB = self.selectTBAreaChannel( rawTB )

        return filterTB



    def getTBByCutOffList(self,cutOffList, minPts, conType, calCode="G2650Local" , getCleanLW=False ,selectionCode="selectFormal"  ):
        """

        :param cutOffList:
        :param minPts:
        :param conType:
        :param calCode:
        :return:
        """
        #
        returnTBList= []
        searchPath="/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/"
        for eachCutoff in cutOffList:

            searchStr= searchPath+"{}dbscanS{}P{}Con{}.fit".format( calCode,  eachCutoff  ,  minPts ,    conType  )
            if getCleanLW:
                searchStr = searchPath + "{}dbscanS{}P{}Con{}*CleanWithLW.fit".format(calCode, eachCutoff, minPts, conType)

            a= sorted( glob.glob(  searchStr  ) ) #sort this in names

            if len(a)==0:
                returnTBList.append(None)

            else:

                tmpTB = Table.read( a[0] )

                #tmpTB=self.selectTB(tmpTB,areaPix=areaPix,conChannel=conChannel) #by pixel
                #tmpTB=self.selectTBByPeak(tmpTB,minDelta=eachCutoff+3 ) #by peak
                tmpTB=self.selectTBByCode(eachCutoff,tmpTB,selectionCode)

                returnTBList.append( tmpTB )


        #remove wrong edges
        returnTBList=doDBSCAN.removeAllEdges(returnTBList)

        return returnTBList



    def getTBNlist(self,TBList):

        Nlist=[]

        for eachTB in TBList:
            if eachTB==None:
                Nlist.append( 0 )
            else:

                Nlist.append(len(eachTB))

        return Nlist



    def selectTBListByPeak(self, TBList, area=4, pixN=16, minDelta= 3 ,cutOff=2 ,minChannel=3  ):
        """
        #the area should be replace with a pattern of 2*2 tile in the projection image
        :return:
        """
        #
        returnTB = []
        for eachTB in TBList:

            filterTB=self.selectTBByPeak(eachTB, area=area, pixN=pixN,minDelta=minDelta, cutOff=cutOff, minChannel=minChannel  )
            returnTB.append(filterTB)

        return returnTB


    def selectTBByPeak(self,TB,  pixN=16, minDelta= 3 ,cutOff=2   ):
        """

        :param TB:
        :param areaPix:
        :param conChannel:
        :return:
        """
        #


        filterTB = TB[TB["pixN"] >= pixN ]

        filterTB = filterTB[filterTB["peak"] >= (minDelta+cutOff)*self.rmsCO12]

        #select by minCHannel

        #filterTB = filterTB[filterTB["allChannel"] >= minChannel  ]

        tmpList=doDBSCAN.removeAllEdges([filterTB])

        return tmpList[0]

    ##########################################################################################
    def selectTBFormal(self,TB,cutOff=2 ,  pixN=16, minDelta= 3 ,hasBeam=True,minChannel=3 ,verbose=False,removeEdge=True ):
        """
        # This is the most strict critera to select, first by pixN, second by miNDelta, which only affect peak, the peak is calculated by cuOff+minDelta,
        # minChannel and has Beam would also be done
        :param TB:
        :param areaPix:
        :param conChannel:
        :return:
        """

        #first, check cuOff, to prevent wrong cutOff

        #first voxel
        filterTB = TB[TB["pixN"] >= pixN ]

        #second peak
        filterTB = filterTB[filterTB["peak"] >= (minDelta+cutOff)*self.rmsCO12]

        #third by beam,
        if hasBeam: # this is
            filterTB = filterTB[filterTB["has22"] >=  0.5  ]

        #select by minCHannel
        filterTB = filterTB[filterTB["allChannel"] >= minChannel  ]


        #reamove edged
        if removeEdge:
            tmpList=doDBSCAN.removeAllEdges([filterTB])

        else:
            return filterTB
            pass

        return tmpList[0]




    def selectTBOld(self,TB,  pixN=16, minDelta= 3 ,cutOff=2    ):
        """

        :param TB:
        :param areaPix:
        :param conChannel:
        :return:
        """
        #




        filterTB = TB[TB["pixN"] >= pixN ]

        filterTB = filterTB[filterTB["peak"] >= (minDelta+cutOff)*self.rmsCO12]

        #select by minCHannel

        tmpList=doDBSCAN.removeAllEdges([filterTB])

        return tmpList[0]


    def selectTBAreaChannel(self,TB,  channelN=3   ):
        """
        ############## ####
        :param TB:
        :param areaPix:
        :param conChannel:
        :return:
        """
        #

        #
        filterTB = TB[TB["has22"] >= 0.5 ]

        filterTB = filterTB[filterTB["allChannel"] >= channelN ]

        #select by minCHannel

        tmpList=doDBSCAN.removeAllEdges([filterTB])

        return tmpList[0]


    def selectTB(self,TB, areaPix=9,conChannel=True   ):
        """

        :param TB:
        :param areaPix:
        :param conChannel:
        :return:
        """
        #

        filterTB = TB[TB["area_exact"] >= areaPix * 0.25]

        # print  len(  filterTB )
        if conChannel:  # at least contain a spectra consecutively
            filterTB = filterTB[filterTB["Nchannel"] >= 0.5]


        #remove edges


        tmpList=doDBSCAN.removeAllEdges([filterTB])

        filterTB = tmpList[0]

        return filterTB

    def getTBListByCon(self,filePath,conType,saveTag, selectionCode = "Formal"  ,   clean=False):
        """
        #selectionCode: Formal (4 criteria), voxAndPeak, areaAndChannel


        :param filePath:
        :param conType:
        :param areaPix is the minimum
        :return:
        """

        doMinList = None

        if conType==1:
            doMinList = self.con1PtsAll

        if conType==2:
            doMinList = self.con2PtsAll

        if conType==3:
            doMinList = self.con3PtsAll

        TBNameList = []
        #getTBList

        for eachMinPt in  doMinList:
            searchStr = filePath +  "{}dbscanS2P{}Con{}.fit".format(saveTag , eachMinPt  ,  conType   )
            a= sorted( glob.glob(  searchStr  ) ) #sort this in names

            TBNameList.append(a[0])

        tbList=[]
        Nlist = []

        for eachTB in TBNameList:
            getTB=   Table.read(eachTB)
            #remove by  area
            cutOff=self.getCutOFF(eachTB)

            filterTB=self.selectTBByCode( cutOff,  getTB,  selectionCode)


            if clean: #remove sources, that are rejected
                #find the corresponding fits file

                searchFITS =  eachTB[0:-4]+".fits"
                saveFITS=eachTB[0:-4]+"_Clean.fits"

                dataCluster,headCluster=myFITS.readFITS( searchFITS )

                cleanCluster = doDBSCAN.cleanLabelFITS(dataCluster,filterTB)
                fits.writeto( saveFITS,   cleanCluster,   headCluster,   overwrite= True  )



            tbList.append(   filterTB     )
            Nlist.append(  len(filterTB)  )
            #produce

        return tbList,Nlist



    def drawNegativeTest(self,areaAndChannel=True):

        drawPath = "/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/"

        areaPix  =  4
        doClean=False
        region1Tag= "negativeTest"
        selectCode1=self.selectAreaAndChannelCode
        TBListCon1_Code1 , NListCon1_Code1 = self.getTBListByCon(drawPath,1, region1Tag , selectionCode= selectCode1   , clean=doClean )
        TBListCon2_Code1 , NListCon2_Code1 = self.getTBListByCon(drawPath,2, region1Tag , selectionCode= selectCode1, clean=doClean )
        TBListCon3_Code1 , NListCon3_Code1 = self.getTBListByCon(drawPath,3, region1Tag , selectionCode= selectCode1 , clean=doClean )


        selectCode2=self.selectVoxAndPeakCode
        TBListCon1_Code2 , NListCon1_Code2 = self.getTBListByCon(drawPath,1, region1Tag , selectionCode= selectCode2 , clean=doClean )
        TBListCon2_Code2 , NListCon2_Code2 = self.getTBListByCon(drawPath,2, region1Tag , selectionCode= selectCode2 , clean=doClean )
        TBListCon3_Code2 , NListCon3_Code2 = self.getTBListByCon(drawPath,3, region1Tag , selectionCode= selectCode2 , clean=doClean )


        selectCode3=self.selectFormalCode
        TBListCon1_Code3 , NListCon1_Code3 = self.getTBListByCon(drawPath,1, region1Tag , selectionCode= selectCode3   , clean=doClean )
        TBListCon2_Code3 , NListCon2_Code3 = self.getTBListByCon(drawPath,2, region1Tag , selectionCode= selectCode3, clean=doClean )
        TBListCon3_Code3 , NListCon3_Code3 = self.getTBListByCon(drawPath,3, region1Tag , selectionCode= selectCode3 , clean=doClean )



        fig = plt.figure(figsize=(19, 6))

        ax1 = fig.add_subplot(1, 3, 2)
        ax2 = fig.add_subplot(1, 3, 1 , sharex=ax1, sharey=ax1)
        ax3 = fig.add_subplot(1, 3, 3,  sharex=ax1, sharey=ax1)

        # fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
        rc('text', usetex=True)

        overAllFontSize = 16

        rc('font', **{'family': 'sans-serif', 'size': overAllFontSize, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
            r'\LARGE',  # force siunitx to use the fonts

        ]

        #plt.rcParams.update({'font.size': 22})

        if 1:
            markersize=5

            ax1.plot(self.con1PtsAll ,  NListCon1_Code1 ,  'D--',color='b' ,lw=1.5, markersize=markersize, label="Connectivity 1" )
            ax1.plot(self.con2PtsAll ,  NListCon2_Code1,  '^--',color='g' ,lw=1.5, markersize=markersize , label="Connectivity 2" )
            ax1.plot(self.con3PtsAll ,  NListCon3_Code1 ,  's--',color='r' ,lw=1.5, markersize=markersize,  markerfacecolor='none' ,label="Connectivity 3" )


            ax2.plot(self.con1PtsAll ,  NListCon1_Code2 ,  'D--',color='b' ,lw=1.5, markersize=markersize )
            ax2.plot(self.con2PtsAll ,  NListCon2_Code2 ,  '^--',color='g' ,lw=1.5, markersize=markersize )
            ax2.plot(self.con3PtsAll ,  NListCon3_Code2 ,  's--',color='r' ,lw=1.5, markersize=markersize,  markerfacecolor='none' )


            ax3.plot(self.con1PtsAll ,  NListCon1_Code3 ,  'D--',color='b' ,lw=1.5, markersize=markersize , label="Connectivity 1" )
            ax3.plot(self.con2PtsAll ,  NListCon2_Code3 ,  '^--',color='g' ,lw=1.5, markersize=markersize , label="Connectivity 2" )
            ax3.plot(self.con3PtsAll ,  NListCon3_Code3 ,  's--',color='r' ,lw=1.5, markersize=markersize,  markerfacecolor='none',label="Connectivity 3" )





            #ax1.plot(self.con2PtsAll  , NListCon2, '^--', color='g' ,lw=1.5)

            #ax1.plot(self.con3PtsAll ,  NListCon3, 's--',color='r' , lw=1.5)

            maxY=10

            ax1.set_xlim(2, 27)
            ax1.set_ylim(-2 , maxY )


            # ax.set_ylim(-1,12)

        ax2.set_ylabel("Number of fake clusters", fontsize=overAllFontSize)
        ax1.set_xlabel("MinPts" , fontsize=overAllFontSize)

        at = AnchoredText('By beam size and velocity channels', loc=4, frameon=False)
        ax1.add_artist(at)



        #ax2.set_ylabel("Number of falsely detected clusters" , fontsize=overAllFontSize)
        ax2.set_xlabel("MinPts" , fontsize=overAllFontSize)
        at = AnchoredText("By voxels and peak", loc=4, frameon=False)
        ax2.add_artist(at)

        ax3.set_xlabel("MinPts" , fontsize=overAllFontSize)
        at = AnchoredText('By voxels, peak, beam size, and velocity channels', loc=4, frameon=False)
        ax3.add_artist(at)


        ax3.legend( loc=1 )
        ax1.tick_params(axis='both', which='major', labelsize=overAllFontSize)
        ax1.tick_params(axis='both', which='minor', labelsize=overAllFontSize)
        ax2.tick_params(axis='both', which='major', labelsize=overAllFontSize)
        ax2.tick_params(axis='both', which='minor', labelsize=overAllFontSize)

        ax3.tick_params(axis='both', which='major', labelsize=overAllFontSize)
        ax3.tick_params(axis='both', which='minor', labelsize=overAllFontSize)


        ax1.set_ylabel("Number of fake clusters", fontsize=overAllFontSize)
        ax3.set_ylabel("Number of fake clusters", fontsize=overAllFontSize)



        fig.tight_layout()


        plt.savefig("dbscanNegativeTest.png", bbox_inches='tight',dpi=600)
        plt.savefig("dbscanNegativeTest.pdf", bbox_inches='tight' )


    def drawEmptyTest(self,areaAndChannel=True):

        drawPath = "/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/"

        areaPix  =  4
        doClean=False
        region1Tag= "testEmpty1"
        TBListCon1_region1, NListCon1_region1 = self.getTBListByCon(drawPath,1, region1Tag, areaAndChannel=areaAndChannel, clean=doClean )
        TBListCon2_region1, NListCon2_region1 = self.getTBListByCon(drawPath,2,region1Tag , areaAndChannel=areaAndChannel,  clean=doClean )
        TBListCon3_region1, NListCon3_region1 = self.getTBListByCon(drawPath,3, region1Tag ,  areaAndChannel=areaAndChannel,  clean=doClean )


        region2Tag= "testEmpty2"
        TBListCon1_region2, NListCon1_region2 = self.getTBListByCon(drawPath,1, region2Tag,  areaAndChannel=areaAndChannel,  clean=doClean )
        TBListCon2_region2, NListCon2_region2 = self.getTBListByCon(drawPath,2, region2Tag,   areaAndChannel=areaAndChannel,  clean=doClean )
        TBListCon3_region2, NListCon3_region2 = self.getTBListByCon(drawPath,3, region2Tag,  areaAndChannel=areaAndChannel,  clean=doClean )




        fig = plt.figure(figsize=(12, 6))
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)

        # fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
        rc('text', usetex=True)

        overAllFontSize = 12

        rc('font', **{'family': 'sans-serif', 'size': overAllFontSize, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
            r'\LARGE',  # force siunitx to use the fonts

        ]

        #plt.rcParams.update({'font.size': 22})

        if 1:
            markersize=5

            ax1.plot(self.con1PtsAll ,  NListCon1_region1,  'D--',color='b' ,lw=1.5, markersize=markersize,label="Connectivity 1" )
            ax2.plot(self.con1PtsAll ,  NListCon1_region2 ,  'D--',color='b' ,lw=1.5, markersize=markersize )

            ax1.plot(self.con2PtsAll ,  NListCon2_region1,  '^--',color='g' ,lw=1.5, markersize=markersize ,label="Connectivity 2" )
            ax2.plot(self.con2PtsAll ,  NListCon2_region2 ,  '^--',color='g' ,lw=1.5, markersize=markersize )


            ax1.plot(self.con3PtsAll ,  NListCon3_region1,  's--',color='r' ,lw=1.5, markersize=markersize,  markerfacecolor='none' ,label="Connectivity 3" )
            ax2.plot(self.con3PtsAll ,  NListCon3_region2 ,  's--',color='r' ,lw=1.5, markersize=markersize,  markerfacecolor='none' )


            #ax1.plot(self.con2PtsAll  , NListCon2, '^--', color='g' ,lw=1.5)

            #ax1.plot(self.con3PtsAll ,  NListCon3, 's--',color='r' , lw=1.5)

            ax1.set_xlim(2, 27)
            ax1.set_ylim(-1, 10 )

            ax2.set_xlim(2, 27)
            ax2.set_ylim(-1, 10)
            # ax.set_ylim(-1,12)

        ax1.set_ylabel("Number of falsely detected clusters", fontsize=overAllFontSize)
        ax1.set_xlabel("MinPts" , fontsize=overAllFontSize)

        at = AnchoredText('Region 1', loc=5, frameon=False)
        ax1.add_artist(at)




        ax2.set_ylabel("Number of falsely detected clusters" , fontsize=overAllFontSize)
        ax2.set_xlabel("MinPts" , fontsize=overAllFontSize)
        at = AnchoredText('Region 2', loc=5, frameon=False)
        ax2.add_artist(at)



        ax1.legend()
        ax1.tick_params(axis='both', which='major', labelsize=overAllFontSize)
        ax1.tick_params(axis='both', which='minor', labelsize=overAllFontSize)
        ax2.tick_params(axis='both', which='major', labelsize=overAllFontSize)
        ax2.tick_params(axis='both', which='minor', labelsize=overAllFontSize)


        plt.savefig("dbscanEmptyTest.png", bbox_inches='tight',dpi=600)
        plt.savefig("dbscanEmptyTest.pdf", bbox_inches='tight' )




    def drawG2650Local(self):
        """

        :return:
        """
        drawPath = "/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/"



        areaPix  =  9
        doClean=False
        region1Tag= "G2650Local"
        TBListCon1_region1, NListCon1_region1 = self.getTBListByCon(drawPath,1, region1Tag, areaPix=areaPix , clean=doClean )
        TBListCon2_region1, NListCon2_region1 = self.getTBListByCon(drawPath,2,region1Tag ,areaPix=areaPix , clean=doClean )
        TBListCon3_region1, NListCon3_region1 = self.getTBListByCon(drawPath,3, region1Tag ,areaPix=areaPix , clean=doClean )


        print len(TBListCon1_region1)




    def getCatalog(self,rawCOFITS, savePath, fitsName, saveMarker ):

        rawDBSCANTBFile = doDBSCAN.getCatFromLabelArray(rawCOFITS, savePath+fitsName, doDBSCAN.TBModel,  saveMarker= saveMarker )


    def pipeLineG2650(self,startIndex=1,redo=False, keepFITSFile=False):

        saveTag = "G2650Local"
        self.saveTag = saveTag

        for eachCutoff in self.cutoffList[startIndex:]:

            for minPts in  self.con1PtsAll :
                self.pureDBSCAN(self.G2650Local, eachCutoff , MinPts=minPts, saveTag=saveTag,  connectivity=1, inputRMS=None,redo=redo,keepFITSFile=keepFITSFile)


            for minPts in  self.con2PtsAll :
                self.pureDBSCAN(self.G2650Local , eachCutoff , MinPts=minPts, saveTag=saveTag, connectivity=2, inputRMS=None,redo=redo ,keepFITSFile=keepFITSFile)


            for minPts in   self.con3PtsAll :
                self.pureDBSCAN(self.G2650Local , eachCutoff , MinPts=minPts, saveTag=saveTag, connectivity=3, inputRMS=None,redo=redo ,keepFITSFile=keepFITSFile)





    def pipeLineG2650ByCut(self, cutOff, redo=False,keepFITSFile = False):

        #
        #self.saveTag = saveTag


        for eachCutoff in [cutOff]:

            for minPts in  self.con1PtsAll :
                self.pureDBSCAN(self.G2650Local, eachCutoff , MinPts=minPts, saveTag= self.saveTag,  connectivity=1, inputRMS=None, redo=redo,keepFITSFile=keepFITSFile)


            for minPts in  self.con2PtsAll :
                self.pureDBSCAN(self.G2650Local , eachCutoff , MinPts=minPts, saveTag= self.saveTag, connectivity=2, inputRMS=None ,redo=redo,keepFITSFile=keepFITSFile)


            for minPts in   self.con3PtsAll :
                self.pureDBSCAN(self.G2650Local , eachCutoff , MinPts=minPts, saveTag= self.saveTag, connectivity=3, inputRMS=None ,redo=redo,keepFITSFile=keepFITSFile)



    def plotNumberSingle(self,ax,minPts,conType):
        """
        according to minPts get the Nlist, and we need a color set for minPints, 3-27,
        the three algrithms will be replaced with three connection types
        :param ax:
        :param minPts:
        :param conType:
        :return:
        """
        tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)



        tbNlist = self.getTBNlist(tbList)
        ax.plot(self.cutoffInKelvin , tbNlist, 'o-',   color= self.mMinPts.to_rgba(minPts) ,  lw=1, markersize=3)
        #self.getColorByMinPts(minPts)

    def setLabelFontSize(self,axList):
        """
        to increase the font size of the labels
        :param axList:
        :return:
        """

        ########




    def plotNumber(self):
        """
        Draw cloud number
        :return:
        """


        #looks like we have many minPts,

        con1Str = "Connectivity 1"
        con2Str = "Connectivity 2"
        con3Str = "Connectivity 3"
        scimesStr = "SCIMES"


        fig = plt.figure(figsize=(28, 7))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 22, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]



        axCon1 = fig.add_subplot(1, 4, 1)
        for p in self.con1G2650:
            self.plotNumberSingle( axCon1, p, 1 )



        axCon2 = fig.add_subplot(1, 4, 2,sharex=axCon1 ,sharey=axCon1 )
        for p in self.con2G2650:
            self.plotNumberSingle( axCon2, p, 2 )
        #axCon2.set_yticks([])
        plt.setp(axCon2.get_yticklabels(), visible=False)

        axCon3 = fig.add_subplot(1, 4, 3,sharex=axCon1  ,sharey=axCon1 )
        for p in self.con3G2650:
            self.plotNumberSingle( axCon3, p, 3 )



        #draw SCIMES
        #
        axCon4 = fig.add_subplot(1, 4, 4, sharex=axCon1  ,sharey=axCon1 )
        tbList, sigmaList = self.getSCIMESTBList()

        tbNlist = self.getTBNlist(tbList)
        axCon4.plot( np.asarray( sigmaList )*0.5, tbNlist, 'o-',   color= 'black' ,  lw=1, markersize=3)

        ###

        plt.setp(axCon4.get_yticklabels(), visible=False)



        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(axCon3 )
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        self.mMinPts._A=[]
        cb = plt.colorbar( self.mMinPts , cax=cax1)

        cax1.set_ylabel("MinPts")

        plt.setp(axCon3.get_yticklabels(), visible=False)

        axCon1.set_xlabel(r"CO cutoff (K)")
        axCon1.set_ylabel(r"Total number of molecular clouds")

        axCon2.set_xlabel(r"CO cutoff (K)")
        axCon3.set_xlabel(r"CO cutoff (K)")
        axCon4.set_xlabel(r"CO cutoff (K)")



        #label connectivity type
        at1 = AnchoredText(con1Str, loc=1, frameon=False)
        axCon1.add_artist(at1)

        at2 = AnchoredText(con2Str, loc=1, frameon=False)
        axCon2.add_artist(at2)

        at3 = AnchoredText(con3Str, loc=1, frameon=False)
        axCon3.add_artist(at3)

        at4 = AnchoredText(scimesStr, loc=1, frameon=False)
        axCon4.add_artist( at4 )

        axCon3.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
        axCon1.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
        axCon2.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
        axCon4.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))


        #
        #set up the color for minPts

        fig.tight_layout()

        plt.savefig("compareParaNumber.pdf", bbox_inches='tight')
        plt.savefig("compareParaNumber.png", bbox_inches='tight', dpi=300)

    def plotMass(self):

        """
        Plot the mass distributuion according the distance and Vlsr relationship
        :return:
        """


        fig = plt.figure(figsize=(25, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        axCon1 = fig.add_subplot(1, 4, 1)

        #draw area distribution

        self.PerFormAllSingleG2650(axCon1,self.plotMassSingle,1 )



        axCon1.set_yscale('log')
        axCon1.set_xscale('log')

        axCon2 = fig.add_subplot(1, 4, 2,sharex=axCon1 ,sharey=axCon1 )
        self.PerFormAllSingleG2650(axCon2,self.plotMassSingle,2 )

        axCon3 = fig.add_subplot(1, 4, 3,sharex=axCon1  ,sharey=axCon1 )
        self.PerFormAllSingleG2650(axCon3,self.plotMassSingle,3 )

        #for scimes
        axCon4 = fig.add_subplot(1, 4, 4,sharex=axCon1  ,sharey=axCon1 )
        self.plotMassSingle(axCon4, 4, 1,  plotSCIMES=True )

        #self.PerFormAllSingleG2650(axCon4)


        #draw color bar
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(axCon3)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        self.mMinPts._A=[]
        cb = plt.colorbar( self.mMinPts , cax=cax1)

        cax1.set_ylabel("MinPts")

        massXLabel=  r"Mass ($M_\odot$)"
        axCon1.set_xlabel( massXLabel )
        axCon2.set_xlabel( massXLabel )
        axCon3.set_xlabel( massXLabel )
        axCon4.set_xlabel( massXLabel )



        #label connectivity type
        at1 = AnchoredText(self.con1Str, loc=1, frameon=False)
        axCon1.add_artist(at1)

        at2 = AnchoredText(self.con2Str, loc=1, frameon=False)
        axCon2.add_artist(at2)

        at3 = AnchoredText(self.con3Str, loc=1, frameon=False)
        axCon3.add_artist(at3)

        at4 = AnchoredText( self.SCIMESStr , loc=1, frameon=False)
        axCon4.add_artist( at4 )



        # axSCI.set_ylabel(r"Bin number of trunks ")
        axCon1.set_ylabel(r"Number of clusters")

        self.drawCompleteMass([axCon1,axCon2,axCon3,axCon4])

        fig.tight_layout()
        plt.savefig("compareParaMass.pdf", bbox_inches='tight')
        plt.savefig("compareParaMass.png", bbox_inches='tight', dpi=300)



    def plotMassSingle(self,ax,minPts,conType, plotSCIMES=False):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:
        """
        #the area is actually the mass

        #areaEdges = np.linspace(0.25 / 3600., 150, 10000)

        areaEdges = np.linspace(0.1 , 10000, 10000)

        areaCenter = doDBSCAN.getEdgeCenter(areaEdges)
        #first, get cloud catalog
        #tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)

        if plotSCIMES:

            tbList,_ = self.getSCIMESTBList()
            plotColor="black"

        else:

            tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)
            plotColor= self.mMinPts.to_rgba(minPts)


        if 1:#two extreme cut off

            drawTB = tbList[0]
            massSigma2= doDBSCAN.getMass(drawTB)
            binN, binEdges =  np.histogram( massSigma2 , bins=areaEdges)

            ## calculate the alpha, is it 2.35?

            #self.getMassAlphaList([drawTB],completeMass=10)

            ax.plot(areaCenter[binN > 0], binN[binN > 0], color= plotColor , linestyle='-',  markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )

            drawTB = tbList[-1]
            massSigma2 = doDBSCAN.getMass(drawTB)
            binN, binEdges = np.histogram(massSigma2, bins=areaEdges)
            ax.plot(areaCenter[binN > 0], binN[binN > 0], color= plotColor , linestyle='--',
                    markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )


    def plotMassAlpha(self): #
        """
        #get and plot the alpha distribution
        :return:
        """

        fig = plt.figure(figsize=(24, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axCon1 = fig.add_subplot(1, 4, 1)
        axCon2 = fig.add_subplot(1, 4, 2, sharex=axCon1, sharey=axCon1)
        axCon3 = fig.add_subplot(1, 4, 3, sharex=axCon1, sharey=axCon1)
        axCon4 = fig.add_subplot(1, 4, 4, sharex=axCon1, sharey=axCon1)

        axList = [axCon1, axCon2, axCon3]

        self.PerFormAllSingleG2650(axCon1, self.plotMassAlphaSingle, 1)

        self.PerFormAllSingleG2650(axCon2, self.plotMassAlphaSingle, 2)

        self.PerFormAllSingleG2650(axCon3, self.plotMassAlphaSingle, 3)
        self.plotMassAlphaSingle(axCon4,1,1,plotSCIMES=True)
        axCon1.set_ylabel(r"$\alpha$ (mass)")

        axCon1.set_xlabel(r"CO cutoff (K)")
        axCon2.set_xlabel(r"CO cutoff (K)")
        axCon3.set_xlabel(r"CO cutoff (K)")
        axCon4.set_xlabel(r"CO cutoff (K)")


        #################

        self.labelContype(axList)
        self.labelSCIMES(axCon4)
        self.drawColorBar(axCon3)

        #self.drawHorizontalLine(axList, self.averageFluxAlpha ) #1.97+/0.06
        #self.drawHorizontalLine([axCon4], self.averageFluxAlpha ) #1.97+/0.06

        self.setCutoffTicks(axList+[axCon4])

        fig.tight_layout()
        plt.savefig("compareParaAlphaMass.pdf", bbox_inches='tight')
        plt.savefig("compareParaAlphaMass.png", bbox_inches='tight', dpi=300)


    def plotMassAlphaSingle(self,ax,minPts,conType, plotSCIMES=False ):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:
        """

        #alphaDendro, alphaDendroError = self.getAlphaList(tb16Den)
        #
        if plotSCIMES:
            tbList,_ = self.getSCIMESTBList()

            plotColor="black"

        else:
            tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)

            plotColor= self.mMinPts.to_rgba(minPts)

        #get tb list
        #tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)
        alphaMass, alphaMassError,completeMass= doDBSCAN.getMassAlphaList(tbList, self.cutoffList)

        ebDBSCAN = ax.errorbar(self.cutoffInKelvin , alphaMass , yerr=alphaMassError , c=  plotColor , marker='^', linestyle="-",    capsize=3, elinewidth=1.0, lw=1,   markerfacecolor='none')

        ebDBSCAN[-1][0].set_linestyle(':')







    def plotArea(self):
        """
        #plot the area, to avoide crowding, only odd MinPts and largest the smallest cutoff were drawn
        :return:
        """

        #only odd MinPts were drawn



        fig = plt.figure(figsize=(25, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        axCon1 = fig.add_subplot(1, 4, 1)

        #draw area distribution

        self.PerFormAllSingleG2650(axCon1,self.plotAreaSingle,1 )



        axCon1.set_yscale('log')
        axCon1.set_xscale('log')

        axCon2 = fig.add_subplot(1, 4, 2,sharex=axCon1 ,sharey=axCon1 )
        self.PerFormAllSingleG2650(axCon2,self.plotAreaSingle,2 )

        axCon3 = fig.add_subplot(1, 4, 3,sharex=axCon1  ,sharey=axCon1 )
        self.PerFormAllSingleG2650(axCon3,self.plotAreaSingle,3 )

        #for scimes
        axCon4 = fig.add_subplot(1, 4, 4,sharex=axCon1  ,sharey=axCon1 )
        #self.PerFormAllSingleG2650(axCon3,self.plotAreaSingle,3 )

        self.plotAreaSCIMES(axCon4)


        #draw color bar
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(axCon3)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        self.mMinPts._A=[]
        cb = plt.colorbar( self.mMinPts , cax=cax1)

        cax1.set_ylabel("MinPts")

        axCon1.set_xlabel(r"Angular area (deg$^2$)")
        axCon2.set_xlabel(r"Angular area (deg$^2$)")
        axCon3.set_xlabel(r"Angular area (deg$^2$)")
        axCon4.set_xlabel(r"Angular area (deg$^2$)")



        #label connectivity type
        at1 = AnchoredText(self.con1Str, loc=1, frameon=False)
        axCon1.add_artist(at1)

        at2 = AnchoredText(self.con2Str, loc=1, frameon=False)
        axCon2.add_artist(at2)

        at3 = AnchoredText(self.con3Str, loc=1, frameon=False)
        axCon3.add_artist(at3)

        at4 = AnchoredText( self.SCIMESStr , loc=1, frameon=False)
        axCon4.add_artist( at4 )



        # axSCI.set_ylabel(r"Bin number of trunks ")
        axCon1.set_ylabel(r"Number of clusters")

        self.drawCompleteArea([axCon1,axCon2,axCon3,axCon4])

        fig.tight_layout()
        plt.savefig("compareParaArea.pdf", bbox_inches='tight')
        plt.savefig("compareParaArea.png", bbox_inches='tight', dpi=300)


    #

    def plotAreaPhysical(self):
        """
        #plot the area, to avoide crowding, only odd MinPts and largest the smallest cutoff were drawn
        :return:
        """

        #only odd MinPts were drawn



        fig = plt.figure(figsize=(25, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        axCon1 = fig.add_subplot(1, 4, 1)

        #draw area distribution

        self.PerFormAllSingleG2650(axCon1,self.plotAreaPhysicalSingle,1 )



        axCon1.set_yscale('log')
        axCon1.set_xscale('log')

        axCon2 = fig.add_subplot(1, 4, 2,sharex=axCon1 ,sharey=axCon1 )
        self.PerFormAllSingleG2650(axCon2,self.plotAreaPhysicalSingle,2 )

        axCon3 = fig.add_subplot(1, 4, 3,sharex=axCon1  ,sharey=axCon1 )
        self.PerFormAllSingleG2650(axCon3,self.plotAreaPhysicalSingle,3 )

        #for scimes
        axCon4 = fig.add_subplot(1, 4, 4,sharex=axCon1  ,sharey=axCon1 )
        #self.PerFormAllSingleG2650(axCon3,self.plotAreaSingle,3 )

        self.plotAreaPhysicalSCIMES(axCon4)



        #draw color bar
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(axCon3)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        self.mMinPts._A=[]
        cb = plt.colorbar( self.mMinPts , cax=cax1)

        cax1.set_ylabel("MinPts")
        xLabelStr = r"Physical area ($\rm pc^{2}$)"
        axCon1.set_xlabel( xLabelStr )
        axCon2.set_xlabel(  xLabelStr )
        axCon3.set_xlabel(  xLabelStr )
        axCon4.set_xlabel(  xLabelStr )



        #label connectivity type
        at1 = AnchoredText(self.con1Str, loc=1, frameon=False)
        axCon1.add_artist(at1)

        at2 = AnchoredText(self.con2Str, loc=1, frameon=False)
        axCon2.add_artist(at2)

        at3 = AnchoredText(self.con3Str, loc=1, frameon=False)
        axCon3.add_artist(at3)

        at4 = AnchoredText( self.SCIMESStr , loc=1, frameon=False)
        axCon4.add_artist( at4 )



        # axSCI.set_ylabel(r"Bin number of trunks ")
        axCon1.set_ylabel(r"Number of clusters")

        self.drawCompleteAreaPhysical([axCon1,axCon2,axCon3,axCon4])

        fig.tight_layout()
        plt.savefig("physicalAreaDistribute.pdf", bbox_inches='tight')
        plt.savefig("physicalAreaDistribute.png", bbox_inches='tight', dpi=300)




    def drawColorBar(self,axCon3):

        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(axCon3)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        self.mMinPts._A=[]
        cb = plt.colorbar( self.mMinPts , cax=cax1)

        cax1.set_ylabel("MinPts")



    def drawCompleteArea(self,axList):

        compoleteArea = 4 * (1500. / 250) ** 2 * 0.25 / 3600.  # 0.0225

        for eachAx in axList:
            eachAx.plot([compoleteArea, compoleteArea], [1, 10000], '--', color='black', lw=1)

    def drawCompleteAreaPhysical(self,axList):

        length= 1500*np.deg2rad(0.5/60)
        compoleteArea =  length**2 * 4    #4 pixels
        for eachAx in axList:
            eachAx.plot([compoleteArea, compoleteArea], [1, 10000], '--', color='black', lw=1)


    def plotAreaPhysicalSCIMES(self, ax ):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:

        """

        physicalEdges = np.linspace(0, 100, 1000)  # square pc^2
        physicalCenter = doDBSCAN.getEdgeCenter(physicalEdges)



        #first, get cloud catalog
        tbList, sigmaList = self.getSCIMESTBList() #(self.cutoffList,minPts,conType)





        #sigmaLowestTB =  tbList[0]

        #binN , binEdges  = np.histogram(sigmaLowestTB["area_exact"] / 3600., bins=areaEdges)
        #ax.plot(areaCenter[binN > 0], binN[binN > 0], 'o-', markersize=1, lw=0.8,   alpha=0.5, color=  self.mMinPts.to_rgba(minPts)  )

        if 1:#two extreme cut off

            drawTB = tbList[0]
            #binN, binEdges = np.histogram(drawTB["area_exact"] / 3600., bins=areaEdges)
            #ax.plot(areaCenter[binN > 0], binN[binN > 0], color= 'black', linestyle='-',  markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )

            realAreaSigma2 = doDBSCAN.getRealArea(drawTB)
            binNSigma2, binEdgesSigma2 = np.histogram(realAreaSigma2, bins=physicalEdges)
            ax.plot(physicalCenter[binNSigma2 > 0], binNSigma2[binNSigma2 > 0], color= 'black' ,  linestyle='-', markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )



            drawTB = tbList[-1]
            realAreaSigma7 = doDBSCAN.getRealArea(drawTB)
            binNSigma7, binEdgesSigma7 = np.histogram(realAreaSigma7, bins=physicalEdges)
            ax.plot(physicalCenter[binNSigma7 > 0], binNSigma7[binNSigma7 > 0], color= 'black' ,  linestyle='--', markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full',  alpha=0.6, )

            #drawTB = tbList[-1]
            #binN, binEdges = np.histogram(drawTB["area_exact"] / 3600., bins=areaEdges)
            #ax.plot(areaCenter[binN > 0], binN[binN > 0], color= 'black' , linestyle='--', markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )





    def plotAreaSCIMES(self, ax ):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:

        """

        areaEdges = np.linspace(0.25 / 3600., 150, 10000)
        areaCenter = doDBSCAN.getEdgeCenter(areaEdges)
        #first, get cloud catalog
        tbList, sigmaList = self.getSCIMESTBList() #(self.cutoffList,minPts,conType)


        if 1:#two extreme cut off

            drawTB = tbList[0]
            binN, binEdges = np.histogram(drawTB["area_exact"] / 3600., bins=areaEdges)
            ax.plot(areaCenter[binN > 0], binN[binN > 0], color= 'black', linestyle='-',  markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )



            drawTB = tbList[-1]
            binN, binEdges = np.histogram(drawTB["area_exact"] / 3600., bins=areaEdges)
            ax.plot(areaCenter[binN > 0], binN[binN > 0], color= 'black' , linestyle='--', markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )




    def plotAreaPhysicalSingle(self,ax,minPts,conType ):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:
        """

        physicalEdges = np.linspace(0, 100, 1000)  # square pc^2
        physicalCenter = doDBSCAN.getEdgeCenter(physicalEdges)






        tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)


        #plot sigma 1
        drawTB = tbList[0]
        realAreaSigma2 = doDBSCAN.getRealArea(drawTB)
        binNSigma2, binEdgesSigma2 = np.histogram(realAreaSigma2, bins=physicalEdges)
        ax.plot(physicalCenter[binNSigma2 > 0], binNSigma2[binNSigma2 > 0], color=self.mMinPts.to_rgba(minPts), linestyle='-',  markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )



        #plot sigma 7

        drawTB = tbList[-1]
        realAreaSigma7 = doDBSCAN.getRealArea(drawTB)
        binNSigma7, binEdgesSigma7 = np.histogram(realAreaSigma7, bins=physicalEdges)
        ax.plot(physicalCenter[binNSigma7 > 0], binNSigma7[binNSigma7 > 0], color=self.mMinPts.to_rgba(minPts),   linestyle='--',  markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )



    def plotAreaSingle(self,ax,minPts,conType ):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:
        """

        areaEdges = np.linspace(0.25 / 3600., 150, 10000)
        areaCenter = doDBSCAN.getEdgeCenter(areaEdges)
        #first, get cloud catalog
        tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)


        #sigmaLowestTB =  tbList[0]

        #binN , binEdges  = np.histogram(sigmaLowestTB["area_exact"] / 3600., bins=areaEdges)
        #ax.plot(areaCenter[binN > 0], binN[binN > 0], 'o-', markersize=1, lw=0.8,   alpha=0.5, color=  self.mMinPts.to_rgba(minPts)  )

        if 0:#all cutoffs
            for i in range( len(tbList) ):

                drawTB=tbList[i]

                cutOff=self.cutoffList[i]

                binN , binEdges  = np.histogram(drawTB["area_exact"] / 3600., bins=areaEdges)

                ax.plot(areaCenter[binN > 0], binN[binN > 0], color=  self.mMinPts.to_rgba(minPts), linestyle='-', markersize=cutOff, lw=0.6, marker='o' , markeredgewidth=0.5, fillstyle='none', alpha=0.6,   )


        if 1:#two extreme cut off

            drawTB = tbList[0]

            cutOff = self.cutoffList[0]

            binN, binEdges = np.histogram(drawTB["area_exact"] / 3600., bins=areaEdges)

            ax.plot(areaCenter[binN > 0], binN[binN > 0], color=self.mMinPts.to_rgba(minPts), linestyle='-',  markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )



            drawTB = tbList[-1]

            #cutOff = self.cutoffList[-5]
            print -1
            binN, binEdges = np.histogram(drawTB["area_exact"] / 3600., bins=areaEdges)
            ax.plot(areaCenter[binN > 0], binN[binN > 0], color=self.mMinPts.to_rgba(minPts), linestyle='--',  markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )



    def plotFlux(self):
        """
        #plot the distribution of flux
        :return:
        """
        #plot the distribution of flux


        fig = plt.figure(figsize=(24, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        #need to draw complete position


        #print "Complete, ", completeFlux, self.rms, sigmaListDen[i]
        #ax.plot([completeFlux, completeFlux], [2, 3000], '--', markersize=1, lw=0.8, alpha=0.5, color=clrs[i])
        ############################################
        axCon1 = fig.add_subplot(1, 4, 1)

        ###########################
        axCon2 = fig.add_subplot(1, 4, 2, sharex=axCon1 ,sharey=axCon1 )

        #####################################
        axCon3 = fig.add_subplot(1, 4, 3, sharex=axCon1  ,sharey=axCon1 )
        axCon4 = fig.add_subplot(1, 4, 4, sharex=axCon1  ,sharey=axCon1 )
        self.plotFluxSingle(axCon4,1,1, plotSCIMES=True )

        self.PerFormAllSingleG2650(axCon1,self.plotFluxSingle,1 )
        self.PerFormAllSingleG2650(axCon2,self.plotFluxSingle,2 )
        self.PerFormAllSingleG2650(axCon3,self.plotFluxSingle,3 )

        at4 = AnchoredText( self.SCIMESStr , loc=1, frameon=False)
        axCon4.add_artist( at4 )

        ######################################

        axList=[axCon1,axCon2,axCon3]

        self.drawCompleteFlux(axList)
        self.drawCompleteFlux( [axCon4] )

        self.labelContype(axList)

        #draw color bar
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(axCon3)
        cax1 = divider.append_axes("right", size="5%", pad=0.05)
        self.mMinPts._A=[]
        cb = plt.colorbar( self.mMinPts , cax=cax1)

        cax1.set_ylabel("MinPts")

        ####
        axCon1.set_ylabel(r"Number of clusters")

        strXlabel = r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_{\rm A}$)"

        axCon1.set_xlabel( strXlabel )
        axCon2.set_xlabel( strXlabel )
        axCon3.set_xlabel( strXlabel )
        axCon4.set_xlabel( strXlabel )


        #label connectivity type

        axCon1.set_yscale('log')
        axCon1.set_xscale('log')


        fig.tight_layout()
        plt.savefig("compareParaFlux.pdf", bbox_inches='tight')
        plt.savefig("compareParaFlux.png", bbox_inches='tight', dpi=300)


    def labelContype(self,axList):

        axCon1,axCon2,axCon3=axList


        at1 = AnchoredText(self.con1Str, loc=1, frameon=False)
        axCon1.add_artist(at1)

        at2 = AnchoredText(self.con2Str, loc=1, frameon=False)
        axCon2.add_artist(at2)

        at3 = AnchoredText(self.con3Str, loc=1, frameon=False)
        axCon3.add_artist(at3)

    def labelSCIMES(self,ax):

        at = AnchoredText(self.SCIMESStr, loc=1, frameon=False)
        ax.add_artist(at)



    def drawCompleteFlux(self,axList):

        completeFluxCut2 = 144 * self.rmsCO12 * 0.2 * self.cutoffList[0] * 3  # K km/s, the last 3  is the two channels
        completeFluxCut7 = 144 * self.rmsCO12 * 0.2 * self.cutoffList[-1] * 3  # K km/s, the last 3  is the two channels

        for eachAx in axList:
            eachAx.plot([completeFluxCut2, completeFluxCut2], [2, 10000], '-', lw=1, alpha=0.5, color='k')
            eachAx.plot([completeFluxCut7, completeFluxCut7], [2, 10000], '--', lw=1, alpha=0.5, color='k')


    def drawCompleteMass(self,axList):

        completeFluxCut2 = 144 * self.rmsCO12 * 0.2 * self.cutoffList[0] * 3  # K km/s, the last 3  is the two channels

        completeMassCut2 = doDBSCAN.calmassByXfactor(completeFluxCut2,1500)

        completeFluxCut7 = 144 * self.rmsCO12 * 0.2 * self.cutoffList[-1] * 3  # K km/s, the last 3  is the two channels
        completeMassCut7 = doDBSCAN.calmassByXfactor(completeFluxCut7,1500)

        for eachAx in axList:
            eachAx.plot([completeMassCut2, completeMassCut2], [1, 1000], '-', lw=1, alpha=0.5, color='k')
            eachAx.plot([completeMassCut7, completeMassCut7], [1, 1000], '--', lw=1, alpha=0.5, color='k')





    def plotFluxSingle(self,ax,minPts,conType,plotSCIMES=False):

        """

        :return:
        """
        #first get tables

        #first, get cloud catalog

        if plotSCIMES:
            tbList,_ = self.getSCIMESTBList()

            plotColor="black"
        else:
            tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)

            plotColor= self.mMinPts.to_rgba(minPts)



        fluxEdges = np.linspace(8, 1e5, 1000)
        fluxCenter = doDBSCAN.getEdgeCenter( fluxEdges )



        if 1:#two extreme cut off

            drawTB = tbList[0]

            cutOff = self.cutoffList[0]

            sum = drawTB["sum"] * 0.2  # K km/s
            binN, binEdges  = np.histogram(sum, bins=fluxEdges)

            ax.plot(fluxCenter[binN > 0], binN[binN > 0], color=plotColor, linestyle='-',
                    markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )



            drawTB = tbList[-1]

            cutOff = self.cutoffList[0]

            sum = drawTB["sum"] * 0.2  # K km/s
            binN, binEdges  = np.histogram(sum, bins=fluxEdges)

            ax.plot(fluxCenter[binN > 0], binN[binN > 0], color=plotColor , linestyle='--',
                    markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, )



    def setCutoffTicks(self,axList ):
        """

        :param axList:
        :return:
        """


        for eachAx in axList:

            eachAx.set_xticks(   np.arange(1, 4, step=0.5 )   )




    def plotTotalFlux(self):
        """

        #draw a comparison to total flux, comparing three algrothms

        :return:
        """

        #first

        fig = plt.figure(figsize=(26, 6.5 ))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axCon1 = fig.add_subplot(1, 4, 1)
        axCon2 = fig.add_subplot(1, 4, 2,sharex=axCon1 ,sharey=axCon1 )
        axCon3 = fig.add_subplot(1, 4, 3,sharex=axCon1  ,sharey=axCon1 )
        axCon4 = fig.add_subplot(1, 4, 4, sharex=axCon1, sharey=axCon1)

        axList=[axCon1,axCon2,axCon3]
        ################# Draw con1

        self.PerFormAllSingleG2650(axCon1,self.plotTotalFluxSingle,1 )

        self.plotTotalFluxSingle( axCon1, 4, 1, plotSCIMES=False, totalFluxCon1=None)


        #draw total flux

        totalFluxListCon1=doAllDBSCAN.getTotalFluxByMask( doAllDBSCAN.fitsPath+ "G2650LocaldbscanS2P4Con1_CO_Masked.fits" )
        axCon1.plot(self.cutoffInKelvin , totalFluxListCon1, 'o--',   color= 'k' ,  lw=1, markersize=3,label="Total flux above the cutoff")
        axCon1.legend(loc =3 )
        print "Fitting line of Connectivity 1"
        self.fitFluxLine( self.cutoffInKelvin,  totalFluxListCon1  )

        ################# Draw con2

        self.PerFormAllSingleG2650(axCon2,self.plotTotalFluxSingle,2 )



        totalFluxListCon2=doAllDBSCAN.getTotalFluxByMask( doAllDBSCAN.fitsPath+ "G2650LocaldbscanS2P8Con2_CO_Masked.fits" )
        axCon2.plot(self.cutoffInKelvin , totalFluxListCon2, 'o--', color= 'k' ,   lw=1, markersize=3  ,label="Total flux above the cutoff")

        axCon2.legend(loc =3 )
        print "Fitting line of Connectivity 2"
        self.fitFluxLine( self.cutoffInKelvin,  totalFluxListCon2  )

        ################# Draw con3Z
        self.PerFormAllSingleG2650(axCon3,self.plotTotalFluxSingle,3 )

        totalFluxListCon3 =doAllDBSCAN.getTotalFluxByMask( doAllDBSCAN.fitsPath+ "G2650LocaldbscanS2P11Con3_CO_Masked.fits" )
        axCon3.plot(self.cutoffInKelvin , totalFluxListCon3 , 'o--',   color= 'k' ,  lw=1, markersize=3 ,label="Total flux above the cutoff")

        axCon3.legend(loc =3 )

        print "Fitting line of Connectivity 3"
        self.fitFluxLine( self.cutoffInKelvin,  totalFluxListCon3  )

        ################# Draw scimes
        self.plotTotalFluxSingle( axCon4,1,1,plotSCIMES=True, totalFluxCon1=totalFluxListCon1  )

        #totalFluxListCon3 =doAllDBSCAN.getTotalFluxByMask(   "/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/ClusterAsgn_2_16Ve20_CO_Masked.fits"   )
        #axCon4.plot(self.cutoffInKelvin , totalFluxListCon3 , 'o--',   color= 'k' ,  lw=1, markersize=3 ,label="Total flux above the cutoff")

        #axCon3.legend(loc =3 )


        at4 = AnchoredText( self.SCIMESStr , loc=1, frameon=False)
        axCon4.add_artist( at4 )



        ########
        axCon1.set_ylabel(r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_\mathrm{A}$)")
        axCon1.set_xlabel(r"CO cutoff (K)")
        axCon2.set_xlabel(r"CO cutoff (K)")
        axCon3.set_xlabel(r"CO cutoff (K)")
        axCon4.set_xlabel(r"CO cutoff (K)")

        ####
        self.setCutoffTicks( axList  +[axCon4] )
        self.labelContype(axList)

        self.drawColorBar(axCon3)

        fig.tight_layout()
        plt.savefig("compareParaTotalFlux.pdf", bbox_inches='tight')
        plt.savefig("compareParaTotalFlux.png", bbox_inches='tight', dpi=300)





    def fitFluxLine(self,cutoffInKelvin,totalFluxList):
        """
        :return:
        """
        #convert sigma list to Kelvin

        z = np.polyfit( cutoffInKelvin ,  totalFluxList, 1)
        p = np.poly1d(z)

        #print totalFluxList, "????????????????????????????"
        print z
        #print p(0.0) / p(3), p(0.0) / p(2)
        print p(1.5) / p(0), p(1) / p(0)



    def getCutOFF(self,tbFileName,breakPoint=None):

        filePath,fileName=os.path.split(tbFileName)

        if breakPoint==None:

            _,getCutoff=fileName.split("dbscanS")

        else:


            _,getCutoff=fileName.split( breakPoint )



        cutSigma, _ =getCutoff.split("P")

        return float(cutSigma)



    def getMinPts(self,tbFileName):

        filePath,fileName=os.path.split(tbFileName)

        _,getminPts= fileName.split("dbscanS")

        _, getminPts =getminPts.split("P")

        getminPts,  _  =getminPts.split("Con")

        return int(getminPts)


    def getConType(self,tbFileName):

        filePath,fileName=os.path.split(tbFileName)

        _,conTypeStr  =fileName.split("Con")

        if "_CleanWithLW" in conTypeStr:
            conType, _ = conTypeStr.split("_CleanWithLW")
            return int(conType )
        else:

            conType, _ = conTypeStr.split(".fit")
            return int(conType )



    def getSaveTag(self,tbFileName):

        filePath, fileName = os.path.split(tbFileName)

        saveTag, _ = fileName.split("dbscanS")

        return saveTag




    def getMaskByLabel(self,FITSfile, tbFile , onlyClean=False ,onlySelect=False,inputCutOff=None,reProduceFITS=True,removeEdge=True,minDelta=3):
        """

        :param FITSfile:
        :param tbFile:
        :return:
        """

        rawTB=Table.read(tbFile)

        #find cutoff automatically

        #split the name of the table
        if inputCutOff==None:
            cutOff= self.getCutOFF( tbFile  )
        else:
            cutOff=inputCutOff

        #the parameters are final

        filterTB=self.selectTBFormal(rawTB,cutOff=cutOff, pixN=16,minDelta=minDelta,hasBeam=True,minChannel=3 ,removeEdge=removeEdge)

        if onlySelect:
            #saveFITS
            return filterTB


        saveFITS = FITSfile[0:-5] + "_Clean.fits"
        maskFITS = FITSfile[0:-5] + "_CO_Masked.fits"

        ## reprodce fits files

        #####
        preExist=os.path.isfile(FITSfile)
        if not preExist and reProduceFITS:
            #reproduct the fits
            cutOff=self.getCutOFF(FITSfile)
            minPts=self.getMinPts(FITSfile)
            conType=self.getConType(FITSfile)
            saveTag=self.getSaveTag(FITSfile)

            self.pureDBSCAN( self.G2650Local, cutOff, MinPts=minPts, saveTag=saveTag, connectivity=conType , inputRMS=None, redo=True,  keepFITSFile=True, onlyGetFITS=True)

        ############






        dataCluster, headCluster = myFITS.readFITS(FITSfile)
        cleanCluster = doDBSCAN.cleanLabelFITS(dataCluster, filterTB)

        #write down clean fits
        fits.writeto(saveFITS, cleanCluster, headCluster, overwrite=True)

        if onlyClean:
            print saveFITS
            print "Program existing due to the set of only cleaning..."
            return



        print "The cutOff is ", cutOff

        print "The number of selected clouds, ", len(filterTB)



        dataCO,headCO=myFITS.readFITS(self.G2650Local)
        noiseLabel=np.min( cleanCluster[0] )
        dataCO[cleanCluster==noiseLabel ] = 0

        fits.writeto(maskFITS, dataCO, headCO, overwrite=True)



    def PerFormAllSingleG2650(self,ax,singleFunction,conType, plot3=False):


        if plot3:
            if conType==1:
                conList=self.con1G2650

            if conType==2:
                conList=self.con2G2650

            if conType==3:
                conList=self.con3G2650


            singleFunction(ax, np.min(conList ) , conType)
            singleFunction(ax, np.max(conList ) , conType)
            singleFunction(ax, self.getMiddleMinPts(conList ) , conType)




        else:

            if conType==1:
                conList=self.con1G2650

            if conType==2:
                conList=self.con2G2650

            if conType==3:
                conList=self.con3G2650

            for p in conList:
                if p%2 == 0:
                    continue

                singleFunction(ax,p, conType)



    def getTotalFluxList(self,tbList):

        fluxList= []

        for eachTB in tbList:

            sumFlux=np.sum(eachTB["sum"])*0.2 # to  K km/s

            fluxList.append(sumFlux)

        return fluxList





    def plotTotalFluxSingle(self,ax,minPts,conType,plotSCIMES=False , totalFluxCon1=None ):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:
        """
        #get all tables

        if plotSCIMES:
            tbList,_=self.getSCIMESTBList()
            plotColor="black"

        else:
            tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)
            plotColor= self.mMinPts.to_rgba(minPts)

        fluxList=self.getTotalFluxList(tbList)

        if totalFluxCon1!=None:
            print "SCIMES ratio to the total flux"
            print 100*np.asarray( fluxList )/totalFluxCon1[0]



        ax.plot(self.cutoffInKelvin , fluxList, 'o-',   color= plotColor  ,  lw=1, markersize=3)


    def getTotalFluxByMask(self,maskedCOFITS):

        dataCOMask,headCOmask=myFITS.readFITS(maskedCOFITS)

        totalFluxList=[]

        for eachCut in self.cutoffInKelvin:

            dataCOMask[ dataCOMask<eachCut]=0
            totalFluxSum=np.sum(dataCOMask)*0.2
            totalFluxList.append(totalFluxSum )

        return totalFluxList



    def getEquivalentLineWidth(self,fitsFile, TBFile,saveSpectral=False, saveSpectralPath="", reProduceFITS=False ,keepFITS=False ):

        """
        get line width and save the table file according to fitsFile and TBFile
        :param fitsFile:
        :param TBFile:
        :return:
        """

        saveTBAs=fitsFile[:-5]+"_CleanWithLW.fit"
        regionName= os.path.split(fitsFile)[1][:-5]

        #####
        preExist=os.path.isfile(fitsFile)
        if not preExist and reProduceFITS:
            #reproduct the fits
            cutOff=self.getCutOFF(fitsFile)
            minPts=self.getMinPts(fitsFile)
            conType=self.getConType(fitsFile)
            saveTag=self.getSaveTag(fitsFile)

            self.pureDBSCAN( self.G2650Local, cutOff, MinPts=minPts, saveTag=saveTag, connectivity=conType , inputRMS=None, redo=True,  keepFITSFile=True, onlyGetFITS=True)




        #first read label and CO data
        dataCluster,headCluster=myFITS.readFITS(fitsFile)
        dataCO,headCO= myFITS.readFITS(self.G2650Local)



        #read TBa
        rawTB=Table.read(TBFile)
        filterTB= self.selectTBFormal(rawTB)  # remove unrelated sources  #self.selectTB(rawTB)


        #based on filter table, calculate line width sanve save the file

        #create 1 dimensional index for data cluster


        #
        #
        cubeCO = SpectralCube.read(self.G2650Local)
        # #
        vAxis = cubeCO.spectral_axis
        # print vAxis.real.to_value()

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
            pass

        i = 0
        for eachR in filterTB:

            testID = int(eachR["_idx"])

            testIndices = doDBSCAN.getIndices(Z0, Y0, X0, values1D, testID)
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

        #save the table
        filterTB.write(saveTBAs, overwrite=True)
        if not keepFITS and reProduceFITS and not preExist:
            os.remove(fitsFile)


    def getFITSAndTBFileList(self,cutOff,conType):
        """

        :param cutOff:
        :param conType:
        :return:
        """

        fitsList=[]

        tbList=[]

        minPtsList=self.conTypeG2650[conType] #

        for eachP in minPtsList:

            fitsName=self.fitsPath+"G2650LocaldbscanS{}P{}Con{}.fits".format( cutOff,  eachP  , conType )
            tableName =  self.fitsPath+"G2650LocaldbscanS{}P{}Con{}.fit".format( cutOff,  eachP  , conType )


            fitsList.append( fitsName)
            tbList.append( tableName )


        return fitsList,tbList


    def plotPeakSingle(self,ax,minPts,conType,  plotSCIMES= False ):
        """
        plot Peak distribution
        #Is this meaning full?

        :return:
        """
        #only plot largest and lowest cutoff

        #get tables

        if plotSCIMES:

            tbList,_ = self.getSCIMESTBList()
            plotColor="black"

        else:
            tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)
            plotColor= self.mMinPts.to_rgba(minPts)

        velEdges = np.linspace(0, 40, 400)
        velCenter = doDBSCAN.getEdgeCenter(velEdges)


        #draw table1
        tableSigma2 =  tbList[0]

        vDisperse = tableSigma2["peak"]

        #peakV = np.max(vDisperse)
        #meanV = np.mean(vDisperse)
        #medianV = np.median(vDisperse)
        #labelSigma2=  r"{:.1f}$\sigma$, P8".format( 2)
        #suffix = "\nMean: {:.1f}, Median: {:.1f}".format(meanV, medianV)
        #print np.min(vDisperse), np.max(vDisperse)

        binN , binEdges  = np.histogram(vDisperse, bins=velEdges)
        #stepa = ax.step(velCenter, binN, lw=1.0, label=labelSigma2 + suffix)
        stepa = ax.step(velCenter, binN, lw= 0.8 ,  linestyle='-', color= plotColor   )


        #stepa = ax.plot(velCenter, binN, lw= 0.8 ,  linestyle='-', color=self.mMinPts.to_rgba(minPts)  )
        #ax.plot( velCenter  , binN, 'o-',   color= self.mMinPts.to_rgba(minPts) , linestyle='-',  lw=0.8, markersize=1)


        tableSigma3 =  tbList[-1]
        vDisperse = tableSigma3["peak"]
        binN , binEdges  = np.histogram(vDisperse, bins=velEdges)

        stepa = ax.step(velCenter, binN, lw=0.8,  linestyle='--' ,  color=  plotColor   )
        #ax.plot( velCenter  , binN, 'o-',   color= self.mMinPts.to_rgba(minPts) ,  linestyle='--',  lw=0.8, markersize=1)


    def numberWithComma(self,axList):


        for eachAx in axList:
            eachAx.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))


    def plotPeak(self):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:
        """

        fig = plt.figure(figsize=(26, 6.3))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 19, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axCon1 = fig.add_subplot(1, 4, 1)
        axCon2 = fig.add_subplot(1, 4, 2, sharex=axCon1  , sharey=axCon1 )
        axCon3 = fig.add_subplot(1, 4, 3, sharex=axCon1  , sharey=axCon1 )
        axCon4 = fig.add_subplot(1, 4, 4, sharex=axCon1  , sharey=axCon1 )

        axList=[axCon1,axCon2,axCon3]
        ###############




        #self.PerFormAllSingleG2650(axCon1,self.plotPeakSingle,1 )
        self.plotPeakSingle(axCon1,  np.min(self.con1G2650),1 )
        self.plotPeakSingle(axCon1, self.getMiddleMinPts(self.con1G2650),1 )
        self.plotPeakSingle(axCon1,  np.max(self.con1G2650),1 )

        self.plotPeakSingle(axCon2,  np.min(self.con2G2650),2 )
        self.plotPeakSingle(axCon2, self.getMiddleMinPts(self.con2G2650),2 )
        self.plotPeakSingle(axCon2,  np.max(self.con2G2650),2 )

        self.plotPeakSingle(axCon3,  np.min(self.con3G2650),3 )
        self.plotPeakSingle(axCon3, self.getMiddleMinPts(self.con3G2650),3 )
        self.plotPeakSingle(axCon3,  np.max(self.con3G2650),3 )

        self.plotPeakSingle(axCon4,  1 ,3, plotSCIMES=True )



        axCon1.set_xlim(0, 10)
        ########
        axCon1.set_ylabel(r"Number of clusters")
        axCon1.set_xlabel(r"Peak brightness temperature (K)" )
        axCon2.set_xlabel( r"Peak brightness temperature (K)" )
        axCon3.set_xlabel( r"Peak brightness temperature (K)" )
        axCon4.set_xlabel( r"Peak brightness temperature (K)" )



        axCon4.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
        axCon3.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
        axCon1.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
        axCon2.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
        ####

        self.labelContype(axList)
        self.labelSCIMES(axCon4 )
        self.drawColorBar(axCon3)


        fig.tight_layout()
        plt.savefig("peakDistribute.pdf", bbox_inches='tight')
        plt.savefig("peakDistribute.png", bbox_inches='tight', dpi=300)

    def calLineWidthPipeLine(self,  cutOff ,reProduceFITS=False, keepFITS=False ):
        """

        :return:
        """



        for conType in [1,2,3]:

            fitsList,tbList=self.getFITSAndTBFileList(cutOff,conType)


            for i in range(len(tbList)):

                fitsName=fitsList[i]
                tableName= tbList[i ]


                if  fitsName[:-5]  == tableName[:-4]:
                    print "Processing,",  fitsName
                    doAllDBSCAN.getEquivalentLineWidth( fitsName ,  tableName ,reProduceFITS=reProduceFITS, keepFITS=keepFITS )



    def plotLineWidthSingle(self,ax, minPts,conType,plotSCIMES=False):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:
        """

        #get all tables
        if plotSCIMES:
            tbList,_ = self.getSCIMESTBList()
            plotColor = "black"
        else:
            tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType, getCleanLW=True)
            plotColor =  self.mMinPts.to_rgba(minPts) 



        velEdges = np.linspace(0, 15, 300)
        velCenter = doDBSCAN.getEdgeCenter(velEdges)
        #################################
        tbSigma2=tbList[0]
        vDisperse = tbSigma2["lineWidth"]
        #maxV = np.max(vDisperse)
        #meanV = np.mean(vDisperse)
        #medianV = np.median(vDisperse)

        binN , binEdges  = np.histogram(vDisperse, bins=velEdges)

        #peakV = areaCenter[binN8.argmax()]

        #suffix = "\nPeak: {:.2f}, Mean: {:.2f}, Median: {:.2f}".format(peakV, meanV, medianV)

        ax.step(velCenter, binN, lw= 0.8 ,  linestyle='-', color=plotColor )

        ###############
        tbSigma7=tbList[-1]
        vDisperse = tbSigma7["lineWidth"]

        binN , binEdges  = np.histogram(vDisperse, bins=velEdges)
        ax.step(velCenter, binN, lw= 0.8 ,  linestyle='--', color= plotColor  )



    def plotVrmsSingle(self,ax, minPts,conType, plotSCIMES=False ):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:
        """

        #get all tables

        #get all tables
        if plotSCIMES:
            tbList, _  = self.getSCIMESTBList()
            plotColor = "black"
        else:
            tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType, getCleanLW=True)
            plotColor =  self.mMinPts.to_rgba(minPts)


        velEdges = np.linspace(0, 15, 300)
        velCenter = doDBSCAN.getEdgeCenter(velEdges)
        #################################
        tbSigma2=tbList[0]
        vDisperse = tbSigma2["v_rms"]
        #maxV = np.max(vDisperse)
        #meanV = np.mean(vDisperse)
        #medianV = np.median(vDisperse)

        binN , binEdges  = np.histogram(vDisperse, bins=velEdges)

        #peakV = areaCenter[binN8.argmax()]

        #suffix = "\nPeak: {:.2f}, Mean: {:.2f}, Median: {:.2f}".format(peakV, meanV, medianV)

        ax.step(velCenter, binN, lw= 0.8 ,  linestyle='-', color= plotColor   )
        #ax.plot(velCenter, binN, lw= 0.8 ,marker='o', markersize=1,    linestyle='-', color=self.mMinPts.to_rgba(minPts)  )

        ###############
        tbSigma7=tbList[-1]
        vDisperse = tbSigma7["v_rms"]

        binN , binEdges  = np.histogram(vDisperse, bins=velEdges)

        ax.step(velCenter, binN, lw= 0.8 ,  linestyle='--', color= plotColor  )
        #ax.plot(velCenter, binN, lw= 0.8 ,marker='o', markersize=1,   linestyle='--', color=self.mMinPts.to_rgba(minPts)  )





    def plotLineWidth(self):

        """
        :return:
        """

        fig = plt.figure(figsize=(26, 12))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 20, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axCon1 = fig.add_subplot(2, 4, 1)
        axCon2 = fig.add_subplot(2, 4, 2, sharex=axCon1 ,sharey=axCon1 )
        axCon3 = fig.add_subplot(2, 4, 3, sharex=axCon1  ,sharey=axCon1 )
        axCon4 = fig.add_subplot(2, 4, 4, sharex=axCon1  ,sharey=axCon1 )


        axList1=[axCon1,axCon2,axCon3]




        ########
        self.PerFormAllSingleG2650(axCon1,self.plotLineWidthSingle,1 ,plot3=True)
        self.PerFormAllSingleG2650(axCon2,self.plotLineWidthSingle,2 ,plot3=True)
        self.PerFormAllSingleG2650(axCon3,self.plotLineWidthSingle,3 ,plot3=True)

        self.plotLineWidthSingle( axCon4, 1, 1, plotSCIMES=True )

        ########## Vrms
        axCon5 = fig.add_subplot(2, 4, 5, sharex=axCon1 ,sharey=axCon1 )
        axCon6 = fig.add_subplot(2, 4, 6, sharex=axCon1 ,sharey=axCon1 )
        axCon7 = fig.add_subplot(2, 4, 7, sharex=axCon1  ,sharey=axCon1 )
        axCon8 = fig.add_subplot(2, 4, 8, sharex=axCon1  ,sharey=axCon1 )

        self.plotVrmsSingle( axCon8,1,1, plotSCIMES=True )



        axList2=[axCon5,axCon6,axCon7]

        #

        #axSCIMES.set_xlabel(r"Equivalent linewidth ($\rm km\ s$$^{-1}$)")
        axCon5.set_ylabel(r"Number of clusters")
        axCon1.set_ylabel(r"Number of clusters")
        
        axCon5.set_xlabel(r"Velocity dispersion ($\rm km\ s$$^{-1}$)")
        axCon6.set_xlabel(r"Velocity dispersion ($\rm km\ s$$^{-1}$)")
        axCon7.set_xlabel(r"Velocity dispersion ($\rm km\ s$$^{-1}$)")
        axCon8.set_xlabel(r"Velocity dispersion ($\rm km\ s$$^{-1}$)")



        #self.plotVrmsSingle(axCon4,3,1)

        self.PerFormAllSingleG2650(axCon5,self.plotVrmsSingle,1 ,plot3=True)
        self.PerFormAllSingleG2650(axCon6,self.plotVrmsSingle,2 ,plot3=True)
        self.PerFormAllSingleG2650(axCon7,self.plotVrmsSingle,3  ,plot3=True)

        axCon5.set_xlim([0, 3])
        axCon5.set_ylim([-500, 3000])


        self.labelContype(axList1)
        self.labelContype(axList2)

        self.labelSCIMES(axCon4)
        self.labelSCIMES(axCon8)



        self.drawColorBar(axCon3)
        self.drawColorBar(axCon7)



        for eachAx in axList1+[axCon4]:

            at = AnchoredText("Eqivalent width", loc=5, frameon=False )
            eachAx.add_artist(at)

        for eachAx in axList2+[axCon8]:

            at = AnchoredText("the second moment", loc=5, frameon=False )
            eachAx.add_artist(at)




        fig.tight_layout()
        plt.savefig("velDistribute.pdf", bbox_inches='tight')
        plt.savefig("velDistribute.png", bbox_inches='tight', dpi=300)


    def plotAreaAlpha(self):
        """
        #get and plot the alpha distribution
        :return:
        """



        fig = plt.figure(figsize=(24, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axCon1 = fig.add_subplot(1, 4, 1)
        axCon2 = fig.add_subplot(1, 4, 2,sharex=axCon1 ,sharey=axCon1 )
        axCon3 = fig.add_subplot(1, 4, 3,sharex=axCon1  ,sharey=axCon1 )
        axCon4 = fig.add_subplot(1, 4, 4,sharex=axCon1  ,sharey=axCon1 )

        
        axList=[axCon1,axCon2,axCon3]



        self.PerFormAllSingleG2650(axCon1,self.plotAreaAlphaSingle,1 )
        self.PerFormAllSingleG2650(axCon2,self.plotAreaAlphaSingle,2 )
        self.PerFormAllSingleG2650(axCon3,self.plotAreaAlphaSingle,3 )
        self.plotAreaAlphaSingle(axCon4, 1,1, plotSCIMES=True)


        axCon1.set_ylabel(r"$\alpha$ (angular area)")

        axCon1.set_xlabel(r"CO cutoff (K)")
        axCon2.set_xlabel(r"CO cutoff (K)")
        axCon3.set_xlabel(r"CO cutoff (K)")
        axCon4.set_xlabel(r"CO cutoff (K)")


        #################

        self.labelContype(axList)
        self.labelSCIMES(axCon4)
        self.drawColorBar(axCon3)
        self.drawHorizontalLine(axList, 1.97 ) #1.97+/0.06
        self.drawHorizontalLine([axCon4], 1.97 ) #1.97+/0.06
        self.setCutoffTicks(axList+[axCon4])

        fig.tight_layout()
        plt.savefig("compareParaAlphaArea.pdf", bbox_inches='tight')
        plt.savefig("compareParaAlphaArea.png", bbox_inches='tight', dpi=300)

    def plotAreaAlphaPhysical(self):
        """
        #plot the power-law index, calculated with physcial
        :return:
        """

        fig = plt.figure(figsize=(24, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axCon1 = fig.add_subplot(1, 4, 1)
        axCon2 = fig.add_subplot(1, 4, 2, sharex=axCon1, sharey=axCon1)
        axCon3 = fig.add_subplot(1, 4, 3, sharex=axCon1, sharey=axCon1)
        axCon4 = fig.add_subplot(1, 4, 4, sharex=axCon1, sharey=axCon1)

        axList = [axCon1, axCon2, axCon3]

        self.PerFormAllSingleG2650(axCon1, self.plotAreaAlphaPhysicalSingle, 1)
        self.PerFormAllSingleG2650(axCon2, self.plotAreaAlphaPhysicalSingle, 2)
        self.PerFormAllSingleG2650(axCon3, self.plotAreaAlphaPhysicalSingle, 3)
        self.plotAreaAlphaPhysicalSingle(axCon4, 1, 1, plotSCIMES=True)

        axCon1.set_ylabel(r"$\alpha$ (angular area)")

        axCon1.set_xlabel(r"CO cutoff (K)")
        axCon2.set_xlabel(r"CO cutoff (K)")
        axCon3.set_xlabel(r"CO cutoff (K)")
        axCon4.set_xlabel(r"CO cutoff (K)")

        #################

        self.labelContype(axList)
        self.labelSCIMES(axCon4)
        self.drawColorBar(axCon3)
        self.drawHorizontalLine(axList, 1.97)  # 1.97+/0.06
        self.drawHorizontalLine([axCon4], 1.97)  # 1.97+/0.06
        self.setCutoffTicks(axList+[axCon4])

        fig.tight_layout()
        plt.savefig("compareParaAlphaAreaPhysical.pdf", bbox_inches='tight')
        plt.savefig("compareParaAlphaAreaPhysical.png", bbox_inches='tight', dpi=300)





    def plotAreaAlphaBeq1(self):
        """
        #get and plot the alpha distribution
        :return:
        """

        fig = plt.figure(figsize=(24, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axCon1 = fig.add_subplot(1, 4, 1)
        axCon2 = fig.add_subplot(1, 4, 2, sharex=axCon1, sharey=axCon1)
        axCon3 = fig.add_subplot(1, 4, 3, sharex=axCon1, sharey=axCon1)
        axCon4 = fig.add_subplot(1, 4, 4, sharex=axCon1, sharey=axCon1)

        axList = [axCon1, axCon2, axCon3]

        self.PerFormAllSingleG2650(axCon1, self.plotAreaAlphaSingleBeq1, 1)
        self.PerFormAllSingleG2650(axCon2, self.plotAreaAlphaSingleBeq1, 2)
        self.PerFormAllSingleG2650(axCon3, self.plotAreaAlphaSingleBeq1, 3)
        self.plotAreaAlphaSingleBeq1(axCon4, 1, 1, plotSCIMES=True)

        axCon1.set_ylabel(r"$\alpha$ (angular area)")

        axCon1.set_xlabel(r"CO cutoff (K)")
        axCon2.set_xlabel(r"CO cutoff (K)")
        axCon3.set_xlabel(r"CO cutoff (K)")
        axCon4.set_xlabel(r"CO cutoff (K)")

        #################

        #self.labelContype(axList)

        ##########
        at1 = AnchoredText(self.con1Str+ " ($|b|>1^{\circ}$)" , loc=1, frameon=False)
        axCon1.add_artist(at1)

        at2 = AnchoredText(self.con2Str+ " ($|b|>1^{\circ}$)" , loc=1, frameon=False)
        axCon2.add_artist(at2)

        at3 = AnchoredText(self.con3Str+ " ($|b|>1^{\circ}$)" , loc=1, frameon=False)
        axCon3.add_artist(at3)

        at = AnchoredText(self.SCIMESStr+" ($|b|>1^{\circ}$)", loc=1, frameon=False)
        axCon4.add_artist(at)


        ###############

        #self.labelSCIMES(axCon4)
        self.drawColorBar(axCon3)
        self.drawHorizontalLine(axList, 1.97)  # 1.97+/0.06
        self.drawHorizontalLine([axCon4], 1.97)  # 1.97+/0.06
        self.setCutoffTicks(axList+[axCon4])

        fig.tight_layout()
        plt.savefig("compareParaAlphaAreaBeq1.pdf", bbox_inches='tight')
        plt.savefig("compareParaAlphaAreaBeq1.png", bbox_inches='tight', dpi=300)





    def drawHorizontalLine(self,axList,yValue):
        """

        :param axList:
        :param yValue:
        :return:
        """
        ####
        ####
        ####
        for eachAx in axList:
            eachAx.axhline(y=yValue , color='black', linestyle='--', lw=1, label="")


    def plotAreaAlphaSingleBeq1(self,ax,minPts,conType,plotSCIMES=False):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:
        """

        #alphaDendro, alphaDendroError = self.getAlphaList(tb16Den)
        #

        #get tb list

        if plotSCIMES:
            tbList,_ = self.getSCIMESTBList()

            plotColor="black"
        else:
            tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)

            plotColor= self.mMinPts.to_rgba(minPts)


        #in this step select only molecular clouds larger than 1 degree

        tbList= self.selectBeq1List( tbList )

        alphaArea, alphaAreaError = doDBSCAN.getAlphaList(tbList)

        ebDBSCAN = ax.errorbar(self.cutoffInKelvin , alphaArea , yerr=alphaAreaError , c= plotColor  , marker='^', linestyle="-",    capsize=3, elinewidth=1.0, lw=1,   markerfacecolor='none')

        ebDBSCAN[-1][0].set_linestyle(':')



    def selectBeq1List(self,TBList):

        newTBList=[]

        for eachTB in TBList:
            hgTB= eachTB[ abs(   eachTB["y_cen"] )>1 ]

            newTBList.append( hgTB )



        return newTBList



    def plotAreaAlphaSingle(self,ax,minPts,conType,plotSCIMES=False):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:
        """

        #alphaDendro, alphaDendroError = self.getAlphaList(tb16Den)
        #

        #get tb list

        if plotSCIMES:
            tbList,_ = self.getSCIMESTBList()

            plotColor="black"
        else:
            tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)

            plotColor= self.mMinPts.to_rgba(minPts)


        alphaArea, alphaAreaError = doDBSCAN.getAlphaList(tbList)

        ebDBSCAN = ax.errorbar(self.cutoffInKelvin , alphaArea , yerr=alphaAreaError , c= plotColor  , marker='^', linestyle="-",    capsize=3, elinewidth=1.0, lw=1,   markerfacecolor='none')

        ebDBSCAN[-1][0].set_linestyle(':')



    def plotAreaAlphaPhysicalSingle(self,ax,minPts,conType,plotSCIMES=False):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:
        """

        #alphaDendro, alphaDendroError = self.getAlphaList(tb16Den)
        #

        #get tb list

        if plotSCIMES:

            tbList,_ = self.getSCIMESTBList()
            plotColor="black"

        else:

            tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)
            plotColor= self.mMinPts.to_rgba(minPts)

        ########################################################


        alphaArea, alphaAreaError = self.getPhysicalAlphaList(tbList)
        ebDBSCAN = ax.errorbar(self.cutoffInKelvin , alphaArea , yerr=alphaAreaError , c= plotColor  , marker='^', linestyle="-",    capsize=3, elinewidth=1.0, lw=1,   markerfacecolor='none')
        ebDBSCAN[-1][0].set_linestyle(':')



    def getPhysicalAlphaList(self,TBList):
        """

        :param TBList:
        :return:
        """

        #####

        alphaList=[]
        alphaErrorList = []

        #physicalEdges = np.linspace(0, 100, 1000)  # square pc^2
        #physicalCenter = self.getEdgeCenter(physicalEdges)
        length= 1500*np.deg2rad(0.5/60)
        compoleteAreaPhysical =  length**2 * 4   #4 pixels

        for eachTB in TBList:
            realArea = doDBSCAN.getRealArea(eachTB)
           # binN, binEdges = np.histogram(realArea, bins=physicalEdges)

            # calculate alpha

            meanA, stdA = doDBSCAN.getAlphaWithMCMC(realArea, minArea=compoleteAreaPhysical, maxArea=None, physicalArea=True)

            alphaList.append( meanA )

            alphaErrorList.append( stdA )

        return alphaList, alphaErrorList




    def getMassAlphaList(self,TBList,completeMass= 5 ):
        """

        :param TBList:
        :return:
        """

        #####

        alphaList=[]
        alphaErrorList = []

        #physicalEdges = np.linspace(0, 100, 1000)  # square pc^2
        #physicalCenter = self.getEdgeCenter(physicalEdges)
 
        #compoleteMass=   5   #4 pixels

        for eachTB in TBList:
            massArray = doDBSCAN.getMass(eachTB)
           # binN, binEdges = np.histogram(realArea, bins=physicalEdges)

            # calculate alpha

            #meanA, stdA = doDBSCAN.getAlphaWithMCMC(massArray, minArea=completeMass, maxArea=None, physicalArea=True)
            meanA, stdA,completeMassList = doDBSCAN.getMassAlphaList(massArray  )

            alphaList.append( meanA )

            alphaErrorList.append( stdA )
            print completeMassList
        return alphaList, alphaErrorList




    def plotFluxAlpha(self):
        """
        #get and plot the alpha distribution
        :return:
        """

        fig = plt.figure(figsize=(24, 6))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 18, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axCon1 = fig.add_subplot(1, 4, 1)
        axCon2 = fig.add_subplot(1, 4, 2, sharex=axCon1, sharey=axCon1)
        axCon3 = fig.add_subplot(1, 4, 3, sharex=axCon1, sharey=axCon1)
        axCon4 = fig.add_subplot(1, 4, 4, sharex=axCon1, sharey=axCon1)

        axList = [axCon1, axCon2, axCon3]

        self.PerFormAllSingleG2650(axCon1, self.plotFluxAlphaSingle, 1)

        self.PerFormAllSingleG2650(axCon2, self.plotFluxAlphaSingle, 2)

        self.PerFormAllSingleG2650(axCon3, self.plotFluxAlphaSingle, 3)
        self.plotFluxAlphaSingle(axCon4,1,1,plotSCIMES=True)
        axCon1.set_ylabel(r"$\alpha$ (flux)")

        axCon1.set_xlabel(r"CO cutoff (K)")
        axCon2.set_xlabel(r"CO cutoff (K)")
        axCon3.set_xlabel(r"CO cutoff (K)")
        axCon4.set_xlabel(r"CO cutoff (K)")

        #################

        self.labelContype(axList)
        self.labelSCIMES(axCon4)
        self.drawColorBar(axCon3)

        self.drawHorizontalLine(axList, self.averageFluxAlpha ) #1.97+/0.06
        self.drawHorizontalLine([axCon4], self.averageFluxAlpha ) #1.97+/0.06

        self.setCutoffTicks(axList+[axCon4])

        fig.tight_layout()
        plt.savefig("compareParaAlphaFlux.pdf", bbox_inches='tight')
        plt.savefig("compareParaAlphaFlux.png", bbox_inches='tight', dpi=300)



    def plotFluxAlphaSingle(self,ax,minPts,conType, plotSCIMES=False ):
        """

        :param ax:
        :param minPts:
        :param conType:
        :return:
        """

        #alphaDendro, alphaDendroError = self.getAlphaList(tb16Den)
        #
        if plotSCIMES:
            tbList,_ = self.getSCIMESTBList()

            plotColor="black"

        else:
            tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)

            plotColor= self.mMinPts.to_rgba(minPts)

        #get tb list
        #tbList = self.getTBByCutOffList(self.cutoffList,minPts,conType)
        alphaFlux, alphaFluxError = doDBSCAN.getFluxAlphaList(tbList,self.cutoffList)

        ebDBSCAN = ax.errorbar(self.cutoffInKelvin , alphaFlux , yerr=alphaFluxError , c=  plotColor , marker='^', linestyle="-",    capsize=3, elinewidth=1.0, lw=1,   markerfacecolor='none')

        ebDBSCAN[-1][0].set_linestyle(':')




    def checkCriteria(self):
        """
        Test criteria to find if the criteria is reasonalble?
        :return:
        """

        ###

        fig = plt.figure(figsize=(20, 8))
        ax = fig.add_subplot(1, 2, 1) # previous
        ax2 = fig.add_subplot(1, 2, 2,sharex=ax,sharey=ax)

        # fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)

        overAllFontSize = 20
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': overAllFontSize, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        savePath = "/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/"

        tbTest = Table.read(savePath + "G2650LocaldbscanS7.0P3Con1.fit")

        #tbTest = tbTest[tbTest["area_exact"] < 1000 * 0.25]

        ax.scatter(tbTest["area_exact"] / 0.25, tbTest["allChannel"], s=3, c='b', label=r"7$\sigma$ cutoff (All clusters)  ")

        ax.set_ylabel(r"Velocity channels", fontsize=overAllFontSize)

        at = AnchoredText("All clusters", loc=5, frameon=False, prop={"color": "k", "size": overAllFontSize})
       # ax.add_artist(at)

        ax.set_xlim(-10,300)

        ax.set_ylim(-5,30)

        ##############################################################

        drawTB=self.selectTBByPeak(tbTest)

        ax2.scatter(drawTB["area_exact"] / 0.25, drawTB["allChannel"], s=3, c='b', label=r"7$\sigma$ cutoff (Clusters last submit)")

        ax2.set_xlabel(r"Area (pixels)", fontsize=overAllFontSize)


        at = AnchoredText("Clusters last submit", loc=5, frameon=False, prop={"color": "k", "size": overAllFontSize})
        #ax2.add_artist(at)
        ###################

        ax.tick_params(axis='both', which='major', labelsize=overAllFontSize)
        ax.tick_params(axis='both', which='minor', labelsize=overAllFontSize)
        ax2.tick_params(axis='both', which='major', labelsize=overAllFontSize)
        ax2.tick_params(axis='both', which='minor', labelsize=overAllFontSize)

        ax.axhline(y=3, color='r', linestyle='--', lw=1, label='3 velocity channels')

        ax.axvline(x=4, linewidth=1, color='g', linestyle='--', label="Beam size (4 pixels)")

        ax2.axhline(y=3, color='r', linestyle='--', lw=1, label='3 velocity channels')

        ax2.axvline(x=4, linewidth=1, color='g', linestyle='--', label="Beam size (4 pixels)")




        ax.legend()
        ax2.legend()
        fig.tight_layout()
        plt.savefig("AreaLineWidth.png", bbox_inches='tight', dpi=600)



    def checkCriteriaWithHist(self,old=False,cutOff=3.0):
        """
        :return:

        """
        from matplotlib.ticker import NullFormatter
        nullfmt = NullFormatter()
        fig = plt.figure(figsize=(10, 8))
        #ax = fig.add_subplot(1, 1, 1) # previous

        # fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)

        overAllFontSize = 16
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': overAllFontSize, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        # definitions for the axes
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left + width + 0.02

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.2]
        rect_histy = [left_h, bottom, 0.2, height]

        # start with a rectangular Figure

        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)

        # no labels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)

        ###################draw scatter



        savePath = "/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/"
        #tbTest = Table.read(savePath + "G2650LocaldbscanS7.0P3Con1.fit")
        tbTest = Table.read(savePath + "G2650LocaldbscanS{}P4Con1.fit".format(cutOff))
        #tbTest = Table.read(savePath + "G2650LocaldbscanS7.0P3Con1.fit")

        rawN = len( tbTest )

        #drawTB=self.selectTBByPeak(tbTest, minDelta=3,cutOff=2 , minChannel=1,pixN=8)



        if old:
            drawTB=self.selectTBOld(tbTest, minDelta=3,cutOff=2 , pixN=16)

            x = drawTB["area_exact"] / 0.25
            y = drawTB["allChannel"]

            peakSigma=int( cutOff+3 )

            axScatter.scatter(x,y , s=3, c='b', label=r"1.5 K cutoff, min\_npix$\geq$16, Peak$\geq$ 3 K" ) #.format( int(cutOff)  ) )
        else:
            #drawTB=tbTest #self.selectTBAreaChannel(tbTest )
            drawTB= self.selectTBAreaChannel(tbTest )

            x = drawTB["area_exact"] / 0.25
            y = drawTB["allChannel"]
            axScatter.scatter(x,y , s=3, c='b', label=r"1.5 K cutoff") #.format(int(cutOff)))



        ###make text
        axScatter.text(  100,1.5 ,"3 channels",color='red')
        axScatter.text(  -1.5,27 ,r"Beam size (2$\times$2 region) ",rotation=90,color="g")


        totalN=len(drawTB)


        axScatter.set_xlabel(r"Area (pixels)", fontsize=overAllFontSize)
        axScatter.set_ylabel(r"Velocity channels", fontsize=overAllFontSize)



        x0=-10
        x1=300
        y0=-5
        y1=30
        binwidth=1
        axScatter.set_xlim(-10,150)
        axScatter.set_ylim(-2,30)



        ##############

        binsx = np.arange( x0 +0.5, x1 + binwidth, binwidth)
        binsy = np.arange( y0+0.5 , y1 + binwidth, binwidth)
        axHistx.hist(x, bins=binsx)
        axHisty.hist(y, bins=binsy, orientation='horizontal')

        axHistx.set_xlim(axScatter.get_xlim())
        axHisty.set_ylim(axScatter.get_ylim())

        at = AnchoredText("{} of {} clusters".format( totalN , rawN ), loc=1, frameon=False, prop={"color": "k", "size": overAllFontSize})
        axHistx.add_artist(at)
        ############

        axScatter.axhline(y=3, color='r', linestyle='--', lw=1, label="")

        axScatter.axvline(x=4, linewidth=1, color='g', linestyle='--', label="")

        axHisty.axhline(y=3, color='r', linestyle='--', lw=1, label='')

        axHistx.axvline(x=4, linewidth=1, color='g', linestyle='--', label="")

        ###########################


        axScatter.legend(loc=4)
        fig.tight_layout()
        if old:
            plt.savefig("AreaLineWidthHistFigure1.png", bbox_inches='tight', dpi=600)
            plt.savefig("AreaLineWidthHistFigure1.pdf", bbox_inches='tight' )

        else:
            plt.savefig("AreaLineWidthHistFigure2.png", bbox_inches='tight', dpi=600)
            plt.savefig("AreaLineWidthHistFigure2.pdf", bbox_inches='tight' )


    def getProjection(self,labelFITS, id,savePath="/home/qzyan/WORK/myDownloads/MWISPcloud/areaCheck/" ):
        """
        #get the projection fits
        :param labelFITS:
        :param id:
        :param savePath:
        :return:
        """
        dataLabel,headLabel=myFITS.readFITS( labelFITS )

        cloudIndex = np.where(dataLabel==id)

        ZZ,YY,XX= cloudIndex

        zeroProjection = np.zeros_like(dataLabel[0])
        zeroProjection[tuple([YY,XX])] = 1



        fits.writeto(savePath+"CLoud{}Area.fits".format(id),zeroProjection,header=headLabel,overwrite=True)



    #

    def compareSCIMES(self,channelNumber=63):
        """
        This function is used to compare many cases of SCIMES parameter, to justify the correct use of parameters
        :return:
        """
        #need to draw fits

        # drawLrange,drawBrange=gaiaDis.box(38.9496172,0.1091115,3263.632 ,2991.628 ,4.9024796e-06)
        # region 2
        # drawLrange,drawBrange=gaiaDis.box(44.8346022,0.8519293,2361.857 ,2141.563 ,4.9024796e-06)
        # drawLrange, drawBrange = gaiaDis.box(42.8611667, 0.1834138, 3122.226, 2901.286, 4.9024796e-06)



        #drawLrange, drawBrange = gaiaDis.box(38.0855621, -1.5002922, 18000.000, 14400.000, 0)



        #drawLrange, drawBrange = gaiaDis.box(47.0438486, -0.9243936, 18000.000 ,14400.000 , 0)

        #drawLrange, drawBrange = gaiaDis.box(45.6405916, -1.5156708, 27000.000,21600.000, 0)
        drawLrange, drawBrange = gaiaDis.box(33.1774400,1.4085280,27000.000 ,21600.000 ,0) #region 2, index 83

        if 1:
            drawLrange, drawBrange = gaiaDis.box(31.8877408,0.9472144,36000.000 ,28800.000 ,0)  # region 2, index 83
            channelNumber = 85



        print drawLrange
        print drawBrange
        # vRange=[1.8,9.8] #km/s
        vRange = [-6, 6]  # km/s

        # first crop fits

        rawCO = "G2650Local30.fits"

        # rawCO="G2650CO12MaskedCO.fits"

        # labelDendroFITS = "G2650minV3minP16_TrunkAsignMask0.fits"
        # labelDBSCANFITS = "DBCLEAN3.0_16Label.fits"
        # labelSCIMESFITS = "./scimesG2650/ClusterAsgn_3_16Ve20.fits"
        fitsSavePath="/home/qzyan/WORK/dataDisk/MWISPcloud/"
        scimesDefaultFITS = fitsSavePath+"ClusterAsgn_scimesDefault_minV3P16.fits"


        scimesVe15 = fitsSavePath+"ClusterAsgn_3_16Ve15.fits"
        scimesVe17 = fitsSavePath+"ClusterAsgn_3_16Ve17.fits"
        scimesVe20 = fitsSavePath+"ClusterAsgn_3_16Ve20.fits"
        scimesVe23 = fitsSavePath+"ClusterAsgn_3_16Ve23.fits"
        #scimesVe25 = fitsSavePath+"ClusterAsgn_3_16Ve25.fits"
        scimesVe25 = fitsSavePath+  "ClusterAsgn_3_16Ve20AndVoLu.fits"

        

        drawLabelFITSList =[ scimesDefaultFITS, scimesVe17 ,   scimesVe20 ,  scimesVe23, scimesVe25  ]

        dataList,headList=[] , []


        for eachLabel in drawLabelFITSList:
            labelData,labelHead= myFITS.readFITS(eachLabel)

            dataList.append( labelData[channelNumber] )
            headList.append( labelHead )



        WCSCrop = WCS(headList[0])

        lIndexRange, bIndexRange = doDBSCAN.convertLBrangeToIndexRange(WCSCrop, drawLrange, drawBrange)

        # get velocity range

        _, _, channelV = WCSCrop.wcs_pix2world(0, 0, channelNumber, 0)
        channelV = channelV / 1000.

        overAllFontSize= 14

        fig = plt.figure(1, figsize=(15,8))
        rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'], "size": overAllFontSize })
        # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

        rc('text', usetex=True)
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axCO = pywcsgrid2.subplot(231, header=WCSCrop)

        dataCO, headCO = myFITS.readFITS(rawCO)

        intData = dataCO[channelNumber]  # np.sum(dataCO ,axis=0 ) *0.2

        cmapCO = plt.cm.bone

        cmapCO.set_bad('black')

        axCO.imshow(np.sqrt(intData), origin='lower', cmap=cmapCO, vmin=0, vmax=3, interpolation='none')

        at = AnchoredText(
            r"$^{{12}}\mathrm{{CO}}~(J=1\rightarrow0)$, $V_{{\rm LSR}}={}\ \rm km\ s^{{-1}}$".format(channelV), loc=3,
            frameon=False,
            prop={"color": "w", "size": 12})
        axCO.add_artist(at)


        #mark the cloud, which could tell which criteria is the best
        #circle(115.09103,653.12465,57.046886)
        circleDraw=plt.Circle((115.09103,653.12465), 57.046886, color='r', fill=False,lw=0.5)
        axCO.add_artist(circleDraw)
        axCO.set_ticklabel_type("absdeg", "absdeg")
        axCO.axis[:].major_ticks.set_color("w")
        axCO.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
        axCO.set_ylabel(r"Galactic Latitude ($^{\circ}$)")




        ###dendrogram
        ax2 = pywcsgrid2.subplot(232, header=WCSCrop, sharex=axCO, sharey=axCO)
        ax3 = pywcsgrid2.subplot(233, header=WCSCrop, sharex=axCO, sharey=axCO)
        ax4 = pywcsgrid2.subplot(234, header=WCSCrop, sharex=axCO, sharey=axCO)
        ax5 = pywcsgrid2.subplot(235, header=WCSCrop, sharex=axCO, sharey=axCO)
        ax6 = pywcsgrid2.subplot(236, header=WCSCrop, sharex=axCO, sharey=axCO)


        axList=[ax2, ax3, ax4, ax5, ax6]

        #markStr=[ "Volume and luminosity criteria",r"$\sigma_s$ = {} km s$^{{-1}}$".format(17), r"$\sigma_s$ = {} km s$^{{-1}}$".format(20),\
                  #r"$\sigma_s$ = {} km s$^{{-1}}$".format(23),r"$\sigma_s$ = {} km s$^{{-1}}$".format(25) ]


        markStr=[ "Volume and luminosity",r"$\sigma_s$ = {} km s$^{{-1}}$".format(17), r"$\sigma_s$ = {} km s$^{{-1}}$".format(20),\
                  r"$\sigma_s$ = {} km s$^{{-1}}$".format(23), r"Volume, luminosity, and $\sigma_s$ = {} km s$^{{-1}}$".format(20)  ]

        for i in range( len(axList) ):


            axSCIMES=axList[i]
            doDBSCAN.showLabels(axSCIMES,  dataList[i]  )


            at = AnchoredText(markStr[i], loc=3, frameon=False, prop={"color": "w", "size": overAllFontSize})
            axSCIMES.add_artist(at)
            axSCIMES.set_ticklabel_type("absdeg", "absdeg")
            axSCIMES.axis[:].major_ticks.set_color("w")

            axSCIMES.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
            axSCIMES.set_ylabel(r"Galactic Latitude ($^{\circ}$)")



        axCO.set_xlim(lIndexRange)
        axCO.set_ylim(bIndexRange)


        import matplotlib.patheffects as path_effects
        text = axCO["gal"].text(46.263, -01.3303, "LDN 673", color='Red', alpha=0.8, fontsize=overAllFontSize,  horizontalalignment='center', verticalalignment='center')
        text.set_path_effects([path_effects.Stroke(linewidth=0.8, foreground="black"), path_effects.Normal()])
        fig.tight_layout()

        #plt.savefig("compareSCIMESPara{}.pdf".format(channelNumber), bbox_inches="tight")
        #plt.savefig("compareSCIMESPara{}.png".format(channelNumber), bbox_inches="tight", dpi=600)

        plt.savefig("compareSCIMESPara.pdf" , bbox_inches="tight")
        plt.savefig("compareSCIMESPara.png" , bbox_inches="tight", dpi=600)


    def getChannelDataHeadByVel(self,fitsFile,  testVelocity ):
        """

        :param fitsFile:
        :param Vrange:
        :return:
        """



        tmpExtractFITS= self.tmpPath+"ChannelExtract.fits"

        doFITS.cropFITS(fitsFile, outFITS= tmpExtractFITS , Vrange=[testVelocity,  testVelocity ], overWrite=True)

        data,head=myFITS.readFITS( tmpExtractFITS )
        if len(data.shape) == 3:
            data=data[0]

        return data,head

    def drawHDBSCAN(self,testVelocity=10):
        """

        Test the result of HDBSCAN
        #HDBSCAN,
        :return:
        """

        #drawLrange, drawBrange = gaiaDis.box(45.6405916, -1.5156708, 27000.000,21600.000, 0)
        #
        drawLrange, drawBrange = gaiaDis.box(37.5000000,0.0000000,18000.000 ,14400.000 ,0)



        print drawLrange
        print drawBrange
        # vRange=[1.8,9.8] #km/s
        vRange = [-6, 6]  # km/s

        # first crop fits

        rawCO = "G2650Local30.fits"


        DBSCANFITS="/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/G2650LocaldbscanS3.0P4Con1_Clean.fits"

        HDBSCANFITS="/home/qzyan/WORK/myDownloads/MWISPcloud/testHdbscan/NormalTesthdbscanProduceLocalS3P16_Clean.fits"

        #first moment according to the veloictyh
        dataCO,headCO=self.getChannelDataHeadByVel(rawCO, testVelocity )
        labelDBSCAN,headCO=self.getChannelDataHeadByVel(DBSCANFITS, testVelocity )
        labelHDBSCAN,headCO=self.getChannelDataHeadByVel(HDBSCANFITS, testVelocity )



        dataList=[labelHDBSCAN, labelDBSCAN  ]

        WCSCrop = WCS( headCO )

        lIndexRange, bIndexRange = doDBSCAN.convertLBrangeToIndexRange(WCSCrop, drawLrange, drawBrange)

        # get velocity range

        #_, _, channelV = WCSCrop.wcs_pix2world(0, 0, channelNumber, 0)
        channelV = testVelocity

        fig = plt.figure(1, figsize=(15 , 4.5))
        rc('font', **{'family': 'sans-serif', 'serif': ['Helvetica'], "size": 15})
        # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

        rc('text', usetex=True)
        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        axCO = pywcsgrid2.subplot(131, header=WCSCrop)




        cmapCO = plt.cm.bone

        cmapCO.set_bad('black')

        axCO.imshow(np.sqrt(dataCO), origin='lower', cmap=cmapCO, vmin=0, vmax=3, interpolation='none')

        at = AnchoredText(
            r"$^{{12}}\mathrm{{CO}}~(J=1\rightarrow0)$, $V_{{\rm LSR}}={}\ \rm km\ s^{{-1}}$".format(channelV), loc=3,
            frameon=False,
            prop={"color": "w", "size": 12})
        axCO.add_artist(at)


        #mark the cloud, which could tell which criteria is the best
        #circle(115.09103,653.12465,57.046886)
        #circleDraw=plt.Circle((115.09103,653.12465), 57.046886, color='r', fill=False,lw=0.5)
        #axCO.add_artist(circleDraw)
        axCO.set_ticklabel_type("absdeg", "absdeg")
        axCO.axis[:].major_ticks.set_color("w")
        axCO.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
        axCO.set_ylabel(r"Galactic Latitude ($^{\circ}$)")


        ###dendrogram
        ax2 = pywcsgrid2.subplot(132, header=WCSCrop, sharex=axCO, sharey=axCO)
        ax3 = pywcsgrid2.subplot(133, header=WCSCrop, sharex=axCO, sharey=axCO)

        axList=[ax2, ax3 ]

        markStr=[ "HDBSCAN (3$\sigma$ cutoff)" , r"DBSCAN (3$\sigma$ cutoff, MinPts=4, connectivity=1)"  ]

        for i in range( len(axList) ):


            axSCIMES=axList[i]
            doDBSCAN.showLabels(axSCIMES,  dataList[i]  )

            at = AnchoredText(markStr[i], loc=3, frameon=False, prop={"color": "w", "size": 12})
            axSCIMES.add_artist(at)
            axSCIMES.set_ticklabel_type("absdeg", "absdeg")
            axSCIMES.axis[:].major_ticks.set_color("w")

            axSCIMES.set_xlabel(r"Galactic Longitude ($^{\circ}$)")
            axSCIMES.set_ylabel(r"Galactic Latitude ($^{\circ}$)")


        axCO.set_xlim(lIndexRange)
        axCO.set_ylim(bIndexRange)


        fig.tight_layout()

        #plt.savefig("compareSCIMESPara{}.pdf".format(channelNumber), bbox_inches="tight")
        #plt.savefig("compareSCIMESPara{}.png".format(channelNumber), bbox_inches="tight", dpi=600)

        plt.savefig("compareHDBSCAN.pdf" , bbox_inches="tight")
        plt.savefig("compareHDBSCAN.png" , bbox_inches="tight", dpi=600)



    def getSCIMESTBList(self):

        # dendroSigmaList=[2,2.5 , 3, 3.5, 4,4.5,5, 5.5, 6,6.5,7]

        dendroSigmaList = [ 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7 ]
        path = "./scimesG2650/"

        #
        tbList=[]
        for sigmas in dendroSigmaList:

            tb16File = path + "ClusterAsgn_{}_16Ve20_CleanWithLW.fit".format(sigmas )
            #print tb16File
            #print os.path.isfile(tb16File)

            #print len(Table.read( tb16File ) )
            tbList.append( Table.read( tb16File ) )



        return tbList, dendroSigmaList




    def getSCIMESFITSList(self):

        # dendroSigmaList=[2,2.5 , 3, 3.5, 4,4.5,5, 5.5, 6,6.5,7]

        dendroSigmaList = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7]
        path = "./scimesG2650/"

        #
        fitsList = []
        for sigmas in dendroSigmaList:

            fits16File = path + "ClusterAsgn_{}_16Ve20.fits".format(sigmas )
            #print os.path.isfile(fits16File)

            fitsList.append( fits16File )

        return fitsList



    def produceSCIMESCat(self,fitsList):
        """

        :param fitsList:
        :return:
        """


        for eachFITS in fitsList:

            saveMarker=    eachFITS[:-5]
            rawDBSCANTBFile = doDBSCAN.getCatFromLabelArray(self.G2650Local, eachFITS, doDBSCAN.TBModel,  saveMarker= saveMarker )
            doAllDBSCAN.getEquivalentLineWidth(eachFITS, rawDBSCANTBFile, reProduceFITS=False, keepFITS=True)



    def pipLineSCIMES(self):
        """
        get scimes fits and clean
        :return:
        """

        #get fits list

        fitsList= self.getSCIMESFITSList(  )
        #self.produceSCIMESCat( fitsList )

    def printTBInfo(self,eachTB):
        """

        :param eachTB:
        :return:
        """

        #print np.min( eachTB["area_exact"] ), np.max(  eachTB["area_exact"] )/3600.

        #print np.min( eachTB["area_exact"] ), np.max(  eachTB["area_exact"] )/3600.

        print  np.max(  eachTB["sum"] )/np.sum(  eachTB["sum"]  )





    def getPeakAndVrmsRange(self):

        """

        :return:
        """

        #
        for eachP in self.con1G2650:

            tbList1 = self.getTBByCutOffList(self.cutoffList, eachP, 1, getCleanLW=True)
            for eachTB in  tbList1:
 
                self.printTBInfo(eachTB)
                #print "Velocity range: ",  np.min(eachTB["lineWidth"]), np.max(eachTB["lineWidth"])
                ##print "Peak range: ",  np.min(eachTB["peak"]), np.max(eachTB["peak"])

        for eachP in self.con2G2650:

            tbList1 = self.getTBByCutOffList(self.cutoffList, eachP, 2, getCleanLW=True)
            for eachTB in  tbList1:
                self.printTBInfo(eachTB)
                #print "Velocity range: ",  np.min(eachTB["lineWidth"]), np.max(eachTB["lineWidth"])
                #print "Peak range: ",  np.min(eachTB["peak"]), np.max(eachTB["peak"])


        for eachP in self.con3G2650:

            tbList1 = self.getTBByCutOffList(self.cutoffList, eachP, 3, getCleanLW=True)
            for eachTB in  tbList1:
                self.printTBInfo(eachTB)
                #print "Velocity range: ",  np.min(eachTB["lineWidth"]), np.max(eachTB["lineWidth"])
                #print "Peak range: ",  np.min(eachTB["peak"]), np.max(eachTB["peak"])


        #tbList2 = self.getTBByCutOffList(self.cutoffList, minPts, 2)
        #tbList3 = self.getTBByCutOffList(self.cutoffList, minPts, 3)


    def getAverageAreaAlpha(self):
        """

        :return:
        """
        ####
        ####

        finalListAlpha= []
        finalListAlphaError = []

        for eachP in self.con1G2650: ##########

            tbList1 = self.getTBByCutOffList(self.cutoffList, eachP, 1 )
            alphaArea, alphaAreaError = doDBSCAN.getAlphaList(tbList1)
            finalListAlpha=finalListAlpha+ alphaArea
            finalListAlphaError = finalListAlphaError +  alphaAreaError


        for eachP in self.con2G2650: ##########

            tbList1 = self.getTBByCutOffList(self.cutoffList, eachP, 2 )
            alphaArea, alphaAreaError = doDBSCAN.getAlphaList(tbList1)
            finalListAlpha=finalListAlpha+ alphaArea
            finalListAlphaError = finalListAlphaError +  alphaAreaError

        for eachP in self.con3G2650: ##########
            tbList1 = self.getTBByCutOffList(self.cutoffList, eachP, 3 )

            alphaArea, alphaAreaError = doDBSCAN.getAlphaList(tbList1)

            finalListAlpha=finalListAlpha+ alphaArea
            finalListAlphaError = finalListAlphaError +  alphaAreaError

        #calculate  average finalListAlpha, and finalListAlphaError
        ################
        print "Average error, ", np.mean(finalListAlphaError)

        errorAlpha = np.mean(finalListAlphaError) ** 2 + np.std(finalListAlpha, ddof=1) ** 2
        errorAlpha = np.sqrt( errorAlpha )

        alphaMean  = np.mean( finalListAlpha )

        print "The mean area alpha is {:.2f}, area error is {:.2f}".format(alphaMean, errorAlpha)





    def getAverageFluxAlpha(self):
        """

        :return:
        """
        ########
        ########



        finalListAlpha= []
        finalListAlphaError = []

        for eachP in self.con1G2650: ##########

            tbList1 = self.getTBByCutOffList(self.cutoffList, eachP, 1 )
            alphaArea, alphaAreaError = doDBSCAN.getFluxAlphaList(tbList1, self.cutoffList)
            finalListAlpha=finalListAlpha+ alphaArea
            finalListAlphaError = finalListAlphaError +  alphaAreaError


        for eachP in self.con2G2650: ##########

            tbList1 = self.getTBByCutOffList(self.cutoffList, eachP, 2 )
            alphaArea, alphaAreaError = doDBSCAN.getFluxAlphaList(tbList1, self.cutoffList)
            finalListAlpha=finalListAlpha+ alphaArea
            finalListAlphaError = finalListAlphaError +  alphaAreaError

        for eachP in self.con3G2650: ##########
            tbList1 = self.getTBByCutOffList(self.cutoffList, eachP, 3 )

            alphaArea, alphaAreaError = doDBSCAN.getFluxAlphaList(tbList1, self.cutoffList)

            finalListAlpha=finalListAlpha+ alphaArea
            finalListAlphaError = finalListAlphaError +  alphaAreaError

        #calculate  average finalListAlpha, and finalListAlphaError
        ################
        print "Average flux alpha error, ", np.mean(finalListAlphaError)

        errorAlpha = np.mean(finalListAlphaError) ** 2 + np.std(finalListAlpha, ddof=1) ** 2
        errorAlpha = np.sqrt( errorAlpha )

        alphaMean  = np.mean( finalListAlpha )

        print "The mean flux alpha is {:.2f}, error is {:.2f}".format(alphaMean, errorAlpha)

    def getPublicCatalog(self, catalogPath, rawCatalog,savePath=None, saveTag= "HarvardDatabase_" ):
        """

        :param TBFileName:
        :return:
        """

        # tb=Table.read(TBFileName)
        # better use
        #rawCatalog,
        #rawCatalog = "minV2minP8dendroMannualCat_LineWidth.fit"
        # rawCatalog="minV2minP8_dendroCatTrunk.fit"

        if savePath==None:

            saveFileName = catalogPath+ saveTag + rawCatalog #  "cloudCatalogDendrgram.fit"
        else:
            saveFileName = savePath+ saveTag + rawCatalog #  "cloudCatalogDendrgram.fit"


        tbRaw = Table.read(catalogPath+rawCatalog)

        #print len(tbRaw)
        ##tbRaw = self.removeWrongEdges(tbRaw)

        #print len(tbRaw)

        publishCat = Table()
        tbRaw["name"] = tbRaw["_idx"].astype(str)
        # add cloudName to the row
        # print tbRaw
        for eachR in tbRaw:
            l = np.float(eachR["x_cen"])
            b = np.float(eachR["y_cen"])
            v = np.float(eachR["v_cen"])  # to km/s
            # print l,b,v

            cloudName = doDBSCAN.getCloudNameByLB(l, b)

            eachR["name"] = cloudName

        #print tbRaw.colnames
        # addName
        publishCat.add_column(tbRaw["name"])

        # add l, b, v, values
        colL = tbRaw["x_cen"]
        colL.name = "l"
        colL.unit = u.deg
        publishCat.add_column(colL)

        colb = tbRaw["y_cen"]
        colb.name = "b"
        colb.unit = u.deg
        publishCat.add_column(colb)

        colv = tbRaw["v_cen"]
        colv.name = "Vlsr"
        colv.unit = u.km / u.s
        publishCat.add_column(colv)

        #######################################

        # lbv, dispersion

        colLrms = tbRaw["l_rms"]
        colLrms.name = "l_sigma"
        colLrms.unit = u.deg
        publishCat.add_column(colLrms)

        colBrms = tbRaw["b_rms"]
        colBrms.name = "b_sigma"
        colBrms.unit = u.deg
        publishCat.add_column(colBrms)


        colVVrms = tbRaw["v_rms"]
        colVVrms.name = "v_sigma"
        colVVrms.unit = u.km / u.s
        publishCat.add_column(colVVrms)



        colVrms = tbRaw["lineWidth"]
        # colBrms.name = "b_sigma"
        colVrms.unit = u.km / u.s
        publishCat.add_column(colVrms)

        colVoxN = tbRaw["pixN"].astype(int)
        colVoxN.name = "voxelN"
        colVoxN.unit = None
        publishCat.add_column(colVoxN)

        colPeak = tbRaw["peak"]
        # colBrms.name = "b_sigma"
        colPeak.unit = u.K
        publishCat.add_column(colPeak)

        ##############
        colArea = tbRaw["area_exact"]
        # colBrms.name = "b_sigma"
        # colArea.unit = u.K
        publishCat.add_column(colArea)
        ##############
        colChannel = tbRaw["allChannel"]
        colChannel.name = "channelNumber"
        colChannel.unit=None
        publishCat.add_column(colChannel)


        ##############
        colBeam = tbRaw["has22"]
        colBeam.name = "hasBeam"
        colBeam.unit=None
        publishCat.add_column(colBeam)


        #########################################

        colFlux = tbRaw["sum"] * 0.2 * 0.25
        colFlux.name = "flux"
        colFlux.unit = u.K * u.km / u.s * u.arcmin ** 2  # *u.def_unit("Omega A")
        publishCat.add_column(colFlux)

        publishCat.write(saveFileName, overwrite=True)




    def getAllPublishCatalog(self):
        
        publishCatalogPath="/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/publishCatalog/"
        catalogPath="/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/"

        searchStr =catalogPath + "G2650*S*P*Con*CleanWithLW.fit"

        if 1:#DBSCAN catalog
            allDBSCANCatalog = glob.glob(searchStr)


            for eachCat in allDBSCANCatalog:



                _,catname= os.path.split(eachCat)
                
                if "P3Con1" in catname:
                    continue
                if "P3Con2" in catname:
                    continue
                if "P4Con2" in catname:
                    continue
                if "P5Con2" in catname:
                    continue
                if "P6Con2" in catname:
                    continue
                if "P7Con2" in catname:
                    continue
                ####
                    
                if "P3Con3" in catname:
                    continue
                if "P4Con3" in catname:
                    continue
                if "P5Con3" in catname:
                    continue
                if "P6Con3" in catname:
                    continue
                if "P7Con3" in catname:
                    continue
                if "P8Con3" in catname:
                    continue
                if "P9Con3" in catname:
                    continue
                if "P10Con3" in catname:
                    continue

                    
                    
                self.getPublicCatalog(  catalogPath, catname,savePath=publishCatalogPath)

        if 1: ##### produce SCIMES catalog
            pass
            sciMESCatPath="/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/"
            searchStr =sciMESCatPath + "C*CleanWithLW.fit"

            allSCIMES = glob.glob(searchStr)


            for eachCat in allSCIMES:
                _,catname= os.path.split(eachCat)
                self.getPublicCatalog(  sciMESCatPath, catname,savePath=publishCatalogPath,saveTag="HarvardDatabase_SCIMES_")


        #first get all lw catalog


    def UMMCDBSCANPipeline(self):
        """

        #process Ursa Major data with DBSCAN

        :return:
        """
        #savePath
        fitsPath =  "/home/qzyan/WORK/myDownloads/MWISPcloud/UMMCTest/"
        rawUMMCCO12="/home/qzyan/WORK/projects/NewUrsaMajorPaper/OriginalFITS/UMMC_CO12_Crop.fits"

        preFix="UMMC_DBSCAN_dbscan"

        # get ALl table

        if 1:#produce all Table

            #first get all fits

            allUMMCfits = glob.glob( fitsPath+preFix+"S3*.fits" )
            for eachFITS in allUMMCfits:
                _,fitsName=os.path.split( eachFITS )
                print fitsName
                #if fitsName !="UMMC_DBSCAN_dbscanS2P4Con1.fits":

                    #continue
                self.getCatalog(rawUMMCCO12, fitsPath, fitsName ,  eachFITS[:-5]  )


    def checkTheAngularAreaPowerLaw(self,TBName):
        """

        :return:
        """
        ##
        colAreaAccurate="area_accurate"
        tb=Table.read(TBName)
        print len(tb)
        self.rmsCO12=0.16
        pureTB=self.selectTBFormal(tb,cutOff=3, removeEdge=False, minDelta= 3 )

        print pureTB.colnames
        #pureTB=pureTB[pureTB["x_cen"]>145]
        #pureTB=pureTB[pureTB["y_cen"]<40]


        print len( pureTB ),"?????????????????Total number of molecular clouds"

        compliteUMMC= 0.001  #0.25*4/3600.  #0.003  #
        #compliteUMMC= 0.25*4/3600.  #0.003  #

        alphaArea, alphaAreaError = doDBSCAN.getAlphaList( [pureTB]  ,minArea= compliteUMMC, colName= colAreaAccurate )

        #draw the power law
        areaEdges = np.linspace(0.25 / 3600., 10, 10000)
        areaCenter = doDBSCAN.getEdgeCenter(areaEdges)

        drawTB= pureTB
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(1, 1, 1)
        # fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 13, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]

        binN, binEdges = np.histogram(drawTB[colAreaAccurate] / 3600., bins=areaEdges)

        ax.plot(areaCenter[binN > 0], binN[binN > 0], color= "black", linestyle='-', markersize=1, lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, label="Ursa Major" )


        ######### draw first quadrant as comparison

        TBLocal=Table.read("/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/G2650LocaldbscanS3.0P4Con1_CleanWithLW.fit")
        binN, binEdges = np.histogram(TBLocal["area_exact"] / 3600., bins=areaEdges)
        ax.plot(areaCenter[binN > 0], binN[binN > 0], color= "black", linestyle='--', markersize=1,     lw=0.6, marker='o', markeredgewidth=0.5, fillstyle='full', alpha=0.6, label="G2650" )


        ax.plot( [compliteUMMC,compliteUMMC],[1,1000]  )
        ax.legend()
        ax.set_yscale('log')
        ax.set_xscale('log')
        plt.savefig("ummcAngularArea.png", dpi=600, bbox_inches='tight')



    def cleanFITS(self,savePath, prefix ):
        """
        GetFITS and TB, make selection and clean the fits
        :param savePath:
        :param prefix:
        :return:
        """

        fitsList= allUMMCfits = glob.glob( fitsPath+prefix+"*.fits" )

        #getMaskByLabel(self, FITSfile, tbFile, onlyClean=False, onlySelect=False, inputCutOff=None, reProduceFITS=True):



    def cloudNumberGrid(self,TBFile,dl=1,db=1):

        """
        find the number cloud per square degree in l-b-grid
        :param TBFile:
        :return:
        """
        ######


        TB=Table.read(TBFile)

        lGrid=np.arange(26,50+dl,dl)
        bGrid=np.arange(-5,5+db, db )

        lCell=zip(lGrid[0:-1]  ,  lGrid[1: ]   )
        bCell=zip(bGrid[0:-1]  ,  bGrid[1: ]   )

        ############
        ##
        NList=[]
        Lcenter=[]
        Bcenter=[]


        for eachLC in lCell:
            for eachBC in bCell:
                ######

                l0,l1=eachLC
                b0,b1=eachBC



                select1 = np.logical_and(TB["x_cen"] >= l0, TB["x_cen"] < l1 )

                select2 = np.logical_and(TB["y_cen"] >= b0, TB["y_cen"] < b1 )

                selectFinal=np.logical_and( select1,select2 )


                subTB=TB[selectFinal]
                NList.append(len(subTB))

                Lcenter.append( (l0+l1)/2. )
                Bcenter.append( (b0 + b1 )/2. )

                #

        return NList,Lcenter,Bcenter


    def drawNVariation(self):
        """
        :return:

        """

        checkTBFile="/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/G2650LocaldbscanS2P4Con1_CleanWithLW.fit"

        N,L,B=self.cloudNumberGrid( checkTBFile )


        print N

        fig = plt.figure(figsize=(12, 6))
        axL = fig.add_subplot(1, 2, 1 )
        axB = fig.add_subplot(1, 2, 2 ,sharey=axL)

        # fig, axs = plt.subplots(nrows=1, ncols=2,  figsize=(12,6),sharex=True)
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 16, 'serif': ['Helvetica']})

        mpl.rcParams['text.latex.preamble'] = [
            r'\usepackage{tgheros}',  # helvetica font
            r'\usepackage{sansmath}',  # math-font matching  helvetica
            r'\sansmath'  # actually tell tex to use it!
            r'\usepackage{siunitx}',  # micro symbols
            r'\sisetup{detect-all}',  # force siunitx to use the fonts
        ]


        axL.scatter(L, N ,s=10 ,c='b')

        axL.set_xlabel(r"$l$")
        axL.set_ylabel(r"Number per square degree")

        axB.scatter(B,N ,s=10,c='b')
        axB.set_xlabel(r"$b$")

        plt.savefig("numberDensityLB.png", bbox_inches='tight',dpi=300)



    def ZZZ(self):
            pass



################


doAllDBSCAN=allDBSCAN()

spectralSavePath = "/home/qzyan/WORK/myDownloads/MWISPcloud/spectraSave/"


if 0: # G13015
    rawFITS="/home/qzyan/WORK/dataDisk/G140/G130150merge12.fits"
    #

    #doAllDBSCAN.pureDBSCAN(rawFITS, 3 , MinPts=4, saveTag="G130150", connectivity=1, inputRMS=None, redo=True, keepFITSFile=True)

    #doAllDBSCAN.getMaskByLabel(doAllDBSCAN.fitsPath+"G130150dbscanS3P4Con1.fits", doAllDBSCAN.fitsPath+"G130150dbscanS3P4Con1.fit",  onlyClean=True )

    TB=doAllDBSCAN.getMaskByLabel(doAllDBSCAN.fitsPath+"G130150dbscanS3P4Con1.fits", doAllDBSCAN.fitsPath+"G130150dbscanS3P4Con1.fit", onlySelect=True )

    TB.write(doAllDBSCAN.fitsPath+"G130150dbscanS3P4Con1_clean.fit", overwrite=True)

    sys.exit()
####################################################################################
if 0:#
    #doAllDBSCAN.checkCriteria()
    doAllDBSCAN.checkCriteriaWithHist( old = True )
    doAllDBSCAN.checkCriteriaWithHist(  )

    sys.exit()




#Average flux alpha error,  0.02146005509641873
#The mean flux alpha is 1.79, error is 0.03
if 0: #test get TBList #

    pass
    #doAllDBSCAN.plotNumber()
    #doAllDBSCAN.plotPeak()

    #doAllDBSCAN.plotLineWidth()

    #doAllDBSCAN.plotArea()
    #doAllDBSCAN.plotFlux()

    #doAllDBSCAN.plotTotalFlux()

    #doAllDBSCAN.plotAreaAlpha()

    #doAllDBSCAN.plotAreaPhysical()
    #doAllDBSCAN.plotAreaAlphaPhysical()


    #doAllDBSCAN.plotFluxAlpha()

    #doAllDBSCAN.plotAreaAlphaBeq1()

    #doDBSCAN.drawCheckCloudsOneChannel()
    #doAllDBSCAN.compareSCIMES(channelNumber= 83 )

    #doAllDBSCAN.checkCriteriaWithHist( old = True )
    #doAllDBSCAN.checkCriteriaWithHist(  )

    #doAllDBSCAN.drawHDBSCAN()

    ####################################
    #doAllDBSCAN.plotMass()

    #doAllDBSCAN.plotMassAlpha()

    #doAllDBSCAN.drawNVariation()


    sys.exit()





if 0:
    #doAllDBSCAN.saveTag = "testAll"
    doAllDBSCAN.saveTag = "negativeTest"
    doAllDBSCAN.drawNegativeTest( )
    sys.exit()









if 0:#test area_exact and area_accurate
    rawG2650FITS="G2650Local30.fits"
    dbSCANLabelPath="/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/"

    fitsName="G2650LocaldbscanS2P5Con2.fits"
    doAllDBSCAN.getCatalog(rawG2650FITS, dbSCANLabelPath, fitsName, dbSCANLabelPath + fitsName[:-5])

    sys.exit()


if 0: #check power-law UMMC

    #doAllDBSCAN.UMMCDBSCANPipeline()

    TbFile="/home/qzyan/WORK/myDownloads/MWISPcloud/UMMCTest/UMMC_DBSCAN_dbscanS3P4Con1.fit"
    doAllDBSCAN.checkTheAngularAreaPowerLaw( TbFile )

    #doAllDBSCAN.getMaskByLabel("/home/qzyan/WORK/myDownloads/MWISPcloud/UMMCTest/UMMC_DBSCAN_dbscanS3P4Con1.fits","/home/qzyan/WORK/myDownloads/MWISPcloud/UMMCTest/UMMC_DBSCAN_dbscanS3P4Con1.fit", onlyClean=True,removeEdge=False )


    sys.exit()
##################











#########################

if 0:
    

    doAllDBSCAN.getAllPublishCatalog()
    sys.exit()



if 0: #get clean fits of HDBSCAN
    #

    hdbscanLabel="/home/qzyan/WORK/myDownloads/MWISPcloud/testHdbscan/NormalTesthdbscanProduceLocalS3P16.fits"
    tbFile="/home/qzyan/WORK/myDownloads/MWISPcloud/testHdbscan/NorMalTestCatalog.fit"
    doAllDBSCAN.getMaskByLabel( hdbscanLabel ,  tbFile ,  onlyClean=True, inputCutOff= 3 )
    sys.exit()


if 0: #get clean fits of HDBSCAN
    #

    hdbscanLabel="/home/qzyan/WORK/myDownloads/MWISPcloud/testHdbscan/ExtendTesthdbscanProduceLocalS3P16.fits"
    tbFile="/home/qzyan/WORK/myDownloads/MWISPcloud/testHdbscan/ExtendTestCatalog.fit"
    doAllDBSCAN.getMaskByLabel( hdbscanLabel ,  tbFile ,  onlyClean=True, inputCutOff= 3 )
    sys.exit()



if 0:
    doAllDBSCAN.drawHDBSCAN()
    sys.exit()




if 0:
    #doAllDBSCAN.getAverageFluxAlpha()
    #doAllDBSCAN.getAverageAreaAlpha()
    doAllDBSCAN.getPeakAndVrmsRange()
    sys.exit()



if 0:# out put spectral save, and draw examples

    ##
    spectralLabelFile="/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/G2650LocaldbscanS2P4Con1.fits"
    spectralTBFile="/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/G2650LocaldbscanS2P4Con1.fit"

    #doAllDBSCAN.getEquivalentLineWidth(  spectralLabelFile, spectralTBFile, saveSpectral=True, saveSpectralPath=spectralSavePath)

    doDBSCAN.getExamineSpectral( savePath=spectralSavePath )
    sys.exit()














if 0: #produce  produce masked fits, used to
    #doAllDBSCAN.getMaskByLabel(doAllDBSCAN.fitsPath+"G2650LocaldbscanS3P4Con1.fits", doAllDBSCAN.fitsPath+"G2650LocaldbscanS3P4Con1.fit",  onlyClean=False )
    
    scimesPath= "/home/qzyan/WORK/myDownloads/MWISPcloud/scimesG2650/"

    doAllDBSCAN.getMaskByLabel( scimesPath + "ClusterAsgn_2_16Ve20.fits", scimesPath  + "ClusterAsgn_2_16Ve20.fit",  onlyClean=False ,inputCutOff =2 )

    #doAllDBSCAN.getMaskByLabel(doAllDBSCAN.fitsPath+"G2650LocaldbscanS2P4Con1.fits", doAllDBSCAN.fitsPath+"G2650LocaldbscanS2P4Con1.fit",  onlyClean=False )
    #doAllDBSCAN.getMaskByLabel(doAllDBSCAN.fitsPath+"G2650LocaldbscanS2P8Con2.fits", doAllDBSCAN.fitsPath+"G2650LocaldbscanS2P8Con2.fit",  onlyClean=False )
    #doAllDBSCAN.getMaskByLabel(doAllDBSCAN.fitsPath+"G2650LocaldbscanS2P11Con3.fits", doAllDBSCAN.fitsPath+"G2650LocaldbscanS2P11Con3.fit",  onlyClean=False )

    sys.exit()



if 0: #process SCIMES
    doAllDBSCAN.getSCIMESTBList()

    #doAllDBSCAN.pipLineSCIMES()
    sys.exit()








if  0:  #test getting line width

    #calculate lineWidth, only for formal selected clouds, which should not too many
    for eachCutOff in  doAllDBSCAN.cutoffList:
        if eachCutOff==2:
            continue
        doAllDBSCAN.calLineWidthPipeLine( eachCutOff,reProduceFITS=True,keepFITS=False )
    sys.exit()








if 0: #produce  produce masked fits, used to
    doAllDBSCAN.getMaskByLabel(doAllDBSCAN.fitsPath+"G2650LocaldbscanS3.0P4Con1.fits", doAllDBSCAN.fitsPath+"G2650LocaldbscanS3.0P4Con1.fit",  onlyClean=False )
    sys.exit()









if 0: #produce clean fits for cloud comparison

    doAllDBSCAN.getMaskByLabel(doAllDBSCAN.fitsPath+"G2650LocaldbscanS2P4Con1.fits", doAllDBSCAN.fitsPath+"G2650LocaldbscanS2P4Con1.fit",  onlyClean=True )
    doAllDBSCAN.getMaskByLabel(doAllDBSCAN.fitsPath+"G2650LocaldbscanS2P8Con2.fits", doAllDBSCAN.fitsPath+"G2650LocaldbscanS2P8Con2.fit",  onlyClean=True )
    doAllDBSCAN.getMaskByLabel(doAllDBSCAN.fitsPath+"G2650LocaldbscanS2P11Con3.fits", doAllDBSCAN.fitsPath+"G2650LocaldbscanS2P11Con3.fit",  onlyClean=True )

    sys.exit()






if 0:  #test getting line width
    #doAllDBSCAN.pipeLineG2650ByCut( 6.5  , redo=True,keepFITSFile=False)
    doAllDBSCAN.pipeLineG2650ByCut( 7.0  , redo=True,keepFITSFile=False)

    sys.exit()



if 0:

    #doAllDBSCAN.pureDBSCAN( "negativeG2650Local.fits", 2, MinPts=6 , saveTag="negativeTest",  connectivity=2,inputRMS=None ,redo=True ,keepFITSFile=True, onlyGetFITS=True )
    #dbscanLabelFITS=doAllDBSCAN.pureDBSCAN( "negativeG2650Local.fits", 2, MinPts=6 , saveTag="negativeTest",  connectivity=3,inputRMS=None ,redo=True ,keepFITSFile=True, onlyGetFITS=True )
    #print dbscanLabelFITS

    tbResult=doAllDBSCAN.getMaskByLabel(doAllDBSCAN.fitsPath+"negativeTestdbscanS2P6Con2.fits", doAllDBSCAN.fitsPath+"negativeTestdbscanS2P6Con3.fit", onlySelect=False ,onlyClean=True)

    sys.exit()






if 0: #test negative cases for DBSCAN
    #the reason for using negative values is to testes all noises, to see if the parameters get too many bad sources
    doAllDBSCAN.G2650Local="negativeG2650Local.fits"
    doAllDBSCAN.saveTag="negativeTest"
    doAllDBSCAN.pipeLineG2650ByCut( 2  , redo=True, keepFITSFile=False)
    sys.exit()




if 0:  #test getting line width
    doAllDBSCAN.pipeLineG2650ByCut( 4.5  , redo=True,keepFITSFile=False)
    doAllDBSCAN.pipeLineG2650ByCut( 5.0  , redo=True,keepFITSFile=False)

    sys.exit()





if 0:
    for eachID in [ 79, 234,524,347,22,750  ]:
        doAllDBSCAN.getProjection("/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/testEmpty2dbscanS2P3Con1.fits",   eachID)
    sys.exit()



if 0:

    tbName="testEmpty2dbscanS2P3Con1.fit"
    tb=Table.read(doAllDBSCAN.fitsPath+tbName)
    tb=tb[tb["has22"]>0.5]

    for eachC in tb:
        id=eachC["_idx"]

        doAllDBSCAN.getProjection("/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/testEmpty2dbscanS2P3Con1.fits", id )
    sys.exit()










if 0:
    pass
    #totalFluxList=doAllDBSCAN.getTotalFluxByMask( doAllDBSCAN.fitsPath+ "G2650LocaldbscanS2P3Con1_CO_Masked.fits" )
    #print totalFluxList







if 0: #G2650 pipeline
    pass
    #doAllDBSCAN.pipeLineG2650ByCut(6.5)

    #doAllDBSCAN.pipeLineG2650ByCut(7.0, keepFITSFile=True,redo=True)



if 0: #G2650 pipeline
    doAllDBSCAN.pipeLineG2650( redo=False )



if 0:# test renew

    savePath="/home/qzyan/WORK/myDownloads/MWISPcloud/newDBSCAN/testAll/"

    doAllDBSCAN.getCatalog(doAllDBSCAN.G2650Local,  savePath ,  "G2650LocaldbscanS2P3Con1.fits",  savePath+"G2650LocaldbscanS2P3Con1NewTB",   )



if 0: # draw #
    #
    #
    doAllDBSCAN.drawG2650Local()







if 0:  # test  local
    pass

    doAllDBSCAN = allDBSCAN()


    sys.exit()






if 0: #test empty fits
    doAllDBSCAN=allDBSCAN()

    saveTag="emptyTest"

    doAllDBSCAN.DBSCANTest(doAllDBSCAN.emptyFITS,saveTag, minAreaPix=1, MinPts=8, connecttivity=3    )



if 0: # test all possible cases for DBSCAN

    doAllDBSCAN = allDBSCAN()

    if 1:
        saveTag="testEmpty2"
        rawCOFITS="G2650LocalEmpty2.fits"

        for minPts in  doAllDBSCAN.con1PtsAll:
            doAllDBSCAN.pureDBSCAN(  rawCOFITS , 2, MinPts=minPts, saveTag=saveTag,  connectivity=1,inputRMS=None  )


        for minPts in  doAllDBSCAN.con2PtsAll:
            doAllDBSCAN.pureDBSCAN( rawCOFITS , 2, MinPts=minPts, saveTag=saveTag,  connectivity=2,inputRMS=None  )

        for minPts in  doAllDBSCAN.con3PtsAll:
            doAllDBSCAN.pureDBSCAN(  rawCOFITS , 2, MinPts=minPts, saveTag=saveTag,  connectivity=3,inputRMS=None  )


    if 1:
        saveTag="testEmpty1"
        rawCOFITS="G2650LocalEmpty1.fits"

        for minPts in  doAllDBSCAN.con1PtsAll:
            doAllDBSCAN.pureDBSCAN(  rawCOFITS , 2, MinPts=minPts, saveTag=saveTag,  connectivity=1,inputRMS=None  )


        for minPts in  doAllDBSCAN.con2PtsAll:
            doAllDBSCAN.pureDBSCAN( rawCOFITS , 2, MinPts=minPts, saveTag=saveTag,  connectivity=2,inputRMS=None  )

        for minPts in  doAllDBSCAN.con3PtsAll:
            doAllDBSCAN.pureDBSCAN(  rawCOFITS , 2, MinPts=minPts, saveTag=saveTag,  connectivity=3,inputRMS=None  )



    sys.exit()

