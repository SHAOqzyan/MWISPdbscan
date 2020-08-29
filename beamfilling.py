#!/usr/bin/python
######################################
#Author Qing-Zeng Yan
#A script of to deal with beam filling factors
#at present, only considering


import numpy as np
import matplotlib as mpl
mpl.use('agg')
import scipy.odr.odrpack as odrpack
import radio_beam
from spectral_cube import SpectralCube
from astropy import units as u

from mwispDBSCAN import MWISPDBSCAN

import matplotlib.pyplot as plt
from myPYTHON import *
from astropy.io import fits
import glob
import os
import sys
import seaborn as sns
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from scipy import optimize
from progressbar import *
from astropy.table import Table,vstack
import gc
import scipy


import matplotlib.cm as cm
from scipy.stats import norm
doFITS=myFITS()
def ffFunctionCutTest(x,a ,b,c,  mu,sigma ) :


    fi = 0.5 * (1 + scipy.special.erf((x - mu) * sigma/np.sqrt(2) )    )

    #return (a * x**2 + b*x+e) * (1 - fi) + (c *x+d) * fi

    #return (a * x**2 - b*x+c) * (1 - fi)  #+ d * fi

    return (a * (x-b)**2+  c) * (1 - fi)  #+ d * fi





class myBFF(object):
    ##################### smooth factors

    meanSNRcol="meanSNR"

    #rawCOFITS=None
    #labelFITS=None




    ##
    ##########


    ffMWISPCol = "fillingFactorMWISP"
    ffMWISPErrorCol = "fillingFactorMWISPError"

    ffCfACol = "fillingFactorCfa"
    ffCfAErrorCol = "fillingFactorErrorCfa"

    ffPeakMWISPCol = "ffPeakMWISP"  # do not record the parameters of fitting with peak, fit if needed
    ffPeakMWISPErrorCol = "ffPeakMWISPError"

    aCol = "para_a"
    bCol = "para_b"
    cCol = "para_c"

    aErrorCol = "error_a"
    bErrorCol = "error_b"
    cErrorCol = "error_c"

    ##################### cutoff factors

    cutFFffMWISPCol = "cutFFfillingFactorMWISP"
    cutFFffMWISPErrorCol = "cutFFfillingFactorMWISPError"

    cutFFffCfACol = "cutFFfillingFactorCfa"
    cutFFffCfAErrorCol = "cutFFfillingFactorErrorCfa"

    cutFFaCol = "cutFFpara_a"
    cutFFbCol = "cutFFpara_b"
    cutFFcCol = "cutFFpara_c"

    cutFFaErrorCol = "cutFFerror_a"
    cutFFbErrorCol = "cutFFerror_b"
    cutFFcErrorCol = "cutFFerror_c"

    cutFFmuCol = "cutFFpara_mu"  # mean of Gausian
    cutFFmuErrorCol = "cutFFerror_mu"  # mean of Gausian

    cutFFsigmaCol = "cutFFpara_sigma"  # sigma of Gausian
    cutFFsigmaErrorCol = "cutFFerror_sigma"  # sigma of Gausian



    ######################### defination
    cov00Col = "cov00"
    cov01Col = "cov01"
    cov02Col = "cov02"

    cov10Col = "cov10"
    cov11Col = "cov11"
    cov12Col = "cov12"

    cov20Col = "cov20"
    cov21Col = "cov21"
    cov22Col = "cov22"
    #####################cutcov
    covCut00Col = "covCut00"
    covCut01Col = "covCut01"
    covCut02Col = "covCut02"
    covCut03Col = "covCut03"
    covCut04Col = "covCut04"

    covCut10Col = "covCut10"
    covCut11Col = "covCut11"
    covCut12Col = "covCut12"
    covCut13Col = "covCut13"
    covCut14Col = "covCut14"

    covCut20Col = "covCut20"
    covCut21Col = "covCut21"
    covCut22Col = "covCut22"
    covCut23Col = "covCut23"
    covCut24Col = "covCut24"

    covCut30Col = "covCut30"
    covCut31Col = "covCut31"
    covCut32Col = "covCut32"
    covCut33Col = "covCut33"
    covCut34Col = "covCut34"

    covCut40Col = "covCut40"
    covCut41Col = "covCut41"
    covCut42Col = "covCut42"
    covCut43Col = "covCut43"
    covCut44Col = "covCut44"

    ##############################
    cutoffFactors =  np.arange(2.0, 20.2 ,  0.2  ) #   from 2 sigma to 10 simg

    meanRMS=None
    velres=None


    idCol="_idx"

    def __init__(self):
        pass


    def quadraticDy(self,x,para):
        a, b, c, mu, sigma = para
        fi = 0.5 * (1 + scipy.special.erf((x - mu) * sigma / np.sqrt(2)))

        dya = (1-fi) *(x-b)**2
        dyb = (1-fi) * ( 2*a*(b-x)  )
        dyc = (1-fi)
        dymu = -(a*(x-b)**2 + c )*1./np.sqrt(np.pi) *np.exp( -(x-mu)**2 *sigma**2/2.  )*sigma/np.sqrt(2)*(-1)

        dysigma = -(a*(x-b)**2 + c )*1./np.sqrt(np.pi) *np.exp( -(x-mu)**2 *sigma**2/2.   )*(x-mu)/np.sqrt(2)


        return np.asarray( [dya, dyb, dyc, dymu, dysigma ]  )
    def getFFandErrorCutoff(self,  para,pcov,targetSNR=2,dimension=5):
        """
        get the flux at SNR 2
        :param para:
        :param pcov:
        :param targetSNR:
        :return:
        """

        if len(para) ==5: #consider quadratic function

            a,b,c,mu,sigma=para

            fx = ffFunctionCutTest(targetSNR,*para)
            f0 =  ffFunctionCutTest(0,*para)

            dfx=self.quadraticDy(targetSNR,para)
            df0= self.quadraticDy(0,para)

            varX = np.matmul(dfx, pcov)
            varX = np.matmul(varX, dfx)

            var0 = np.matmul(df0, pcov)
            var0 = np.matmul(var0, df0)

            varFF = varX / f0 ** 2 + (-fx / f0 ** 2) ** 2 * var0


            return fx / f0, np.sqrt(varFF)


        if len(para) == 3:#Gaussian model
            a,b,c=para
            maxPoint= max( [0,c])

            fx = ffFunctionCutSG(targetSNR,*para)
            f0 =  ffFunctionCutSG(maxPoint,*para)

            dfx=self.SGDy(targetSNR,para)


            df0= self.SGDy(maxPoint,para)

            varX = np.matmul(dfx, pcov)
            varX = np.matmul(varX, dfx)

            var0 = np.matmul(df0, pcov)
            var0 = np.matmul(var0, df0)

            varFF = varX / f0 ** 2 + (-fx / f0 ** 2) ** 2 * var0


            return fx / f0, np.sqrt(varFF)


    def getPcovCutOff(self,eachRow, dimension=5):
        """

        :param eachRow:
        :param dimension:
        :return:
        """
        pcov=np.zeros( (dimension,dimension) )

        ###

        pcov[0, 0] = eachRow[self.covCut00Col]
        pcov[0, 1] = eachRow[self.covCut01Col]
        pcov[0, 2] = eachRow[self.covCut02Col]

        pcov[1, 0] = eachRow[self.covCut10Col]
        pcov[1, 1] = eachRow[self.covCut11Col]
        pcov[1, 2] = eachRow[self.covCut12Col]

        pcov[2, 0] = eachRow[self.covCut20Col]
        pcov[2, 1] = eachRow[self.covCut21Col]
        pcov[2, 2] = eachRow[self.covCut22Col]

        if dimension==3:
            return pcov


        if dimension!=5:
            return None

        pcov[0, 3] = eachRow[self.covCut03Col]
        pcov[0, 4] = eachRow[self.covCut04Col]

        pcov[1, 3] = eachRow[self.covCut13Col]
        pcov[1, 4] = eachRow[self.covCut14Col]

        pcov[2, 3] = eachRow[self.covCut23Col]
        pcov[2, 4] = eachRow[self.covCut24Col]



        pcov[3, 0] = eachRow[self.covCut30Col]
        pcov[3, 1] = eachRow[self.covCut31Col]
        pcov[3, 2] = eachRow[self.covCut32Col]
        pcov[3, 3] = eachRow[self.covCut33Col]
        pcov[3, 4] = eachRow[self.covCut34Col]


        pcov[4, 0] = eachRow[self.covCut40Col]
        pcov[4, 1] = eachRow[self.covCut41Col]
        pcov[4, 2] = eachRow[self.covCut42Col]
        pcov[4, 3] = eachRow[self.covCut43Col]
        pcov[4, 4] = eachRow[self.covCut44Col]

        return pcov




    def savePcovCutoff(self,eachRow,pcov):


        #with three or five


        eachRow[self.covCut00Col] = pcov[0, 0]
        eachRow[self.covCut01Col] = pcov[0, 1]
        eachRow[self.covCut02Col] = pcov[0, 2]

        eachRow[self.covCut10Col] = pcov[1, 0]
        eachRow[self.covCut11Col] = pcov[1, 1]
        eachRow[self.covCut12Col] = pcov[1, 2]


        eachRow[self.covCut20Col] = pcov[2, 0]
        eachRow[self.covCut21Col] = pcov[2, 1]
        eachRow[self.covCut22Col] = pcov[2, 2]




        if pcov.shape[0]==5: #for five paramters


            eachRow[self.covCut03Col] = pcov[0, 3]
            eachRow[self.covCut04Col] = pcov[0, 4]


            eachRow[self.covCut13Col] = pcov[1, 3]
            eachRow[self.covCut14Col] = pcov[1, 4]


            eachRow[self.covCut23Col] = pcov[2, 3]
            eachRow[self.covCut24Col] = pcov[2, 4]
            ###########

            eachRow[self.covCut30Col] = pcov[3, 0]
            eachRow[self.covCut31Col] = pcov[3, 1]
            eachRow[self.covCut32Col] = pcov[3, 2]
            eachRow[self.covCut33Col] = pcov[3, 3]
            eachRow[self.covCut34Col] = pcov[3, 4]


            #################
            eachRow[self.covCut40Col] = pcov[4, 0]
            eachRow[self.covCut41Col] = pcov[4, 1]
            eachRow[self.covCut42Col] = pcov[4, 2]
            eachRow[self.covCut43Col] = pcov[4, 3]
            eachRow[self.covCut44Col] = pcov[4, 4]


    def getIndices(self,labelSets, choseID ,return2D=False):
        Z0, Y0, X0, values1D =labelSets
        cloudIndices = np.where(values1D == choseID)

        cX0 = X0[cloudIndices]
        cY0 = Y0[cloudIndices]
        cZ0 = Z0[cloudIndices]

        if not return2D:

            return tuple([cZ0, cY0, cX0])
        else:
            return tuple([cZ0, cY0, cX0]), tuple([  cY0, cX0])



    def getLabelSet(self,labelData):
        """

        :param labelData:
        :return:
        """
        noiseMask= np.min(   labelData[0] )

        clusterIndex1D = np.where(labelData > noiseMask )
        clusterValue1D = labelData[clusterIndex1D]

        Z0, Y0, X0 = clusterIndex1D
        labelSets=[Z0, Y0, X0, clusterValue1D ]

        return  labelSets


    def addCutFFColnames(self, ffTB):
        cleanTB=ffTB
        ffTB[self.cutFFffMWISPCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFffMWISPCol].unit = ""

        ffTB[self.cutFFffCfACol] = cleanTB["peak"] * 0
        ffTB[self.cutFFffCfACol].unit = ""

        ffTB[self.cutFFaCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFaCol].unit = ""

        ffTB[self.cutFFmuCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFmuCol].unit = ""

        ffTB[self.cutFFmuErrorCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFmuErrorCol].unit = ""

        ffTB[self.cutFFsigmaCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFsigmaCol].unit = ""

        ffTB[self.cutFFsigmaErrorCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFsigmaErrorCol].unit = ""

        ###############################

        ffTB[self.cutFFbCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFcCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFaErrorCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFbErrorCol] = cleanTB["peak"] * 0
        ffTB[self.cutFFcErrorCol] = cleanTB["peak"] * 0

        ffTB[self.cutFFbCol].unit = ""
        ffTB[self.cutFFcCol].unit = ""
        ffTB[self.cutFFaErrorCol].unit = ""
        ffTB[self.cutFFbErrorCol].unit = ""
        ffTB[self.cutFFcErrorCol].unit = ""

        #############
        ffTB[self.covCut00Col] = cleanTB["peak"] * 0
        ffTB[self.covCut01Col] = cleanTB["peak"] * 0
        ffTB[self.covCut02Col] = cleanTB["peak"] * 0
        ffTB[self.covCut03Col] = cleanTB["peak"] * 0
        ffTB[self.covCut04Col] = cleanTB["peak"] * 0

        ffTB[self.covCut10Col] = cleanTB["peak"] * 0
        ffTB[self.covCut11Col] = cleanTB["peak"] * 0
        ffTB[self.covCut12Col] = cleanTB["peak"] * 0
        ffTB[self.covCut13Col] = cleanTB["peak"] * 0
        ffTB[self.covCut14Col] = cleanTB["peak"] * 0

        ffTB[self.covCut20Col] = cleanTB["peak"] * 0
        ffTB[self.covCut21Col] = cleanTB["peak"] * 0
        ffTB[self.covCut22Col] = cleanTB["peak"] * 0
        ffTB[self.covCut23Col] = cleanTB["peak"] * 0
        ffTB[self.covCut24Col] = cleanTB["peak"] * 0

        ffTB[self.covCut30Col] = cleanTB["peak"] * 0
        ffTB[self.covCut31Col] = cleanTB["peak"] * 0
        ffTB[self.covCut32Col] = cleanTB["peak"] * 0
        ffTB[self.covCut33Col] = cleanTB["peak"] * 0
        ffTB[self.covCut34Col] = cleanTB["peak"] * 0

        ffTB[self.covCut40Col] = cleanTB["peak"] * 0
        ffTB[self.covCut41Col] = cleanTB["peak"] * 0
        ffTB[self.covCut42Col] = cleanTB["peak"] * 0
        ffTB[self.covCut43Col] = cleanTB["peak"] * 0
        ffTB[self.covCut44Col] = cleanTB["peak"] * 0
        ffTB[self.meanSNRcol] =cleanTB["peak"]*0

    def getFluxColNameCutoff(self,cutoffFactor):
        """
        The pix name is is to record the tao pix Number of the flux, which is used to estimate the error of the total flux
        :param smFactor:
        :return:
        """
        return "fluxCut{:.1f}".format( cutoffFactor ), "pixCut{:.1f}".format( cutoffFactor )



    def addCutOffFluxColnames(self, cleanTB):
        """
        adding columns to record cutoff flux
        :param cleanTB:
        :return:
        """
        for eachCut in self.cutoffFactors:
            colName, colPixName = self.getFluxColNameCutoff(eachCut)

            cleanTB[colName] = cleanTB["peak"] * 0
            cleanTB[colPixName] = cleanTB["peak"] * 0
            cleanTB[colPixName].unit = ""



    def getSmoothFluxColCutoff(self, rawCOFITS, labelFITS,TB  ,   rmsFITS = None,meanrms=0.5  ):
        """

        get the change of flux according to the cutoff sigmas, remember, this has to be cut according to the rmsData
        need rmsFITS to calculate

        :return:
        """

        ####
        print "Extracting cutoff flux from ", rawCOFITS

        dataLabel , headLabel = doFITS.readFITS( labelFITS )
        dataRaw , headRaw = doFITS.readFITS( rawCOFITS )
        Nz,Ny,Nx = dataRaw.shape
        if rmsFITS is not None:
            rmsData,rmsHead=doFITS.readFITS( rmsFITS )
        else:

            rmsData=np.zeros( (Ny,Nx)  )+meanrms
        #TB=Table.read(TBFile)
        #self.addCutOffFluxColnames(TB)
        velResolution = headRaw["CDELT3"]/1000.

        self.meanRMS=np.nanmean(rmsData)
        self.velres= velResolution






        labelSets= self.getLabelSet(dataLabel) #[Z0, Y0, X0, clusterValue1D ]

        #colName,colPixName=self.getFluxColNameNoiseChange(cutoffFactor)

        widgets = ['Extracting cutoff flux:', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
                   FileTransferSpeed()]  # see docs for other options
        pbar = ProgressBar(widgets=widgets, maxval=len(TB))
        pbar.start()

        TB.sort("area_exact")
        TB.reverse()
        i=0
        for eachRow in TB:
            i=i+1
            #gc.collect()

            ID = eachRow["_idx"]

            cloudIndex,cloudIndex2D = self.getIndices(labelSets, ID,return2D=True)  # np.where( cleanDataSM1==ID )
            coValues=   dataRaw[cloudIndex]
            rmsValues= rmsData[cloudIndex2D]


            for eachCutoff in self.cutoffFactors:

                selectByCut=   coValues>=rmsValues*eachCutoff
                selectCoValues= coValues[  selectByCut ]

                #print selectCoValues

                colName, colPixName = self.getFluxColNameCutoff(eachCutoff)

                eachRow[colName]= np.sum( selectCoValues ,dtype=float )* velResolution
                eachRow[colPixName] = len(selectCoValues)
                #print ID,  eachRow[colPixName] ,  eachRow[colName]

            eachRow[self.meanSNRcol]= np.mean(  coValues/rmsValues  )
            pbar.update(i)

        pbar.finish()
        #TB.write("cutoffFlux_"+TBFile , overwrite=True)

        TB=self.calculateFillingFactorCutoff( TB,  drawFigure= False ,dim=5 )

        #return TB

        #TB.write("cutoffFF_"+TBFile , overwrite=True)

    def getFluxList(self, row):
        """

        :param row:
        :return:
        """


        fluxList = []
        fluxErrorList = []



        for eachCutoff in self.cutoffFactors:
            colName, colPixName = self.getFluxColNameCutoff(eachCutoff)
            fluxList.append(row[colName])  # the flux is already K km/s
            totalVox = row[colPixName]



            totalVox = max([1, totalVox])  # assign 1 pixel error to 0 flux

            eRRor = np.sqrt(totalVox) * self.meanRMS * self.velres

            fluxErrorList.append(eRRor)

        return np.asarray(fluxList), np.asarray(fluxErrorList)

    def calculateFillingFactorCutoff(self,TB,drawFigure=False,inputID=None,saveTag="cutoffBFF", dim=5):

        """
        :param TB:
        :return:
        """
        #TB=Table.read(TBFile)
        #saveName="cutFF_"+TBFile

        #TB=self.addCutFFColnames(TB)
        cutoffArray= self.cutoffFactors #* self.getMeanRMS()



        #add progress bar
        widgets = ['Calculating cutoff filling factors: ', Percentage(), ' ', Bar(marker='0', left='[', right=']'), ' ', ETA(), ' ',
                   FileTransferSpeed()]  # see docs for other options
        pbar = ProgressBar(widgets=widgets, maxval=len(TB))
        pbar.start()

        i=0


        if inputID!=None: #replace TB

            #select TB by ID
            TB=TB[  TB[self.idCol]==inputID  ]




        for eachRow in TB:
            i=i+1
            pbar.update(i)
            fluxList,fluxError = self.getFluxList(eachRow  )

            #only use position values for cutoff


            ID=eachRow[ self.idCol ]
            if inputID!=None and ID==inputID: #for debug

                wmsipFilling, cfaFilling, fittingParaAndError,pcov = self.getFillingFactorAndDrawCutoff(cutoffArray, fluxList,   fluxError, calID=ID,   drawFigure=True,drawCode = saveTag, individuals= True,   dim= dim )
                print "filling factor", wmsipFilling
                print "Parameters", fittingParaAndError[0]
                print "Para error", fittingParaAndError[1]

                return

            #crop

            wmsipFilling, cfaFilling,  fittingParaAndError ,pcov =self.getFillingFactorAndDrawCutoff(cutoffArray,fluxList,fluxError,calID=ID,drawFigure=drawFigure, drawCode =saveTag, individuals= False  ,dim= dim )




            para, paraError = fittingParaAndError
            eachRow[self.cutFFffMWISPCol]=wmsipFilling #will be recalculate
            eachRow[self.cutFFffCfACol]=cfaFilling #0



            eachRow[self.cutFFaCol]=para[0]
            eachRow[self.cutFFbCol]=para[1]
            eachRow[self.cutFFcCol] = para[2]
            eachRow[self.cutFFmuCol] = para[3]
            eachRow[self.cutFFsigmaCol] = para[4]

            eachRow[self.cutFFaErrorCol]=paraError[0]
            eachRow[self.cutFFbErrorCol]=paraError[1]
            eachRow[self.cutFFcErrorCol]=paraError[2]
            eachRow[self.cutFFmuErrorCol]=paraError[3]
            eachRow[self.cutFFsigmaErrorCol ]=paraError[4]




            self.savePcovCutoff(eachRow,pcov )



        pbar.finish()
        #TB.write( saveName , overwrite=True )
        return TB


    def fittingCutOFF(self , beamList,fluxList,fluxError ,dim=5  ):
        x=np.asarray( beamList )
        y = np.asarray( fluxList )
        yError= np.asarray( fluxError )





            #y= y [2:   ]
        #x   =   x[ 2:     ]  #( x[0:-1] +  x[1: ] )/2.
        #yError =  yError[ 2:     ]

        params1, paramas_covariance1 = optimize.curve_fit( ffFunctionCutTest , x, y, sigma=yError, absolute_sigma=True, p0=         [  np.max(y)  /10  ,    5    ,   np.max(y)/10    ,   7  , 1  ]  )

        return params1, paramas_covariance1, ffFunctionCutTest

    def getFillingFactorAndDrawCutoff(self,beamList,fluxList,fluxError,calID,drawFigure=False,drawCode="",testPoly=False,individuals=False, dim=5  ):
        """

        :param beamList:
        :param fluxList:
        :return:
        """
        x = np.asarray(beamList)

        y = np.asarray(fluxList)

        yError = np.asarray(fluxError)

        #yError=yError/np.max(y)
        #y=y/np.max(y)

        try:

            params, paramas_covariance, useFunction = self.fittingCutOFF(x, y, yError, dim=dim)

            #print calID,params

            errors = np.sqrt(np.diag(paramas_covariance))

        except:
            #params, paramas_covariance, useFunction = self.fittingCutOFF(x, y, yError, dim=dim)

            # params, paramas_covariance = optimize.curve_fit(useFunction, x, y, sigma=yError, absolute_sigma=True,
            # p0=p0My)

            return 0, 0, [[0, 0, 0,0,0], [0, 0, 0,0,0]], np.asarray([[0, 0, 0,0,0], [0, 0, 0,0,0], [0, 0, 0,0,0], [0, 0, 0,0,0], [0, 0, 0,0,0]])



        # errors = np.sqrt(np.diag(paramas_covariance))

        fittingParaAndError = [params, errors]


        mwispFF, mwispFFError = self.getFFandErrorCutoff( params, paramas_covariance,  targetSNR=2 ,dimension= len(params)  )

        # print mwispFF,wmsipFilling,"Equal??Equal??Equal??Equal??Equal??Equal??"

        if not drawFigure:
            return mwispFF, 0, fittingParaAndError, paramas_covariance

        # drawTheFigure
        fig = plt.figure(figsize=(10, 8))
        rc('text', usetex=True)
        rc('font', **{'family': 'sans-serif', 'size': 22, 'serif': ['Helvetica']})

        if 1:
            mpl.rcParams['text.latex.preamble'] = [
                r'\usepackage{tgheros}',  # helvetica font
                r'\usepackage{sansmath}',  # math-font matching  helvetica
                r'\sansmath',  # actually tell tex to use it!
                r'\usepackage{siunitx}',  # micro symbols
                r'\sisetup{detect-all}',  # force siunitx to use the fonts
            ]

        axFitting = fig.add_subplot(1, 1, 1)

        axFitting.errorbar(x, y, yerr=yError, markersize=2, fmt='o', c='red', capsize=0.5, elinewidth=1, lw=1)

        fittingX = np.arange(0, np.max(x), 0.01)


        if useFunction== ffFunctionCutSG:
            label="single Gaussian "
        else:
            label="Double Gaussian "
        if useFunction == ffFunctionCutTest:
            label= "$\mathit f={:.3f}\pm{:.3f}$".format(mwispFF ,mwispFFError  ) ##"Quadratic with Gaussian CDF "


        ffs = axFitting.plot(fittingX, useFunction(fittingX, *params), color='blue', lw=1.0 ,label=label     )


        axFitting.axvline(x=0, ls="--", color='black', lw=1)



        axFitting.legend(loc=1, handlelength=1, fontsize=20)

        axFitting.set_xlabel(r"Cutoff ($\sigma$)")
        # axFitting.set_ylabel("Flux (arcmin)")
        axFitting.set_ylabel(r"Flux ($\rm K\ km\ s$$^{-1}$ $\mathrm{\Omega}_\mathrm{A}$)")
        ###########################
        saveTagFigue = "{}_factorFitting_ID{}{}".format(self.calCode, calID, drawCode)




        if calID != 0:

            if individuals and not testPoly:
                saveTagFigue = saveTagFigue + "_indivi"
                plt.savefig(self.figurePath + "{}.pdf".format(saveTagFigue), bbox_inches='tight')

            # if calID!=102:

            # plt.savefig(self.figurePath+"{}.png".format( saveTagFigue ), bbox_inches='tight', dpi=100)
            # plt.savefig(self.figurePath+"{}.pdf".format( saveTagFigue ), bbox_inches='tight' )
            if testPoly:
                saveTagFigue = saveTagFigue + "_indivi"
                plt.savefig(self.figurePath + "{}Poly.pdf".format(saveTagFigue), bbox_inches='tight')

            if not testPoly and not individuals:
                plt.savefig(self.figurePath + "{}.png".format(saveTagFigue), bbox_inches='tight', dpi=100)



        else:
            plt.savefig(self.paperFigurePath + "{}.png".format(saveTagFigue), bbox_inches='tight', dpi=100)
            plt.savefig(self.paperFigurePath + "{}.pdf".format(saveTagFigue), bbox_inches='tight')

        plt.close(fig)
        gc.collect()

        return mwispFF, 0, fittingParaAndError, paramas_covariance


    def getFillingAndErrorCutoffByTB(self,  ffTB, targetSNR=2 ,dim=5):

        """
        :param ffTB:
        :return:
        """


        ffList = []
        ffListError = []

        for eachRow in ffTB:

            # first care dimensioan 5, with quadratic
            x=  targetSNR #arcmin
            a=eachRow[self.cutFFaCol]
            b=eachRow[self.cutFFbCol]
            c=eachRow[self.cutFFcCol]

            para = np.asarray([a, b, c ])

            if dim==5:

                mu=eachRow[self.cutFFmuCol]
                sigma=eachRow[self.cutFFsigmaCol ]


                para=np.asarray( [a,b,c,mu,sigma] )


            pcov=self.getPcovCutOff(eachRow, dimension=dim)

            f,fe=self.getFFandErrorCutoff(para,pcov,x,dimension= len(para) )





            ffList.append( f )

            ffListError.append( fe )


        return np.asarray(ffList),np.asarray( ffListError )


    def addMWISPFFerrorCutoff(self, TB ,dim=5 ):
        """

        :return:
        """

        mwispFF,mwispFFError= self.getFillingAndErrorCutoffByTB(  TB  , targetSNR=2,dim= dim )

        TB[self.cutFFffMWISPCol] = mwispFF

        TB[self.cutFFffMWISPErrorCol] = mwispFFError





    def getBFFcutoff(self, rawCOFITS, labelFITS,  tbFile, rmsFITS=None,saveTBname=None ):
        """
        append cutoff factors to the tbsanve saveFITS
        :param rawCOFITS:
        :param labelFITS:
        :param tbFile:
        :return:
        """
        if saveTBname is None:
            saveTBname=tbFile[0:-4]+"_BFF.fit"

        #step1, add colnames
        cleanTB=Table.read(tbFile)

        self.addCutOffFluxColnames(cleanTB)
        self.addCutFFColnames(cleanTB)
        ##get flux

        self.getSmoothFluxColCutoff(rawCOFITS,labelFITS,cleanTB,rmsFITS=rmsFITS)


        self.addMWISPFFerrorCutoff(cleanTB,dim=5)
        cleanTB.write( saveTBname )

        #
        #step2, get flux

    def recalCulateBFF(self,tbName,velResolution=0.2,  rmsFITS=None,meanrms =0.5):


        if rmsFITS is not None:
            rmsData,rmsHead=doFITS.readFITS( rmsFITS )
        else:

            rmsData=np.zeros( (Ny,Nx)  )+meanrms
        #TB=Table.read(TBFile)
        #self.addCutOffFluxColnames(TB)


        self.meanRMS=np.nanmean(rmsData)
        self.velres= velResolution

        print self.velres,"The value of velocity resolutuion  "


        TB=Table.read(tbName)

        TB=self.calculateFillingFactorCutoff( TB,  drawFigure= False ,dim=5 )

        self.addMWISPFFerrorCutoff(TB)
        TB.write( tbName,overwrite=True )
        print TB

    def ZZZ(self):
        pass


doFF=myBFF()

doFF.getBFFcutoff(  "/T620/ysu/t40160/mosaic_U02.fits", "mosaic_U02dbscanS2P4Con1_Clean.fits",  "mosaic_U02dbscanS2P4Con1_Clean.fit",   rmsFITS="/T620/ysu/t40160/mosaic_U02_rms.fits" )
doFF.getBFFcutoff(  "/T620/ysu/t40160/mosaic_L02.fits", "mosaic_L02dbscanS2P4Con1_Clean.fits", "mosaic_L02dbscanS2P4Con1_Clean.fit",rmsFITS="/T620/ysu/t40160/mosaic_L02_rms.fits" )
doFF.getBFFcutoff(  "/T620/ysu/t40160/mosaic_U05_40160.fits", "mosaic_U05_40160dbscanS2P4Con1_Clean.fits", "mosaic_U05_40160dbscanS2P4Con1_Clean.fit",rmsFITS="/T620/ysu/t40160/mosaic_U05_40160_rms.fits" )
doFF.getBFFcutoff(  "/T620/ysu/t40160/mosaic_L05_40160.fits", "mosaic_L05_40160dbscanS2P4Con1_Clean.fits", "mosaic_L05_40160dbscanS2P4Con1_Clean.fit",rmsFITS="/T620/ysu/t40160/mosaic_L05_40160_rms.fits" )



#doFF.recalCulateBFF("mosaic_U02dbscanS2P4Con1_Clean_BFF.fit", rmsFITS="/T620/ysu/t40160/mosaic_U02_rms.fits",velResolution=0.2 )
#doFF.recalCulateBFF("mosaic_L02dbscanS2P4Con1_Clean_BFF.fit", rmsFITS="/T620/ysu/t40160/mosaic_L02_rms.fits",velResolution=0.2 )

#doFF.recalCulateBFF("mosaic_U05_40160dbscanS2P4Con1_Clean_BFF.fit", rmsFITS="/T620/ysu/t40160/mosaic_U05_40160_rms.fits",velResolution=0.5 )
#doFF.recalCulateBFF("mosaic_L05_40160dbscanS2P4Con1_Clean_BFF.fit", rmsFITS="/T620/ysu/t40160/mosaic_L05_40160_rms.fits",velResolution=0.5 )

#doFF.recalCulateBFF("mosaic_U02dbscanS2P4Con1_Clean_BFF.fit","Q1Sub.fits" )

#for testing

#doFF=myBFF()
#doFF.getBFFcutoff(  "Q1Sub.fits", "Q1SubdbscanS2P4Con1_Clean.fits", "Q1SubdbscanS2P4Con1_Clean.fit" )