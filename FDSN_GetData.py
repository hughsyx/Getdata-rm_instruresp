 #encoding=utf8
import glob
from obspy.core import UTCDateTime, Stream
import obspy.signal.invsim
import os
import sys
import numpy
import scipy.signal
import obspy.clients.fdsn
import scipy.interpolate
from obspy.io.xseed import Parser
from pylab import *
import warnings
warnings.filterwarnings("ignore")
client = obspy.clients.fdsn.Client()
##### AVAILABLE DATACENTERS 
# GFZ     http://geofon.gfz-potsdam.de
# IRIS    http://service.iris.edu
# SCEDC   http://service.SCEDC.org
# NERIES  http://www.seismicportal.eu
# ORFEUS  http://www.orfeus-eu.org
# RESIF   http://ws.resif.fr
# USGS    http://comcat.cr.usgs.gov
# USP     http://sismo.iag.usp.br
#clientIRIS = obspy.fdsn.Client("IRIS")
clientSCEDC = obspy.clients.fdsn.Client("SCEDC")
# .... if you need others

def getData(arg_list):

    #######################
    ##### PARAMETERS ######
    #######################
    DictOfStation =  {'ARV':'CI','JTH'} #{'station1': 'network1', 'station2': 'network'...}
    DayStart = 1
    YearStart = 2015
    DayEnd = 2
    YearEnd = 2015
    ListOfComponents = ['BHZ'] # you can put a list here ...
    Direc = 'test'
    NewFreq = 5 # resampling at ...
    f_prefilt=(0.005, 0.01, 2.0, 2.25) # prefilter before deconvolution ...
    sizeMin= 250000 # REMOVE FILES SMALLER THAN THAT (Bytes).. in case the 
    FormatOut = 'SAC'
    ###################################
    ###################################
    ###################################

    NumberOfSubList = 1 #int(arg_list[0])    ## if you want to make Sublists ...
    SubList = 0 #int(arg_list[1])
    SubDictOfStation ={}
    ListOfStation = makeListOfSubListFromList(DictOfStation.keys(), NumberOfSubList)[SubList-1]
    ListOfNetwork = makeListOfSubListFromList(DictOfStation.values(), NumberOfSubList)[SubList-1]
    SubDictOfStation = dict(zip(ListOfStation, ListOfNetwork)) 
    DirecFig  = Direc + os.sep + 'FIG'

    for Year in range(YearStart,YearEnd+1):

        tempDayStart = 1
        tempDayEnd = 365

        if Year == YearStart:
            tempDayStart = DayStart

        if Year == YearEnd:
            tempDayEnd = DayEnd

        if IsLeapYear(Year) and DayEnd == 365:
            DayEnd = 366

        for i in range(tempDayStart,tempDayEnd+1):

            #################################################
            ####    date
            #################################################         
            if i < 10:
                sti = "00" + str(i)
            elif i < 100:
                sti = "0" + str(i)
            else:
                sti = str(i)
         
            stT1 = str(Year) + sti + "T000000.0"
            t1 = UTCDateTime(stT1)
            stT2 = str(Year) + sti + "T235959.999999"
            t2 = UTCDateTime(stT2)
            print stT1
        
            for station, network in SubDictOfStation.items():
                #print '\t' + station + ' ' + network
                for compo in ListOfComponents:

                    DirecData = Direc + os.sep + compo

                    #################################################
                    ####    directories
                    #################################################
                    
                    nameDir    = DirecData + os.sep + str(Year) + os.sep + network + '_' + station
                    nameDirFIG = DirecFig + os.sep + str(Year) + os.sep + network + '_' + station
                    nameDirOut_nodata = DirecData +  os.sep + 'NO_DATA' + os.sep + network + '_' + station
                    nameFileOut_nodata = nameDirOut_nodata + os.sep + network + '.' + station + '.' + compo + '.' + sti + '.' + str(Year)
                    
                    if (glob.glob(nameDir + os.sep + network + '.' + station + '.*.' + compo + '.' + sti + '*.SAC') == []) and (glob.glob(nameFileOut_nodata + '.NO_DATA') == []) and (glob.glob(nameFileOut_nodata + '.AlmostNO_DATA') == []):
                        
                    #################################################
                    ####    DownloadFromClient
                    #################################################
                    
                        try:
                            DownloadFromClient("SCEDC",network, station,compo, t1, t2, sti, nameDir, nameDirFIG, NewFreq, f_prefilt, FormatOut)
                        except:
                            Nodata(nameDirOut_nodata,network,station,compo,sti,Year) ### if you want to try HH in case that BH is not working
                            pass                                                     ### just comment that an uncoment the following lines
                                #if compo[0:2]=="BH":
                                #   compo2 = "HH" + compo[2:3]
                                #try:
                                #    DownloadFromClient("SCEDC",network, station,compo2, t1, t2, sti, nameDir, nameDirFIG, NewFreq, f_prefilt, FormatOut)
                                #except:
                                #   try:                        
                                #       DownloadFromClient("SCEDC",network, station,compo2, t1, t2, sti, nameDir, nameDirFIG, NewFreq, f_prefilt, FormatOut)
                                #    except:
                                #       Nodata(nameDirOut_nodata,network,station,compo,sti,Year)
                                #       pass

                        CheckOutput(nameDir, nameDirOut_nodata,FormatOut,sizeMin,network,station,compo,sti,Year)

#################################################
####    FUNCTIONS
#################################################


def DownloadFromClient(CLIENT,network, station,compo, t1, t2, sti, nameDir, nameDirFIG, NewFreq, f_prefilt, FormatOut):
    if CLIENT=="SCEDC":
        st = clientSCEDC.get_waveforms(network, station,'A0',compo, t1, t2, attach_response=True)
        print('   SCEDC --->')
    MakeDir(nameDir)
#st = DeconvDecimate(st,NewFreq,f_prefilt)
    nameFileStream = nameDir + os.sep + network + '.' + station + '.' + st[0].stats.location + '.' + compo + '.' + sti
    print('         ' + network + '.' + station + '.' + st[0].stats.location + '.' + compo + '.' + sti)
    st.write(nameFileStream + '.SAC', format = FormatOut)


def DeconvDecimate(st,NewFreq,f_prefilt):
    for tr in st:
        tr.detrend(type='constant')
        tr.detrend(type='linear')
        tr.data=glitchCorrectionWithFactorStd(tr.data, 10)
        tr=MY_makeNewFrequenceTrace(tr, NewFreq)
        tr.remove_response(output="VEL",taper = True, taper_fraction=0.05,pre_filt= f_prefilt ,water_level = 60.0)
    return st


def MY_makeNewFrequenceTrace(st, NewFrequence):#destroy the trace,
    old_freq=st.stats.sampling_rate
    rateFreq = float(st.stats.sampling_rate)/float(NewFrequence)
    if st.stats.sampling_rate==int(st.stats.sampling_rate) and rateFreq==int(rateFreq):
        if int(rateFreq)>1:
            st.filter('lowpass_Cheby_2',freq=float(NewFrequence))
            st.decimate(int(rateFreq),no_filter=True)
    else:
        rateFreq_closest=round(rateFreq)
        temp_freq=float(NewFrequence)*float(rateFreq_closest)
        st.data=makeInterpolationNumpy(st.data, float(st.stats.sampling_rate), temp_freq)
        st.stats.sampling_rate=temp_freq
        new_rateFreq = int(float(st.stats.sampling_rate)/float(NewFrequence))
        st.filter('lowpass_Cheby_2',freq=float(NewFrequence))
        st.decimate(new_rateFreq,no_filter=True)
    return st


def makeInterpolationNumpy(Trace, Frequence, NewFrequence):
    VectorPeriodTrace = numpy.arange(0,len(Trace)/float(Frequence),1.0/float(Frequence))
    newVectorPeriodTrace = numpy.arange(0,int(len(Trace)/Frequence),1.0/NewFrequence)
    return numpy.interp(newVectorPeriodTrace,VectorPeriodTrace,Trace)


def makeListOfSubListFromList(List, numberSubList=1):
    if numberSubList < 1:
        numberSubList = 1
    LenShortSubList = int(len(List)/numberSubList)
    LenLongSubList = int(len(List)/numberSubList)+1
    LenTotalLongSublist = (len(List)-LenShortSubList*numberSubList)*LenLongSubList
    ListOfSubList = []
    for indice in range(0, LenTotalLongSublist, LenLongSubList):
        ListOfSubList.append(List[indice:indice+LenLongSubList])
    if not LenShortSubList == 0:
        for indice in range(LenTotalLongSublist, len(List), LenShortSubList):
            ListOfSubList.append(List[indice:indice+LenShortSubList])
    return ListOfSubList


def glitchCorrectionWithFactorStd(Trace, FactorTestStd, NumberOfStd = 1, FactorReplaceWithStd = 0):#Trace is destroyed
    for i in range(NumberOfStd):
        arrayReplace = numpy.ones(len(Trace), dtype ='float')*numpy.std(Trace)*FactorReplaceWithStd
        arrayReplace *= numpy.sign(Trace)
        Trace = scipy.where(scipy.absolute(Trace)>FactorTestStd*numpy.std(Trace), arrayReplace, Trace)
    return Trace


def CheckOutput(nameDir, nameDirOut_nodata,FormatOut,sizeMin,network,station,compo,sti,Year):
    nameFileStream = nameDir + os.sep + network + '.' + station + '*' + compo + '.' + sti
    if glob.glob(nameFileStream + '*')!=[]:
        for filetest in glob.glob(nameFileStream + '*'):
            if os.path.getsize(filetest) <= sizeMin:
                os.remove(filetest)
                #print('         ' + filetest + ' is too small !')
        if glob.glob(nameFileStream + '*')==[]:
            Smalldata(nameDirOut_nodata,network,station,compo,sti,Year)
    try:
        os.rmdir(nameDir)
    except:
        pass


def Nodata(nameDirOut_nodata,network,station,compo,sti,Year):
    MakeDir(nameDirOut_nodata)
    nameFileOut_nodata = nameDirOut_nodata + os.sep + network + '.' + station + '.' + compo + '.' + sti + '.' + str(Year)
    nameFileStream = network + '.' + station + '.' + compo + '.' + sti
    print('      ! NODATA ! ' + nameFileStream)
    touch(nameFileOut_nodata + '.NO_DATA')


def Smalldata(nameDirOut_nodata,network,station,compo,sti,Year):
    MakeDir(nameDirOut_nodata)
    nameFileOut_nodata = nameDirOut_nodata + os.sep + network + '.' + station + '.' + compo + '.' + sti + '.' + str(Year)
    nameFileStream = network + '.' + station + '.' + compo + '.' + sti
    print('      ! AlmostNODATA ! ' + nameFileStream)
    touch(nameFileOut_nodata + '.AlmostNO_DATA')


def touch(fname, times=None):
    with file(fname, 'a'):
        os.utime(fname, times)

def MakeDir(nameDir):
    try:
        os.makedirs(nameDir)
    except:
        pass


def IsLeapYear(year):
   if((year % 4) == 0):
      if((year % 100) == 0):
         if( (year % 400) == 0):
            return 1
         else:
            return 0
      else:
         return 1
   return 0


if __name__ == "__main__":
    getData(sys.argv[1:])



