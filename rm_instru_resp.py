#encoding=utf8
#!/usr/bin/python

import os
import glob
import subprocess
import shutil
import pdb
def MakeDir(nameDir):
    try:
        os.makedirs(nameDir)
    except:
        pass

def findPZ(path):
    list = os.listdir(path)
    for i in range(1,len(list)):
        if list[i].find('PZ')>0:
            return list[i]


def trans(net,sta,chan):
    year = '2014'
    net_sta = net+'_'+sta
    #path2 = 'new/'+net_sta
    os.putenv("SAC_DISPLAY_COPYRIGHT","0")
    p = subprocess.Popen(['sac'],stdin=subprocess.PIPE)
    trans = "transfer from polezero s SAC_PZs_ALL to VEL freqlimits 0.005 0.01 9.5 9.75 \n"
    s = ""
    s += "echo on\n"
    file_dir =  chan + '/' + year + '/' + net_sta
    MakeDir(file_dir)
    for filename in glob.glob(file_dir+'/*.SAC'):
        s += "r %s \n" % (filename)
        s += "rglitches\n"
        #s += "taper type cosine width 0.015\n"
        s += "rmean; rtr; taper \n"
        s += trans
        s += "interpolate delta 0.05 \n"
        s += "w %s\n" % ('new'+'/'+filename)
        #s += "cd ..\n"
    
    s += "quit \n"
    p.communicate(s.encode())


STAList = ['BPH01','BPH02','BPH03','BPH04','BPH05','BPH06','BPH07','BPH08','BPH09','BPH10','BPH11','BPH12','BPH13'] #'PAS','USC','GR2','SOT'
for Sta in STAList:
    for chan in ['BHE','BHN','BHZ']:
        trans('PY',Sta,chan)




