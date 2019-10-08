import os,glob
import sys
import random
import numpy as np


FileList = glob.glob("/ion/pCT_data/raw_data/CPC_2016_08_13/HIT_Dec/"+sys.argv[1]+"*.dat")
logList = glob.glob("/ion/pCT_data/raw_data/CPC_2016_08_13/HIT_Dec/"+sys.argv[1]+"*.log")
FileList.sort()
logList.sort()
print FileList
print logList

for i in range(0,len(FileList)):
    print FileList[i],logList[i]

    outputDir = "/ion/home/hect/Data/CPC16/Processed/"+sys.argv[1]+str(i)
    os.system("mkdir -p "+outputDir)

    conf_temp = open("/ion/home/hect/Preprocessing/PreProcLenny/pCT_config_temp.txt")
    conf = open("/ion/home/hect/Preprocessing/PreProcLenny/pCT_config.txt","w")
    for lines in conf_temp.readlines():
	if(lines.count("$OUTPUTDIR")):
	    newLine = "outputDir = "+outputDir+" # Leave as '.' only if you want all the output to end up in the current directory"
	    conf.write(newLine+"\n")
	elif(lines.count("$LOGFILE")):
	    newLine = "log = "+logList[i]+" # Specify a path to the DAQ log file, to be used in continuous scans to find the starting angle, NULL if stepped scan"
	    conf.write(newLine+"\n")
	    print newLine	
	else:
	    conf.write(lines)
    os.system("cp /ion/home/hect/Preprocessing/PreProcLenny/pCT_config.txt "+outputDir)	

    preproc = "/ion/home/hect/Preprocessing/PreProcLenny/bin/pCT_Preprocessing "+FileList[i]
    sub_temp=open("/ion/home/hect/Preprocessing/PreProcLenny/submit_PreProc_temp.sh")
    sub=open("/ion/home/hect/Preprocessing/PreProcLenny/submit_PreProc.sh","w")

    for lines in sub_temp.readlines():           
	if(lines.count("$COMMANDLINE")):
            sub.write(preproc+" &\n")
            print preproc
        else:
            sub.write(lines)

    print "-----------------------------------------------------"
    sub_temp.close()
    sub.close()
    os.system("cp /ion/home/hect/Preprocessing/PreProcLenny/submit_PreProc.sh "+outputDir)
    os.system("qsub /ion/home/hect/Preprocessing/PreProcLenny/submit_PreProc.sh")
    
