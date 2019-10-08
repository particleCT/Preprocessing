import os,sys,glob

FileList = glob.glob("/ion/pCT_data/raw_data/CDH17/Linepair*.dat") #Water_003*
print FileList
outputDir = "/ion/home/hect/CDH17/Linepair/"
os.system("mkdir -p "+outputDir) 
os.chdir("outputDir")
for filename in FileList:
    command = "/ion/home/hect/Preprocessing/PreProcLenny/bin/pCT_Preprocessing "+filename
    print command
    os.system(command)
    print "--------------------"


