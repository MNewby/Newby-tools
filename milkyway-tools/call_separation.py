import subprocess as sp
import glob as glob

pfolder = "/home/newbym2/milkyway_parameters_sansSgr"
sfolder = "/home/newbym2/milkyway_stars_sansSgr"
ofolder = "/home/newbym2/milkyway_out_sansSgr"
params = glob.glob(pfolder+"/*sanSgr.txt")

for p in params:
    wedge = p[-24:-22]
    search = p[-13:-11]
    s = sfolder+"/stars-"+wedge+"-sansSgr.txt"
    o = ofolder+"/out-sansSgr-"+wedge+"-"+search+".txt"
    command = "./milkyway_separation -i -a "+p+" -s "+s+" -o "+o
    print command
    sts = sp.call(command, shell=True)
print "# --- Done"
#p-09-3s-free_de_sanSgr.txt
