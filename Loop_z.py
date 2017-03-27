import numpy as np
from optparse import OptionParser
import os
import time

parser = OptionParser()
parser.add_option("-z", "--zmin", type="float", default=0.0, help="filter [%default]")
parser.add_option("-Z", "--zmax", type="float", default=1.2, help="filter [%default]")
parser.add_option("-N", "--nevts", type="int", default=10, help="filter [%default]")
parser.add_option("-m", "--model", type="string", default='salt2-extended', help="filter [%default]")
parser.add_option("-v", "--version", type="string", default='', help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default='WFD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="int", default=309, help="filter [%default]")
parser.add_option("-o", "--opsimrun", type="string", default='minion_1016', help="filter [%default]")
parser.add_option("-r", "--runtype", type="string", default='Simulation', help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=-1, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default='Ia', help="filter [%default]")
parser.add_option("-n","--nrolling", type="int", default='3', help="filter [%default]")
parser.add_option("-p","--percent_merge", type="int", default='80', help="filter [%default]")
parser.add_option("-d", "--dbFile", type="string",default='None', help="dbFile to process [%default]")

opts, args = parser.parse_args()


def cmd(zmin,zmax,nevts,fieldname,fieldid,season,runtype,sntype,nrolling,percent,dbfile):
    cmd='python Ana_Metrics.py --zmin '+str(zmin)+' --zmax '+str(zmax)+' --nevts '+str(nevts)+' --fieldname '+fieldname+' --fieldid '+str(fieldid)+' --season '+str(season)+' --runtype '+runtype+' --sntype '+sntype+' --dbFile '+dbfile
    if runtype.count('Rolling'):
        cmd+=' --nrolling '+str(nrolling)+' --percent_merge '+str(percent)
    return cmd

if opts.runtype.count('Simulation') > 0 or opts.runtype.count('Rolling') > 0:

    cwd = os.getcwd()
    dirScript= cwd + "/scripts"

    if not os.path.isdir(dirScript) :
        os.makedirs(dirScript)
        
    dirLog = cwd + "/logs"
    if not os.path.isdir(dirLog) :
        os.makedirs(dirLog)    

    zstep=0.1
    
    diffz=opts.zmax-opts.zmin

    common_info=opts.runtype+'_'+opts.sntype+'_'+opts.fieldname+'_'+str(opts.fieldid)

    print 'hello ?',zstep,diffz
    if diffz >= 0.99*zstep:
        print 'pass'    
        zrange=np.arange(opts.zmin,opts.zmax,zstep)
        for i in range(len(zrange)-1):
           
            zmin=(zrange[i] if zrange[i]>0. else zrange[i]+0.01)
            zmax=zrange[i+1]
            mycm=cmd(zmin,zmax,opts.nevts,opts.fieldname,opts.fieldid,opts.season,opts.runtype,opts.sntype,opts.nrolling,opts.percent_merge,opts.dbFile)
            print mycm

            name_id=common_info+'_'+str(zmin)+'_'+str(zmax)+'_'+str(opts.nevts)+'_season_'+str(opts.season)
            log = dirLog + '/'+name_id+'.log'
            qsub = "qsub -P P_lsst -l sps=1,ct=120000,h_vmem=16G -j y -o "+ log + " <<EOF"
            scriptName = dirScript+'/'+name_id+'.sh'

            script = open(scriptName,"w")
            script.write(qsub + "\n")
            script.write("#!/usr/local/bin/bash\n")
            script.write(" cd " + cwd + "\n")
            script.write("bash " + "\n")
            #script.write("ls /usr/local/grid/emi-3/WN/SL6_64/3.10.0-1.2/usr/lib64/" + "\n")
            script.write(" source setups.sh\n")
            script.write(mycm+" \n")
            script.write("EOF" + "\n")
            script.close()
    #os.system("chmod +x " + scriptName)
            os.system("sh "+scriptName)
            #os.system(mycm)
            time.sleep(3)
        
        zmin=zrange[len(zrange)-1]
        zmax=zrange[len(zrange)-1]+zstep
        mycm=cmd(zmin,zmax,opts.nevts,opts.fieldname,opts.fieldid,opts.season,opts.runtype,opts.sntype,opts.nrolling,opts.percent_merge,opts.dbFile)
        print mycm
        name_id=common_info+'_'+str(zmin)+'_'+str(zmax)+'_'+str(opts.nevts)+'_season_'+str(opts.season)
        log = dirLog + '/'+name_id+'.log'
        qsub = "qsub -P P_lsst  -l sps=1,ct=160000,h_vmem=16G -j y -o "+ log + " <<EOF"
        scriptName = dirScript+'/'+name_id+'.sh'
        
        script = open(scriptName,"w")
        script.write(qsub + "\n")
        script.write("#!/usr/local/bin/bash\n")
        script.write(" cd " + cwd + "\n")
        script.write("bash " + "\n")
        #script.write("ls /usr/local/grid/emi-3/WN/SL6_64/3.10.0-1.2/usr/lib64/" + "\n")
        script.write(" source setups.sh\n")
        script.write(mycm+" \n")
        script.write("EOF" + "\n")
        script.close()
        os.system("sh "+scriptName)

        #os.system(mycm)

if opts.runtype == 'Observation':
    
    cmd='python Ana_Metrics.py --fieldname '+opts.fieldname+' --fieldid '+str(opts.fieldid)+' --runtype '+opts.runtype
    os.system(cmd)
