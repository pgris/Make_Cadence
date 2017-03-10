import matplotlib.pyplot as plt
import lsst.sims.maf.metricBundles as metricBundles
#import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.db as db
import lsst.sims.maf.utils as utils
from MyMetric import AnaMetric
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-z", "--zmin", type="float", default=0.01, help="filter [%default]")
parser.add_option("-Z", "--zmax", type="float", default=0.5, help="filter [%default]")
parser.add_option("-N", "--nevts", type="int", default=10, help="filter [%default]")
parser.add_option("-m", "--model", type="string", default='salt2-extended', help="filter [%default]")
parser.add_option("-v", "--version", type="string", default='1.0', help="filter [%default]")
parser.add_option("-f", "--fieldname", type="string", default='DD', help="filter [%default]")
parser.add_option("-i", "--fieldid", type="float", default=290, help="filter [%default]")
parser.add_option("-o", "--opsimrun", type="string", default='minion_1016', help="filter [%default]")
parser.add_option("-r", "--runtype", type="string", default='Simulation', help="filter [%default]")
parser.add_option("-s", "--season", type="int", default=-1, help="filter [%default]")
parser.add_option("-t", "--sntype", type="string", default='Ia', help="filter [%default]")
parser.add_option("-n","--nrolling", type="int", default='3', help="filter [%default]")
parser.add_option("-p","--percent_merge", type="int", default='80', help="filter [%default]")

opts, args = parser.parse_args()

outDir ='Test'

#dbFile = '/data/pgris/sims_operation/Run_OpSim/enigma_1189_sqlite.db'
dbFile = '/sps/lsst/data/dev/pgris/sims_operations/DB_Files/'+opts.opsimrun+'_sqlite.db'
#dbFile = '/data/pgris/sims_operation/Run_OpSim/clrlsstsrv_1068_sqlite.db'
opsimdb = utils.connectOpsimDb(dbFile)
resultsDb = db.ResultsDb(outDir=outDir)

propinfo, proptags = opsimdb.fetchPropInfo()
print 'hello',proptags,propinfo

#field='DD'
#numID=744
#field='WFD'


#metric=metrics.SupernovaMetric(m5Col='fiveSigmaDepth', redshift=0.1, resolution=5.,Nbetween=20)
metric=AnaMetric(m5Col='fiveSigmaDepth',zmin=opts.zmin,zmax=opts.zmax,Nevts=opts.nevts,model=opts.model,version=opts.version,fieldname=opts.fieldname,fieldID=opts.fieldid,opsimrun=opts.opsimrun,runtype=opts.runtype,season=opts.season,sntype=opts.sntype,nrolling=opts.nrolling,percent_merge=opts.percent_merge)
#slicer = slicers.HealpixSlicer(nside=256)


slicer=slicers.OpsimFieldSlicer()
sqlconstraint = utils.createSQLWhere(opts.fieldname, proptags)
mb = metricBundles.MetricBundle(metric, slicer, sqlconstraint)

mbD = {0:mb}


mbg =  metricBundles.MetricBundleGroup(mbD, opsimdb,
                                      outDir=outDir, resultsDb=resultsDb)

mbg.runAll()


#mbg.plotAll(closefigs=False)
#plt.show()
