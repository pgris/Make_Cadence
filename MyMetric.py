import numpy as np
from lsst.sims.maf.metrics import BaseMetric

from lsst.sims.photUtils import SignalToNoise
from lsst.sims.photUtils import PhotometricParameters
import sncosmo
from astropy.table import vstack,Table
import astropy.units as u
import matplotlib.pyplot as plt
from lsst.sims.photUtils.EBV import EBVbase
from lsst.sims.photUtils import Bandpass
import math
import time

from SN_Rates import SN_Rate
from cosmology import *
from SN_Object import SN_Object
from Throughputs import Throughputs
#import pickle as pkl
import cPickle as pkl
import glob
import os


class AnaMetric(BaseMetric):
    """
    Measure how many time series meet a given time and filter distribution requirement.
    """
    """
    def __init__(self, metricName='AnaMetric',
                 mjdCol='expMJD', filterCol='filter', m5Col='fiveSigmaDepth',
                 units='', redshift=0.,
                 Tmin = -20., Tmax = 60., Nbetween=7, Nfilt={'u':0,'g':10,'r':10,'i':10,'z':5,'y':5}, Nfilt_obs=5,Tless = -5., Nless=1,
                 Tmore = 30., Nmore=1, peakGap=2., snrCut=10., singleDepthLimit=23.,
                 resolution=5., badval=-666,
                 uniqueBlocks=False, **kwargs):
    """
    def __init__(self, metricName='AnaMetric',mjdCol='expMJD', filterCol='filter', m5Col='fiveSigmaDepth',units='', badval=-666,uniqueBlocks=False,zmin=0.01,zmax=0.5,Nevts=10, model='salt2-extended',version='1.0',fieldname='DD',fieldID=290,opsimrun='minion_1016',snrate='Flat',runtype='Simulation',season=-1,sntype='Ia',**kwargs):
    
       """
        
        Tmin = the minimum day to consider the SN.
        Tmax = the maximum to consider.
        Nbetween = the number of observations to demand between Tmin and Tmax
        Nfilt = number of unique filters that must observe the SN above the snrCut
        Tless = minimum time to consider 'near peak'
        Tmore = max time to consider 'near peak'
        Nless = number of observations to demand before Tless
        Nmore = number of observations to demand after Tmore
        peakGap = maximum gap alowed between observations in the 'near peak' time
        snrCut = require snr above this limit when counting Nfilt XXX-not yet implemented
        singleDepthLimit = require observations in Nfilt different filters to be this
        deep near the peak.  This is a rough approximation for the Science Book
        requirements for a SNR cut.  Ideally, one would import a time-variable SN SED,
        redshift it, and make filter-keyed dictionary of interpolation objects so the
        magnitude of the SN could be calculated at each observation and then use the m5col
        to compute a SNR.
        resolution = time step (days) to consider when calculating observing windows
        uniqueBlocks = should the code count the number of unique sequences that meet
        the requirements (True), or should all sequences that meet the conditions
        be counted (False).

        The filter centers are shifted to the SN restframe and only observations
        with filters between 300 < lam_rest < 900 nm are included

        In the science book, the metric demands Nfilt observations above a SNR cut.
        Here, we demand Nfilt observations near the peak with a given singleDepthLimt.
        """
       self.mjdCol = mjdCol
       self.m5Col = m5Col
       self.filterCol = filterCol
       self.dateCol = 'expDate'
       self.fieldRA='fieldRA'
       self.fieldDec='fieldDec'
       self.fieldID='fieldID'
       self.ditheredRA='ditheredRA'
       self.ditheredDec='ditheredDec'
       self.visitTime='visitExpTime'
       self.finSeeing='finSeeing'
       self.rawSeeing='rawSeeing'
       self.moonPhase='moonPhase'
       self.airmass='airmass'
       self.filtSkyBrightness='filtSkyBrightness'

       self.zmin=zmin
       self.zmax=zmax
       self.Nevts=Nevts
       self.snrate=snrate
       self.model=model
       self.version=version
       self.fieldName=fieldname
       self.fieldID_ref=int(fieldID)
       self.outputdir='Sim_'+opsimrun
       self.runtype=runtype
       self.season=season
       self.sntype=sntype

       self.time_begin=time.time()
       if not os.path.exists(self.outputdir):
           os.makedirs(self.outputdir)

       #super(AnaMetric, self).__init__(col=[self.mjdCol, self.m5Col, self.filterCol, self.dateCol,self.fieldRA,self.fieldDec, self.ditheredRA,self.ditheredDec,self.visitTime,self.finSeeing,self.rawSeeing,self.moonPhase,self.airmass,self.filtSkyBrightness,self.fieldID],
       super(AnaMetric, self).__init__(col=[self.mjdCol, self.m5Col, self.filterCol, self.dateCol,self.fieldRA,self.fieldDec, self.ditheredRA,self.ditheredDec,self.visitTime,self.rawSeeing,self.moonPhase,self.airmass,self.filtSkyBrightness,self.fieldID],
                                              metricName=metricName, units=units, badval=badval,
                                              **kwargs)
       """
       self.redshift = redshift
       self.Tmin = Tmin
       self.Tmax = Tmax
       self.Nbetween = Nbetween
       self.Nfilt = Nfilt
       self.Nfilt_obs = Nfilt_obs
       self.Tless = Tless
       self.Nless = Nless
       self.Tmore = Tmore
       self.Nmore = Nmore
       self.peakGap = peakGap
       self.snrCut = snrCut
       self.resolution = resolution
       """
       self.uniqueBlocks = uniqueBlocks
       self.filterNames = np.array(['u','g','r','i','z','y'])
        # Set rough values for the filter effective wavelengths.
       self.singleDepthLimit = -1
       """
       self.Cm = {'u':23.60,
                  'g':24.57,
                  'r':24.57,
                  'i':24.47,
                  'z':24.19,
                  'y':23.45,
                  'y1':23.74}
       
       self.kAtm = {'u':0.48,
                    'g':0.21,
                    'r':0.10,
                    'i':0.07,
                    'z':0.06,
                    'y':0.06}
       """

       #Values taken from sims_operations/tools/opsimUtils.py
       """
       self.Cm = {'u': 22.74,
                  'g': 24.38,
                  'r': 24.43,
                  'i': 24.30,
                  'z': 24.15,
                  'y': 23.70}

       self.dCm_infinity = {'u': 0.75,
                            'g': 0.19,
                            'r': 0.10,
                            'i': 0.07,
                            'z': 0.05,
                            'y': 0.04}

       """                  
       #Values taken from sims_operations/tools/schema_tools/gen_output.py 

       self.kAtm = {'u': 0.50,
                    'g': 0.21,
                    'r': 0.13,
                    'i': 0.10,
                    'z': 0.07,
                    'y': 0.18} 
       

       self.msky = {'u': 22.95,
               'g': 22.24,
               'r': 21.20,
               'i': 20.47,
               'z': 19.60,
               'y': 18.63} 

       self.Cm = {'u':22.94,
             'g':24.46,
             'r':24.48,
             'i':24.34,
             'z':24.18,
             'y':23.73}

       self.dCm_infinity = {'u':0.56,
                       'g':0.12,
                       'r':0.06,
                       'i':0.05,
                       'z':0.03,
                       'y':0.02}
      

       #FWHM_500 = seeing at 500 nm
       # FWHM_Sys_Zenith = sqrt(telSeeing**2 + opticalDesSeeing**2 + cameraSeeing**2)
       # Filter_Wavelength_Correction = (500 nm / Filter_Effective_Wavelength)**0.3
       # Airmass_Correction = airmass**0.6
       # FWHM_Sys = FWHM_Sys_Zenith * Airmass_Correction
       # FWHM_Atm = FWHM_500 * Filter_Wavelength_Correction * Airmass_Correction
       # FWHM_Eff = scaleToNeff * sqrt(FWHM_Sys**2 + atmNeffFactor * FWHM_Atm**2)
       # FWHM_Eff is the value in ObsHistory.finSeeing for the observations filter
       #
       # Units = unitless, Format = float, no default
       #

       self.telSeeing = 0.250 # design goal
       self.opticalDesSeeing = 0.08
       self.cameraSeeing = 0.300
       # Scaling factors for above seeing calculation
       self.scaleToNeff = 1.16
       self.atmNeffFactor = 1.04
       self.FWHM_Sys_Zenith = np.sqrt(self.telSeeing**2 + self.opticalDesSeeing**2 + self.cameraSeeing**2)

       self.filterWave = {'u': 367.0, 'g': 482.5, 'r': 622.2, 'i': 754.5, 'z': 869.1, 'y': 971.0}

    """     
    def run(self, dataSlice, slicePoint=None):
         print 'fieldIds',dataSlice[self.fieldID][0]
    """ 

    def run(self, dataSlice, slicePoint=None):

        if self.runtype=='Observation':
            self.run_observation(dataSlice,slicePoint)
        else:
            if self.runtype.count('Simulation') >0:
                if self.runtype.count('Rolling') >0:
                    self.rolling=True
                else:
                    self.rolling=False
                self.run_simulation(dataSlice,slicePoint, self.sntype) 
       

    def run_observation(self, dataSlice, slicePoint=None):
             
        dataSlice = dataSlice[np.where(dataSlice[self.fieldID]==self.fieldID_ref)]
        if len(dataSlice) > 0:
            dataSlice.sort(order=self.mjdCol)
            
            #print 'hello',dataSlice[self.fieldRA],len(dataSlice[self.fieldRA]),set(dataSlice[self.fieldRA]),len(set(dataSlice[self.fieldRA]))
            ra_field=list(set(dataSlice[self.fieldRA]))[0]
            dec_field=list(set(dataSlice[self.fieldDec]))[0]

            lsstmwebv = EBVbase()
            ebvofMW = lsstmwebv.calculateEbv(
            equatorialCoordinates=np.array([[ra_field], [dec_field]]))[0]

            dictout={}

            dictout['ebvofMW']=ebvofMW

            dictout['dataSlice']=dataSlice

            outdir=self.outputdir.replace('Sim','Obs')
            #name_for_pkl=outdir+'/Observations_'+self.fieldName+'_'+str(int(self.fieldID_ref))+'_'+str(ra_field)+'_'+str(dec_field)
            name_for_pkl=outdir+'/Observations_'+self.fieldName+'_'+str(int(self.fieldID_ref))
       
            pkl_file = open(name_for_pkl+'.pkl','wb')
            
            pkl.dump(dictout, pkl_file)
            
            pkl_file.close()
        """
        for filtre in self.filterNames:
            obs_filtre=dataSlice[np.where(dataSlice['filter']==filtre)]
            for obs in obs_filtre:
                print obs['expMJD'],obs['fiveSigmaDepth']
            
        """

    def run_simulation(self, dataSlice, slicePoint=None, sn_type='Ia'):
       
          
        if dataSlice[self.fieldID][0]==self.fieldID_ref:

            addit=''
            if self.rolling==True:
                addit='Rolling_'


            TSeason_min=-1
            TSeason_max=-1

            if self.season > -1:
                sfile=open('Seasons/Seasons_'+addit+self.fieldName+'_'+str(self.fieldID_ref)+'.txt', 'r')
                for line in sfile.readlines():
                    if int(line.split(' ')[0]) == self.season:
                        TSeason_min=float(line.split(' ')[1])
                        TSeason_max=float(line.split(' ')[2])
                        break

                if TSeason_min==-1 and TSeason_max==-1:
                    print 'Big problem : season not found - Stop'
                    return{}


            #range around the max
            lowrange=-30
            highrange=50

            # Read throughputs
            
            transmission=Throughputs()

            # Register LSST band pass (system) in sncosmo

            for filtre in self.filterNames:
                band=sncosmo.Bandpass(transmission.lsst_system[filtre].wavelen, transmission.lsst_system[filtre].sb, name='LSST::'+filtre,wave_unit=u.nm)
                sncosmo.registry.register(band, force=True)
        
            zmin=self.zmin
            zmax=self.zmax
            zstep=0.001
            
             #Load x1 and c asymmetric distributions
            # values from Scolnic & Kessler may 2016 arXiv:1603.01559v2
            x1_mean=0.964
            sig_m_x1=1.467
            sig_p_x1=0.235
            c_mean=-0.099
            sig_m_c=0.003
            sig_p_c=0.119

            if zmax <= 0.1:
               x1_mean=0.419
               sig_m_x1=3.024
               sig_p_x1=0.742
               c_mean=-0.069
               sig_m_c=0.003
               sig_p_c=0.148 

            x1_vals, x1_weights = self.gauss_asym_distrib(x1_mean,sig_m_x1,sig_p_x1)
            c_vals, c_weights = self.gauss_asym_distrib(c_mean,sig_m_c,sig_p_c)
           
            #fieldID_ref= np.array([309])
            #fieldID_ref= np.array([290])
            #print "attention",type(dataSlice)
            """
            goodFilters_f = np.in1d(dataSlice['filter'],self.filterNames)
            goodFilters = np.in1d(goodFilters_f[self.fieldID],fieldID_ref)
            dataSlice = dataSlice[goodFilters]
            """

            dataSlice = dataSlice[np.where(dataSlice[self.fieldID]==self.fieldID_ref)]
            if self.rolling == True:
                dataSlice = pkl.load(open('Rolling_Cadence_'+ self.fieldName+'_'+str(self.fieldID_ref)+'.pkl','rb'))

            if dataSlice.size == 0:
            #print 'Data slice sizee',dataSlice.size
                return (self.badval, self.badval,self.badval)

            #sort dataslice depending on date

            dataSlice.sort(order=self.mjdCol)
   
            """
            print 'dataslice',dataSlice['fieldRA'],dataSlice['fieldDec'],dataSlice['expMJD']
            time = dataSlice[self.mjdCol]-dataSlice[self.mjdCol].min()
            """

            """
            for data in dataSlice: 
            myfilter=data['filter']
            print 'LSST::%s %f %f %f -1 %f %f %f %f %f %f' % (data['filter'],data[self.fieldID],data['expMJD'], data[self.visitTime], data[self.finSeeing],data[self.moonPhase],data[self.filtSkyBrightness],self.kAtm[myfilter],data[self.airmass],data[self.m5Col])
            """

        
            m_begin_date=dataSlice['expMJD'].min()
            m_end_date=dataSlice['expMJD'].max()

            if self.season > -1:
                m_begin_date=TSeason_min
                m_end_date=TSeason_max
 
            #print 'hello',m_begin_date,m_end_date
            rate_SN=SN_Rate(GeneralCosmo(0.27, 0.73, -1.,0),zmin,zmax,zstep,m_end_date-m_begin_date,3.6,self.snrate)
            
            N_sn=int(np.sum(rate_SN.hz_weight))
            print 'About to generate',N_sn,'supernovae'
            
            weights=rate_SN.hz_weight/np.sum(rate_SN.hz_weight)


            #print 'time',dataSlice['expMJD']  

            outdict={}

            ra_field=dataSlice[self.fieldRA][0]
            dec_field=dataSlice[self.fieldDec][0]

            N_sn=self.Nevts

            #name_for_pkl='SuperNova_'+self.fieldName+'_'+str(self.fieldID_ref)+'_'+str(ra_field)+'_'+str(dec_field)+'_'+str(zmin)+'_'+str(zmax)+'_'+self.model
           
            name_for_pkl='SuperNova_'+addit+sn_type+'_'+self.fieldName+'_'+str(self.fieldID_ref)+'_'+str(zmin)+'_'+str(zmax)

            """
            if self.version != '':
                name_for_pkl+='_'+self.version
            """

            name_for_pkl+='_'+str(N_sn)

            if self.season > -1:
                name_for_pkl+='_season_'+str(self.season)

            name_for_pkl=self.outputdir+'/'+name_for_pkl
            num=len(glob.glob(name_for_pkl+'*.pkl'))
            pkl_file = open(name_for_pkl+'_'+str(num)+'.pkl','wb')

            #print 'hello',zmin,zmax
            """
            self.dict_out={}
            self.dict_out['t0']=[]
            self.dict_out['X1']=[]
            self.dict_out['c']=[]
            self.dict_out['wave']=[]
            self.dict_out['fluxes']=[]
            """

            #N_sn=1
            for i in range(0,N_sn):
            
                outdict={}

                T0 = np.random.uniform(m_begin_date,m_end_date)
               
                c=np.random.choice(c_vals,1,p=c_weights)[0]
                x1=np.random.choice(x1_vals,1,p=x1_weights)[0]

                
                 #print 'alors c x1',c,X1,c_weights
                """
                z=0.04250000000000000
                T0= 51595.826132848888
                #x0      = 0.00252041986763165
                X1= 1.2478457417636148
                c = 0.0059916789377777759
                """

                timelow = T0+ lowrange
                timehigh = T0 + highrange
                
                #T0=0.5*(m_begin_date+m_end_date)
            
                z=np.random.choice(rate_SN.hz,1,p=weights)[0]

                """
                T0=0.5*(TSeason_min+TSeason_max)
                c=0
                x1=0
                z=0.15
                """

                outdict['t0']=T0
                outdict['c']=c
                outdict['x1']=x1
                outdict['z']=z
                outdict['ra']=ra_field
                outdict['dec']=dec_field
                #print 'alors c x1',c,X1,z,T0
                """
                outdict['status']='unkown'
                outdict['Obs_LC']=None
                outdict['Obs_LC_for_fit']=None
                outdict['sncosmo_res']=None
                outdict['sncosmo_fitted']={}
                """
                outdict['status']='unkown'
                outdict['fit']=None
                outdict['mbsim']=-999.
                outdict['observations']=None

            #print 'redshift',z

                phase_of_first_point = (m_begin_date-T0) / (1+z)
                phase_of_last_point = (m_end_date-T0) / (1+z)

                if phase_of_first_point > -5. or phase_of_last_point < 20.: #days
                    #print '[KILL]',T0,z,phase_of_first_point,phase_of_last_point
                    outdict['status']='killed'
                    pkl.dump(outdict, pkl_file)
                    continue

                observations=dataSlice[np.where(np.logical_and(dataSlice['expMJD']>timelow,dataSlice['expMJD']<timehigh))]
                if len(observations) > 0:
                    ra=observations[self.fieldRA][0]
                    dec=observations[self.fieldDec][0]
                    #print 'This is T0',T0,c,x1
                    self.SN=SN_Object(ra=np.rad2deg(ra),dec=np.rad2deg(dec),z=z,t0=T0,c=c,x1=x1,model=self.model,version=self.version,sn_type=sn_type)
                    outdict['sn_type']=self.SN.sn_type
                    outdict['sn_model']=self.SN.model
                    outdict['sn_version']=self.SN.version
                    self.thetime=time.time()
                    #print 'before fitting',self.thetime
                    dict_fit,mbsim,myobs=self.Simulate_and_Fit_LC(observations,transmission,zmin,zmax)
                    #print 'after fitting',time.time()-self.thetime
                    outdict['fit']=dict_fit
                    outdict['mbsim']=mbsim
                    outdict['observations']=myobs
                    outdict['status']='try_fit'
                    #print 'dumping in pkl',i
                    pkl.dump(outdict, pkl_file)

                else:
                #print 'No obs'
                    outdict['status']='No obs in [T0'+str(int(lowrange))+';T0+'+str(int(highrange))+']'
                    pkl.dump(outdict, pkl_file)


                """
                #print 'fit status',stat
                    if res == None:
                    #print 'not fitted'
                        outdict['status']='Not fitted'
                        pkl.dump(outdict, pkl_file)
                        continue
                    else:
                    #print 'resultat',res
                    #print 'model',fitted_model.param_names
                        outdict['sncosmo_res']=res
                        for i,par in enumerate(fitted_model.param_names):
                            outdict['sncosmo_fitted'][par]=fitted_model.parameters[i]
                        outdict['status']='fitted'
                        pkl.dump(outdict, pkl_file)
                    #break
                """
            

            pkl_file.close()

        #print 'Finished',(time.time()-self.time_begin)/60.
        return{}
   
    def Simulate_and_Fit_LC(self,observations,transmission,zmin,zmax):

        #print 'time',observations['expMJD'],observations['filter']
 
        #print 'simulate and fit'
        ra=observations[self.fieldRA][0]
        dec=observations[self.fieldDec][0]

        
        if self.SN.sn_type=='Ia':
            mbsim=self.SN.SN._source.peakmag('bessellb','vega')
        else:
            mbsim=-1


        #This will be the data for sncosmo fitting
        table_for_fit={}
        table_for_fit['error_calc'] = Table(names=('time','flux','fluxerr','band','zp','zpsys'), dtype=('f8', 'f8','f8','S7','f4','S4'))
        table_for_fit['error_coadd_calc'] = Table(names=('time','flux','fluxerr','band','zp','zpsys'), dtype=('f8', 'f8','f8','S7','f4','S4'))
        table_for_fit['error_opsim'] = Table(names=('time','flux','fluxerr','band','zp','zpsys'), dtype=('f8', 'f8','f8','S7','f4','S4'))
        table_for_fit['error_coadd_opsim'] = Table(names=('time','flux','fluxerr','band','zp','zpsys'), dtype=('f8', 'f8','f8','S7','f4','S4'))


        mytype=[('obsHistID',np.int), ('filtSkyBrightness', np.float), ('airmass', np.float), ('moonPhase', np.float), ('fieldRA', np.float), ('fieldDec', np.float), ('visitExpTime', np.float), ('expDate', np.int), ('filter',np.dtype('a15')), ('fieldID', np.int), ('fiveSigmaDepth', np.float), ('ditheredDec', np.float), ('expMJD', np.float), ('ditheredRA',np.float), ('rawSeeing', np.float),('flux',np.float),('err_flux',np.float),('err_flux_opsim',np.float),('finSeeing',np.float),('katm_opsim',np.float),('katm_calc',np.float),('m5_calc',np.float),('Tb',np.float),('Sigmab',np.float),('Cm',np.float),('dCm',np.float),('mag_SN',np.float),('snr_m5_through',np.float),('snr_m5_opsim',np.float),('gamma_through',np.float),('gamma_opsim',np.float),('snr_SED',np.float)]

        myobservations=np.zeros((60,1),dtype=mytype)

        #print 'Nobservations',len(observations)
        nobs=-1
        for filtre in self.filterNames:
            obs_filtre=observations[np.where(observations['filter']==filtre)]
            #print 'ehehe',obs_filtre
                              
            for obs in obs_filtre:

                nobs+=1

                if len(myobservations) <= nobs:
                    myobservations=np.resize(myobservations,(len(myobservations)+100,1))

                for name in observations.dtype.names:
                    myobservations[name][nobs]=obs[name]

                #print 'time uu',obs['expMJD']
                seeing=obs['rawSeeing']
                #seeing=obs['finSeeing']
                time_obs=obs['expMJD']
                m5_opsim=obs['fiveSigmaDepth']

                #print 'getting SED'
                sed_SN=self.SN.get_SED(time_obs)
                #print 'got SED',sed_SN.wavelen,sed_SN.flambda,obs['expMJD']

                

                """
                outf = open('SN_'+str(time)+'.dat', 'wb')
                for i,wave in enumerate(sn.SEDfromSNcosmo.wavelen):
                    print >> outf,wave,sn.SEDfromSNcosmo.flambda[i]
                outf.close()
                """
                #print 'loading transmission airmass'
                transmission.Load_Atmosphere(obs['airmass'])
                flux_SN=sed_SN.calcFlux(bandpass=transmission.lsst_atmos_aerosol[filtre])
                #print 'this is my flux',flux_SN
                #flux_SN=sed_SN.calcFlux(bandpass=transmission.lsst_system[filtre]) / 3631.0

                
                myup=transmission.darksky.calcInteg(transmission.lsst_system[filtre])
               
                """
                wavelen, sb = transmission.lsst_system[filtre].multiplyThroughputs(transmission.lsst_atmos[filtre].wavelen, transmission.lsst_atmos[filtre].sb)
                lsst_total= Bandpass(wavelen=wavelen, sb=sb)
                """
                Tb=self.Calc_Integ(transmission.lsst_atmos[filtre])
                Sigmab=self.Calc_Integ(transmission.lsst_system[filtre])
                katm=-2.5*np.log10(Tb/Sigmab)
                
                mbsky_through=-2.5*np.log10(myup/(3631.*Sigmab))
                
                #print 'there mbsky',filtre,mbsky_through,obs['filtSkyBrightness'],katm,self.kAtm[filtre],Tb,Sigmab,obs['airmass']

                Filter_Wavelength_Correction = np.power(500.0 / self.filterWave[filtre], 0.3)
                Airmass_Correction = math.pow(obs['airmass'],0.6)
                FWHM_Sys = self.FWHM_Sys_Zenith * Airmass_Correction
                FWHM_Atm = seeing * Filter_Wavelength_Correction * Airmass_Correction
                finSeeing = self.scaleToNeff * math.sqrt(np.power(FWHM_Sys,2) + self.atmNeffFactor * np.power(FWHM_Atm,2))


                #print 'hello pal',filtre,finSeeing,obs['visitExpTime']
                Tscale = obs['visitExpTime']/ 30.0 * np.power(10.0, -0.4*(obs['filtSkyBrightness'] - self.msky[filtre]))
                dCm = self.dCm_infinity[filtre] - 1.25*np.log10(1 + np.power(10.,0.8*self.dCm_infinity[filtre]- 1.)/Tscale)

                m5_recalc=dCm+self.Cm[filtre]+0.5*(obs['filtSkyBrightness']-21.)+2.5*np.log10(0.7/finSeeing)-self.kAtm[filtre]*(obs['airmass']-1.)+1.25*np.log10(obs['visitExpTime']/30.)

                myobservations['Cm'][nobs]=self.Cm[filtre]
                myobservations['dCm'][nobs]=dCm
                myobservations['finSeeing'][nobs]=finSeeing
                myobservations['Tb'][nobs]=Tb
                myobservations['Sigmab'][nobs]=Sigmab
                myobservations['katm_calc'][nobs]=katm
                myobservations['katm_opsim'][nobs]=self.kAtm[filtre]

                if flux_SN >0:
                 
                    #print 'positive flux',flux_SN
                    mag_SN=-2.5 * np.log10(flux_SN / 3631.0)
                    
                    FWHMeff = SignalToNoise.FWHMgeom2FWHMeff(finSeeing)
                    #FWHMeff = SignalToNoise.FWHMgeom2FWHMeff(seeing)
                    photParams = PhotometricParameters()
                    snr_SN= SignalToNoise.calcSNR_sed(sed_SN,transmission.lsst_atmos_aerosol[filtre], transmission.darksky, transmission.lsst_system[filtre], 
                                                      photParams, FWHMeff=FWHMeff, verbose=False)
                    m5_calc=SignalToNoise.calcM5(transmission.darksky,transmission.lsst_atmos_aerosol[filtre],transmission.lsst_system[filtre],photParams=photParams,FWHMeff=FWHMeff)
                    snr_m5_through,gamma_through=SignalToNoise.calcSNR_m5(mag_SN,transmission.lsst_atmos_aerosol[filtre],m5_calc,photParams)
                    snr_m5_opsim,gamma_opsim=SignalToNoise.calcSNR_m5(mag_SN,transmission.lsst_atmos_aerosol[filtre],m5_opsim,photParams)

                    #print 'm5 diff',filtre,m5_calc,m5_opsim,m5_calc/m5_opsim,m5_recalc,(m5_opsim/m5_recalc)
                    err_flux_SN=flux_SN/snr_SN
                    err_flux_SN_opsim=flux_SN/snr_m5_opsim

                    myobservations['mag_SN'][nobs]=mag_SN
                    myobservations['flux'][nobs]=flux_SN
                    myobservations['err_flux'][nobs]=err_flux_SN
                    myobservations['err_flux_opsim'][nobs]=err_flux_SN_opsim
                    myobservations['m5_calc'][nobs]=m5_calc
                    myobservations['snr_m5_through'][nobs]=snr_m5_through
                    myobservations['snr_m5_opsim'][nobs]=snr_m5_opsim
                    myobservations['gamma_through'][nobs]=gamma_through
                    myobservations['gamma_opsim'][nobs]=gamma_opsim
                    myobservations['snr_SED'][nobs]=snr_SN
                  
                    #print 'SNR',flux_SN,flux_SN/err_flux_SN,flux_SN/err_flux_SN_opsim
                    #if flux_SN/err_flux_SN >=5:
                    table_for_fit['error_calc'].add_row((time_obs,flux_SN,err_flux_SN,'LSST::'+filtre,25,'ab'))
                    #if flux_SN/err_flux_SN_opsim >=5.:
                    table_for_fit['error_opsim'].add_row((time_obs,flux_SN,err_flux_SN_opsim,'LSST::'+filtre,25,'ab'))
                    #print 'Getting fluxes and errors',time.time()-self.thetime,filtre,nobs
                else:
                    err_flux_SN=-999.
                    err_flux_SN_opsim=-999.
                    myobservations['mag_SN'][nobs]=-999
                    myobservations['flux'][nobs]=flux_SN
                    myobservations['err_flux'][nobs]=-999.
                    myobservations['err_flux_opsim'][nobs]=-999.
                    myobservations['m5_calc'][nobs]=-999. 
                    myobservations['snr_m5_through'][nobs]=-999
                    myobservations['snr_m5_opsim'][nobs]=-999
                    myobservations['gamma_through'][nobs]=-999
                    myobservations['gamma_opsim'][nobs]=-999
                    myobservations['snr_SED'][nobs]=-999

                #print 'flux SN',flux_SN,err_flux_SN,mag_SN,snr_SN,snr_m5_through,snr_m5_opsim

                #t.add_row((time,flux_SN,err_flux_SN,'LSST::'+filtre,0.,'vega'))
               
                #break
           

        #t = Table([list(data['time']), list(data['band']),data['flux'],data['fluxerr'],data['flux_aero'],data['fluxerr_aero'],data['fluxerr_new'],data['zp'],data['zpsys']], names=('time','band','flux','fluxerr','flux_aero','fluxerr_aero','fluxerr_new','zp','zpsys'), meta={'name': 'first table'})

        #model.set(z=0.5)
        #print SN.SN

        myobservations=np.resize(myobservations,(nobs+1,1))
        #print 'there obs',myobservations
        #print 'Getting coadds',time.time()-self.thetime
        for band in ['u','g','r','i','z','y']:
            sela=table_for_fit['error_calc'][np.where(table_for_fit['error_calc']['band']=='LSST::'+band)]
            #sela=sela[np.where(np.logical_and(sela['flux']/sela['fluxerr']>5.,sela['flux']>0.))]
            selb=table_for_fit['error_opsim'][np.where(table_for_fit['error_opsim']['band']=='LSST::'+band)]
            #selb=selb[np.where(np.logical_and(selb['flux']/selb['fluxerr']>5.,selb['flux']>0.))]

            table_for_fit['error_coadd_calc']=vstack([table_for_fit['error_coadd_calc'],self.Get_coadd(sela)])
            table_for_fit['error_coadd_opsim']=vstack([table_for_fit['error_coadd_opsim'],self.Get_coadd(selb)])

        #print 'There we go fitting',time.time()-self.thetime
        dict_fit={}
        #for val in ['error_calc','error_coadd_calc','error_opsim','error_coadd_opsim']:
        for val in ['error_coadd_calc','error_coadd_opsim']:
            dict_fit[val]={}
            dict_fit[val]['sncosmo_fitted']={}
            dict_fit[val]['table_for_fit']=table_for_fit[val]
            #print 'fit',val,time.time()-self.thetime
            res,fitted_model,mbfit,fit_status=self.Fit_SN(table_for_fit[val],zmin,zmax)
            if res is not None:
                dict_fit[val]['sncosmo_res']=res
            #self.dict_fit[val]['fitted_model']=fitted_model
                for i,par in enumerate(fitted_model.param_names):
                    dict_fit[val]['sncosmo_fitted'][par]=fitted_model.parameters[i]
                dict_fit[val]['mbfit']=mbfit
            dict_fit[val]['fit_status']=fit_status
       
        return dict_fit,mbsim,myobservations 


    def Fit_SN(self,t,zmin,zmax):


        select=t[np.where(np.logical_and(t['flux']/t['fluxerr']>5.,t['flux']>0.))]
        #selectb=select[np.where(select['band']=='LSST::g')]
        #print 'oh boy',len(select),len(selectb)
        #if len(select)>=1 and len(selectb)>=1:
        if len(select)>=5:
            #print 'data to be fitted',select
            try:
                
                
                #print 'before fit here'
                #print 'SN parameters',self.SN.SN
                
                #print 'fitting',select
                
                z_sim=self.SN.z
                #print 'hello z',z_sim
                #print 'fit it',val,time.time()-self.thetime
                res, fitted_model = sncosmo.fit_lc(select, self.SN.SN_fit_model,['z', 't0', 'x0', 'x1', 'c'],bounds={'z':(z_sim-0.01, z_sim+0.01)})
                #res, fitted_model = sncosmo.fit_lc(select, self.SN.SN,['t0', 'x0', 'x1', 'c'])
    

                #print 'after fit',res.keys()
                #print res.keys()
                
                """
                print 'after fit'
                print res['parameters'],res['errors']
                """

                mbfit=fitted_model._source.peakmag('bessellb','vega')
                #print 'oooo test',-2.5*np.log10(res['parameters'][2])+10.635,fitted_model.bandmag('bessellb','vega',res['parameters'][1]),mbsim,mbfit,mbsim-mbfit

                """
                sncosmo.plot_lc(t, model=fitted_model,color='k',pulls=False)
                
                plt.show()
                """
                return res,fitted_model,mbfit,'ok'
            
            except (RuntimeError, TypeError, NameError):
                #print 'crashed'
                return None,None,-1,'crash'

        else:
            return None,None,-1,'Noobs'
            


    def Calc_Integ(self,bandpass):
        resu=0.
        dlam=0
        for i,wave in enumerate(bandpass.wavelen):
            if i < len(bandpass.wavelen)-1:
                dlam=bandpass.wavelen[i+1]-wave
                resu+=dlam*bandpass.sb[i]/wave
            #resu+=dlam*bandpass.sb[i]

        return resu    

    def gauss_asym_distrib(self,mean,sigma_minus,sigma_plus):
        xmin=mean-5.*sigma_minus
        xmax=mean+5.*sigma_plus

        pas=1.e-4

        nsteps=int((xmax-xmin)/pas)

        xvals=[]
        weights=[]

        for i in range(nsteps):
            x=xmin+float(i)*pas
            if x < mean:
                res=np.exp(-np.power(x-mean,2)/(2*np.power(sigma_minus,2.)))
            else:
                res=np.exp(-np.power(x-mean,2)/(2*np.power(sigma_plus,2.)))

            xvals.append(x)
            weights.append(res)


        return xvals,weights/np.sum(weights)


    def Get_coadd(self,table_for_fit):
    #names_ref=('time','flux','fluxerr','band','zp','zpsys')
        dtype_ref=('f8', 'f8','f8','S7','f4','S4')
        names_ref=table_for_fit.colnames
    #dtype_ref=table_for_fit.dtype
        
    #print 'before filtering',table_for_fit
        out_table=Table(names=names_ref,dtype=dtype_ref)
        dict_for_coadd={}
        if len(table_for_fit) > 0:
            inum=0
            dict_for_coadd[inum]=Table(names=names_ref,dtype=dtype_ref)
        #print 'timediff',24.*60.*60.*(table_for_fit['time']-table_for_fit['time'][0])
                                
            iloop=0
        #print 'blablabla',dict_for_coadd[inum]
            dict_for_coadd[inum].add_row(table_for_fit[iloop])
                                
            if len(table_for_fit) > 1:
                while iloop < len(table_for_fit)-1:   
                    diff_time_sec=24.*60.*60.*(table_for_fit['time'][iloop+1]-table_for_fit['time'][iloop])
                #print 'alors ???',diff_time_sec,inum
                    if diff_time_sec > 40.:
                        inum+=1
                        dict_for_coadd[inum]=Table(names=names_ref,dtype=dtype_ref)
                
                    dict_for_coadd[inum].add_row(table_for_fit[iloop+1])
                    
                    iloop+=1
        #print 'thedict',dict_for_coadd

    #print 'coadding'
        for key,vals in dict_for_coadd.items():
        #print key,vals
            #out_table.add_row((np.mean(vals['time']),np.mean(vals['flux']),np.sqrt(np.sum(vals['fluxerr']*vals['fluxerr']))/np.sqrt(float(len(vals))),vals['band'][0],vals['zp'][0],vals['zpsys'][0]))
            mean_pond=np.sum(vals['flux']*(np.power(vals['fluxerr'],-2))/np.sum(np.power(vals['fluxerr'],-2)))
            sigsum=1./np.sqrt(np.sum(np.power(vals['fluxerr'],-2)))
            out_table.add_row((np.mean(vals['time']),mean_pond,sigsum,vals['band'][0],vals['zp'][0],vals['zpsys'][0]))


    #print 'after filtering',out_table

        return out_table

        

    def reduceMedianMaxGap(self,data):
        #The median maximum gap near the peak of the light curve
        """
        print 'hello pal',np.median(data['maxGap'])
        print 'passed'
        result = np.median(data['maxGap'])
        if np.isnan(result):
            result = self.badval
        return result
        """
        return -1
    """
    def reduceNsequences(self,data):
        #The number of sequences that met the requirements
        return data['result']
    def reduceMedianNobs(self,data):
        #Median number of observations covering the entire light curve
        result = np.median(data['Nobs'])
        if np.isnan(result):
            result = self.badval
        return result

    def reduceMedianNobs_peak(self,data):
        #Median number of observations covering the entire light curve
        result = np.median(data['Nobs_peak'])
        if np.isnan(result):
            result = self.badval
        return result
    """
