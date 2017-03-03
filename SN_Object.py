import numpy as np
import sncosmo
from lsst.sims.photUtils.EBV import EBVbase
from lsst.sims.photUtils import Sed
from astropy import (cosmology, units as u, constants as const)
from cosmology import GeneralCosmo
from astropy.cosmology import FlatLambdaCDM
import astropy.constants as const
import astropy.units as u
from lsst.sims.photUtils import Bandpass
from sncosmo.io import read_griddata_ascii, read_griddata_fits
from scipy.interpolate import (InterpolatedUnivariateSpline as Spline1d,
                               RectBivariateSpline as Spline2d,
                               splmake, spleval)

class SN_Object:
    def __init__(self,ra,dec,z,t0,c=0,x1=0,peakAbsMagBesselB=-19.3000,model='salt2-extended',version='1.0',sn_type='Ia',mwdust=True):
        self.radeg=ra
        self.decdeg=dec
        self.z=z
        self.t0=t0
        self.c=c
        self.x1=x1
        self.peakAbsMagBesselB=peakAbsMagBesselB
        self.model=model
        self.version=version
        self.sn_type=sn_type
        self.id_SED=''
        self.cosmology=cosmology.WMAP9
        astropy_cosmo=FlatLambdaCDM(H0= self.cosmology.H0, Om0=self.cosmology.Om0)
        self.lumidist=astropy_cosmo.luminosity_distance(self.z).value*1.e6
        
        #print 'sntype',self.sn_type
        dust = sncosmo.OD94Dust()
        self.lsstmwebv = EBVbase()
        self.ebvofMW = self.lsstmwebv.calculateEbv(
            equatorialCoordinates=np.array([[np.radians(self.radeg)], [np.radians(self.decdeg)]]))[0]
        
        if self.sn_type== 'Ia':

            if version == '':
                source=sncosmo.get_source(self.model)
            else:
                source=sncosmo.get_source(self.model,version=self.version)
                      
            self.SN=sncosmo.Model(source=source,effects=[dust, dust],
                              effect_names=['host', 'mw'],
                              effect_frames=['rest', 'obs'])
            
            #self.SN=sncosmo.Model(source=self.source)
            
        #self.SN=sncosmo.Model(source=self.source)
       
        #print 'blabla',self.SN.minwave(),self.SN.maxwave()
            lowrange = -30.
            highrange=50.

            self.SN.set(z=self.z)
            self.SN.set(t0=self.t0)
            self.SN.set(c=self.c)
            self.SN.set(x1=self.x1)

            self.SN.set_source_peakabsmag(self.peakAbsMagBesselB, 'bessellB', 'vega',cosmo=self.cosmology)    

            #self.SN.set(mwebv=self.ebvofMW)
           
        else:
            self.Fill_nonIa_ID()
            resu=self.Get_nonIa_Template_ID()
            #id_SED='SDSS-018109'
            #id_SED='SDSS-018457'
            """
            if id_SED != None:
                self.id_SED=id_SED
                print 'tagged',self.id_SED
                #self.SED_all_phases=self.Load_SED('NON1A/'+id_SED+'.SED')
                phase, wave, values = read_griddata_ascii('NON1A/'+id_SED+'.SED')
                print 'min max',np.min(wave),np.max(wave)
                self.SED_template = Spline2d(phase, wave, values, kx=2, ky=2)
            """
            
            
            self.model=thename=resu[0]
            self.version=resu[1]
        
            
            #print 'the choice',resu[0],resu[1]
            """
            thename='snana-2007kw' 
            theversion='1.0'
            """

            if thename != None:
            
                #print 'Getting the source'
                source=sncosmo.get_source(self.model,self.version)
                #print 'Getting the model'
                self.SN=sncosmo.Model(source=source,effects=[dust, dust],
                                      effect_names=['host', 'mw'],
                                      effect_frames=['rest', 'obs'])

                #print 'Setting z and T0'
                self.SN.set(z=self.z)
                self.SN.set(t0=self.t0)
                #print 'SN here',self.z,self.t0,thename,theversion
                self.SN.set_source_peakabsmag(-18., 'bessellB', 'vega',cosmo=self.cosmology)
        #MW extinction
        if mwdust:
            self.SN.set(mwebv=self.ebvofMW)  

        self.model_for_fit()

    def get_SED(self,time):
    
        if self.sn_type== 'Ia':

            if self.model == 'salt2' and self.version=='2.4':
                model_min=2000.
                model_max=9200.0
                wave_min=model_min*(1.+self.z)
                wave_max=model_max*(1.+self.z)
            
            
            if self.model == 'salt2-extended':
                model_min=300.
                model_max=180000.
                wave_min=3000.
                wave_max=11501.
            
        
            wave= np.arange(wave_min,wave_max,1.)
            self.fluxes=10.*self.SN.flux(time,wave)
            self.wavelength=wave/10.
          
            SED_time = Sed(wavelen=self.wavelength, flambda=self.fluxes)
        
            """
            ax, bx = self.SEDfromSNcosmo.setupCCMab()
            self.SEDfromSNcosmo.addCCMDust(a_x=ax, b_x=bx, ebv=self.ebvofMW)
            """

        else:
            wave_min=3000.
            wave_max=11501.            
            wave= np.arange(wave_min,wave_max,1.)

            #print 'getting the flux'
            self.fluxes=10.*self.SN.flux(time,wave)
            self.wavelength=wave/10.
          
            #SED_time = Sed(wavelen=self.wavelength, flambda=self.fluxes/np.power(self.lumidist,2.))
                           #/np.power(self.lumidist,2.))
            
            SED_time = Sed(wavelen=self.wavelength, flambda=self.fluxes)
            """
            a=1./(1.+self.z)
            phase=(time-self.t0)*a
            restwave=wave*a
            flambda=self.SED_template(phase, restwave)[0]
            #print 'alors',len(wave),len(flambda),flambda
           
            SED_time=Sed(wavelen=wave/10., flambda=10.*a*flambda /np.power(self.lumidist,2.))
            """

        return SED_time

    def Load_SED(self, filename):

        sfile=open(filename, 'r')
        mytype=[('phase', np.float),('wavelength', np.float),('flambda', np.float)]
        tab_sed=np.zeros((60,1),dtype=[type for type in mytype])

        inum=-1
        for line in sfile.readlines():
            if line[0]!='#':     
                thesplb=line.split(' ')
                thespl=[spl for spl in thesplb if spl!='']
                #print 'hello pal',thespl
                phase=float(thespl[0])
                wave=float(thespl[1])
                flambda=float(thespl[2].strip())
                inum+=1
                if len(tab_sed) <= inum:
                    tab_sed=np.resize(tab_sed,(len(tab_sed)+100,1))

                tab_sed['phase'][inum]=phase
                tab_sed['wavelength'][inum]=wave
                tab_sed['flambda'][inum]=flambda

        return np.resize(tab_sed,(inum+1,1))


    def model_for_fit(self):

        #print 'there man, to fit'
        dust = sncosmo.OD94Dust() 
        self.fit_model='salt2-extended'
        self.fit_version='1.0'
        
        source=sncosmo.get_source(self.fit_model,version=self.fit_version)
        
        self.SN_fit_model=sncosmo.Model(source=source,effects=[dust, dust],
                              effect_names=['host', 'mw'],
                              effect_frames=['rest', 'obs'])
        
        self.SN_fit_model.set(z=self.z)
        self.SN_fit_model.set_source_peakabsmag(self.peakAbsMagBesselB, 'bessellB', 'vega',cosmo=cosmology.WMAP9)
        
        self.lsstmwebv = EBVbase()
        self.ebvofMW = self.lsstmwebv.calculateEbv(
            equatorialCoordinates=np.array([[np.radians(self.radeg)], [np.radians(self.decdeg)]]))[0]
        
        self.SN_fit_model.set(mwebv=self.ebvofMW)

    def Fill_nonIa_ID(self):

        self.dict_nonIa={}
        #filename='NON1A/NON1A.LIST'
        filename='NON1A.list'
        sfile=open(filename, 'r')
        
        for line in sfile.readlines():
            #print 'hello line',line,line.count('NONIA:')
            #if line.count('NON1A:') > 0:
                
            thesplb=line.split(' ')
            thespl=[spl for spl in thesplb if spl!='']
            thename=thespl[0].strip()
            thetype=thespl[2].strip()
            theversion=thespl[1].strip()
            if self.dict_nonIa.has_key(thetype):
                self.dict_nonIa[thetype].append((thename,theversion))
            else:
                self.dict_nonIa[thetype]=[]
                self.dict_nonIa[thetype].append((thename,theversion))

    def Fill_nonIa_ID_mine(self):
        
        self.dict_nonIa={}
        filename='NON1A/NON1A.LIST'
        #filename='NON1A.list'
        sfile=open(filename, 'r')
        
        for line in sfile.readlines():
            #print 'hello line',line,line.count('NONIA:')
            if line.count('NON1A:') > 0:
                
                thesplb=line.split(' ')
                thespl=[spl for spl in thesplb if spl!='']
                thename=thespl[2].strip()
                thetype=thespl[1].strip()
                if self.dict_nonIa.has_key(thetype):
                    self.dict_nonIa[thetype].append(thename)
                else:
                    self.dict_nonIa[thetype]=[]
                    self.dict_nonIa[thetype].append(thename)

    def Get_nonIa_Template_ID(self):
        
        if not self.dict_nonIa.has_key(self.sn_type):
            print 'problem here - SN type not valid'
            print 'Possible SN types',self.dict_nonIa.keys()
            return None
        else:

            choose=np.random.randint(0,len(self.dict_nonIa[self.sn_type]))
            #print 'result',choose,self.dict_nonIa[self.sn_type][choose]
            return self.dict_nonIa[self.sn_type][choose]
        
    def Get_SED_Restframe(self, Sed_time):
        
        a=1./(1.+self.z)
        bandpass_besselb=Bandpass(wavelen=sncosmo.get_bandpass('bessellB').wave, sb=sncosmo.get_bandpass('bessellB').trans)
        print 'before',Sed_time.wavelen,Sed_time.flambda
        #print 'there we go',Sed_time.wavelen,Sed_time.flambda
        SED_rest=Sed(wavelen=Sed_time.wavelen*a*10., flambda=Sed_time.flambda*np.power(self.lumidist,2.)/a/10./HC_ERG_AA)
        print 'hello',Sed_time.wavelen*a*10,Sed_time.flambda*np.power(self.lumidist,2.)/a/10./HC_ERG_AA
        print 'heelp',SED_rest.wavelen,SED_rest.flambda

        SED_new=Sed(wavelen=SED_rest.wavelen/a, flambda=a*SED_rest.flambda/np.power(self.lumidist,2.))
        #print 'ici',SED_new.wavelen,SED_new.flambda
        #estimate the flux in the B band bessellb 
        
        flux_B_rest=SED_rest.calcFlux(bandpass=bandpass_besselb)
        flux_B_new=SED_new.calcFlux(bandpass=bandpass_besselb)
        #now the magnitudes (apparent and absolute)
        

        vega_SED=Sed()

        vega_SED.readSED_flambda('vega.txt')
        flux_vega=vega_SED.calcFlux(bandpass=bandpass_besselb)

        mag_B_rest=-2.5*np.log10(flux_B_rest/flux_vega)
        mag_B_new=-2.5*np.log10(flux_B_new/3631.)

        print 'hello',len(vega_SED.wavelen),len(vega_SED.flambda),len(bandpass_besselb.sb),flux_vega,flux_B_rest

        bwave=bandpass_besselb.wavelen
        btrans=bandpass_besselb.sb
        vega_wave=vega_SED.wavelen
        vega_flambda=vega_SED.flambda

        mask = ((vega_wave > bwave[0]) & (vega_wave < bwave[-1]))
        d = vega_wave[mask]
        f = vega_flambda[mask]

        trans = np.interp(d, bwave, btrans)
        binw = np.gradient(d)
        ftot = np.sum(f * trans * binw)
        
        print 'vega flux',flux_vega,ftot
        return SED_rest,mag_B_rest,mag_B_new
