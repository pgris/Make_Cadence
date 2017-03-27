import os
#matplotlib inline
import matplotlib.pyplot as plt
from lsst.sims.photUtils import Bandpass
from lsst.sims.photUtils import Sed
import numpy as np

class Throughputs:
    def __init__(self,through_dir='LSST_THROUGHPUTS_BASELINE',atmos_dir='THROUGHPUTS_DIR',aerosol=True):

        self.throughputsDir = os.getenv(through_dir)
        self.atmosDir = os.path.join(os.getenv(atmos_dir), 'atmos')
        
        self.filterlist = ('u', 'g', 'r', 'i', 'z', 'y')
        self.filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}

        self.lsst_std = {}
        self.lsst_system = {}
        self.lsst_detector = {}
        self.lsst_atmos={}
        self.lsst_atmos_aerosol={}
        self.airmass=-1.
        self.aerosol=aerosol
        self.Load_System()
        self.Load_DarkSky()
        self.Load_Atmosphere()
        self.lsst_telescope={}
        self.Load_Telescope()

    def Load_System(self):
        
        for f in self.filterlist:
            self.lsst_std[f] = Bandpass()
            #self.lsst_std[f].readThroughput(os.path.join(self.throughputsDir, 'total_'+f+'.dat'))
            self.lsst_system[f] = Bandpass()
            self.lsst_system[f].readThroughputList(['detector.dat', 'lens1.dat', 'lens2.dat', 'lens3.dat', 
                                                    'm1.dat', 'm2.dat', 'm3.dat', 'filter_'+f+'.dat'], 
                                                   rootDir=self.throughputsDir)

    def Load_Telescope(self):
        
        toload=['detector', 'lens1', 'lens2', 'lens3', 'm1', 'm2', 'm3']
        for filtre in self.filterlist:
            toload.append('filter_'+filtre)

        for system in toload:
            self.lsst_telescope[system]=Bandpass()
            self.lsst_telescope[system].readThroughputList([system+'.dat'], 
                                                   rootDir=self.throughputsDir)
            


    def Load_DarkSky(self):

        self.darksky = Sed()
        self.darksky.readSED_flambda(os.path.join(self.throughputsDir, 'darksky.dat'))
        
    def Load_Atmosphere(self, airmass=1.2):

        self.airmass=airmass
        atmosphere = Bandpass()
        atmosphere.readThroughput(os.path.join(self.atmosDir, 'atmos_%d.dat' %(self.airmass*10)))
        self.atmos= Bandpass(wavelen=atmosphere.wavelen, sb=atmosphere.sb)
       

        for f in self.filterlist:
            wavelen, sb = self.lsst_system[f].multiplyThroughputs(atmosphere.wavelen, atmosphere.sb)
            self.lsst_atmos[f]= Bandpass(wavelen=wavelen, sb=sb)

        if self.aerosol:
            atmosphere_aero = Bandpass()
            atmosphere_aero.readThroughput(os.path.join(self.atmosDir, 'atmos_%d_aerosol.dat' %(self.airmass*10)))
            self.atmos_aerosol= Bandpass(wavelen=atmosphere_aero.wavelen, sb=atmosphere_aero.sb)

            for f in self.filterlist:
                wavelen, sb = self.lsst_system[f].multiplyThroughputs(atmosphere_aero.wavelen, atmosphere_aero.sb)
                self.lsst_atmos_aerosol[f]= Bandpass(wavelen=wavelen, sb=sb)


    def Plot_Throughputs(self):

        #colors=['b','g','r','m','c',[0.8,0,0]]
        #self.filtercolors = {'u':'b', 'g':'c', 'r':'g', 'i':'y', 'z':'r', 'y':'m'}
        style=[',',',',',',',']
        for i,band in enumerate(['u','g','r','i','z','y']):
            #plt.plot(self.lsst_std[band].wavelen,self.lsst_std[band].sb,linestyle='-',color=self.filtercolors[band], label='%s - std' %(band))
            plt.plot(self.lsst_system[band].wavelen,self.lsst_system[band].sb,linestyle='--',color=self.filtercolors[band], label='%s - syst' %(band))
            plt.plot(self.lsst_atmos[band].wavelen,self.lsst_atmos[band].sb,linestyle='-.',color=self.filtercolors[band], label='%s - syst+atm' %(band))
            if len(self.lsst_atmos_aerosol) > 0:
                plt.plot(self.lsst_atmos_aerosol[band].wavelen,self.lsst_atmos_aerosol[band].sb,linestyle='-',color=self.filtercolors[band], label='%s - syst+atm+aero' %(band))

        plt.plot(self.atmos.wavelen, self.atmos.sb, 'k:', label='X =%.1f atmos' %(self.airmass),linestyle='-')
        if len(self.lsst_atmos_aerosol) > 0:
            plt.plot(self.atmos_aerosol.wavelen, self.atmos_aerosol.sb, 'k:', label='X =%.1f atm+aero' %(self.airmass),linestyle='--')
        plt.legend(loc=(0.85, 0.1), fontsize='smaller', fancybox=True, numpoints=1)
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Sb (0-1)')
        plt.title('System throughput')
        #plt.show()

    def Plot_DarkSky(self):

        plt.plot(self.darksky.wavelen,self.darksky.flambda,'k:',linestyle='-')
        plt.xlabel('Wavelength (nm)')                                                                                   
        plt.ylabel('flambda (ergs/cm$^2$/s/nm)')                                                                                          
        plt.title('Dark Sky SED')                                                                                  
        #plt.show()
 
    def Plot_Throughputs_Spectrum(self,wavelength,fluxes,z):
        
        fig, ax1 = plt.subplots()
        style=[',',',',',',',']

        for i,band in enumerate(['u','g','r','i','z','y']):
            print 'Mean wave',band,np.mean(self.lsst_std[band].wavelen)
            #plt.plot(self.lsst_std[band].wavelen,self.lsst_std[band].sb,linestyle='-',color=self.filtercolors[band], label='%s - std' %(band))
            #plt.plot(self.lsst_system[band].wavelen,self.lsst_system[band].sb,linestyle='--',color=self.filtercolors[band], label='%s - system' %(band))
            plt.plot(self.lsst_system[band].wavelen,self.lsst_system[band].sb,linestyle='-',color=self.filtercolors[band], label='%s - system' %(band))
            #plt.plot(self.lsst_atmos[band].wavelen,self.lsst_atmos[band].sb,linestyle='-.',color=self.filtercolors[band], label='%s - system+atmos' %(band))
            #ax1.plot(self.lsst_atmos_aerosol[band].wavelen,self.lsst_atmos_aerosol[band].sb,linestyle='-',color=self.filtercolors[band], label='%s - system+atmos+aerosol' %(band))

        
        
        #plt.plot(self.atmos.wavelen, self.atmos.sb, 'k:', label='X =%.1f atmos' %(self.airmass),linestyle='-')
        #plt.plot(self.atmos_aerosol.wavelen, self.atmos_aerosol.sb, 'k:', label='X =%.1f atmos+aero' %(self.airmass),linestyle='--')
        #plt.legend(loc=(0.85, 0.1), fontsize='smaller', fancybox=True, numpoints=1)
        ax1.set_xlabel('Wavelength (nm)')
        ax1.set_ylabel('Sb (0-1)')

        ax2 = ax1.twinx()
        ax2.plot(wavelength,fluxes, color='lightgrey',marker='.',markersize=0.02,alpha=0.5)
        ax2.set_ylabel('Flux (Unites Arbitraires)', color='k')

        #z=0.5

        wave_shifted=[]
        for val in wavelength:
            wave_shifted.append(val*(1+z))
        ax2.plot(wave_shifted,fluxes, color='lightgreen',marker='.',markersize=0.02,alpha=0.5)

        #plt.title('System throughput')

        ax1.set_xlim([300,1150])
        #plt.show()
