from math import pow
import numpy as np

class SN_Rate:

    def __init__(self,cosmology,zmin,zmax,zstep,duration,area,rate_name='Ripoche'):

      self.cosmology=cosmology
      self.zmin=zmin
      self.zmax=zmax
      self.rate_name=rate_name
      self.hz=[]
      self.hz_weight=[]
      self.area=area
      self.duration=duration

      self.FillPdfHists(zstep,1)

    def FillPdfHists(self,zstep,sf):
        print 'Filling histogram for redshift',self.zmin,self.zmax,zstep
        nbins = int((self.zmax-self.zmin) / zstep)
        for i in range(0,nbins):
            dzmin = self.zmin+i*zstep
            dzmax = dzmin+zstep
            weight=sf*self.NumberOfSupernovae(dzmin, dzmax)
            self.hz.append(0.5*(dzmin+dzmax))
            self.hz_weight.append(weight)
            #print '++>',weight,'  dzmin, dzmax=',dzmin,dzmax 
        #print '==>',np.sum(self.hz_weight)  
        
    def NumberOfSupernovae(self,zmin,zmax):
        c = 3.e5 #km/s
        H0 = 70
        #steradians2SquareDegrees=(180.*180./M_PI/M_PI)
        steradians2SquareDegrees=(180.*180./np.pi/np.pi)
        nstep = int(np.ceil(zmax-zmin)/0.02)
        zstep = (zmax-zmin)/nstep
        total = 0
        #print 'nsuper', steradians2SquareDegrees,nstep,zstep
        for istep in range(0,nstep):
            zmin_step = zmin+istep*zstep
            zmax_step = zmin_step+zstep
            vol = (self.cosmology.Volume(zmax_step) - self.cosmology.Volume(zmin_step))/steradians2SquareDegrees
            zmean_step = 0.5*(zmax_step+zmin_step)

            if self.rate_name == 'Ripoche':
                rate = self.ripoche_rate(zmean_step)

            if self.rate_name == 'Flat':
                rate = self.flat_rate(zmean_step)

      #if self.use_mannucci_rate :
        #  rate = mannucci_rate(zmean_step)          
            if self.rate_name=='perret':
                rate = self.perrett_rate(zmean_step)
        
            nsn= rate*vol*pow(c/H0,3.)*self.area*self.duration/(1+zmean_step)
            #print 'nsn',nsn,rate,vol,self.area,self.duration,zmean_step
            if self.rate_name == 'Flat':
                nsn=rate*self.area*self.duration/(1+zmean_step)
            total += nsn

        return total


    def flat_rate(self,z):
        return (1.+z)

    def ripoche_rate(self,z):
        rate = 1.53e-4*0.343
  # h70^3 Mpc-3 y-1, with h70 = H_0/70 km/s/Mpc at z= 0.5*/
        exposant = 2.14 #see Ripoche at moriond 2008
        # since it was not measured beyond 1, choose 
  #increasing rate up to z=1, static afterwards (this mimicks more or less Mannucci et al 2007)
  #zinc = (z<1) ? z : 1
        zinc= z if z<1 else 1
  #print 'Ripoche rate',z,rate*pow((1+zinc)/1.5,exposant)
        return rate*pow((1+zinc)/1.5,exposant)

    def perrett_rate(z):
        rate = 0.17E-4
        exposant = 2.11
  #since it was not measured beyond 1, choose 
  #increasing rate up to z=1, static afterwards (this mimicks more or less Mannucci et al 2007)
        zinc= z if z<1 else 1
        return rate*pow((1+zinc),exposant)

"""
def mannucci_rate(z)
{
  static General1DFunction *rate=NULL;
  if (!rate)
    {
      string rateFileName = LocateFile("rates/mannucci.rate");
      rate = new General1DFunction(rateFileName);
    }
  return rate->Value(z);
}
"""
