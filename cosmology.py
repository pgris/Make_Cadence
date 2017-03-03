import numpy as np

class Cosmology:
    def __init__(self,Omegak):
        self.omegak=Omegak

        self.sqrt_abs_omegak = np.sqrt(np.absolute(Omegak)) 
        
        if Omegak > 0.0001:
            self.Sin = np.sinh
            self.Cos = np.cosh
        else:
            if Omegak < -0.0001:
                self.Sin = np.sin
                self.Cos = np.cos
            else:
                self.Sin = self.ident
                self.Cos = self.one
                self.sqrt_abs_omegak = 1.
                
        self.omegar=0

    def ident(self,x):
        return x

    def one(self,x):
        return 1

    def integ_hinv(self,zmin,zmax,nstep):
        return -1

    def E(self,z):
        return -1

    def Dm(self, z, Nsteps=100): # Proper motion distance
        return self.Sin(self.sqrt_abs_omegak*self.integ_hinv(0,z,Nsteps))/self.sqrt_abs_omegak

    def Dl(self,z,Nsteps=100): #the Luminosity distance (in units of c/H0)
        return (1+z)*self.Dm(z,Nsteps)
    
    def Da(self,z,Nsteps=100): #the Angular distance (in units of c/H0)
        return self.Dm(z,Nsteps)/(1+z)
        
    def dVdz(self,z):
        # still to be checked against CPT 92 curves 
        # checked "experimentally" however that this routine is the derivative of Volume(z) 
        integral = self.integ_hinv(0,z)
        dm = self.Sin(self.sqrt_abs_omegak*integral)/self.sqrt_abs_omegak
        ddmdz = self.Cos(self.sqrt_abs_omegak*integral) / E(z)
        return dm*dm*ddmdz/np.sqrt(1+self.omegak*pow(dm,2.))
  
#from 0 to z. it is in fact H_0^3 V
    def Volume(self,z):
        integral = self.integ_hinv(0,z)
        #print 'integral',integral
        dm = self.Sin(self.sqrt_abs_omegak*integral)/self.sqrt_abs_omegak
        #print 'dm',dm
        x = self.omegak*pow(dm,2.)
        # I checked experimentaly the continuity over omegak = 0 */
        if np.absolute(x) > 1e-3:
            #use exact (however unstable) expression (CPT eq 27.)
            return 0.5*(dm*np.sqrt(1+x) - integral)/self.omegak
        else:
            # use taylor expansion
            return pow(dm,3.)*(1 - 0.3*x + 9.*pow(x,2.)/56.- 5.*pow(x,3.)/48.)/3.
"""
> a:=int(x^2/sqrt(1+ok*x^2),x=0..z);
                           2 1/2   3/2        1/2              2 1/2
                z (1 + ok z )    ok    - ln(ok    z + (1 + ok z )   ) ok
       a := 1/2 --------------------------------------------------------
                                           5/2
                                         ok
 
(Equivalent to Carrol Press Turner eq 27)

> series(a,ok=0);
                3         5            7   2          9   3       4
           1/3 z  - 1/10 z  ok + 3/56 z  ok  - 5/144 z  ok  + O(ok )
 
>                                                                               
"""

class GeneralCosmo(Cosmology):
    def __init__(self,OmegaM,OmegaX,W0,W1):
        Cosmology.__init__(self,1.-OmegaM-OmegaX)
        self.omegam = OmegaM 
        self.omegax = OmegaX 
        self.w0 = W0 
        self.w1 = W1

    def w(self,z):
        return self.w0+self.w1*z

    def E(self,z):
        #cannot use easily the routine above, because of the change of variable
        return np.sqrt(self.omegam*pow(1.+z,3.) + self.omegax*
	      np.exp(3.*(self.w1*z + (1.+ self.w0-self.w1)*np.log(1.+z)))
	      + pow(1+z,2.)*(self.omegak + pow(1+z,2.)*self.omegar))
    
    def integ_hinv(self, ZMin,ZMax,NStep=100):
        #change z to u=1/sqrt(1+z)
        umin = 1./np.sqrt(1+ZMax)
        umax = 1./np.sqrt(1+ZMin)
        step = (umax-umin)/NStep
        u = umin+0.5*step
        #print 'umin,...',umin,umax,step,u
        thesum =0
        for i in range(0,NStep):
            u2 = pow(u,2.)
            fact_ox = np.exp(3*(-2*np.log(u)*(self.w0-self.w1)+self.w1*(1./u2-1)))
            u3hinv = np.sqrt(self.omegam + self.omegax*fact_ox + self.omegak*u2 + self.omegar/u2)
            thesum += 1./u3hinv;
            u+=step
            #print 'there pal',i,u2,fact_ox,u3hinv,thesum
        return 2.*thesum*step
