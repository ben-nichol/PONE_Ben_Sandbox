import numpy as np
import os

class MuonGenerator () :

    def ReadFluxTable(self) :

        infile = open(os.getenv('PONESRCDIR')+"/data/fluxtable_1460.csv","r")
        lines = infile.readlines()
        #first 4 lines are header
        Zenithline = 5

        zenithbinlines = lines[Zenithline].split(" ",1000)
        self.zenithbins = list()
        self.energybins = list()
        self.flux = list()
        for i in range(1,len(zenithbinlines)) :
            if zenithbinlines[i] == '' :
                continue
            if zenithbinlines[i] == ' ' :
                continue
            try :
                _flux = float(zenithbinlines[i]) 
                self.zenithbins.append(max(_flux*(np.pi/180.),0.0))
                self.flux.append(list())
            except :
                #do nothing
                k=0

        self.maxzenith = self.zenithbins[-1]
        #print(self.maxzenith*np.pi/180.)
        self.de = 0.1
        self.minzenith = self.zenithbins[0]
        self.dz = ((self.maxzenith-self.minzenith)/(len(self.zenithbins)-1))

        #for i in range(len(self.zenithbins)) :
        #    print(str(self.zenithbins[i]*np.pi/180.)+" "+str(float(i)*self.dz))

        for l in range(Zenithline+1,len(lines)) :
            splitline = lines[l].split(" ",1000)
            self.energybins.append(np.log10(float(splitline[0]))-3.0)

            for i in range(1,len(splitline)) :
                try :
                    _flux = float(splitline[i])
                    self.flux[i-1].append(max(_flux*np.cos(self.zenithbins[i-1]),0.0))
                    #self.flux[-i].append(max(_flux*np.cos(self.zenithbins[i-1]),0.0))
                except :
                    k=0

        self.minenergy = self.energybins[0] #convert MeV to GeV
        self.de = (self.energybins[-1]-self.energybins[0])/(len(self.energybins)-1)
        
        #print("max energy "+str(self.minenergy+self.de*len(self.flux[0])))

    def GetArea():
        return np.pi*(self.maxinjectionR**2.0)

    def ComputeAcceptanceRange(self,injectR) :

        if injectR <= self.cylR :
            return 0.0, np.pi

        maxzenith = np.arctan((injectR-self.cylR)/self.cylH)
        maxtheta = np.arctan(self.cylR/(injectR-self.cylR))

        return maxzenith, maxtheta

    def GetZenithCDF(self,injectR) :

        self.zenithcdf = list()
        minZenith, maxtheta = self.ComputeAcceptanceRange(injectR)

        ZenithWeight = list()
        for z in range(int(minZenith/dz),int(np.pi/dz)) :
            rate = 0.0
            zenith = dz*z
            for e in range(len(self.flux[z])) :
                rate += self.flux[z][e]*self.GetEnergyde(e)
            ZenithWeight.append(rate*np.sin(self.minzenith+(0.5+float(z))*self.dz))
        totalflux = sum(ZenithWeight)

        z=0
        cdf = 0.0
        self.cdf_inv = list()
        for i in range(self.d_cdf) :
            p = float(0.5+i)/self.d_cdf -cdf
            while Zenithweight[z]/totalflux <= p :
                p-= ZenithWeight[z]/totalflux
                cdf += ZenithWeight[z]/totalflux
                z+=1
            self.zenithcdf.append(self.minzenith+(float(z)+(p/(ZenithWeight[z]/totalflux)))*self.dz)

    def MakeEnergyCDFs(self) :
        self.energyinvers_cdfs = list()

        for z in range(len(self.flux)) :
            eflux = self.flux[z].copy()
            for e in range(len(eflux)):
                energy = self.minenergy + e*self.de
                if energy < self.mingenenergy:
                    eflux[e] = 0.0
                if energy > self.maxgenenergy :
                    eflux[e] = 0.0
                eflux[e] *= self.GetEnergyde(e)

            totalflux = sum(eflux)
            #print(totalflux)
            e=0
            cdf = 0.0
            inverse_cdf = list()
            for i in range(self.d_cdf) :
                p = float(0.5+i)/self.d_cdf -cdf
                while eflux[e]/totalflux <= p :
                    p-= eflux[e]/totalflux
                    cdf += eflux[e]/totalflux
                    e+=1
                inverse_cdf.append(self.minenergy+(e+(p/(eflux[e]/totalflux)))*self.de)
            self.energyinvers_cdfs.append(inverse_cdf)

    def GetEnergyde(self,e) :
        return 10.**(self.minenergy + self.de*float(e+1)) - 10.**(self.minenergy + self.de*float(e))

    def MakeZenithCDF(self) :
        self.FullZenithCDF = list()

        for z in range(len(self.flux)) :
            _flux = 0.0
            for e in range(len(self.flux[z])) :
                _flux += self.flux[z][e]*self.GetEnergyde(e)*np.sin(self.minzenith+(0.5+float(z))*self.dz)
            self.FullZenithCDF.append(_flux)
        totalFlux = sum(self.FullZenithCDF)

        #for z in range(len(self.FullZenithCDF)) :
        #    print("zenith = "+str(z*self.dz)+" flux = "+str(self.FullZenithCDF[z]))

        self.FullZenithCDF[0] /= totalFlux

        for z in range(1,len(self.FullZenithCDF)) :
            self.FullZenithCDF[z] = self.FullZenithCDF[z-1]+self.FullZenithCDF[z]/totalFlux


    def GetZenithBin(self,zenith) :
        return max(0,min(int((zenith),len(self.flux)-1))

    def GetZenithCDF(self,radius):
        minZenith, maxtheta = self.ComputeAcceptanceRange(radius)
        thiscdf = self.FullZenithCDF.copy()
        #print(str(radius)+ " " +str(minZenith) + " " + str(self.dz*(len(thiscdf)-1))+" "+str(self.maxzenith))
        for z in range(len(thiscdf)) :
            if z*self.dz < minZenith :
                thiscdf[z] = 0.0

        totalflux = sum(thiscdf)

        for z in range(len(thiscdf)) :
            thiscdf[z] /= totalflux
        cdf = 0.0
        p = 0.0
        z = 0
        inverscdf = list()
        for i in range(self.d_cdf) :
            p = (0.5+float(i))/self.d_cdf -cdf
            while thiscdf[z] <= p:
                p-= thiscdf[z]
                cdf += thiscdf[z]
                z+=1
            _z = (float(z)+(p/thiscdf[z]))*self.dz
            inverscdf.append(_z)
        #print(len(thiscdf))
        #print(inverscdf)
        return inverscdf,cdf

    def GetZenithCDFvalue(self, rbin, maxzenith) :

        return self.Zenithcfd[rbin][GetZenithBin(maxZenith)]

    def ComputeRadialWeights(self) :

        RadialWeight = list()
        RadialWeightCFD = list()
        for i in range(int(self.maxinjectionR)) :
            radius = float(i)+0.5
            rate = 0.0
            minZenith, maxtheta = self.ComputeAcceptanceRange(radius)
         #   print("radius = "+str(radius)+" minZenith = "+str(minZenith)+" maxtheta = "+str(maxtheta))
            for z in range(len(self.flux)):
                zenith =  (float(z)+0.5)*self.dz
                if zenith < minZenith :
                    continue
                for e in range(len(self.flux[z])) :
                    rate += self.flux[z][e]*(np.sin(zenith)*dz)*(2.0*maxtheta*radius)*self.GetEnergyde(e)

            RadialWeight.append(rate)

        totalflux = sum(RadialWeight)
        maxflux = max(RadialWeight)

        #for i in range(len(RadialWeight)) :
        #    print("radius "+str(i)+" totalfluxfrac "+str(RadialWeight[i]/maxflux))

        r=0
        cdf = 0.0
        self.cdf_inv = list()
        self.Zenithcdf_inv = list()
        self.Zenithcdf = list()
        for i in range(self.d_cdf) :
            p = (0.5+float(i))/self.d_cdf -cdf
            while RadialWeight[r]/totalflux <= p :
                p-= RadialWeight[r]/totalflux
                cdf += RadialWeight[r]/totalflux
                r+=1
            radius = float(r)+(p/(RadialWeight[r]/totalflux))
            self.cdf_inv.append(radius)
            zenithinvcdf, zenithcdf = self.GetZenithCDF(radius)
            self.Zenithcdf_inv.append(zenithinvcdf)
            self.Zenithcdf.append(zenithcdf)

    def __init__(self,_cylR,_minenergy,_maxenergy,seed=0):

        self.cylR = _cylR
        self.cylH = 1400.0
        self.d_cdf = int(1e4)

        self.mingenenergy = float(_minenergy)
        self.maxgenenergy = float(_maxenergy)
        self.maxinjectionR = 7000.
        np.random.seed(seed)

        self.ReadFluxTable()
        self.MakeEnergyCDFs()
        self.MakeZenithCDF()
        self.ComputeRadialWeights()

    def GenerateRandomEvent(self) :
        
        radius = 10000.
        while radius > 6000. :
            random = np.random.rand()
            rbin = random*self.d_cdf
            rbin_int = int(rbin)
            radius = self.cdf_inv[rbin_int] 
            if rbin_int < len(self.cdf_inv)-1 :
                radius += (rbin-rbin_int)*(self.cdf_inv[rbin_int+1]-self.cdf_inv[rbin_int])

        random = np.random.rand()
        injecttheta = 2.0*np.pi*random
    
        random = np.random.rand()
        zbin = int(random*self.d_cdf)
        #print(zbin) 
        zenith = self.Zenithcdf_inv[rbin_int][zbin]
    
        random = np.random.rand()
        ebin = int(random*self.d_cdf)
        z = self.GetZenithBin(zenith)
        
        energyexp = self.energyinvers_cdfs[z][ebin]
        energy = 10.0**(energyexp)

        #print("energy = "+str(energy)+ " minincdf = "+str(self.energyinvers_cdfs[z][0]))

        minzenith, maxtheta = self.ComputeAcceptanceRange(radius)
        zenbin = self.GetZenithBin(minzenith)
        zenbin_extra = minzenith/self.dz - zenbin
        p = self.FullZenithCDF[zenbin]
        if zenbin < len(self.FullZenithCDF)-1 and abs(minzenith-zenbin*self.dz) > 0.0 :
            p_high = self.FullZenithCDF[zenbin+1]
            p = p+(p_high-p)*(zenbin_extra)
        #if radius > 7000. :
        #print("radius = "+str(radius))
        #print(self.FullZenithCDF[zenbin])
        #print(self.FullZenithCDF[zenbin+1])
        #print(zenbin)
        #print(zenith)
        #print(minzenith)
        #print(p)
        if not p < 1.0 :
            #print("p = "+str(p))
            p=0.999999
        weight = (np.pi/maxtheta)*(1./(1. - p))
        #weight = p
        random = np.random.rand()
        theta = (-1.0+2.0*random)*maxtheta

        efluxbin = int(min(max(0,(energyexp-self.minenergy)/self.de),len(self.flux[0])-1))
        zfluxbin = int(min(max(0,zenith/self.dz),len(self.flux)-1))
        flux = self.flux[zfluxbin][efluxbin]

        vertex = [radius*np.cos(injecttheta),radius*np.sin(injecttheta),self.cylH/2.0]
        direction = [-np.cos(injecttheta)*np.sin(zenith),-np.sin(injecttheta)*np.sin(zenith),-np.cos(zenith)]
        direction_rotated = [direction[0]*np.cos(theta) -direction[1]*np.sin(theta),
                             direction[0]*np.sin(theta) + direction[1]*np.cos(theta),
                             direction[2]]

        return vertex, direction_rotated, energy, weight, flux
