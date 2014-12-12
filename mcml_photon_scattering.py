
import numpy as np

WEIGHT = 1e-4       # critical weight for roulette

class Medium:
    """Medium class - optical medium class defining the optical properties
        Class instance variables:
            n - refractive index
            mua - absorption coefficient. [1/cm]
            mus - scattering coefficient. [1/cm]
            g - anisotropy
        Methods:
            
    """
    
    def __init__(self, mediumName = 'dermis'):
        
        # initialize medium optical properties
        if mediumName.lower() == 'TYPE_II_EPIDERMIS'.lower():
            self.n = 1.3
            self.mua = 5.0
            self.mus = 200.0
            self.g = 0.70
        elif mediumName.lower() == 'DERMIS'.lower():
            self.n = 1.4
            self.mua = 0.26
            self.mus = 137.0
            self.g = 0.90
        else:
            self.n = 1.4
            self.mua = 0.26
            self.mus = 137.0
            self.g = 0.90     
        

class LayerStruct:
    """LayerStruct class - multi-layered structure
        Class instance variables:
            nIn - refractive index of the incidence medium
            nOut - refractive index of the exit medium
            numLayers - number of layers
            layer - list of layer objects
            layerThickness - layer thickness in [cm]
            layerZ - layer depth z coordinates, top and bottom [cm]
            cosCrit - ciritical angle cosines of each layer, top and bottom
        Methods:
            
    """
    
    def __init__(self, layerName = 'BARE_DERMIS'):
        if layerName.lower() == 'TYPE_II_SKIN'.lower():
            self.nIn = 1.0      # incidence medium index
            self.nOut = 1.4     # exit medium index
            self.numLayers = 2  # number of layers
            self.layer = [Medium('TYPE_II_EPIDERMIS'), Medium('DERMIS')]
            self.layerThickness = [0.006, 0.3]  # layer thickness in [cm]
        elif layerName.lower() == 'BARE_DERMIS'.lower():
            self.nIn = 1.0
            self.nOut = 1.4
            self.numLayers = 1
            self.layer = [Medium('DERMIS')]
            self.layerThickness = [0.3]     
        else:
            self.nIn = 1.0
            self.nOut = 1.4
            self.numLayers = 1
            self.layer = [Medium('DERMIS')]
            self.layerThickness = [0.3]
        self.layerZ = []
        self.cosCrit = []
        z = 0   # incidence first medium z coordinate [cm]
        # find the z depth coordinates and cosine critical angles for each
        #   layer
        for i in range(self.numLayers):
            self.layerZ.append([z, z+self.layerThickness[i]])
            z = self.layerZ[-1][1]
            # calculate the critical angle cosines for each layer
            # crticial angle at top interface of the current layer
            n1 = self.layer[i].n
            if i == 0:
                n2 = self.nIn
            else:
                n2 = self.layer[i-1].n
            if n1 > n2:
                cosCrit0 = (1.0 - n2*n2/(n1*n1))**0.5
            else:
                cosCrit0 = 0.0
            # crticial angle at bottom interface of the current layer
            if (i+1) == self.numLayers:
                n2 = self.nOut
            else:
                n2 = self.layer[i+1].n
            if n1 > n2:
                cosCrit1 = (1.0 - n2*n2/(n1*n1))**0.5
            else:
                cosCrit1 = 0.0
            self.cosCrit.append([cosCrit0, cosCrit1])
 

class Photon:
    """Photon class - MCML photon class for Monte Carlo scattering model in
        multilayered turbid media. 
        Class instance variables:
            x = Cartesian coordinate x [cm]
            y = Cartesian coordinate y [cm]
            z = Cartesian coordinate z [cm]
            ux = directional cosine x of a photon
            uy = directional cosine y of a photon
            uz = directional cosine z of a photon
            w - weight
            dead - true if photon is terminated
            layer - index to layer where the photon packet resides
            s - current step size [cm]
            sleft - step size left, dimensionless [-]           
        Methods:
            
    """

    def __init__(self, layerObj = LayerStruct('BARE_DERMIS'), rSpecular = 0.017):
        
        # initialize a photon
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.ux = 0.0
        self.uy = 0.0
        self.uz = 1.0
        self.w = 1.0 - rSpecular
        self.dead = False
        self.layer = 0
        self.s = 0
        self.sleft = 0

        # take care of the case when the first layer is glass
        if (layerObj.layer[0].mua == 0) and (layerObj.layer[0].mus == 0):
            self.layer = 1      # skip to next layer
            self.z = layerObj.layerZ[1][0]  # use z0 from the next layer
        


class ModelInput:
    """ModelInput class - multi-layered photon scattering model input
        Class instance variables:
            numPhotons - number of photons to be traced
            Wth - play roulette if photon weight < Wth
            dz - z grid separation [cm]
            dr - r grid separation [cm]
            da - alpha grid separation [radian]
            nz - array range 0..nz-1
            nr - array range 0..nr-1
            na - array range 0..na-1
            layerObj - medium layer structure class instance
        Methods:
            
    """
    
    def __init__(self, modelName = 'BARE_DERMIS', numPhotonsSet = 1000):
        if modelName.lower() == 'TYPE_II_SKIN'.lower():
            self.layerObj = LayerStruct('TYPE_II_SKIN')
            self.dz = 20e-4
            self.dr = 20e-4
            self.nz = 10
            self.nr = 50
            self.na = 10
        elif modelName.lower() == 'BARE_DERMIS'.lower():
            self.layerObj = LayerStruct('BARE_DERMIS')
            self.dz = 100e-4
            self.dr = 100e-4
            self.nz = 30
            self.nr = 50
            self.na = 10
        else:
            self.layerObj = LayerStruct('BARE_DERMIS')
            self.dz = 100e-4
            self.dr = 100e-4
            self.nz = 30
            self.nr = 50
            self.na = 10
        self.numPhotons = numPhotonsSet
        self.Wth = WEIGHT
        self.da = 0.5*np.pi/self.na


class MCMLModel(ModelInput):
    """MCMLModel class - multi-layered photon scattering model, inherits from
        ModelInput layer structure setup
        Class instance variables:
            Rsp - specular reflectance [-]
            Rd - total diffuse reflectance [-]
            A - total absorption probability [-]
            Tt - total transmittance [-]
            Rd_ra - 2D distribution of diffuse reflectance [1/(cm2 sr)]
            Rd_r - 1D radial distribution of diffuse reflectance [1/cm2]
            Rd_a - 1D angular distribution of diffuse reflectance [1/sr]
            A_rz - 2D probability density in turbid media over r & z [1/cm3]
            A_z - 1D probability density over z [1/cm]
            A_l - each layer's absorption probability [-]
            Tt_ra - 2D distribution of total transmittance [1/(cm2 sr)]
            Tt_r - 1D radial distribution of transmittance [1/cm2]
            Tt_a - 1D angular distribution of transmittance [1/sr]
        Methods:
            
    """

    def __init__(self, modelName = 'BARE_DERMIS', numPhotonsSet = 1000):
        ModelInput.__init__(self, modelName, numPhotonsSet)
        # initialize the model grid arrays    
        self.Rsp = 0.0
        self.Rd = 0.0
        self.A = 0.0
        self.Tt = 0.0
        self.Rd_ra = np.matrix(np.zeros((self.nr, self.na)))
        self.Rd_r = np.zeros(self.nr)
        self.Rd_a = np.zeros(self.na)
        self.A_rz = np.matrix(np.zeros((self.nr, self.nz)))
        self.A_z = np.zeros(self.nz)
        self.A_l = np.zeros(2 + self.layerObj.numLayers)
        self.Tt_ra = np.matrix(np.zeros((self.nr, self.na)))
        self.Tt_r = np.zeros(self.nr)
        self.Tt_a = np.zeros(self.na)          
            

