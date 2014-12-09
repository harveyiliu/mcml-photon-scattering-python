

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
        if mediumName.lower() == 'TypeIIEpidermis'.lower():
            self.n = 1.3
            self.mua = 5.0
            self.mus = 200.0
            self.g = 0.70
        elif mediumName.lower() == 'dermis'.lower():
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
            n - refractive index
            mua - absorption coefficient. [1/cm]
            mus - scattering coefficient. [1/cm]
            g - anisotropy
        Methods:
            
    """
    
    def __init__(self, layerName = 'TypeIISkin'):
        if layerName.lower() == 'TypeIISkin'.lower():
            self.nIn = 1.0      # incidence medium index
            self.nOut = 1.4     # exit medium index
            self.layer = [Medium('TypeIIEpidermis'), Medium('dermis')]
            self.layerThickness = [0.006, 0.3]  # layer thickness in [cm]
        elif layerName.lower() == 'dermis'.lower():
            self.nIn = 1.0
            self.nOut = 1.4
            self.layer = [Medium('dermis')]
            self.layerThickness = [0.3]     
        else:
            self.nIn = 1.0
            self.nOut = 1.4
            self.layer = [Medium('dermis')]
            self.layerThickness = [0.3]
        self.layerZ = []
        self.cosCrit = []
        z = 0   # incidence first medium z coordinate [cm]
        # find the z depth coordinates and cosine critical angles for each
        #   layer
        for i in range(len(self.layerThickness)):
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
            if (i+1) == len(self.layerThickness):
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

    def __init__(self, rSpecular = 0.017, \
        layerObj = LayerStruct('TypeIISkin')):
        
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
        


   

        
              
            

