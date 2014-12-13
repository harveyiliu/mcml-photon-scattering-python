
import numpy as np

WEIGHT = 1e-4       # critical weight for roulette
CHANCE = 0.1		    # Chance of roulette survival
PARTIALREFLECTION = 0     # 1=split photon, 0=statistical reflection.
COSZERO = 1.0 - 1.0e-12     # cosine of about 1e-6 rad
COS90D = 1.0e-6     # cosine of about 1.57 - 1e-6 rad

class Medium:
    """Medium class - optical medium class defining the optical properties
        Class instance variables:
            n - refractive index
            mua - absorption coefficient. [1/cm]
            mus - scattering coefficient. [1/cm]
            g - anisotropy
        Methods:
            
    """
    
    def __init__(self, mediumName = 'DERMIS'):
        
        # initialize medium optical properties
        if mediumName.lower() == 'AIR'.lower():
            self.n = 1.0
            self.mua = 0.0
            self.mus = 0.0
            self.g = 0.0
        elif mediumName.lower() == 'DERMIS'.lower():
            self.n = 1.4
            self.mua = 0.26
            self.mus = 137.0
            self.g = 0.90
        elif mediumName.lower() == 'TYPE_II_EPIDERMIS'.lower():
            self.n = 1.3
            self.mua = 5.0
            self.mus = 200.0
            self.g = 0.70    
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
        if layerName.lower() == 'BARE_DERMIS'.lower():
            self.numLayers = 1
            self.layer = [Medium('AIR'), Medium('DERMIS'), Medium('DERMIS')]
            self.layerThickness = [0.3]
        elif layerName.lower() == 'TYPE_II_SKIN'.lower():
            self.numLayers = 2  # number of layers
            self.layer = [Medium('AIR'), Medium('TYPE_II_EPIDERMIS'), \
                Medium('DERMIS'), Medium('DERMIS')]
            self.layerThickness = [0.006, 0.3]  # layer thickness [cm]           
        else:
            self.numLayers = 1
            self.layer = [Medium('AIR'), Medium('DERMIS'), Medium('DERMIS')]
            self.layerThickness = [0.3]
        self.layerZ = []
        self.cosCrit = []
        z = 0   # incidence first medium z coordinate [cm]
        self.layerZ.append([0, 0])  # first incidence medium, not used
        self.cosCrit.append([0, 0])  # first incidence medium, not used
        # find the z depth coordinates and cosine critical angles for each
        #   layer
        for i in range(1, self.numLayers+1):
            self.layerZ.append([z, z+self.layerThickness[i-1]])
            z = self.layerZ[-1][1]
            # calculate the critical angle cosines for each layer
            # crticial angle at top interface of the current layer
            n1 = self.layer[i].n
            n2 = self.layer[i-1].n
            if n1 > n2:
                cosCrit0 = (1.0 - n2*n2/(n1*n1))**0.5
            else:
                cosCrit0 = 0.0
            # crticial angle at bottom interface of the current layer
            n2 = self.layer[i+1].n
            if n1 > n2:
                cosCrit1 = (1.0 - n2*n2/(n1*n1))**0.5
            else:
                cosCrit1 = 0.0
            self.cosCrit.append([cosCrit0, cosCrit1])


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
        if modelName.lower() == 'BARE_DERMIS'.lower():
            self.layerObj = LayerStruct('BARE_DERMIS')
            self.dz = 100e-4
            self.dr = 100e-4
            self.nz = 30
            self.nr = 50
            self.na = 10
        elif modelName.lower() == 'TYPE_II_SKIN'.lower():
            self.layerObj = LayerStruct('TYPE_II_SKIN')
            self.dz = 20e-4
            self.dr = 20e-4
            self.nz = 10
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
        # extend the ModelInput base class instance variables
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

    def __init__(self, layerObj = LayerStruct('BARE_DERMIS'), \
        rSpecular = 0.017):
        
        # initialize a photon
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        self.ux = 0.0
        self.uy = 0.0
        self.uz = 1.0
        self.w = 1.0 - rSpecular
        self.dead = False
        self.layer = 1
        self.s = 0
        self.sleft = 0

        # take care of the case when the first layer is glass
        if (layerObj.layer[1].mua == 0.0) and (layerObj.layer[1].mus == 0.0):
            self.layer = 2      # skip to next layer
            self.z = layerObj.layerZ[2][0]  # use z0 from the next layer
        
    
    def run_one_photon(self, model):
        while self.dead == False:
            self.hop_drop_spin(model)


    def hop_drop_spin(self, model):
        layer = self.layer
        if (model.layerObj.layer[layer].mua == 0.0) and \
            (model.layerObj.layer[layer].mus == 0.0):
            # glass layer
            self.hop_in_glass(model)
        else:
            self.hop_drop_spin_in_tissue(model)
        if (self.w < model.Wth) and (self.dead == False):
            self.roulette

    def hop_in_glass(self, model):
        # Move the photon packet in glass layer.
        # Horizontal photons are killed because they will
        # never interact with tissue again.
        if self.uz == 0.0:
            # horizontal photon in glass is killed
            self.dead = True
        else:
            self.step_size_in_glass(model)
            self.hop()
            self.cross_or_not(model)


    def hop_drop_spin_in_tissue(self, model):
        # Set a step size, move the photon, drop some weight, 
        # choose a new photon direction for propagation.
        # When a step size is long enough for the photon to 
        # hit an interface, this step is divided into two steps. 
        # First, move the photon to the boundary free of 
        # absorption or scattering, then decide whether the 
        # photon is reflected or transmitted.
        # Then move the photon in the current or transmission 
        # medium with the unfinished stepsize to interaction 
        # site.  If the unfinished stepsize is still too long, 
        # repeat the above process.
        self.step_size_in_tissue(model)
        if self.hit_boundary(model) == True:
            self.hop()      # move to boundary plane
            self.cross_or_not(model)
        else:
            self.hop()
            self.drop(model);
            self.spin(model.layerObj.layer[self.layer].g)

    def step_size_in_glass(self, model):
        # If uz != 0, return the photon step size in glass, 
        # Otherwise, return 0.
        # The step size is the distance between the current 
        # position and the boundary in the photon direction.
        # Make sure uz !=0 before calling this function.

        # Stepsize to the boundary.
        layer = self.layer
        uz = self.uz
        if uz > 0.0:
            dl_b = (model.layerObj.layerZ[layer][1] - self.z)/uz
        elif uz < 0.0:
            dl_b = (model.layerObj.layerZ[layer][0] - self.z)/uz
        else:
            dl_b = 0.0
  
        self.s = dl_b;


    def step_size_in_tissue(self, model):
        # Pick a step size for a photon packet when it is in 
        # tissue.
        # If the member sleft is zero, make a new step size 
        # with: -log(rnd)/(mua+mus).
        # Otherwise, pick up the leftover in sleft.
        # Layer is the index to layer.
        # In_Ptr is the input parameters.

        layer = self.layer
        mua = model.layerObj.layer[layer].mua;
        mus = model.layerObj.layer[layer].mus;

        if self.sleft == 0.0:      # make a new step
          rnd = np.random.random_sample()
          self.s = -np.log(rnd)/(mua + mus)
        else:                       # take the leftover
	        self.s = self.sleft/(mua + mus)
	        self.sleft = 0.0            

    def hop(self):
        # Move the photon s away in the current layer of medium.
        s = self.s
        self.x += s*self.ux
        self.y += s*self.uy
        self.z += s*self.uz

    def cross_or_not(self, model):
        if self.uz < 0.0:
            self.cross_up_or_not(model)
        else:
            self.cross_dn_or_not(model)

    def cross_up_or_not(self, model):
        # Decide whether the photon will be transmitted or 
        # reflected on the upper boundary (uz<0) of the current 
        # layer.
        # If "layer" is the first layer, the photon packet will 
        # be partially transmitted and partially reflected if 
        # PARTIALREFLECTION is set to 1,
        # or the photon packet will be either transmitted or 
        # reflected determined statistically if PARTIALREFLECTION 
        # is set to 0.
        # Record the transmitted photon weight as reflection.
        # If the "layer" is not the first layer and the photon 
        # packet is transmitted, move the photon to "layer-1".
        # Update the photon parmameters.

        uz = self.uz
        r = 0.0     # reflectance
        layer = self.layer
        ni = model.layerObj.layer[layer].n
        nt = model.layerObj.layer[layer-1].n
  
        # Get r.
        if -uz <= model.layerObj.cosCrit[layer][0]:
            r = 1.0		      # total internal reflection
        else:
            r, uz1 = RFresnel(ni, nt, -uz)
  
        if PARTIALREFLECTION == 1:
            if (layer == 1) and (r < 1.0):      # partially transmitted
                self.uz = -uz1          # transmitted photon
                self.record_R(r, model)
                self.uz = -uz           # reflected photon		
            elif np.random.random_sample() > r:     # transmitted to layer-1
                self.layer -= 1
                self.ux *= ni/nt
                self.uy *= ni/nt
                self.uz = -uz1
            else	:		      		# reflected
                self.uz = -uz
        else:
            if np.random.random_sample() > r:   # transmitted to layer-1
                if layer == 1:
                    self.uz = -uz1
                    self.record_R(0.0, model)
                    self.dead = True
                else:
                    self.layer -= 1
                    self.ux *= ni/nt
                    self.uy *= ni/nt
                    self.uz = -uz1
            else: 						# reflected
                self.uz = -uz


    def cross_dn_or_not(self, model):
        # Decide whether the photon will be transmitted  or be 
        # reflected on the bottom boundary (uz>0) of the current 
        # layer.
        # If the photon is transmitted, move the photon to 
        # "layer+1". If "layer" is the last layer, record the 
        # transmitted weight as transmittance. See comments for 
        # CrossUpOrNot.
        # Update the photon parmameters.
        
        uz = self.uz        # z directional cosine
        r = 0.0             # reflectance
        layer = self.layer
        ni = model.layerObj.layer[layer].n
        nt = model.layerObj.layer[layer+1].n
  
        # Get r
        if uz <= model.layerObj.cosCrit[layer][1]: 
            r=1.0		# total internal reflection
        else:
            r, uz1 = RFresnel(ni, nt, uz)
  
        if PARTIALREFLECTION	 == 1:
            if (layer == model.layerObj.numLayers) and (r < 1.0):
                self.uz = uz1
                self.record_T(r, model)
                self.uz = -uz
            elif np.random.random_sample() > r:     # transmitted to layer+1
                self.layer += 1
                self.ux *= ni/nt
                self.uy *= ni/nt
                self.uz = uz1
            else: 						# reflected
                self.uz = -uz;
        else:
            if np.random.random_sample() > r:	# transmitted to layer+1
                if layer == model.layerObj.numLayers:
                    self.uz = uz1
                    self.record_T(0.0, model)
                    self.dead = True
                else:
                    self.layer += 1
                    self.ux *= ni/nt
                    self.uy *= ni/nt
                    self.uz = uz1
            else: 						# reflected
                self.uz = -uz

    def hit_boundary(self, model):
        # Check if the step will hit the boundary.
        # Return 1 if hit boundary.
        # Return 0 otherwise.
        # If the projected step hits the boundary, the members
        # s and sleft of Photon_Ptr are updated.

        layer = self.layer
        uz = self.uz
  
        # Distance to the boundary.
        if uz > 0.0:
            dl_b = (model.layerObj.layerZ[layer][1] - self.z)/uz	    # dl_b>0
        elif uz < 0.0:
            dl_b = (model.layerObj.layerZ[layer][0] - self.z)/uz    # dl_b>0
  
        if (uz != 0.0) and (self.s > dl_b):
            # not horizontal & crossing
            mut = model.layerObj.layer[layer].mua + \
                model.layerObj.layer[layer].mus
            self.sleft = (self.s - dl_b)*mut
            self.s = dl_b
            hit = True
        else:
            hit = False
  
        return hit


    def drop(self, model):
        # Drop photon weight inside the tissue (not glass).
        # The photon is assumed not dead.
        # The weight drop is dw = w*mua/(mua+mus).
        # The dropped weight is assigned to the absorption array 
        # elements.

        x = self.x
        y = self.y
        layer = self.layer
 
        # compute array indices
        iz = int(self.z/model.dz)
        if iz > (model.nz - 1):
            iz = model.nz - 1
  
        ir = int((x*x + y*y)**0.5/model.dr)
        if ir > (model.nr - 1):
            ir = model.nr - 1
  
        # update photon weight.
        mua = model.layerObj.layer[layer].mua
        mus = model.layerObj.layer[layer].mus
        dwa = self.w * mua/(mua+mus)
        self.w -= dwa
  
        # assign dwa to the absorption array element.
        model.A_rz[ir, iz] += dwa


    def spin(self, g):
        # Choose a new direction for photon propagation by 
        # sampling the polar deflection angle theta and the 
        # azimuthal angle psi.
        # Note:
        # theta: 0 - pi so sin(theta) is always positive 
        # feel free to use sqrt() for cos(theta). 
        # psi:   0 - 2pi 
        # for 0-pi  sin(psi) is + 
        # for pi-2pi sin(psi) is - 

        ux = self.ux
        uy = self.uy
        uz = self.uz

        cost = SpinTheta(g)
        sint = (1.0 - cost*cost)**0.5	
        # sqrt() is faster than sin().

        psi = 2.0*np.pi*np.random.random_sample()       # spin psi 0-2pi
        cosp = np.cos(psi)
        if psi < np.pi:
            sinp = (1.0 - cosp*cosp)**0.5
            # sqrt() is faster than sin().
        else:
            sinp = -(1.0 - cosp*cosp)**0.5	
  
        if np.fabs(uz) > COSZERO:       # normal incident.
            self.ux = sint*cosp
            self.uy = sint*sinp
            self.uz = cost*np.sign(uz)	
            # SIGN() is faster than division.
        else:       # regular incident.
            temp = (1.0 - uz*uz)**0.5
            self.ux = sint*(ux*uz*cosp - uy*sinp)/temp + ux*cost
            self.uy = sint*(uy*uz*cosp + ux*sinp)/temp + uy*cost
            self.uz = -sint*cosp*temp + uz*cost


    def record_R(self, refl, model):
        # Record the photon weight exiting the first layer(uz<0), 
        # no matter whether the layer is glass or not, to the 
        # reflection array.
        # Update the photon weight as well.
        x = self.x
        y = self.y
  
        ir = int((x*x + y*y)**0.5/model.dr)
        irMax = int(model.nr - 1)
        if ir > irMax:
            ir = irMax
  
        ia = int(np.arccos(-self.uz)/model.da) 
        if ia > (model.na - 1):
            ia = model.na - 1
  
        # assign photon to the reflection array element.
        model.Rd_ra[ir, ia] += self.w*(1.0 - refl)
        self.w *= refl



    def record_T(self, refl, model):
        # Record the photon weight exiting the last layer(uz>0), 
        # no matter whether the layer is glass or not, to the 
        # transmittance array.
        # Update the photon weight as well.
  
        x = self.x
        y = self.y
 
        ir = int((x*x + y*y)**0.5/model.dr)
        irMax = int(model.nr - 1)
        if ir > irMax:
            ir = irMax
  
        ia = int(np.arccos(self.uz)/model.da)
        if ia > (model.na - 1):
            ia = model.na - 1
  
        # assign photon to the transmittance array element.
        model.Tt_ra[ir, ia] += self.w*(1.0-refl) 
        self.w *= refl



    def roulette(self):
        # The photon weight is small, and the photon packet tries 
        # to survive a roulette.
        if self.w == 0.0:	
            self.dead = True
        elif np.random.random_sample() < CHANCE:    # survived the roulette.
            self.w /= CHANCE
        else: 
            self.dead = True



def RFresnel(n1, n2, ca1):
    # Compute the Fresnel reflectance.
    # Make sure that the cosine of the incident angle a1
    # is positive, and the case when the angle is greater 
    # than the critical angle is ruled out.
    # Avoid trigonometric function operations as much as
    # possible, because they are computation-intensive.
  
    if n1 == n2:			  	# matched boundary
        ca2 = ca1
        r = 0.0
    elif ca1 > COSZERO:     # normal incident
        ca2 = ca1
        r = (n2-n1)/(n2+n1)
        r *= r
    elif ca1 < COS90D:      # very slant
        ca2 = 0.0
        r = 1.0
    else:           # general	
        # sine of the incident and transmission angles
        sa1 = (1.0 - ca1*ca1)**0.5
        sa2 = n1*sa1/n2
        if sa2 >= 1.0:
            # double check for total internal reflection
            ca2 = 0.0
            r = 1.0
        else:
            # cosines of the sum ap or
            # difference am of the two
            # angles. ap = a1+a2
            # am = a1 - a2     
            ca2 = (1.0 - sa2*sa2)**0.5;     
            cap = ca1*ca2 - sa1*sa2     # c+ = cc - ss
            cam = ca1*ca2 + sa1*sa2     # c- = cc + ss
            sap = sa1*ca2 + ca1*sa2     # s+ = sc + cs
            sam = sa1*ca2 - ca1*sa2     # s- = sc - cs
            r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam) 
    return r, ca2



def SpinTheta(g):
    # Choose (sample) a new theta angle for photon propagation
    # according to the anisotropy.
    # If anisotropy g is 0, then
    # cos(theta) = 2*rand-1.
    # otherwise
    # sample according to the Henyey-Greenstein function.
    # Returns the cosine of the polar deflection angle theta.
  
    if g == 0.0: 
        cost = 2*np.random.random_sample() - 1
    else:
        temp = (1 - g*g)/(1 - g + 2*g*np.random.random_sample())
        cost = (1 + g*g - temp*temp)/(2*g)
        if cost < -1:
            cost = -1.0
        elif cost > 1:
            cost = 1.0
    return cost




