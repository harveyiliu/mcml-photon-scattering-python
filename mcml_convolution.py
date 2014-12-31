
import numpy as np
from scipy import integrate
import mcml_photon_scattering

class Beam:
    """Beam class - incident light beam class
        Parameters to describe a photon beam.
        Pencil: infinitely narrow beam. This is default for the
            beam from the mcml output.
        Flat:	Flat beam with radius R.
        Gaussian:	Gaussian with 1/e2 radius R.
        Others: general beam described by points with interpolation.
        Class instance variables:
            type - incident beam type, FLAT or GAUSSIAN
            P - total beam power/energy [W or J]
            R - beam radius, defined as 1/e^2 for Gaussian beam [cm]
        Methods:
            
    """
    
    def __init__(self, beamName = 'TRIA_HRL'):
        if beamName.lower() == 'TRIA_HRL'.lower():
            self.type = 'FLAT'
            self.P = 20     # total power. [J or W]
            self.R = 0.5    # radius. [cm]
        elif beamName.lower() == 'TRIA_FAN'.lower():
            self.type = 'GAUSSIAN'
            self.P = 0.012      # total power. [J or W]
            self.R = 0.025      # radius. [cm]
        else:
            self.type = 'FLAT'
            self.P = 20     # total power. [J or W]
            self.R = 0.5    # radius. [cm]



class Node:
    """Node class - node link list binary tree class
        Data structures for the binary tree used to store part of
        the integrand evaluation.
        
        Class instance variables:
            x - x grid node position
            y - y grid node position
            left - left node pointer
            right - right node pointer          
        Methods:
            
    """
    
    def __init__(self):
        self.x = None
        self.y = None
        self.left = None
        self.right = None
 
   
def fill_node(x, y):
    l = Node()
    l.x = x
    l.y = y
    return l


def search_node(tree, x):
    l = tree
    found = False
    while (l != None and found == False):
        if x < l.x:
            l = l.left
        elif x > l.x:
            l = l.right
        else:
            found = True
    return l

def insert_node(tree, x, y):
    l1 = None
    l2 = tree
    while l2 != None:
        l1 = l2
        if x < l2.x:
            l2 = l2.left
        else:
            l2 = l2.right

    if l1 == None:		# Empty tree.
        tree = fill_node(x, y)
    elif x < l1.x:
        l1.left = fill_node(x, y)
    else:
        l1.right = fill_node(x, y)



class ConvVar:
    """ConvVar class - convoluation variables class
        A global structure to pass the current coordinate of the
        physical quantities being evaluated and the pointers of the
        input and output parameters to the integration function.
        
        Class instance variables:
            r - r position
            iz - iz index
            ia - ia index
            tree - A tree to store ITheta() & ExpBessI0().           
        Methods:
            
    """
    
    def __init__(self):
        self.r = 0
        self.iz = 0
        self.ia = 0
        self.tree = None      # A tree to store ITheta() & ExpBessI0().


class ConvInput:
    """ConvInput class - beam convolution input class
        Input parameters for each independent run.
        z and r are for the cylindrical coordinate system. [cm]
        a is for the angle alpha between the photon exiting
        direction and the surface normal. [radian]
        The grid line separations in z, r, and alpha
        directions are dz, dr, and da respectively.  The numbers
        of grid lines in z, r, and alpha directions are
        nz, nr, and na respectively.
        The member layerspecs will point to an array of
        structures which store parameters of each layer.
        This array has (number_layers + 2) elements. One
        element is for a layer.
        The layers 0 and (num_layers + 1) are for top ambient
        medium and the bottom ambient medium respectively.
        For convolution, the grid line separations in z, and alpha
        directions are still dz, and da respectively.  The numbers
        of grid lines in z, and alpha directions are still
        nz, and na respectively. However, the grid line separation
        and the number of grid lines in r direction are drc and
        nrc respectively.
        Class instance variables:
            beam - incident beam class instance object
            drc - convolution r grid separation.[cm]
            nrc - convolution array range 0..nrc-1.
            eps - relative error in convolution           
        Methods:
            
    """
    
    def __init__(self, mcmlModel, convName = 'TRIA_HRL'):
        if convName.lower() == 'TRIA_HRL'.lower():
            self.beam = Beam('TRIA_HRL')    # incident beam of finite size
            self.drc = 0.005        # convolution r grid separation.[cm]
            self.nrc = 150          # convolution array range 0..nrc-1.
               
        elif convName.lower() == 'TRIA_FAN'.lower():
            self.beam = Beam('TRIA_FAN')
            self.drc = 0.001
            self.nrc = 100 
        else:
            self.beam = Beam('TRIA_HRL')
            self.drc = 0.005
            self.nrc = 200
        self.mcmlModel = mcmlModel
        self.convVar = ConvVar()
            


class MCMLConv(ConvInput):
    """MCMLConv class - multi-layered photon scattering model beam convolution
        inherits from ConvInput beam setup
        Structures for scored physical quantities
        from mcml and to be convolved for photon
        beams of finite size.  Therefore, "Out"
        here means the output of both mcml and conv.
        The member allocated is used to keep the status
        of the arrays.  It is set to 1 if all the arrays
        are allocated and assigned values.  It is set to
        0 otherwise.
        z and r represent z and r coordinates of the
        cylindrical coordinate system. [cm]
        a is the angle alpha between the photon exiting
        direction and the normal to the surfaces. [radian]
        See comments of the InputStruct.
        See manual for the physcial quantities.
        Class instance variables:
            Rd_rac - convolved data. [J/(cm2 sr)]
            Rd_rc - 1D radial distribution of diffuse reflectance [J/cm2]
            A_rzc - 2D probability density in turbid media over r & z [J/cm3]
            Tt_rac - 2D distribution of total transmittance [J/(cm2 sr)]
            Tt_rc - 1D radial distribution of transmittance [J/cm2]
        Methods:
            
    """

    def __init__(self, mcmlModel, convName = 'TRIA_HRL'):
        # extend the ConvInput base class instance variables
        ConvInput.__init__(self, mcmlModel, convName)
        # initialize the model grid arrays    
        self.Rd_rac = np.matrix(np.zeros((self.nrc, self.mcmlModel.na)))
        self.Rd_rc = np.zeros(self.nrc)
        self.A_rzc = np.matrix(np.zeros((self.nrc, self.mcmlModel.nz)))
        self.Tt_rac = np.matrix(np.zeros((self.nrc, self.mcmlModel.na)))
        self.Tt_rc = np.zeros(self.nrc)
        self.F_rzc = np.matrix(np.zeros((self.nrc, self.mcmlModel.nz)))

    
    def run_conv(self):
        self.conv_Rd_ra()
        self.conv_Rd_r()
        self.conv_A_rz()
        self.conv_Tt_ra()
        self.conv_Tt_r()
        self.conv_A2F()

    
    def conv_Rd_ra(self):
        P = self.beam.P
        R = self.beam.R
        for irc in range(self.nrc):
            rc = (irc + 0.5)*self.drc
            self.convVar.r = rc
            self.convVar.tree = None    	# init the tree
            for ia in range(self.mcmlModel.na):
                self.convVar.ia = ia
                if self.beam.type.lower() == 'FLAT'.lower():
                    self.Rd_rac[irc, ia] = 2*P/(R*R) \
                        *flat_integration(Rd_ra_FG_integrand, self)
                else	:       # Gaussian
                    self.Rd_rac[irc, ia] = 4*P/(R*R) \
                        *Gauss_integration(Rd_ra_FG_integrand, self)


    def conv_Rd_r(self):
        P = self.beam.P
        R = self.beam.R
        for irc in range(self.nrc):
            rc = (irc + 0.5)*self.drc
            self.convVar.r = rc
            if self.beam.type.lower() == 'FLAT'.lower():
                self.Rd_rc[irc] = 2*P/(R*R) \
                    *flat_integration(Rd_r_FG_integrand, self)
            else	:       # Gaussian
                self.Rd_rc[irc] = 4*P/(R*R) \
                    *Gauss_integration(Rd_r_FG_integrand, self)



    def conv_A_rz(self):
        P = self.beam.P
        R = self.beam.R
        for irc in range(self.nrc):
            rc = (irc + 0.5)*self.drc
            self.convVar.r = rc
            self.convVar.tree = None    	# init the tree
            for iz in range(self.mcmlModel.nz):
                self.convVar.iz = iz
                if self.beam.type.lower() == 'FLAT'.lower():
                    self.A_rzc[irc, iz] = 2*P/(R*R) \
                        *flat_integration(A_rz_FG_integrand, self)
                else:       # Gaussian
                    self.A_rzc[irc, iz] = 4*P/(R*R) \
                        *Gauss_integration(A_rz_FG_integrand, self)
    


    def conv_Tt_ra(self):
        P = self.beam.P
        R = self.beam.R
        for irc in range(self.nrc):
            rc = (irc + 0.5)*self.drc
            self.convVar.r = rc
            self.convVar.tree = None        # init the tree
            for ia in range(self.mcmlModel.na):
                self.convVar.ia = ia
                if self.beam.type.lower() == 'FLAT'.lower():
                    self.Tt_rac[irc, ia] = 2*P/(R*R) \
                        *flat_integration(Tt_ra_FG_integrand, self)
                else:       # Gaussian
                    self.Tt_rac[irc, ia] = 4*P/(R*R) \
                        *Gauss_integration(Tt_ra_FG_integrand, self)


    def conv_Tt_r(self):
        P = self.beam.P
        R = self.beam.R
        for irc in range(self.nrc):
            rc = (irc + 0.5)*self.drc
            self.convVar.r = rc
            if self.beam.type.lower() == 'FLAT'.lower():
                self.Tt_rc[irc] = 2*P/(R*R) \
                    *flat_integration(Tt_r_FG_integrand, self)
            else	:       # Gaussian
                self.Tt_rc[irc] = 4*P/(R*R) \
                    *Gauss_integration(Tt_r_FG_integrand, self)

    def conv_A2F(self):
        nz = self.mcmlModel.nz
        for ir in range(self.nrc):
            for iz in range(nz):
                mua = self.mcmlModel.layerObj.layer[self.iz_to_layer(iz)].mua
                if (mua > 0.0):
	                self.F_rzc[ir, iz] = self.A_rzc[ir, iz]/mua     # F in J/cm2


    def iz_to_layer(self, iz):
        i = 1       # index to layer
        numLayers = self.mcmlModel.layerObj.numLayers
        dz = self.mcmlModel.dz
        while ((iz + 0.5)*dz >= self.mcmlModel.layerObj.layerZ[i][1] \
                and i < numLayers):
            i += 1
        return i

    
    def center_half_max_depth(self):
        nz = self.mcmlModel.nz
        dz = self.mcmlModel.dz
        depth = 0
        for iz in range(nz):
            if self.F_rzc[0, iz] <= 0.5*self.F_rzc[0, 0]:
                depth = (iz + 0.5)*dz
                break
        return depth


    def surface_half_max_width(self):
        width = 0
        for irc in range(self.nrc):
            if self.F_rzc[irc, 0] <= 0.5*self.F_rzc[0, 0]:
                width = 2*(irc + 0.5)*self.drc
                break
        return width    
            

      


def I_theta(r, r2, R):
    if R >= (r + r2):
        temp = 1
    elif np.fabs(r - r2) <= R:
        temp = (r*r + r2*r2 - R*R)/(2*r*r2)
        if np.fabs(temp) > 1:
            temp = np.sign(temp)
        temp = np.arccos(temp)/np.pi
    else	:		# R < fabs(r-r2)
        temp = 0
    return temp



def exp_Bess_I0(r, r2, R):
    _RR = 1/(R*R)
    x = 4*r*r2*_RR
    y = 2*(r2*r2 + r*r)*_RR
    expbess = np.exp(-y + x)*np.i0(x)
    return (expbess)


def RT_ra_interp(r2, RT_ra, mcmlConv):
# Interpolate for the arrays Rd_ra[] or Tt_ra[].
    nr = mcmlConv.mcmlModel.nr

    ir2 = r2/mcmlConv.mcmlModel.dr
    ia = mcmlConv.convVar.ia
    if nr < 3:
        RT_at_r2 = RT_ra[0, ia]
    elif ir2 < (nr - 1.5):      	# interpolation
        ir2lo = np.maximum(0, int(ir2 - 0.5))	    # truncation
        RT_lo = RT_ra[ir2lo, ia]
        RT_hi = RT_ra[ir2lo + 1, ia]
        RT_at_r2 = RT_lo + (RT_hi - RT_lo)*(ir2 - ir2lo - 0.5)
    else:			# extrapolation
        ir2lo = nr - 3
        RT_lo = RT_ra[ir2lo, ia]
        RT_hi = RT_ra[ir2lo + 1, ia]
        if RT_lo >= RT_hi:		# Noise test
            RT_at_r2 = RT_lo + (RT_hi - RT_lo)*(ir2 - ir2lo - 0.5)
        else:
            RT_at_r2 = 0.0
    return np.maximum(0, RT_at_r2)




def RT_r_interp(r2, RT_r, mcmlConv):
# Interpolate for the arrays Rd_r[] or Tt_r[].
    nr = mcmlConv.mcmlModel.nr

    ir2 = r2/mcmlConv.mcmlModel.dr
    if nr < 3:
        RT_at_r2 = RT_r[0]
    elif ir2 < (nr - 1.5):      # interpolation
        ir2lo = np.maximum(0, int(ir2 - 0.5))       # truncation
        RT_lo = RT_r[ir2lo]
        RT_hi = RT_r[ir2lo + 1]
        RT_at_r2 = RT_lo + (RT_hi - RT_lo)*(ir2 - ir2lo - 0.5)
    else:           # extrapolation
        ir2lo = nr - 3
        RT_lo = RT_r[ir2lo]
        RT_hi = RT_r[ir2lo + 1]
        if RT_lo >= RT_hi:      # Noise test
            RT_at_r2 = RT_lo + (RT_hi - RT_lo)*(ir2 - ir2lo - 0.5)
        else:
            RT_at_r2 = 0.0
    return np.maximum(0, RT_at_r2)




def A_rz_interp(r2, mcmlConv):
# Interpolate for the arrays A_rz[]
    A_rz = mcmlConv.mcmlModel.A_rz
    nr = mcmlConv.mcmlModel.nr

    ir2 = r2/mcmlConv.mcmlModel.dr
    iz = mcmlConv.convVar.iz
    if nr < 3:
        A_at_r2 = A_rz[0, iz]
    elif ir2 < (nr - 1.5):       # interpolation
        ir2lo = np.maximum(0, int(ir2 - 0.5))   # truncation
        A_lo = A_rz[ir2lo, iz]
        A_hi = A_rz[ir2lo + 1, iz]
        A_at_r2 = A_lo + (A_hi - A_lo)*(ir2 - ir2lo - 0.5)
    else:       # extrapolation
        ir2lo = nr - 3
        A_lo = A_rz[ir2lo, iz]
        A_hi = A_rz[ir2lo + 1, iz]
        if A_lo >= A_hi:       # Noise test
            A_at_r2 = A_lo + (A_hi - A_lo)*(ir2 - ir2lo - 0.5)
        else:
            A_at_r2 = 0.0
    return np.maximum(0, A_at_r2)







def Rd_ra_FG_integrand(r2, mcmlConv):
    # Convolution integrand for either flat or gaussian beams.
    # See comments for A_rzFGIntegrand().
    # r" in the integration.

    RT_ra = mcmlConv.mcmlModel.Rd_ra
    Rd_at_r2 = RT_ra_interp(r2, RT_ra, mcmlConv)
    R = mcmlConv.beam.R
    r = mcmlConv.convVar.r
    tree = mcmlConv.convVar.tree
    link = search_node(tree, r2)
    if link != None:	    # f in tree.
        f = link.y
    else:
        if mcmlConv.beam.type.lower() == 'FLAT'.lower():
            f = I_theta(r, r2, R)
        else:			# Gaussian
            f = exp_Bess_I0(r, r2, R)
        insert_node(tree, r2, f)

    f *= Rd_at_r2*r2
    return f





def Rd_r_FG_integrand(r2, mcmlConv):
# Convolution integrand for either flat or gaussian beams.
# See comments for A_rzFGIntegrand().
# r" in the integration

    RT_r = mcmlConv.mcmlModel.Rd_r
    Rd_at_r2 = RT_r_interp(r2, RT_r, mcmlConv)
    R = mcmlConv.beam.R
    r = mcmlConv.convVar.r
    if mcmlConv.beam.type.lower() == 'FLAT'.lower():
        f = Rd_at_r2*I_theta(r, r2, R)*r2
    else	:       # Gaussian
        f = Rd_at_r2*exp_Bess_I0(r, r2, R)*r2
    return f




def A_rz_FG_integrand(r2, mcmlConv):
# Convolution integrand for either flat or gaussian beams.
# Return the integrand for the convolution integral.
# r2 is the r" in the formula shown in the manual.
# When r2 is in the range of recorded array, interpolation
# is used to evaluate the diffuse reflectance at r2.
# Note that since the last grid elements collect all the
# photon weight that falls beyond the grid system, we should
# avoid using them in the convolution.
# r" in the integration

    A_at_r2 = A_rz_interp(r2, mcmlConv)
    R = mcmlConv.beam.R
    r = mcmlConv.convVar.r
    tree = mcmlConv.convVar.tree
    link = search_node(tree, r2)
    if link != None:        # f in tree
        f = link.y
    else:
        if mcmlConv.beam.type.lower() == 'FLAT'.lower():
            f = I_theta(r, r2, R)
        else:       # Gaussian
            f = exp_Bess_I0(r, r2, R)
        insert_node(tree, r2, f)

    f *= A_at_r2*r2
    return f




def Tt_ra_FG_integrand(r2, mcmlConv):
# Convolution integrand for either flat or gaussian beams.
# See comments for A_rzFGIntegrand().
# r" in the integration.

    TT_ra = mcmlConv.mcmlModel.Tt_ra
    Tt_at_r2 = RT_ra_interp(r2, TT_ra, mcmlConv)
    R = mcmlConv.beam.R
    r = mcmlConv.convVar.r
    tree = mcmlConv.convVar.tree
    link = search_node(tree, r2)
    if link != None:        # f in tree
        f = link.y
    else:
        if mcmlConv.beam.type.lower() == 'FLAT'.lower():
            f = I_theta(r, r2, R)
        else	:       # Gaussian
            f = exp_Bess_I0(r, r2, R)
        insert_node(tree, r2, f)

    f *= Tt_at_r2*r2
    return f




def Tt_r_FG_integrand(r2, mcmlConv):
# Convolution integrand for either flat or gaussian beams.
# See comments for A_rzFGIntegrand().
# r" in the integration

    TT_r = mcmlConv.mcmlModel.Tt_r
    Tt_at_r2 = RT_r_interp(r2, TT_r, mcmlConv)
    R = mcmlConv.beam.R
    r = mcmlConv.convVar.r
    if mcmlConv.beam.type.lower() == 'FLAT'.lower():
        f = Tt_at_r2*I_theta(r, r2, R)*r2
    else	:       # Gaussian
        f = Tt_at_r2*exp_Bess_I0(r, r2, R)*r2
    return f



def flat_integration(func, mcmlConv):
    rc = mcmlConv.convVar.r
    R = mcmlConv.beam.R
    b_max = (mcmlConv.mcmlModel.nr - 0.5)*mcmlConv.mcmlModel.dr
    a = np.maximum(0, rc - R)
    b = np.minimum(b_max, rc + R)

    if (a >= b):
        return 0
    else:
        return integrate.quad(func, a, b, args=(mcmlConv,), \
            epsabs=1.0e-5, epsrel=1.0e-5, limit=500)[0]



def Gauss_integration(func, mcmlConv):
# Used by convolution over Gaussian beam.  Ignore the value
# beyond GAUSSLIMIT radius.

    GAUSSLIMIT = 4
    rc = mcmlConv.convVar.r
    R = mcmlConv.beam.R
    b_max = (mcmlConv.mcmlModel.nr - 0.5)*mcmlConv.mcmlModel.dr
    a = np.maximum(0, (rc - GAUSSLIMIT * R))
    b = np.minimum(b_max, (rc + GAUSSLIMIT * R))

    if (a >= b):
        return 0
    else:
        return integrate.quad(func, a, b, args=(mcmlConv,), \
            epsabs=1.0e-5, epsrel=1.0e-5, limit=500)[0]


            
