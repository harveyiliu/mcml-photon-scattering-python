#!/usr/bin/python
import mcml_photon_scattering as mcml
import mcml_convolution as conv
import time


if __name__ == "__main__":
    """ MCML photon scattering and beam convolution demo
            Input arguments:
                modelName = predefined model name, current set of models are:
                    1. BARE_DERMIS (800-nm) - default
                    2. TYPE_II_SKIN (800-nm)
                    3. CORNEA (1060-nm)
                    4. EYE_ANTERIOR (1060-nm)
                N = number of photons to be used in the Monte-Carlo simulation
                beamType = predefined beam type,
                    1. FLAT
                    2. GAUSSIAN
                P = incident beam power/energy (W or J)
                W = beam width (mm)
        Run command example:
            ./mcml_convolution_demo.py BARE_DERMIS 1000 FLAT 1.0 10.0 &      
    """
    
    import sys
    modelName = sys.argv[1]     # MCML model name
    N = int(sys.argv[2])        # number of photons used for Monte-Carlo
    beamType = sys.argv[3]      # beam type
    P = float(sys.argv[4])      # beam power/energy [W or J]
    W = float(sys.argv[5])      # beam width [mm]
    R = 0.1*W/2                 # beam radius [cm]
    
    timeStart = time.time()
    model = mcml.MCMLModel(modelName)
    model.do_one_run(N)
    model.sum_scale_result()
    print '\nPhoton number: {0:d}'.format(model.numPhotons)
    print 'Specular reflection Rsp (%): {0:.2f}'.format(100*model.Rsp)
    print 'Diffused reflection Rd (%): {0:.2f}'.format(100*model.Rd)
    print 'Absorption A (%): {0:.2f}'.format(100*model.A)
    print 'Transmission Tt (%): {0:.2f}'.format(100*model.Tt)
    
    mcmlConv = conv.MCMLConv(model, beamType, P, R)
    mcmlConv.run_conv()
    halfMaxDepth = mcmlConv.center_half_max_depth()
    halfMaxWidth = mcmlConv.surface_half_max_width()
    runTime = time.time() - timeStart   # total runtime in [sec.]
    print '\nCenter half max depth (mm): {0:.3f}'.format(10*halfMaxDepth)
    print 'Surface half max width (mm): {0:.3f}'.format(10*halfMaxWidth)
    print 'Total run time (min.): {0:.2f}'.format(runTime/60)

    
