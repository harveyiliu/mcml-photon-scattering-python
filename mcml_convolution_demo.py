#!/usr/bin/python
import mcml_photon_scattering as mcml
import mcml_convolution as conv


if __name__ == "__main__":
    """ MCML photon scattering and beam convolution demo
            Input arguments:
                modelName = predefined model name, current set of models are:
                    1. BARE_DERMIS (800-nm) - default
                    2. TYPE_II_SKIN (800-nm)
                    3. CORNEA (1060-nm)
                    4. EYE_ANTERIOR (1060-nm)
                N = number of photons to be used in the Monte-Carlo simulation
                convName = predefined beam convolution system name,
                    1. TRIA_HRL
                    2. TRIA_FAN
        Run command example:
            ./mcml_convolution_demo.py TYPE_II_SKIN 1000 TRIA_HRL &      
    """
    
    import sys
    modelName = sys.argv[1]     # MCML model name
    N = int(sys.argv[2])        # number of photons used for Monte-Carlo
    convName = sys.argv[3]      # convolution system name
    model = mcml.MCMLModel(modelName)
    model.do_one_run(N)
    model.sum_scale_result()
    print '\nPhoton number: {0:d}'.format(model.numPhotons)
    print 'Specular reflection Rsp (%): {0:.2f}'.format(100*model.Rsp)
    print 'Diffused reflection Rd (%): {0:.2f}'.format(100*model.Rd)
    print 'Absorption A (%): {0:.2f}'.format(100*model.A)
    print 'Transmission Tt (%): {0:.2f}'.format(100*model.Tt)
    
    mcmlConv = conv.MCMLConv(model, convName)
    mcmlConv.run_conv()
    halfMaxDepth = mcmlConv.center_half_max_depth()
    halfMaxWidth = mcmlConv.surface_half_max_width()
    print '\nCenter half max depth (mm): {0:.3f}'.format(10*halfMaxDepth)
    print 'Surface half max width (mm): {0:.2f}'.format(10*halfMaxWidth)

    
