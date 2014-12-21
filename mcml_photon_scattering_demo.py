#!/usr/bin/python
import mcml_photon_scattering as mcml


if __name__ == "__main__":
    """ MCML photon scattering modeling demo
            Input arguments:
                modelName = predefined model name, current set of models are:
                    1. BARE_DERMIS (800-nm) - default
                    2. TYPE_II_SKIN (800-nm)
                    3. CORNEA (1060-nm)
                N = number of photons to be used in the Monte-Carlo simulation
        Run command example:
            ./mcml_photon_scattering_demo.py TYPE_II_SKIN 1000 &      
    """
    
    import sys
    modelName = sys.argv[1]     # MCML model name
    N = int(sys.argv[2])        # number of photons used for Monte-Carlo
    model = mcml.MCMLModel(modelName)
    model.do_one_run(N)
    model.sum_scale_result()
    print '\nPhoton number: {0:d}'.format(model.numPhotons)
    print 'Rsp: {0:.3f}'.format(model.Rsp)
    print 'Rd: {0:.3f}'.format(model.Rd)
    print 'A: {0:.3f}'.format(model.A)
    print 'Tt: {0:.3f}'.format(model.Tt)
    
