#!/usr/bin/python
import mcml_photon_scattering as mcml


if __name__ == "__main__":
    """ MCML photon scattering modeling demo
            Input arguments:
                N = number of photons to be used in the Monte-Carlo simulation
    """
    
    import sys
    N = int(sys.argv[1])        # number of photons used for Monte-Carlo
    model = mcml.MCMLModel('TYPE_II_SKIN')
    model.do_one_run(N)
    model.sum_scale_result()
    print '\nRsp: {0:.3f}'.format(model.Rsp)
    print 'Rd: {0:.3f}'.format(model.Rd)
    print 'A: {0:.3f}'.format(model.A)
    print 'Tt: {0:.3f}'.format(model.Tt)
    
