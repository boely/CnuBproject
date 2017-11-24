import os
import os.path
import ConfigParser
import string
import sys

sys.path.append("src/.")
from NeutrinoFlux import *

def usage():
    print "Usage:  %s  <config-file>  <outfile> \n" % os.path.basename
(sys.argv[0])

def main():
    args = sys.argv[1:]
    print len(args)
    if "-h" in args or "--help" in args or len(args) < 2:
        usage()
        sys.exit(2)

    config = ConfigParser.RawConfigParser()
    config.read(args[0])
    sys.stderr.write('Reading config file %s \n'% args[0])

    flux = NeutrinoFlux()
    flux.z_max = config.getfloat("Flux", "z_max")
    flux.alpha = config.getfloat("Flux", "alpha") 
    flux.n     = config.getfloat("Flux", "n") 
    flux.m_n   = config.getfloat("Flux", "m_n") 
    flux.eta0  = config.getfloat("Flux", "eta0") 
    flux.j     = config.getfloat("Flux", "j") 
    flux.Eresolution = config.getfloat("Telescope", "energy_resolution")

    prefix_outfilename = args[1]
    flux.save_to_ascii(prefix_outfilename)

if __name__ == "__main__":
    sys.exit(main())
