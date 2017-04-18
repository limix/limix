from optparse import OptionParser
import scipy as SP

import os
from . import simulator as sim
from .read_utils import readCovarianceMatrixFile
from .read_utils import readBimFile 

def genPhenoCube(sim,Xr,vTotR=4e-3,nCausalR=10,pCommonR=0.8,vTotBg=0.4,pHidd=0.6,pCommon=0.8):
    # region
    nCommonR  = int(SP.around(nCausalR*pCommonR))
    # background
    vCommonBg = pCommon*vTotBg
    # noise
    vTotH   = pHidd*(1-vTotR-vTotBg)
    vTotN   = (1-pHidd)*(1-vTotR-vTotBg)
    vCommonH = pCommon*vTotH
    all_settings = {
        'vTotR':vTotR,'nCommonR':nCommonR,'nCausalR':nCausalR,
        'vTotBg':vTotBg,'vCommonBg':vCommonBg,'pCausalBg':1.,'use_XX':True,
        'vTotH':vTotH,'vCommonH':vCommonH,'nHidden':10,
        'vTotN':vTotN,'vCommonN':0.}

    Y,info = sim.genPheno(Xr,**all_settings)
    return Y,info

def simPheno(options):

    print('importing covariance matrix')
    if options.cfile is None: options.cfile=options.bfile
    XX = readCovarianceMatrixFile(options.cfile,readEig=False)['K']

    print('simulating phenotypes')
    SP.random.seed(options.seed)
    simulator = sim.CSimulator(bfile=options.bfile,XX=XX,P=options.nTraits)
    Xr,region = simulator.getRegion(chrom_i=options.chrom,size=options.windowSize,min_nSNPs=options.nCausalR,pos_min=options.pos_min,pos_max=options.pos_max)
 
    Y,info    = genPhenoCube(simulator,Xr,vTotR=options.vTotR,nCausalR=options.nCausalR,pCommonR=options.pCommonR,vTotBg=options.vTotBg,pHidd=options.pHidden,pCommon=options.pCommon)

    print('exporting pheno file')
    if options.pfile is not None:
        outdir = os.path.split(options.pfile)[0]
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    else:
        identifier = '_seed%d_nTraits%d_wndSize%d_vTotR%.2f_nCausalR%d_pCommonR%.2f_vTotBg%.2f_pHidden%.2f_pCommon%.2f'%(options.seed,options.nTraits,options.windowSize,options.vTotR,options.nCausalR,options.pCommonR,options.vTotBg,options.pHidden,options.pCommon)
        options.pfile = os.path.split(options.bfile)[-1] + '%s'%identifier

    pfile  = options.pfile + '.phe'
    rfile  = options.pfile + '.phe.region'

    SP.savetxt(pfile,Y)
    SP.savetxt(rfile,region)
