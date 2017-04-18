import sys
import h5py
import pdb
import scipy as SP
import scipy.stats as ST
import scipy.linalg as LA
import time as TIME
import copy
import warnings
import os
import csv

def splitGeno(pos,method='slidingWindow',size=5e4,step=None,annotation_file=None,cis=1e4,funct=None,out_file=None):
    """
    split geno into windows and store output in csv file
    Args:
        pos:    genomic position in the format (chrom,pos)
        method:     method used to slit the windows:
                    'slidingWindow':    uses a sliding window
                    'geneWindow':       uses windows centered on genes
        size:       window size used in slidingWindow method
        step:       moving step used in slidingWindow method
        annotation_file:    file containing the annotation file for geneWindow method
        out_file:   output csv file
    """
    assert method in ['slidingWindow','geneWindow'], 'method not known'

    # create folder if does not exists
    out_dir, fname = os.path.split(out_file)
    if (out_dir!='') and (not os.path.exists(out_dir)): os.makedirs(out_dir)

    # calculates windows using the indicated method
    if method=='slidingWindow':
        nWnds,nSnps = splitGenoSlidingWindow(pos,out_file,size=size,step=step)
    elif method=='geneWindow':
        #out = splitGenoGeneWindow(pos,out_file,annotation_file=annotation_file,cis=cis,funct=funct)
        pass

    return nWnds,nSnps

def splitGenoSlidingWindow(pos,out_file,size=5e4,step=None):
    """
    split into windows using a slide criterion
    Args:
        size:       window size
        step:       moving step (default: 0.5*size)
    Returns:
        wnd_i:      number of windows
        nSnps:      vector of per-window number of SNPs
    """
    if step is None:    step = 0.5*size
    chroms = SP.unique(pos[:,0])

    RV = []
    wnd_i = 0
    wnd_file = csv.writer(open(out_file,'w'),delimiter='\t')
    nSnps = [] 
    for chrom_i in chroms:
        Ichrom = pos[:,0]==chrom_i
        idx_chrom_start = SP.where(Ichrom)[0][0]
        pos_chr = pos[Ichrom,1]
        start = pos_chr.min()
        pos_chr_max = pos_chr.max()
        while 1:
            if start>pos_chr_max: break
            end = start+size
            Ir = (pos_chr>=start)*(pos_chr<end)
            _nSnps = Ir.sum()
            if _nSnps>0:
                idx_wnd_start = idx_chrom_start+SP.where(Ir)[0][0]
                nSnps.append(_nSnps)
                line = SP.array([wnd_i,chrom_i,start,end,idx_wnd_start,_nSnps],dtype=int)
                wnd_file.writerow(line)
                wnd_i+=1
            start += step
    nSnps = SP.array(nSnps)
    return wnd_i,nSnps 

def splitGenoGeneWindow(self,annotation_file=None,cis=1e4,funct='protein_coding'):
    """
    split into windows based on genes
    """
     #1. load annotation
    assert annotation_file is not None, 'Splitter:: specify annotation file'
    try:
        f = h5py.File(annotation_file,'r')
        geneID = f['geneID'][:]
        gene_chrom = f['chrom'][:]
        gene_start = f['start'][:]
        gene_end = f['end'][:]
        gene_strand = f['strand'][:]
        gene_function = f['function'][:]
        f.close()
    except:
        print('Splitter:: format annotation file not valid')

    # if funct is not None, it has to be a list
    if funct is not None and funct!=list:   funct=[funct]

    windows = []
    nSnps   = []
    Igene   = []
    #2. calculates windows 
    for gene_i in range(geneID.shape[0]):
        if funct is not None:
            if gene_function[gene_i] not in funct:
                Igene.append(False)
                continue
        wnd = [gene_chrom[gene_i],gene_start[gene_i]-cis,gene_end[gene_i]+cis]
        Ir = (self.chrom==wnd[0])*(self.pos>=wnd[1])*(self.pos<=wnd[2])
        _nSnps = Ir.sum()
        if _nSnps>=minSnps and _nSnps<=maxSnps:
            windows.append(wnd)
            nSnps.append(_nSnps)
            Igene.append(True)
        else:
            Igene.append(False)
    Igene = SP.array(Igene)
    self.info['nSnps'] = SP.array(nSnps)
    self.info['geneID']        = geneID[Igene]
    self.info['gene_start']    = gene_start[Igene]
    self.info['gene_end']      = gene_end[Igene]
    self.info['gene_strand']   = gene_strand[Igene]
    self.info['gene_function'] = gene_function[Igene]
    return SP.array(windows)

if __name__ == "__main__":

    data  = './../data/1000G_chr22/chrom22'
    window_size = 1e4
    precompute_windows(data,size=window_size,plot=True)

