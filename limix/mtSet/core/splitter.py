import sys
sys.path.append('./..')
import h5py
import pdb
import scipy as SP
import scipy.stats as ST
import scipy.linalg as LA
import time as TIME
import copy
import warnings
import os


class Splitter():

    def __init__(self,pos=None,chrom=None,separate_chroms=False):
        """
        Constructor
        Args:
            pos:        position
            chrom:      chromosome
        """
        # assert
        assert pos is not None, 'Slider:: set pos'
        assert chrom is not None, 'Slider:: set chrom'

        self.pos   = pos
        self.chrom = chrom
        # sep chroms
        self.separate_chroms = separate_chroms
        # windows
        self.windows = None
        # additional info for windows
        self.info   = {}
        pass

    def splitGeno(self,method='slidingWindow',size=5e4,step=None,annotation_file=None,cis=1e4,funct=None,minSnps=1.,maxSnps=SP.inf,cache=False,out_dir='./cache',fname=None,rewrite=False):
        """
        split into windows
        Args:
            method:     method used to slit the windows:
                        'slidingWindow':    uses a sliding window
                        'geneWindow':       uses windows centered on genes
            size:       window size used in slidingWindow method
            step:       moving step used in slidingWindow method
            annotation_file:    file containing the annotation file for geneWindow method
            minSnps:    only windows with nSnps>=minSnps are considered
            maxSnps:    only windows with nSnps>=maxSnps are considered
            cache:      if results need to be cashed
            out_dir:    folder used for caching
            fname:      name of the hdf5 file used for caching
            rewrite:    if true, rewrite the existing cache hdf5 file
        """
        self.info= {} # forget everything
        # check if it is necessary to read form file or not
        read_from_file = False
        if cache:
            assert fname is not None, 'Splitter:: specify fname'
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            out_file = os.path.join(out_dir,fname)
            read_from_file = os.path.exists(out_file) and not rewrite

        self.size = size

        if read_from_file:
            # reads from file
            f = h5py.File(out_file,'r')
            self._wnd_pos = f['wnd_pos'][:]
            self._idx_wnd_start = f['idx_wnd_start'][:]
            self._nSnps = f['nSnps'][:]
            f.close()
        else:
            # calculates windows using the indicated method
            if method=='slidingWindow':
                self._splitGenoSlidingWindow(size=size,step=step,minSnps=minSnps,maxSnps=maxSnps)
            elif method=='geneWindow':
                self._splitGenoGeneWindow(annotation_file=annotation_file,cis=cis,funct=funct,minSnps=minSnps,maxSnps=maxSnps)

            if cache:
                f = h5py.File(out_file,'w')
                f.create_dataset('wnd_pos', data=self._wnd_pos)
                f.create_dataset('idx_wnd_start', data=self._idx_wnd_start)
                f.create_dataset('nSnps', data=self._nSnps)
                f.close()

    def _splitGenoSlidingWindow(self,size=5e4,step=None,minSnps=1.,maxSnps=SP.inf):
        """
        split into windows using a slide criterion
        Args:
            size:       window size
            step:       moving step (default: 0.5*size)
            minSnps:    only windows with nSnps>=minSnps are considered
            maxSnps:    only windows with nSnps>=maxSnps are considered
        """
        if step is None:    step = 0.5*size
        chroms  = SP.unique(self.chrom)
        wnd_pos       = []
        idx_wnd_start = []
        nSnps         = []
        wnd_i = 0

        nSnps = []
        for chrom_i in chroms:
            start = 0
            Ichrom = self.chrom==chrom_i
            idx_chrom_start = SP.where(Ichrom)[0][0]
            pos_chr = self.pos[Ichrom]
            pos_chr_max = pos_chr.max()
            while 1:
                if start>pos_chr_max: break
                end = start+size
                Ir = (self.pos>=start)*(self.pos<end)
                _nSnps = Ir.sum()
                if _nSnps>minSnps and _nSnps<maxSnps:
                    wnd_pos.append([chrom_i,start,start+size])
                    nSnps.append(_nSnps)
                    idx_wnd_start.append(idx_chrom_start+SP.where(Ir)[0][0])
                    wnd_i+=1
                start += step
        self._wnd_pos = SP.array(wnd_pos)
        self._idx_wnd_start = SP.array(idx_wnd_start)
        self._nSnps = SP.array(nSnps)

    def _splitGenoGeneWindow(self,annotation_file=None,cis=1e4,funct='protein_coding',minSnps=1.,maxSnps=SP.inf):
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

    @property
    def wnd_pos(self):
        return self._wnd_pos

    @property
    def idx_wnd_start(self):
        return self._idx_wnd_start

    @property
    def nSnps(self):
        return self._nSnps

    @property
    def nWindows(self):
        return self.wnd_pos.shape[0]

