# Copyright(c) 2014, The LIMIX developers (Christoph Lippert, Paolo Francesco Casale, Oliver Stegle)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# modified by Barbara Rakitsch
import pdb
import os
import numpy as SP
import csv

def standardize(snps): 
    mdat = SP.ma.masked_array(snps,SP.isnan(snps))  # mask nan values
    N = (mdat - SP.mean(mdat,axis=0)) / SP.std(mdat,axis=0)  # standardize column-wise to mean 0 and std 1
    return N.filled(0.0)  # assign 0.0 to masked missing values

def readBIM(basefilename,usecols=None):
    """
    helper method for speeding up read BED
    """
    bim = basefilename+ '.bim'
    bim = SP.loadtxt(bim,dtype=bytes,usecols=usecols)
    return bim
    

def readFAM(basefilename,usecols=None):
    """
    helper method for speeding up read FAM
    """
    fam = basefilename+'.fam'
    fam = SP.loadtxt(fam,dtype=bytes,usecols=usecols)
    return fam


def readBED(basefilename, useMAFencoding=False,blocksize = 1, start = 0, nSNPs = SP.inf, startpos = None, endpos = None, order  = 'F',standardizeSNPs=False,ipos = 2,bim=None,fam=None):
    '''
    read [basefilename].bed,[basefilename].bim,[basefilename].fam
    --------------------------------------------------------------------------
    Input:
    basefilename    : string of the basename of [basename].bed, [basename].bim,
                      and [basename].fam
    blocksize       : load blocksize SNPs at a time (default 1)
    start           : index of the first SNP to be loaded from the .bed-file
                      (default 0)
    nSNPs           : load nSNPs from the .bed file (default SP.inf, meaning all)
    startpos        : starting position of the loaded genomic region[chr,bpdist]
    endpos          : end-position of the loaded genomic region     [chr,bpdist]
    order           : memory layout of the returned SNP array (default 'F')
                      'F'   : Fortran-style column-major array (SNP-major)
                      'C'   : C-style row-major array (individual-major)
    standardizeSNPs : bool indeicator if the resulting SNP array is supposed to 
                      be zero-mean and unit-vatiance with mean imputed missing
                      values (default False)
    ipos            : the index of the position index to use (default 2)
                        1 : genomic distance
                        2 : base-pair distance
    useMAFencoding  : if set to one, the minor allele is encoded with 2, the major allele with 0.
                      otherwise, the plink coding is used (default False).
    --------------------------------------------------------------------------
    Output dictionary:
    'rs'     : [S] array rs-numbers
    'pos'    : [S*3] array of positions [chromosome, genetic dist, basepair dist]
    'snps'   : [N*S] array of snp-data
    'iid'    : [N*2] array of family IDs and individual IDs
    --------------------------------------------------------------------------    
    '''
    
    if bim is None: bim = readBIM(basefilename,usecols=(0,1,2,3))
    if fam is None: fam = readFAM(basefilename,usecols=(0,1))

    
    rs = bim[:,1]
    pos = SP.array(bim[:,(0,2,3)],dtype = 'float')


    
    if startpos is not None:
        #pdb.set_trace()
        i_c = pos[:,0]==startpos[0]
        i_largerbp = pos[:,ipos]>=startpos[ipos]
        start = which(i_c * i_largerbp)
        while (start-1 >= 0 and pos[start-1,ipos] == startpos[ipos]):
            start = start -1
        i_c = pos[:,0]==endpos[0]
        i_smallerbp = pos[:,ipos]>=endpos[ipos]
        end = which(i_c * i_smallerbp)
        while (end+1 < pos.shape[0] and pos[end+1,ipos] == endpos[ipos]):
            end = end + 1
        nSNPs = end - start
        if (nSNPs<=0) or (end==0) or (start<=0):
            ret = {
                'pos':SP.zeros((0,3)),
                'rs':SP.zeros((0)),
                'iid':fam,
                'snps':SP.zeros((fam.shape[0],0))
                }
            return ret
        pass
    N = fam.shape[0]
    S = bim.shape[0]
    S_res = min(S,start + nSNPs)
    nSNPs = min(S-start,nSNPs)
    #if startpos is not None:
	#print("start: " + str(start))
	#print("end: " + str(end))
    #print("S_res: " + str(S_res))
    #print("nSNPs: " + str(nSNPs))
    if nSNPs<=0:
        ret = {
            'rs'     :rs[start:start],
            'pos'    :pos[start:start,:],
            #'snps'   :SNPs[0:N,start:start],
            'snps'   :SP.zeros((N,0)),
            'iid'    :fam
            }
        return ret
    SNPs = SP.zeros(((SP.ceil(0.25*N)*4),nSNPs),order=order)
    bed = basefilename + '.bed'
    with open(bed, "rb") as f:
        mode = f.read(2)
        if mode != b'l\x1b':
            raise Exception('No valid binary PED file')
        mode = f.read(1) #\x01 = SNP major \x00 = individual major
        if mode != b'\x01':
            raise Exception('only SNP-major is implemented')
        startbit = SP.ceil(0.25*N)*start+3
        f.seek(int(startbit))
        for blockStart in SP.arange(0,nSNPs,blocksize, dtype=int):
            blockEnd = int(min(S,blockStart+blocksize))
            Sblock = min(nSNPs-blockStart,blocksize)
            nbyte = int(SP.ceil(0.25*N)*Sblock)
            bytes = SP.array(bytearray(f.read(nbyte))).reshape((SP.ceil(0.25*N),Sblock),order='F')
            
            SNPs[3::4,blockStart:blockEnd][bytes>=64]=SP.nan
            SNPs[3::4,blockStart:blockEnd][bytes>=128]=1
            SNPs[3::4,blockStart:blockEnd][bytes>=192]=2
            bytes=SP.mod(bytes,64)
            SNPs[2::4,blockStart:blockEnd][bytes>=16]=SP.nan
            SNPs[2::4,blockStart:blockEnd][bytes>=32]=1
            SNPs[2::4,blockStart:blockEnd][bytes>=48]=2
            bytes=SP.mod(bytes,16)
            SNPs[1::4,blockStart:blockEnd][bytes>=4]=SP.nan
            SNPs[1::4,blockStart:blockEnd][bytes>=8]=1
            SNPs[1::4,blockStart:blockEnd][bytes>=12]=2
            bytes=SP.mod(bytes,4)
            SNPs[0::4,blockStart:blockEnd][bytes>=1]=SP.nan
            SNPs[0::4,blockStart:blockEnd][bytes>=2]=1
            SNPs[0::4,blockStart:blockEnd][bytes>=3]=2
    
    if 0: #the binary format as described in the documentation (seems wrong)
        SNPs[3::4][bytes>=128]=SP.nan
        SNPs[3::4][bytes>=192]=1
        bytes=SP.mod(bytes,128)
        SNPs[3::4][bytes>=64]+=1
        bytes=SP.mod(bytes,64)
        SNPs[2::4][bytes>=32]=SP.nan
        SNPs[2::4][bytes>=48]=1
        bytes=SP.mod(bytes,32)
        SNPs[2::4][bytes>=16]+=1
        bytes=SP.mod(bytes,16)
        SNPs[1::4][bytes>=8]=SP.nan
        SNPs[1::4][bytes>=12]=1
        bytes=SP.mod(bytes,8)
        SNPs[1::4][bytes>=4]+=1
        bytes=SP.mod(bytes,4)
        SNPs[0::4][bytes>=2]=SP.nan
        SNPs[0::4][bytes>=3]=1
        bytes=SP.mod(bytes,2)
        SNPs[0::4][bytes>=1]+=1
    snps = SNPs[0:N,:]

    if useMAFencoding:
        imaf = SP.sum(snps==2,axis=0)>SP.sum(snps==0,axis=0)
        snps[:,imaf] = 2 - snps[:,imaf]
        
    if standardizeSNPs:
        snps = standardize(snps)
    ret = {
            'rs'     :rs[start:S_res],
            'pos'    :pos[start:S_res,:],
            'snps'   :snps,
            'iid'    :fam
            }
    return ret

