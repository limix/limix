import scipy as sp
import h5py
#re for regular expressions
import re
import pdb

class SocialData():

    def __init__(self, task = None, kinship_type = "", subset = None):
        assert task is not None, 'Specify task!'
        self.task=task
        self.kinship_type=kinship_type
        self.subset = subset
        self.load()


    def load(self):

        if self.task == 'HSmice_data_REVISIONS':
            in_file = '/nfs/leia/research/stegle/abaud/HSmice/data/HSmice_wPleth_REVISIONS.hdf5'
            f = h5py.File(in_file,'r')

            self.measures = f['data_bcNcovariates']['cols_measures']['measures'][:]
            self.all_pheno = f['data_bcNcovariates']['array'][:].T
            self.pheno_ID = f['data_bcNcovariates']['rows_subjects']['outbred'][:]
            self.all_pheno = f['data_bcNcovariates']['array'][:].T
            self.pheno_ID = f['data_bcNcovariates']['rows_subjects']['outbred'][:]
            self.all_covs2use = f['data_bcNcovariates']['cols_measures']['covariatesUsed']
            self.all_covs = f['data_bcNcovariates']['array'][:].T
            self.covs_ID = f['data_bcNcovariates']['rows_subjects']['outbred'][:]

            if self.kinship_type == 'old_IBS':
                self.kinship_full = f['kinships']['old_IBS']['array'][:]
                self.kinship_full_ID = f['kinships']['old_IBS']['cols_subjects'][:]
            elif self.kinship_type == 'Andres_kinship':
                self.kinship_full = f['kinships']['Andres_kinship']['array'][:]
                self.kinship_full_ID = f['kinships']['Andres_kinship']['cols_subjects'][:]
            self.cage_full = f['data_bcNcovariates']['rows_subjects']['cage'][:]
            self.cage_full_ID = f['data_bcNcovariates']['rows_subjects']['outbred'][:]
            if self.subset is not None:
                self.subset_IDs = f['subsets'][self.subset][:]
            else:
                self.subset_IDs = None

        elif self.task == 'HSmice_hip_REVISIONS':
            in_file = '/nfs/leia/research/stegle/abaud/HSmice/data/HSmice_wPleth_REVISIONS.hdf5'
            f = h5py.File(in_file,'r')

            self.measures = f['hip_expr']['cols_probes']['probes'][:]
            self.all_pheno = f['hip_expr']['array'][:].T
            self.pheno_ID = f['hip_expr']['rows_subjects'][:]
            self.all_covs = f['data_bcNcovariates']['array'][:].T
            self.covs_ID = f['data_bcNcovariates']['rows_subjects']['outbred'][:]
            if self.kinship_type == 'old_IBS':
                self.kinship_full = f['kinships']['old_IBS']['array'][:]
                self.kinship_full_ID = f['kinships']['old_IBS']['cols_subjects'][:]
            elif self.kinship_type == 'Andres_kinship':
                self.kinship_full = f['kinships']['Andres_kinship']['array'][:]
                self.kinship_full_ID = f['kinships']['Andres_kinship']['cols_subjects'][:]
            self.cage_full = f['data_bcNcovariates']['rows_subjects']['cage'][:]
            self.cage_full_ID = f['data_bcNcovariates']['rows_subjects']['outbred'][:]
            if self.subset is not None:
                self.subset_IDs = f['subsets'][self.subset][:]
            else:
                self.subset_IDs = None


        elif 'HSmice_simulations' in self.task:
            if self.task == 'HSmice_simulations_realized':
                in_file = '/nfs/research/stegle/projects/IGE_outbreds/HSmice/simulations/revisions_HSmice_simulations_realized.hdf5'
            elif self.task == 'HSmice_simulations_all':
                in_file = '/nfs/research/stegle/projects/IGE_outbreds/HSmice/simulations/negcorr_revisions_HSmice_simulations_around_best.hdf5'
            f = h5py.File(in_file,'r')

            self.measures = f['simulations']['cols_measures']['measures'][:]
            self.all_pheno = f['simulations']['array'][:].T
            self.pheno_ID = f['simulations']['rows_subjects']['outbred'][:]
            self.covs_ID = None

            in_file = '/nfs/research/stegle/projects/IGE_outbreds/HSmice/HSmice_wPleth_new.hdf5'
            f = h5py.File(in_file,'r')

            if self.kinship_type == 'Andres_kinship':
                self.kinship_full = f['kinships']['Andres_kinship']['array'][:]
                self.kinship_full_ID = f['kinships']['Andres_kinship']['cols_subjects'][:]
            self.cage_full = f['residualsNcovariates']['rows_subjects']['cage'][:]
            self.cage_full_ID = f['residualsNcovariates']['rows_subjects']['outbred'][:]
            if self.subset is not None:
                self.subset_IDs = f['subsets'][self.subset][:]
            else:
                self.subset_IDs = None

        else:
            print("Nothing done: task unknown!")


    def get_data(self,col):

        if self.task == 'HSmice_data_REVISIONS':
            self.trait=self.measures[col]
            self.pheno = self.all_pheno[:,col]
            covs2use = self.all_covs2use[col].split(',')
            Ic = sp.zeros(self.measures.shape[0],dtype=bool)
            for cov in covs2use:
                Ic = sp.logical_or(Ic,self.measures==cov)
            self.covs = self.all_covs[:,Ic]

        elif self.task == 'HSmice_hip_REVISIONS':
            self.trait=self.measures[col]
            self.pheno = self.all_pheno[:,col]
            covs2use = ['GENDER','group_size']
            Ic = sp.zeros(self.measures.shape[0],dtype=bool)
            for cov in covs2use:
                Ic = sp.logical_or(Ic,self.measures==cov)
            self.covs = self.all_covs[:,Ic]

        elif 'HSmice_simulations' in self.task:
            self.trait = self.measures[col]
            self.pheno = self.all_pheno[:,col]
            self.covs = None

        else:
            print("Nothing done: task unknown!")


        return {'trait' : self.trait,
                'pheno' : self.pheno,
                'pheno_ID' : self.pheno_ID,
                'covs' : self.covs,
                'covs_ID' : self.covs_ID,
                'kinship_type' : self.kinship_type,
                'kinship_full' : self.kinship_full,
                'kinship_full_ID' : self.kinship_full_ID,
                'cage_full' : self.cage_full,
                'cage_full_ID' : self.cage_full_ID,
                'subset_IDs' : self.subset_IDs}
