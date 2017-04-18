import scipy as SP
import pdb

if __name__=='__main__':

    Cg_mtSet1VCNull   = SP.loadtxt('./res_mtSet1VCnull_base_Cg.txt')
    Cn_mtSet1VCNull   = SP.loadtxt('./res_mtSet1VCnull_base_Cn.txt')

    Cg_mtSetNull   = SP.loadtxt('./res_mtSetNull_base_Cg.txt')
    Cn_mtSetNull   = SP.loadtxt('./res_mtSetNull_base_Cn.txt')
    Cg_mtSetNullec = SP.loadtxt('./res_mtSetNull_eigenCache_Cg.txt')
    Cn_mtSetNullec = SP.loadtxt('./res_mtSetNull_eigenCache_Cn.txt')

    Cr_mtSet   = SP.loadtxt('./res_mtSet_base_Cr.txt')
    Cg_mtSet   = SP.loadtxt('./res_mtSet_base_Cg.txt')
    Cn_mtSet   = SP.loadtxt('./res_mtSet_base_Cn.txt')
    Cr_mtSetec = SP.loadtxt('./res_mtSet_eigenCache_Cr.txt')
    Cg_mtSetec = SP.loadtxt('./res_mtSet_eigenCache_Cg.txt')
    Cn_mtSetec = SP.loadtxt('./res_mtSet_eigenCache_Cn.txt')

    # Cg_mtSet1VCNull = 0
    print(((Cg_mtSet1VCNull==0).all()))

    # Eigenvalue chaching works
    print((((Cg_mtSetNull-Cg_mtSetNullec)**2<1e-4).all()))
    print((((Cn_mtSetNull-Cn_mtSetNullec)**2<1e-4).all()))
    print((((Cr_mtSet-Cr_mtSetec)**2<1e-4).all()))
    print((((Cg_mtSet-Cg_mtSetec)**2<1e-4).all()))
    print((((Cn_mtSet-Cn_mtSetec)**2<1e-4).all()))

