def entry_point():
    from optparse import OptionParser
    from ..mtSet.core.iset_utils import calc_emp_pv_eff
    import pandas as pd
    import glob
    import os

    parser = OptionParser()
    parser.add_option("--resdir", dest='resdir', type=str, default='./')
    parser.add_option("--outfile", dest='outfile', type=str, default=None)
    # parser.add_option("--manhattan_plot", dest='manhattan',action="store_true",default=False)
    parser.add_option("--tol", dest='tol', type=float, default=4e-3)
    (options, args) = parser.parse_args()

    resdir = options.resdir
    out_file = options.outfile
    tol = options.tol

    print('.. load permutation results')
    file_name = os.path.join(resdir, '*.iSet.perm')
    files = glob.glob(file_name)
    df0 = pd.DataFrame()
    for _file in files:
        print(_file)
        df0 = df0.append(pd.read_csv(_file, index_col=0))

    print('.. load real results')
    file_name = os.path.join(resdir, '*.iSet.real')
    files = glob.glob(file_name)
    df = pd.DataFrame()
    for _file in files:
        print(_file)
        df = df.append(pd.read_csv(_file, index_col=0))

    # calculate P values for the three tests
    for test in ['mtSet', 'iSet', 'iSet-het']:
        df[test + ' pv'] = calc_emp_pv_eff(df[test + ' LLR'].values,
                                           df0[test + ' LLR0'].values)

    print(('.. saving %s' % out_file + '.res'))
    df.to_csv(out_file + '.res')
