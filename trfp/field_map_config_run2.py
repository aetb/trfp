## paths and dictionaries

path = '/data1/newg2/DataProduction/Offline/ArtTFSDir/v9_52_00_dev2/FieldPlainRootOutput_'

runs = ['2a', '2b', '2c', '2d', '2e', '2f', '2g', '2h', '2i']
prod_runs = ['2c', '2d', '2e', '2f', '2g', '2h', '2i']
subruns = ['2a1', '2b1', '2b2', '2c1', '2c2', '2c3', '2c4', '2d1', '2d2', '2d3',
           '2d4', '2d5', '2d6', '2e1', '2e2', '2f1', '2f2', '2g1', '2h1', '2i1']

root_dict = {'2a':{'fp_df_1':range(6630, 6672+1),
                   'tr_df_1':[6627], 'tr_df_2':[6676],
                   'subrun_df':[24433, 24474]},
             '2b':{'fp_df_1':range(6780, 6837+1), 'fp_df_2':range(6845, 6877+1),
                   'tr_df_1':[6777], 'tr_df_2':[6843], 'tr_df_3':[6880],
                   'subrun_df':[24499, 24648]},
             '2c':{'fp_df_1':range(6883, 6934+1), 'fp_df_2':range(6946, 6985+1),
                   'fp_df_3':range(6992, 7029+1), 'fp_df_4':range(7036, 7065+1),
                   'tr_df_1':[6880], 'tr_df_2':[6937], 'tr_df_3':[6989],
                   'tr_df_4':[7032], 'tr_df_5':[7070],
                   'subrun_df':[24683, 25045]},
             '2d':{'fp_df_1':range(7082, 7104), 'fp_df_2':range(7124, 7149),
                   'fp_df_3':range(7155, 7185), 'fp_df_4':range(7191, 7212),
                   'fp_df_5':range(7218, 7248), 'fp_df_6':range(7256, 7290),
                   'tr_df_1':[7078], 'tr_df_2':[7107],
                   'tr_df_3':[7152], 'tr_df_4':[7188],
                   'tr_df_5':[7215], 'tr_df_6':[7253],
                   'tr_df_7':[7293],
                   'subrun_df':[25894, 26384]},
             '2e':{'fp_df_1':range(7396, 7427+1), 'fp_df_2':range(7435, 7465+1),
                   'tr_df_1':[7392], 'tr_df_2':[7432], 'tr_df_3':[7468],
                   'subrun_df':[26459, 26624]},
             '2f':{'fp_df_1':range(7480, 7511+1), 'fp_df_2':range(7521, 7546+1),
                   'tr_df_1':[7477], 'tr_df_2':[7514], 'tr_df_3':[7549],
                   'subrun_df':[26675, 26804]},
             '2g':{'fp_df_1':range(7611, 7635+1),
                   'tr_df_1':[7608], 'tr_df_2':[7638],
                   'subrun_df':[26996, 27043]},
             '2h':{'fp_df_1':range(7678, 7696+1),
                   'tr_df_1':[7675], 'tr_df_2':[7699],
                   'subrun_df':[27166, 27215]},
             '2i':{'fp_df_1':range(7845, 7873+1),
                   'tr_df_1':[7842], 'tr_df_2':[7876],
                   'subrun_df':[27415, 27439]}
             }

## Defining a dictionary that, given a `fp` key, returns an array of the two bookending `tr` keys
## for passing to the baselines dictionary.

pair_dict = {'2a':{'fp_df_1':['tr_df_1', 'tr_df_2']},             
             '2b':{'fp_df_1':['tr_df_1', 'tr_df_2'], 'fp_df_2':['tr_df_2', 'tr_df_3']},             
             '2c':{'fp_df_1':['tr_df_1', 'tr_df_2'], 'fp_df_2':['tr_df_2', 'tr_df_3'],
                   'fp_df_3':['tr_df_3', 'tr_df_4'], 'fp_df_4':['tr_df_4', 'tr_df_5']},             
             '2d':{'fp_df_1':['tr_df_1', 'tr_df_2'], 'fp_df_2':['tr_df_2', 'tr_df_3'],
                   'fp_df_3':['tr_df_3', 'tr_df_4'], 'fp_df_4':['tr_df_4', 'tr_df_5'],
                   'fp_df_5':['tr_df_5', 'tr_df_6'], 'fp_df_6':['tr_df_6', 'tr_df_7']},             
             '2e':{'fp_df_1':['tr_df_1', 'tr_df_2'], 'fp_df_2':['tr_df_2', 'tr_df_3']},             
             '2f':{'fp_df_1':['tr_df_1', 'tr_df_2'], 'fp_df_2':['tr_df_2', 'tr_df_3']},             
             '2g':{'fp_df_1':['tr_df_1', 'tr_df_2']},             
             '2h':{'fp_df_1':['tr_df_1', 'tr_df_2']},             
             '2i':{'fp_df_1':['tr_df_1', 'tr_df_2']}}

# interp_file_dict = {'2a':'hdf5/2021-04-26_run_2a.h5', '2b':'hdf5/2021-04-26_run_2b.h5',
#                     '2c':'hdf5/2021-04-26_run_2c.h5', '2d':'hdf5/2021-04-26_run_2d.h5',
#                     '2e':'hdf5/2021-04-26_run_2e.h5', '2f':'hdf5/2021-04-26_run_2f.h5',
#                     '2g':'hdf5/2021-04-26_run_2g.h5', '2h':'hdf5/2021-04-26_run_2h.h5',
#                     '2i':'hdf5/2021-04-26_run_2i.h5'}
