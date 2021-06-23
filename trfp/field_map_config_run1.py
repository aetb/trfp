## paths and dictionaries

prefix = 'data1/newg2/DataProduction/Offline/ArtTFSDir/v9_21_05_dev/FieldPlainRootOutput_'

runs = ['1a', '1b', '1c', '1d']
subruns = ['1a1', '1b1', '1b2', '1c1', '1c2', '1c3', '1d2', '1d3', '1d4', '1d5', '1d6']

root_dict = {'1a':{'fp_df_1':range(3959,3995),
                   'tr_df_1':[3956], 'tr_df_2':[3997],
                   'subrun_df':[15921, 15992]},
             '1b':{'fp_df_1':range(4061,4096), 'fp_df_2':range(4101,4122)+range(4123,4136),
                   'tr_df_1':[4058], 'tr_df_2':[4098], 'tr_df_3':[4138],
                   'subrun_df':[16110, 16256]},
             '1c':{'fp_df_1':range(4141,4179), 'fp_df_2':range(4193,4223),
                   'fp_df_3':range(4229,4263), 'fp_df_4':range(4283,4489),
                   'tr_df_1':[4138], 'tr_df_2':[4181], 'tr_df_3':[4189],
                   'tr_df_4':[4226], 'tr_df_5':[4265], 'tr_df_6':[4493],
                   'subrun_df':[16355, 16539]},
             '1d':{#'fp_df_1':range(5000,5050),
                   'fp_df_2':range(5057,5101), 'fp_df_3':range(5120,5155),
                   'fp_df_4':range(5172,5215), 'fp_df_5':range(5220,5257),
                   'fp_df_6':range(5262,5301),
                   #'tr_df_1':[4197],
                   'tr_df_2':[5054], 'tr_df_3':[5103],
                   'tr_df_4':[5117], 'tr_df_5':[5157], 'tr_df_6':[5169],
                   'tr_df_7':[5217], 'tr_df_8':[5259], #'tr_df_9':[5303],
                   'subrun_df':[16908, 17528]}}

lowkick_dict = {'fp_df_1':range(4542,4582),
                'tr_df_1':[4539], 'tr_df_2':[4584],
                'subrun_df':[16669, 16714]}

## Defining a dictionary that, given a `fp` key, returns an array of the two bookending `tr` keys
## for passing to the baselines dictionary.

pair_dict = {'1a':{'fp_df_1':['tr_df_1', 'tr_df_2']},
             '1b':{'fp_df_1':['tr_df_1', 'tr_df_2'], 'fp_df_2':['tr_df_2', 'tr_df_3']},
             '1c':{'fp_df_1':['tr_df_1', 'tr_df_2'], 'fp_df_2':['tr_df_3', 'tr_df_4a'],
                   'fp_df_3':['tr_df_4', 'tr_df_5'], 'fp_df_4':['tr_df_5', 'tr_df_6']},
             '1d':{'fp_df_2':['tr_df_2', 'tr_df_3'], 'fp_df_3':['tr_df_4', 'tr_df_5'],
                   'fp_df_4':['tr_df_6', 'tr_df_7'], 'fp_df_5':['tr_df_7', 'tr_df_8'],
                   'fp_df_6':['tr_df_8', 'tr_df_9']}}

interp_file_dict = {'1a':'hdf5/2020-09-30_run_1a.h5', '1b':'hdf5/2020-09-30_run_1b.h5',
                    '1c':'hdf5/2020-09-30_run_1c.h5', '1d':'hdf5/2020-09-30_run_1d.h5'}
