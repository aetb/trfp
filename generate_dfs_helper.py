import csv,datetime,time
import psycopg2
import sys

import pandas as pd
import numpy as np

def get_subrun_df(start_run, end_run):
    
    dsn  = "dbname=gm2_online_prod user=gm2_reader host=g2db-priv port=5433"
    conn = psycopg2.connect(dsn)
    curr = conn.cursor()
    
    # get times by subrun
    subrun_time_columns = ['run', 'subrun', 'start_time', 'end_time', 'start_gps', 'end_gps']
    sql = "select "+", ".join(subrun_time_columns)+" from gm2dq.subrun_time where run >= %i and run <= %i order by run, subrun" % (start_run, end_run)
    curr.execute(sql)
    conn.commit()
    subrun_time = curr.fetchall()
    subrun_time_df = pd.DataFrame(subrun_time, columns=subrun_time_columns)
    
    # get ctags by subruns
    ctagswithdqc_columns = ['run', 'subrun', 'ctags', 't0val', 'fills']
    sql = "select "+", ".join(ctagswithdqc_columns)+" from gm2dq.ctagswithdqc where run >= %i and run <= %i order by run, subrun" % (start_run, end_run)
    curr.execute(sql)
    conn.commit()
    ctagswithdqc = curr.fetchall()
    ctagswithdqc_df = pd.DataFrame(ctagswithdqc, columns=ctagswithdqc_columns)
    
    # get subrun status database
    subrun_status_columns = ['run', 'subrun', 'quad_condition', 'kicker_condition', 'quad_ok',
                             'ctags_ok', 'losses_ok', 'fillcuts_ok', 'field_ok', 'trolley_period', 'field_period',
                             'ctags_loose_ok', 'quad_loose_ok', 'ctags_repeat_ok', 'losses_repeat_ok', 'fillcuts_repeat_ok']
    sql = "select "+", ".join(subrun_status_columns)+" from gm2dq.subrun_status where run >= %i and run <= %i order by run, subrun" % (start_run, end_run)
    curr.execute(sql)
    conn.commit()
    subrun_status = curr.fetchall()
    subrun_status_df = pd.DataFrame(subrun_status, columns=subrun_status_columns)
    
    # merge times, ctags, status into one subrun dataframe
    subrun_df = pd.merge(subrun_time_df, ctagswithdqc_df)
    subrun_df = pd.merge(subrun_df, subrun_status_df)
    subrun_df['ok'] = subrun_df['quad_ok'] & subrun_df['ctags_ok'] & subrun_df['losses_ok'] & subrun_df['fillcuts_ok'] & subrun_df['field_ok']
    subrun_df['start_time'] = subrun_df['start_time'].astype(np.int64)/1e9 + 5*60*60
    subrun_df['end_time'] = subrun_df['end_time'].astype(np.int64)/1e9 + 5*60*60
    subrun_df['start_gps'] = subrun_df['start_gps'].astype(np.int64)/1.0e9 + 5*60*60
    subrun_df['end_gps'] = subrun_df['end_gps'].astype(np.int64)/1.0e9 + 5*60*60
    
    print subrun_df.shape
    
    return subrun_df