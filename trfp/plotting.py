import matplotlib.pyplot as plt
import pandas as pd

single_column_small = (6.202, 3.833)  #inches
single_column_med = (6.202, 6.202)
single_column_large = (6.202, 7.666)

def plt_unix_time_to_CST(ax):
    ax.locator_params(axis='x', nbins=5)
    xticks = ax.get_xticks()
    ax.set_xticklabels([pd.to_datetime(tm, unit='s').tz_localize('UTC').tz_convert('US/Central').strftime('%Y-%m-%d\n %H:%M:%S %Z')
                          for tm in xticks], rotation=30, fontdict={'size':12, 'family':'serif'})

def plt_set_labels(ax, x_label, y_label, title):
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)

def fig_watermark(fig):
    fig.text(0.5, 0.5, 'BLINDED,\nPRELIMINARY',
         fontsize=24, color='black', rotation=45,
         ha='center', va='center', alpha=0.25)
    