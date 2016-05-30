"""
Filename:     plot_qq.py
Author:       Mitchell Black, mtblack@student.unimelb.edu.au
Description:  Create quantile-quantile plot from n 1-D NetCDF timeseries
"""

import numpy as np
import cdms2
import MV2
import os
import sys
import copy
import argparse
from dateutil.parser import *
import datetime 
import matplotlib.pyplot as plt
from scipy import stats

EU = os.path.expanduser
sys.path.append(EU("/Users/mtblack/Documents/University.dir/PhD.dir/Code.dir/modules.dir/"))
sys.path.append(EU("/home/mblack/PhD_Code/modules.dir/"))

import modules_netcdf as nio
import modules_general as gio
import modules_plot as pio

def main(inargs):
    
    # Check input arguments
    plot_dict = vars(inargs)    
    del inargs
    
    mand_argnames = ['obs_infile','obs_varname','obs_units']
    pio.check_inargs_mand_rpt(plot_dict, args2check = mand_argnames, check_against='obs_infile')

    mand_argnames = ['ens_infile','ens_varname','ens_units','ens_color','ens_label','ens_pctlerange','ens_linestyle']
    pio.check_inargs_mand_rpt(plot_dict, args2check = mand_argnames, check_against='ens_infile')
 
    # Initialise figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Get attributes for plotting
    conf_range = plot_dict['conf_range']
    alpha = plot_dict['alpha'][0]
    ensure_taxis_same = plot_dict['ensure_taxis_same'][0]
    
    # Calculate quantiles for observations
   
    obs_infiles = pio.read_infile_list(plot_dict['obs_infile'])  
    obs_varname = plot_dict['obs_varname'][0]   
    obs_units   = pio.fix_label(plot_dict['obs_units'][0])

    pio.print_plotting_status(0,obs_varname,obs_units)

    obs_tseries = nio.concat_nctseries2array(obs_infiles,obs_varname,obs_units,ensure_taxis_same=ensure_taxis_same)
   
    obs_quantiles = calc_quantiles(obs_tseries)

    # Calculate quantiles for model ensemble(s) 

    for i in range(len(plot_dict['ens_infile'])):
        ens_infiles    = pio.read_infile_list(plot_dict['ens_infile'][i])
        ens_varname    = plot_dict['ens_varname'][i][0]
        ens_units      = pio.fix_label(plot_dict['ens_units'][i][0])
        ens_color      = pio.fix_label(plot_dict['ens_color'][i][0])
        ens_label      = plot_dict['ens_label'][i][0]
        ens_pctlerange = plot_dict['ens_pctlerange'][i][0]
        ens_linestyle  = pio.line_alias2style(plot_dict['ens_linestyle'][i][0])
        mark_ptiles    = plot_dict['mark_ptiles']
        tickmark_size = plot_dict['tickmark_size'][0]

        pio.print_plotting_status(i+1,ens_varname, ens_units)      
        
        ens_tseries = nio.concat_nctseries2array(ens_infiles,ens_varname,ens_units,ensure_taxis_same=ensure_taxis_same)

        ens_quantiles = calc_quantiles(ens_tseries)

        ax.plot(obs_quantiles, ens_quantiles, ens_linestyle, color=ens_color, label=ens_label,lw=1.5, zorder=3)
       
        if ens_pctlerange:
            n_qtles = len(obs_quantiles)
            n_ensmembers = ens_tseries.shape[0]
            
            store = np.zeros((n_ensmembers,n_qtles),'f')
            store[:] = np.nan
            
            for i in range(0,n_ensmembers):
                store[i,:] = calc_quantiles(ens_tseries[i,:])
            
            conf_range_lower = np.zeros(n_qtles, 'f')
            conf_range_upper = np.zeros(n_qtles, 'f')
            conf_range_lower[:] = np.nan
            conf_range_upper[:] = np.nan     
          
            for qtle in range(0,n_qtles):
                conf_range_lower[qtle] = stats.scoreatpercentile(store[:,qtle], conf_range[0]*100)
                conf_range_upper[qtle] = stats.scoreatpercentile(store[:,qtle], conf_range[1]*100)
            
            ax.fill_between(obs_quantiles, conf_range_lower, conf_range_upper, facecolor = ens_color, edgecolor = ens_color, alpha = alpha)
  
        if mark_ptiles:
            ens_range = ens_quantiles[-1] - ens_quantiles[0]
            for p in mark_ptiles:
                obs_p = stats.scoreatpercentile(obs_tseries.flatten(), p)
                ens_p = stats.scoreatpercentile(ens_tseries.flatten(), p)
                ax.plot(obs_p, ens_p,'ko', ms=5, zorder=1)
                ax.text(obs_p+0.03*ens_range, ens_p-0.04*ens_range, str(p), ha='center', va='bottom',zorder=1,fontsize=tickmark_size)
 


    # Add figure display options
   
    lgnd_size = plot_dict['lgnd_size'][0]
    lgnd_loc  = plot_dict['lgnd_loc'][0]
    xylim = plot_dict['xylim']
    xlabel = plot_dict['xlabel'][0]
    ylabel = plot_dict['ylabel'][0]
    title = plot_dict['title'][0]
    xy_mnrticks = plot_dict['xy_mnrticks']
    ofile = plot_dict['ofile'][0]
    
    pio.plot_legend_no_repeats(ax, font_size = lgnd_size, loc = lgnd_loc)
    
    plt.tick_params(labelsize=tickmark_size)

    if xylim:
        assert xylim[0] < xylim[1]
        ax.set_xlim([xylim[0], xylim[1]])
        ax.set_ylim([xylim[0], xylim[1]])
        ax.plot(xylim,xylim,'k-')
    else:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        xylim = ( min(xlim[0],ylim[0]), max(xlim[1],ylim[1]) )
        ax.set_xlim(xylim[0],xylim[1])
        ax.set_xlim(xylim[0],xylim[1])
        ax.plot(xylim,xylim,'k-')
        
    if xlabel:
        ax.set_xlabel(pio.fix_label(xlabel))

    if ylabel:
        ax.set_ylabel(pio.fix_label(ylabel))

    if title:
        ax.set_title(pio.fix_label(title))
    
    if xy_mnrticks:
        pio.set_yaxis_minorticks(ax,xy_mnrticks)
        pio.set_xaxis_minorticks(ax,xy_mnrticks)

    # Save figure to file
    fig.savefig(ofile)

def calc_quantiles(vals):
     quantiles = []
     for i in range(1, 100):
         quantiles.append(stats.scoreatpercentile(vals.flatten(), i))
     return quantiles
 
if __name__ == '__main__':

    extra_info = """
Usage:
    python plot_qq.py -h
    legend options:
       location:   1 upper right, 2 upper left, 3 lower left, 4 lower right, 5 right, 6 center left,
                   7 center right, 8 lower center, 9 upper center, 10 center, None no legend
    label options:
       keys:       ~openp~ ( , ~closep~ ), ~deg~ '$^\circ$', ~hash~ # 
Author:
  Mitchell Black, mtblack@student.unimelb.edu.au
"""
    description='Create quantile-quantile plot'
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create group for required arguments
    parser_required = parser.add_argument_group('required arguments')
    
    # Mandantory input arguments
    parser_required.add_argument("--obs_infile", type=str, nargs='*', help="Input file(s) for observations. Can use wildcard '*.nc' or infile '.list'")
    parser_required.add_argument("--obs_varname", type=str, nargs=1, help="Input variable name")
    parser_required.add_argument("--obs_units", type=str, nargs=1, help="Required units for observations")

    # Mandantory input arguments [repeated as desired]
    parser_required.add_argument("--ens_infile", type=str, nargs='*', action='append', help="Input file(s) for model ensemble. Can use wildcard '*.nc' or infile '.list'")
    parser_required.add_argument("--ens_varname", type=str, nargs=1, action='append', help="Input variable name")
    parser_required.add_argument("--ens_units", type=str, nargs=1, action='append', help="Required units for model data")
    parser_required.add_argument("--ens_pctlerange", type=str, nargs=1, action='append', choices=['True','False'], help="Plot ensemble range?")
    parser_required.add_argument("--ens_color", type=str, nargs=1, action='append', help="Specify color for plotting (inc. hexcode, e.g. ~hash~31a354")
    parser_required.add_argument("--ens_label", type=str, nargs=1, action='append', help="Label for figure legend")
    parser_required.add_argument("--ens_linestyle", type=str, nargs=1, action='append', choices=['solid', 'dashdash', 'dotdash', 'dot'], help="Line style")
   
    # Optional input arguments [not repeated]
    parser.add_argument("--marker_size", type=float, nargs=1, default=[2.0], help="Marker size [default: 2.0]")
    parser.add_argument("--conf_range", type=float, nargs=2, default=(0.05, 0.95), metavar=('MIN', 'MAX'), help="Range for confidence bounds [default: 0.05 0.95]")
    parser.add_argument("--mark_ptiles", type=int, nargs='*', default=None, help="List percentiles to mark on image [default: None]")
    parser.add_argument("--alpha", type=float, nargs=1, default=[0.3], help="Fill shape transparency [default: 0.3]")
    parser.add_argument("--title", type=str, nargs=1, default=[None], help="Figure title. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--ylabel", type=str, nargs=1, default=[None], help="Y-axis label. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--xlabel", type=str, nargs=1, default=[None], help="X-axis label. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--xylim", type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help="Set limit for x- and y-axis bounds [default: auto]")
    parser.add_argument("--xy_mnrticks", type=int, nargs=1, default=None, help="No. minor ticks on x- and y-axis [default: None]")
    parser.add_argument("--lgnd_loc", type=int, nargs=1, default=[2], choices=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], help="Location of the primary figure legend [default: 2]")
    parser.add_argument("--lgnd_size", type=str, nargs=1, default=['medium'], choices=['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], help="Size of the legend text [default: medium]")
    parser.add_argument("--ensure_taxis_same", type=str, nargs=1, default=['True'], choices=['True','False'], help="Ensure time axis are same for files read in? [default: True]")
    parser.add_argument("--tickmark_size", type=float, nargs=1, default=[18.0], help="Tick mark font size") 
    parser.add_argument("--ofile", type=str, nargs=1, default=['plot.png'], help="Name of ouput file [default: None]")

    args = parser.parse_args()
    
    main(args)
