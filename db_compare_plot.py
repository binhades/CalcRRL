#! /home/bliu/anaconda3/bin/python

# Filename: db_fit_spec.py
# Aim: to load the data base and then to fit spectra.

import argparse, sqlite3
import numpy as np
import matplotlib.pyplot as plt
import pyrrl as rrl
#===================================================================
# Change the default behavior of few things
import matplotlib as mpl
mpl.rcParams['mathtext.default'] = "regular"
mpl.rcParams['xtick.labelsize'] = "xx-large"
mpl.rcParams['ytick.labelsize'] = "xx-large"

def comp_plot_fwhm(fwhm_x, fwhm_y,fwhm_corr, fileout):

    xlabel = r'FARLS FWHM (km$\,$s$^{-1}$)'
    ylabel = r'WISE FWHM (km$\,$s$^{-1}$)'
    xlim = [0, 120]
    ylim = [0, 120]

    plt.figure(figsize = [7.8,10])
    plt.subplots_adjust()
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.0)
#---------------------------------------------------------------
    ax0 = plt.subplot(gs[0])
    #ax0.errorbar(xdt1_corr, ydt1_corr, xerr=xerr1, yerr=yerr1, lw=0.5, fmt='o', \
    ax0.errorbar(fwhm_x, fwhm_y,  lw=0.5, fmt='o', \
        ecolor='0.5', mfc='black', mec='black', mew=1,markersize=4, capsize=5)
    #x = np.logspace(-4, 4, 10000)
    x = np.arange(0,120,0.2)
    y = np.arange(0,120,0.2)
    ax0.plot(x, x, color='black', linestyle='--', lw=2)
#    ax0.plot(x, y, color='black', linestyle='--', lw=2)
#    ax0.plot(y, x, color='black', linestyle='--', lw=2)
    ax0.set_ylabel(ylabel, fontsize=18)
    ax0.set_xlim(xlim)
    ax0.set_ylim(ylim)
#---------------------------------------------------------------
    ax1 = plt.subplot(gs[1])
#    ax1.errorbar(xdt1_corr, ydt1a_corr, xerr=xerr1, yerr=yerr1a, lw=0.5, \
    ax1.errorbar(fwhm_x, fwhm_corr, lw=0.5, fmt='o',  markersize=4,ecolor='0.5', mfc='black', mec='black', mew=1,  capsize=5)
    ax1.set_yscale('linear')
    ax1.set_xlabel(xlabel, fontsize=18)
    ax1.yaxis.set_ticks(np.arange(-50, 50, 20))
    ax1.yaxis.set_label_coords(-0.1, 0.5)

    ylabel_2 = r'$\Delta$FWHM (km$\,$s$^{-1}$)'
    ax1.set_ylabel(ylabel_2, fontsize=18)
    ax1.hlines(0, 0, 150, colors='black', linewidth=2, linestyles='dashed')
    ax1.set_xlim(xlim)
    ax1.set_ylim([-50,50])
    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.savefig(fileout, format = 'eps', dpi = 300, bbox_inches='tight')

    return 0

def comp_plot_vlsr(vlsr_x, vlsr_y,vlsr_corr,flux_x=None, fileout='compare_vlsr'):

    xlabel = r'FARLS V$_{LSR}$ (km$\,$s$^{-1}$)'
    ylabel = r'WISE V$_{LSR}$ (km$\,$s$^{-1}$)'
    xlim = [-15, 125]
    ylim = [-15, 125]
    #ylim = [-50, 125]

    if flux_x is not None:
        m_s = flux_x * 5000
    else:
        m_s = 4

    plt.figure(figsize = [10,7.7])
    plt.subplots_adjust()
    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.0)
#---------------------------------------------------------------
    ax0 = plt.subplot(gs[0])
    #ax0.errorbar(vlsr_x, vlsr_y, lw=0.5, fmt='o', ecolor='0.5', mfc='black', mec='black', mew=1,markersize=4, capsize=5)
    ax0.scatter(vlsr_x, vlsr_y, marker='o', color='black',s=m_s)
    #x = np.logspace(-4, 4, 10000)
    x = np.arange(-100,200,0.2)
    y = np.arange(-135,165,0.2)
    ax0.plot(x, x, color='black', linestyle='--', lw=2)
#    ax0.plot(x, y, color='black', linestyle='--', lw=2)
#    ax0.plot(y, x, color='black', linestyle='--', lw=2)
    ax0.set_ylabel(ylabel, fontsize=18)
    ax0.set_xlim(xlim)
    ax0.set_ylim(ylim)
#---------------------------------------------------------------
    ax1 = plt.subplot(gs[1])
    #ax1.errorbar(vlsr_x, vlsr_corr, lw=0.5, fmt='o',  markersize=4,ecolor='0.5', mfc='black', mec='black', mew=1,  capsize=5)
    ax1.scatter(vlsr_x, vlsr_corr, marker='o', color='black',s=m_s)
    ax1.set_yscale('linear')
    ax1.set_xlabel(xlabel, fontsize=18)
    ax1.yaxis.set_ticks(np.arange(-100, 100, 20.0))
    ax1.yaxis.set_label_coords(-0.1, 0.5)

    ylabel_2 = r'$\Delta$V$_{LSR}$ (km$\,$s$^{-1}$)'
    ax1.set_ylabel(ylabel_2, fontsize=18)
#    ax1.hlines(1e2, 1e-1, 1e2, colors='black', linewidth=2, linestyles='dashed')
    ax1.hlines(0, -100, 150, colors='black', linewidth=2, linestyles='dashed')
    ax1.set_xlim(xlim)
    ax1.set_ylim([-55,55])
    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.tight_layout()
    plt.savefig(fileout+'.eps', format = 'eps', dpi = 300)#, bbox_inches='tight')
    plt.savefig(fileout+'.png', format = 'png', dpi = 300)#, bbox_inches='tight')
    plt.show()

    return 0


def main(args):

    db_conn = rrl.db.create_connection(args.file_db)
    con_wise = sqlite3.connect(args.wise_db)
    cur = con_wise.cursor()
    db_cur = db_conn.cursor()

    with db_conn:

        spfith_list = rrl.db.read_from_table(db_conn,'SpFitHa')
        vlsr_fast = []
        flux_fast = []
        vlsr_wise = []

        for i, spfith in enumerate(spfith_list):
            spfith = spfith_list[i]
            select_query0 = """ select Glon,Glat from Source where "GName" = ? """
            select_query1= """ select VLSR from wisehii where "GLong" = ? and "GLat" = ? """
            db_cur.execute(select_query0, (spfith['GName'],))
            raw = db_cur.fetchall()
            l = str(round(raw[0][0],3))
            b = str(round(raw[0][1],3))
            cur.execute(select_query1, (l,b))
            records = cur.fetchall()
            if len(records) > 0:
                vlsr_str = records[0][0]
                flux_fast.append(spfith['flux_peak1'])
                vlsr_fast.append(spfith['velo_lsr1'])
                vlsr_wise.append(float(vlsr_str.split(';')[0]))
                print(spfith['Gname'],spfith['velo_lsr1'],spfith['fwhm1'])

    con_wise.close()
    flux_fast = np.array(flux_fast)
    vlsr_fast = np.array(vlsr_fast)
    vlsr_wise = np.array(vlsr_wise)
    vlsr_corr = vlsr_wise - vlsr_fast

    comp_plot_vlsr(vlsr_fast,vlsr_wise,vlsr_corr,flux_x = flux_fast)
            
    return 0

#----------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_db',  type=str,help='The sqlite database file')
    parser.add_argument('wise_db',  type=str,help='The sqlite database file of WISE HII catalog')
    parser.add_argument('--plot',  action='store_true',help='if set, plot the spectra with fitting')
    parser.add_argument('--write',  action='store_true',help='only works when "--plot" is set, write the plots into files')
    args = parser.parse_args()
    main(args)
