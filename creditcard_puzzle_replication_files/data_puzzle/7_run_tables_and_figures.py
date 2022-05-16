# -*- coding: utf-8 -*-
"""
Created on Fri Jun 03 08:49:28 2016

@author: gmf123
"""

from __future__ import division, print_function, absolute_import
import os
clear = lambda: os.system('cls')
clear()

import numpy as np
import pandas as pd


import matplotlib.pyplot as plt
from matplotlib import cm
plt.rcParams.update({'font.size': 8,
                     'font.family': 'STIXGeneral',
                     'mathtext.fontset': 'stix'})

Color1 = (3.0/255.0,103.0/255.0,166.0/255.0)
Color2 = (242.0/255.0,62.0/255.0,46.0/255.0)
Color3 = (3.0/255.0,166.0/255.0,166.0/255.0)
Color4 = (242.0/255.0,131.0/255.0,68.0/255.0)
Color5 = (242.0/255.0,100.0/255.0,48.0/255.0)
ColorGrey = (65.0/255.0,68.0/255.0,81.0/255.0)

def fig_time_shares(y, name):
            
    fig = plt.figure(frameon=False, dpi=800)   
    fig.set_size_inches(5.5, 3.5)    
    ax = fig.add_subplot(1,1,1)
        
    ax.spines["top"].set_linewidth(0.5) 
    ax.spines["bottom"].set_linewidth(0.5) 
    ax.spines["right"].set_linewidth(0.5)  
    ax.spines["left"].set_linewidth(0.5)  
    ax.tick_params(which='both', color=ColorGrey)
    ax.grid(True, zorder=0, color=ColorGrey)
    
    yearfirst = 1989
    yearlast = 2013
    
    ax.plot(np.arange(yearfirst,yearlast+3,3), y[:,1],
            lw=2, ls='-', marker='.', markersize=5,
            color=Color1,
            label=r'puzzle')
    
    ax.plot(np.arange(yearfirst,yearlast+3,3), y[:,2],
            lw=2.5, ls='--', marker='.', markersize=5,
            color=Color2,
            label=r'borrower')

    ax.plot(np.arange(yearfirst,yearlast+3,3), y[:,3],
            lw=2.5, ls='-.', marker='.', markersize=5,
            color=Color3,
            label=r'saver')

    ax.plot(np.arange(yearfirst,yearlast+3,3), y[:,4],
            lw=2.5, ls=':', marker='.', markersize=5,
            color=Color4,
            label=r'corner')

    ax.set_xlabel(r'year', fontsize=10)
    ax.set_ylabel(r'percent', fontsize=10)
    ax.set_ylim([0, 60])
    ax.set_xlim([1991, 2014])
    ax.set_xticks(np.arange(yearfirst,yearlast+3,3))
    
    legend = ax.legend(loc='upper center', ncol=2, fontsize=8)
    frame = legend.get_frame()        
    frame.set_linewidth(0.4)   
    frame.set_edgecolor(ColorGrey)

    [line.set_zorder(3) for line in ax.lines]
    fig.tight_layout()

    fig.savefig(os.getcwd()  + '\\tables_and_figures\\' + name + '.pdf')
    plt.close('all')
    return

if __name__ == "__main__":
          
    names = ['mean', 'p5', 'p15', 'p25', 'p50', 'p75', 'p85', 'p95']
    percs = [5, 15, 25, 50, 75, 85, 95]
    
    vars_other   = ['income_disp', 'install', 'houses', 'homeeq', 'networth']
    labels_head_other = ['After-tax income (mean)', 'Installment loans (mean)', 
                         'Home value (mean)', 'Home equity (mean)', 'Total net worth (mean)']
    labels_sub_other = ['\\hspace{3mm} 25th percentile', '\\hspace{3mm} 50th percentile', '\\hspace{3mm} 75th percentile']
    names_other = ['mean', 'p25', 'p50', 'p75']
        
    for cutoff in ['500']:

        with open(os.getcwd()  + '\\tables_and_figures\\table_4_1_' + cutoff + '.txt', "w") as text_file:
            
            # a. table 4.1
            headrow = ''
            headrow += " & \\multicolumn{1}{c}{Puzzle}"
            headrow += " & \\multicolumn{1}{c}{Borrower}"
            headrow += " & \\multicolumn{1}{c}{Saver}"
            headrow += " & \\multicolumn{1}{c}{Corner}"
            headrow += " & \\multicolumn{1}{c}{All}"           
            headrow += " \\" + "\\" + "\n"         
            text_file.write(headrow) 
            text_file.write("\\midrule \n") 
            text_file.write("\\addlinespace \n")
            text_file.write(" & \multicolumn{5}{c}{\\textit{percent}}")
            text_file.write(" \\" + "\\ \n")
            text_file.write("\\addlinespace \n\n")        
                
            for i in range(0,1+2*len(names)+len(names)):   
                
                string = ''
                
                if i == 0:  
                    
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '.xls')
                    df = xl.parse("shares")
                    values = df.values
                    
                    string = ' Share '
                    
                    for j in range(0,values.size):
                        string += ' & {:3.1f} '.format(values[0,j])            
                    string += ' & {:3.1f} '.format(100)
                    
                elif i >= 1 and i < 1+len(names):
                        
                    if i == 1:
                        string = ' Credit card debt, $D_t$ (mean)'      
                    else:
                        string = '\\hspace{3mm}' + str(percs[i-2]) + 'th percentile '
                    
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '.xls')
                    df = xl.parse("ccdebt")
                    values = df[names[i-1]].values
                    for j in range(0,values.size):
                        string += ' & {:3.2f} '.format(values[j]/100)
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '.xls')
                    df = xl.parse("ccdebt_tot")
                    values = df[names[i-1]].values       
                    string += ' & {:3.2f} '.format(values[0]/100)
                    
                elif 1+len(names) >= 1 and i < 1+2*len(names):
                        
                    if i == 1+len(names):
                        string = ' Liquid assets, $A_t$ (mean)'  
                    else:
                        string = '\\hspace{3mm}' + str(percs[i-2-len(names)]) + 'th percentile '
                    
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '.xls')
                    df = xl.parse("liqass")
                    values = df[names[i-1-len(names)]].values
                    for j in range(0,values.size):
                        string += ' & {:3.2f} '.format(values[j]/100)
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '.xls')
                    df = xl.parse("liqass_tot")
                    values = df[names[i-1-len(names)]].values       
                    string += ' & {:3.2f} '.format(values[0]/100)
                    
                else:
                        
                    if i == 1+2*len(names):
                        string = 'Liquid net worth, $N_t$ (mean)'  
                    else:
                        string = '\\hspace{3mm}' + str(percs[i-2-2*len(names)]) + 'th percentile '
                    
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '.xls')
                    df = xl.parse("liqnetworth")
                    values = df[names[i-1-2*len(names)]].values
                    for j in range(0,values.size):
                        string += ' & {:3.2f} '.format(values[j]/100)
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '.xls')
                    df = xl.parse("liqnetworth_tot")
                    values = df[names[i-1-2*len(names)]].values       
                    string += ' & {:3.2f} '.format(values[0]/100)
                    
                text_file.write(string + " \\\\ \n")                
                
                if i == 0: 
                    text_file.write("\\addlinespace \n")  
                    text_file.write("\\midrule \n")   
                    text_file.write("\\addlinespace \n")
                    text_file.write(" &  \multicolumn{5}{c}{\\textit{relative to mean quarterly income}}")
                    text_file.write(" \\" + "\\ \n")
                    text_file.write("\\addlinespace \n\n")
                if i in [1, 1+len(names)-1, 1+len(names),1+2*len(names)-1, 1+2*len(names)]:
                    text_file.write("\\addlinespace \n\n")
            
        # b. puzzle shares over time
        xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '.xls')
        df = xl.parse("shares_year")
        values = df.values
        fig_time_shares(values, 'time_data_shares_' + cutoff + '')

        # c. table 4.2          
        with open(os.getcwd()  + '\\tables_and_figures\\table_4_2_' + cutoff + '.txt', "w") as text_file:
                              
            headrow = ''
            headrow += " & \\multicolumn{1}{c}{Puzzle}"
            headrow += " & \\multicolumn{1}{c}{Borrower}"
            headrow += " & \\multicolumn{1}{c}{Saver}"
            headrow += " & \\multicolumn{1}{c}{Corner}"
            headrow += " & \\multicolumn{1}{c}{All}"           
            headrow += " \\" + "\\" + "\n"         
            text_file.write(headrow)
            text_file.write("\\midrule \n") 
            text_file.write("\\addlinespace \n")
            text_file.write(" & \multicolumn{5}{c}{\\textit{relative to mean quarterly income}}")
            text_file.write(" \\" + "\\ \n")
            text_file.write("\\addlinespace \n\n")        
                
            for i in range(0,len(vars_other)):
                for j in range(0,len(names_other)):
                    
                    if j == 0:
                        string = labels_head_other[i]
                    else:
                        string = labels_sub_other[j-1]
                    
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '.xls')
                    df = xl.parse(vars_other[i])
                    values = df[names_other[j]].values
                    for k in range(0,values.size):
                        string += ' & {:3.2f} '.format(values[k]/100)
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '.xls')
                    df = xl.parse(vars_other[i] + '_tot')
                    values = df[names_other[j]].values       
                    string += ' & {:3.2f} '.format(values[0]/100)                                        
                    text_file.write(string + " \\\\ \n")  
                    if j == 0:
                        text_file.write("\\addlinespace \n\n")
                
                if i != len(vars_other)-1:
                    text_file.write("\\addlinespace \n\n")
