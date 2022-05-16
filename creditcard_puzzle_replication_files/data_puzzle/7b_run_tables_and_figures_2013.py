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


if __name__ == "__main__":
          
    names = ['mean', 'p5', 'p15', 'p25', 'p50', 'p75', 'p85', 'p95']
    percs = [5, 15, 25, 50, 75, 85, 95]
    
    vars_other   = ['income_disp', 'install', 'houses', 'homeeq', 'networth']
    labels_head_other = ['After-tax income (mean)', 'Installment loans (mean)', 
                         'Home value (mean)', 'Home equity (mean)', 'Total net worth (mean)']
    labels_sub_other = ['\\hspace{3mm} 25th percentile', '\\hspace{3mm} 50th percentile', '\\hspace{3mm} 75th percentile']
    names_other = ['mean', 'p25', 'p50', 'p75']
        
    for cutoff in ['500']:

        with open(os.getcwd()  + '\\tables_and_figures\\table_4_1_' + cutoff + '_2013.txt', "w") as text_file:
            
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
                    
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '_2013.xls')
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
                    
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '_2013.xls')
                    df = xl.parse("ccdebt")
                    values = df[names[i-1]].values
                    for j in range(0,values.size):
                        string += ' & {:3.2f} '.format(values[j]/100)
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '_2013.xls')
                    df = xl.parse("ccdebt_tot")
                    values = df[names[i-1]].values       
                    string += ' & {:3.2f} '.format(values[0]/100)
                    
                elif 1+len(names) >= 1 and i < 1+2*len(names):
                        
                    if i == 1+len(names):
                        string = ' Liquid assets, $A_t$ (mean)'  
                    else:
                        string = '\\hspace{3mm}' + str(percs[i-2-len(names)]) + 'th percentile '
                    
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '_2013.xls')
                    df = xl.parse("liqass")
                    values = df[names[i-1-len(names)]].values
                    for j in range(0,values.size):
                        string += ' & {:3.2f} '.format(values[j]/100)
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '_2013.xls')
                    df = xl.parse("liqass_tot")
                    values = df[names[i-1-len(names)]].values       
                    string += ' & {:3.2f} '.format(values[0]/100)
                    
                else:
                        
                    if i == 1+2*len(names):
                        string = 'Liquid net worth, $N_t$ (mean)'  
                    else:
                        string = '\\hspace{3mm}' + str(percs[i-2-2*len(names)]) + 'th percentile '
                    
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '_2013.xls')
                    df = xl.parse("liqnetworth")
                    values = df[names[i-1-2*len(names)]].values
                    for j in range(0,values.size):
                        string += ' & {:3.2f} '.format(values[j]/100)
                    xl = pd.ExcelFile(os.getcwd() + '\\output_' + cutoff + '_2013.xls')
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
