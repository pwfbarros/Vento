# -*- coding: utf-8 -*-
"""
Created on Sat Jul 15 09:06:05 2017
@author: pbarros
"""
import matplotlib.pyplot as plt

def graph(fig, ax, x, y1, c = 'b', l = '-', l_width = 1.0, n = None, title = None, lab_x = None, lab_y = None,
          xlim = None, ylim = None, grid = True, minor = True, semilogx = False, m = None, m_size = 1):      
    if title != None: fig.suptitle(title, fontsize=14)
    fig.subplots_adjust(hspace=1.)
    if semilogx == True: ax.semilogx(x, y1, label = n, color = c, linestyle = l, linewidth = l_width, marker = m, markersize = m_size)
    else: ax.plot(x, y1, label = n, color = c, linestyle = l, linewidth = l_width, marker = m, markersize = m_size)
    if xlim != None: ax.set_xlim(xlim)
    if ylim != None: ax.set_ylim(ylim)
    if lab_x != None: ax.set_xlabel(lab_x)
    if lab_y != None: ax.set_ylabel(lab_y)
    if n != None: plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    if grid == True: ax.grid(grid)
    if minor == True: ax.grid(minor, which='minor') 
