import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Adds small UMAP1/UMAP2 labeled axes to bottom right corner of figure
def add_umap_axes(ax,umap_axis_length = 0.1, umap_axis_offset = 0.04, 
              axis_lw = 0.5, umap_text_size = 5):
    
    # Determine appropriate location for axes
    ymin,ymax=ax.get_ylim()
    xmin,xmax=ax.get_xlim()
    ylen=ymax-ymin
    xlen=xmax-xmin
    umap_axis_xmin = xmin+umap_axis_offset*xlen
    umap_axis_ymin = ymin+umap_axis_offset*ylen
    
    # Add umap1 axis
    umap1_axis = plt.Rectangle((umap_axis_xmin, umap_axis_ymin), umap_axis_length*xlen, 0, 
                           linewidth=axis_lw, fill=None) 
    ax.add_patch(umap1_axis)

    # Add umap1 axis label
    umap1_text = plt.Rectangle((umap_axis_xmin, ymin), umap_axis_length*xlen, umap_axis_offset*ylen, 
                               linewidth=0, fill=None) 
    ax.add_patch(umap1_text)
    ax.annotate('UMAP1', (umap_axis_xmin + umap1_text.get_width()/2.0, 
                         ymin + umap1_text.get_height()/2.0), 
                color='black', fontsize=umap_text_size, ha='center', va='center')

    # Add umap2 axis
    umap2_axis = plt.Rectangle((umap_axis_xmin, umap_axis_ymin), 0, umap_axis_length*ylen, 
                               linewidth=axis_lw, fill=None) 
    ax.add_patch(umap2_axis)

    # Add umap2 axis label
    umap2_text = plt.Rectangle((xmin, umap_axis_ymin), umap_axis_offset*xlen, umap_axis_length*ylen,
                               linewidth=0, fill=None) 
    ax.add_patch(umap2_text)
    ax.annotate('UMAP2', (xmin + umap2_text.get_width()/2.0, 
                         umap_axis_ymin + umap2_text.get_height()/2.0), 
                color='black', fontsize=umap_text_size, ha='center', va='center',rotation=90)
    
    # Avoid disrupting overall figure layout
    umap1_axis.set_in_layout(False)
    umap2_axis.set_in_layout(False)
    umap2_text.set_in_layout(False)
    umap2_text.set_in_layout(False)