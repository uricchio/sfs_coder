import sfsplot
import sys
import os

data = sfsplot.PlotData()

data.annot_file = os.path.join(os.getcwd(),'input_files','sfscode_annotation_chr2_134545415-138594750_CDSbuff_neut.txt')
data.pi_file = os.path.join(os.getcwd(),'pi.txt')

data.get_annotations()
data.get_pi(pi_0=0.001)

data.pi_plot_regions()

