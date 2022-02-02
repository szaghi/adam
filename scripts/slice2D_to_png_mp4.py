#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

args_parser = argparse.ArgumentParser()
args_parser.add_argument('-b', type=str, required=True, help="basename of group of files or single filename")
args_parser.add_argument('-v', type=int, default=3, help="number of contour variable")
args_parser.add_argument('-clabel', type=str, default="rho", help="contour variable label")
args_parser.add_argument('-levels', type=float, default=None, nargs=3, help="contour variable levels [vmin, vmax, n]")
args_parser.add_argument('-mp4', type=str, default=None, help="mp4 filename")
args = args_parser.parse_args()
slice_files_pattern = args.b
v = args.v

if args.mp4 is not None:
    os.system('ffmpeg -framerate 10 -pattern_type glob -i "' + slice_files_pattern + '*.png" ' + args.mp4)
else:
    for file_name in os.listdir("./"):
        if file_name.startswith(slice_files_pattern):
            print('plot file ' + file_name)
            base_name = file_name.split('.')[0]
            # load data
            with open(base_name + ".dat") as sf:
              head = [next(sf) for x in range(3)]
            # number of variables saved
            hnv = head[0].strip()
            nv = len(hnv.split('= ')[1].split(' '))
            # number of slice points
            hnijk = head[2].strip()
            ni = int(hnijk.split(',')[0].split('=')[1])
            nj = int(hnijk.split(',')[1].split('=')[1])
            nk = int(hnijk.split(',')[2].split('=')[1])
            slice_o = np.loadtxt(base_name + ".dat", delimiter=" ", skiprows=3)
            slice_r = slice_o.reshape(ni,nj,nk,nv)
            if ni == 1:
                x = slice_r[0,:,0,1]
                y = slice_r[0,0,:,2]
                r = slice_r[0,:,:,v]
                x_label = 'Y'
                y_label = 'Z'
            elif nj == 1:
                x = slice_r[0,0,:,0]
                y = slice_r[:,0,0,2]
                r = slice_r[:,0,:,v]
                x_label = 'X'
                y_label = 'Z'
            elif nk == 1:
                x = slice_r[0,:,0,0]
                y = slice_r[:,0,0,1]
                r = slice_r[:,:,0,v]
                x_label = 'X'
                y_label = 'Y'
            # create contour plot
            fig, ax = plt.subplots(1, 1)
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)
            ax.set_title('2D Slice ' + x_label + y_label)
            if args.levels is not None:
                levels = np.linspace(args.levels[0], args.levels[1], int(args.levels[2]) + 1)
                cf = plt.contourf(x, y, r, levels=levels, cmap='RdGy')
            else :
                cf = plt.contourf(x, y, r, levels=20, cmap='RdGy')
            cbar = plt.colorbar(cf)
            cbar.ax.set_ylabel(args.clabel, labelpad=10, rotation=270)
            plt.savefig(base_name + '.png')
            plt.close()
