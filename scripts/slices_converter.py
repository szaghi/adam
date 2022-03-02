#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

def _subparser_basename():
    """
    Construct a cli subparser with the basename group of arguments.

    Returns
    -------
    parser : argparse.ArgumentParser()
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser_group = parser.add_argument_group('basename')
    parser_group.add_argument('-b', type=str, required=True, help="basename of group of files or single filename")
    return parser


def _subparser_mat2dat():
    """
    Construct a cli subparser with the mat2dat group of arguments.

    Returns
    -------
    parser : argparse.ArgumentParser()
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser_group = parser.add_argument_group('mat2dat')
    parser_group.add_argument('-vnames', type=str, default="'VARIABLES = "+'"x" "y" "z" "rho" "rhou" "rhov" "rhow" "rhoe" "f"'+"'",
                              help="variables names")
    return parser


def _subparser_dat2png():
    """
    Construct a cli subparser with the dat2png group of arguments.

    Returns
    -------
    parser : argparse.ArgumentParser()
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser_group = parser.add_argument_group('dat2png')
    parser_group.add_argument('-v', type=int, default=3, help="number of contour variable")
    parser_group.add_argument('-clabel', type=str, default="rho", help="contour variable label")
    parser_group.add_argument('-levels', type=float, default=None, nargs=3, help="contour variable levels [vmin, vmax, n]")
    return parser


def _subparser_png2mp4():
    """
    Construct a cli subparser with the dat2png group of arguments.

    Returns
    -------
    parser : argparse.ArgumentParser()
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser_group = parser.add_argument_group('png2mp4')
    parser_group.add_argument('-mp4', type=str, default=None, help="mp4 filename")
    return parser


def _parser_mat2dat(clisubparsers):
    """
    Construct the mat2dat cli parser.

    Parameters
    ----------
    clisubparsers : argparse subparser object
    """
    basename = _subparser_basename()
    mat2dat = _subparser_mat2dat()
    mat2dat_parser = clisubparsers.add_parser('mat2dat',
                                              help='convert slices from .mat to .dat files',
                                              parents=[basename, mat2dat])
    mat2dat_parser.set_defaults(which='mat2dat')
    return


def _parser_dat2png(clisubparsers):
    """
    Construct the dat2png cli parser.

    Parameters
    ----------
    clisubparsers : argparse subparser object
    """
    basename = _subparser_basename()
    dat2png = _subparser_dat2png()
    dat2png_parser = clisubparsers.add_parser('dat2png',
                                              help='convert slices from .dat to .png files',
                                              parents=[basename, dat2png])
    dat2png_parser.set_defaults(which='dat2png')
    return


def _parser_png2mp4(clisubparsers):
    """
    Construct the png2mp4 cli parser.

    Parameters
    ----------
    clisubparsers : argparse subparser object
    """
    basename = _subparser_basename()
    png2mp4 = _subparser_png2mp4()
    png2mp4_parser = clisubparsers.add_parser('png2mp4',
                                              help='convert slices from .png to .mp4 video',
                                              parents=[basename, png2mp4])
    png2mp4_parser.set_defaults(which='png2mp4')
    return


def _parser_mat2mp4(clisubparsers):
    """
    Construct the mat2mp4 cli parser.

    Parameters
    ----------
    clisubparsers : argparse subparser object
    """
    basename = _subparser_basename()
    mat2dat = _subparser_mat2dat()
    dat2png = _subparser_dat2png()
    png2mp4 = _subparser_png2mp4()
    mat2mp4_parser = clisubparsers.add_parser('mat2mp4',
                                              help='convert slices from .mat to .mp4 files',
                                              parents=[basename, mat2dat, dat2png, png2mp4])
    mat2mp4_parser.set_defaults(which='mat2mp4')
    return


def cli_parser():
    """
    Create the slices_converter.py Command Line Interface (CLI).

    Returns
    -------
    parser : argparse.ArgumentParser()
    """
    cliparser = argparse.ArgumentParser(prog="slices_converter.py",
                                        description="convert ADAM slices files to .dat/.png/.mp4 files",
                                        formatter_class=argparse.RawDescriptionHelpFormatter,
                                        epilog="Examples" +
                                               "\n  slices_converter.py mat2dat -b slice_01" +
                                               "\n  slices_converter.py mat2dat -b slice_01 -vnames VARIABLES = "+'"x" "y" "z" "rho"' +
                                               "\n  slices_converter.py dat2png -b slice_01 -v 3 -levels 1 12 32" +
                                               "\n  slices_converter.py png2mp4 -b slice_01 -mp4 slice_01.mp4" +
                                               "\n  slices_converter.py mat2mp4 -b slice_01 -mp4 slice_01.mp4" +
                                               "\nFor more detailed commands help use" +
                                               "\n  slices_converter.py mat2dat -h,--help" +
                                               "\n  slices_converter.py dat2png -h,--help" +
                                               "\n  slices_converter.py png2mp4 -h,--help" +
                                               "\n  slices_converter.py mat2mp4 -h,--help")
    clisubparsers = cliparser.add_subparsers(title='Commands', description='Valid commands')
    _parser_mat2dat(clisubparsers)
    _parser_dat2png(clisubparsers)
    _parser_png2mp4(clisubparsers)
    _parser_mat2mp4(clisubparsers)
    return cliparser

def is_file_to_convert(src_file, dst_file):
    """
    Check if file is to be converted or is already converted (mat2dat, dat2png).

    Parameters
    ----------
    src_file : str, source file to be converted
    dst_file : str, destination file converted

    Returns
    -------
    is_to_convert : bool, true or false
    """
    is_to_convert = False
    if os.path.exists(dst_file):
       if os.path.getmtime(src_file) >= os.path.getmtime(dst_file):
          is_to_convert = True
    else:
       is_to_convert = True
    return is_to_convert


def run_mat2dat(slice_files_pattern, vnames):
    """
    Run mat2dat command.

    Parameters
    ----------
    slice_files_pattern : str
    vnames : str
    """
    for file_name in os.listdir("./"):
       if file_name.startswith(slice_files_pattern) and file_name.endswith(".mat"):
          base_name = file_name.split('.')[0]
          if is_file_to_convert(file_name, base_name + '.dat'):
             print('convert file ' + file_name)
             os.system('ascot -i ' + file_name + ' -v ' + vnames + ' -o ' + base_name + '.dat')
          else:
             print('file ' + file_name + ' already converted')


def run_dat2png(slice_files_pattern, v, clabel, levels):
    """
    Run dat2png command.

    Parameters
    ----------
    slice_files_pattern : str
    v : int
    levels : list
    """
    for file_name in os.listdir("./"):
       if file_name.startswith(slice_files_pattern) and file_name.endswith(".dat"):
          base_name = file_name.split('.')[0]
          if is_file_to_convert(file_name, base_name + '.png'):
             print('convert file ' + file_name)
             # load data
             with open(file_name) as sf:
                head = [next(sf) for x in range(2)]
             # number of variables saved
             hnv = head[0].strip()
             nv = len(hnv.split('= ')[1].split(' '))
             # number of slice points
             hnijk = head[1].strip()
             ni = int(hnijk.split(',')[1].split('=')[1])
             nj = int(hnijk.split(',')[2].split('=')[1])
             nk = int(hnijk.split(',')[3].split('=')[1])
             slice_o = np.loadtxt(file_name, delimiter=" ", skiprows=2)
             slice_r = slice_o.reshape(ni,nj,nk,nv)
             slice_type = '3D'
             if   ni == 1 and nj == 1:
                slice_type = '1D'
                x = slice_r[0,0,:,2]
                r = slice_r[0,0,:,v]
                x_label = 'Z'
                y_label = clabel
             elif ni == 1 and nk == 1:
                slice_type = '1D'
                x = slice_r[0,:,0,1]
                r = slice_r[0,:,0,v]
                x_label = 'Y'
                y_label = clabel
             elif nj == 1 and nk == 1:
                slice_type = '1D'
                x = slice_r[:,0,0,0]
                r = slice_r[:,0,0,v]
                x_label = 'X'
                y_label = clabel
             elif ni == 1:
                slice_type = '2D'
                x = slice_r[0,0,:,1]
                y = slice_r[0,:,0,2]
                r = slice_r[0,:,:,v]
                x_label = 'Y'
                y_label = 'Z'
             elif nj == 1:
                slice_type = '2D'
                x = slice_r[0,0,:,0]
                y = slice_r[:,0,0,2]
                r = slice_r[:,0,:,v]
                x_label = 'X'
                y_label = 'Z'
             elif nk == 1:
                slice_type = '2D'
                x = slice_r[0,:,0,0]
                y = slice_r[:,0,0,1]
                r = slice_r[:,:,0,v]
                x_label = 'X'
                y_label = 'Y'
             if slice_type == '2D':
                # create contour plot
                fig, ax = plt.subplots(1, 1)
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_label)
                ax.set_title('2D Slice ' + x_label + y_label)
                if levels is not None:
                   clevels = np.linspace(levels[0], levels[1], int(levels[2]) + 1)
                   cf = plt.contourf(x, y, r, levels=clevels, cmap='RdGy')
                else :
                   cf = plt.contourf(x, y, r, levels=20, cmap='RdGy')
                cbar = plt.colorbar(cf)
                cbar.ax.set_ylabel(clabel, labelpad=10, rotation=270)
                plt.savefig(base_name + '.png')
                plt.close()
             elif slice_type == '1D':
                # create xy plot
                fig, ax = plt.subplots(1, 1)
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_label)
                ax.set_title('1D Slice ' + x_label)
                if levels is not None:
                   ax.set_ylim([levels[0], levels[1]])
                plt.grid()
                xy = plt.plot(x, r)
                plt.savefig(base_name + '.png')
                plt.close()
          else:
             print('file ' + file_name + ' already converted')


def run_png2mp4(slice_files_pattern, mp4):
    """
    Run png2mp4 command.

    Parameters
    ----------
    slice_files_pattern : str
    mp4 : str
    """
    print('generate MP4 video')
    os.system('ffmpeg -framerate 10 -pattern_type glob -i "' + slice_files_pattern + '*.png" ' + mp4)


if __name__ == '__main__':
    cliparser = cli_parser()
    cliargs = cliparser.parse_args()
    if cliargs.which == 'mat2dat':
        run_mat2dat(slice_files_pattern=cliargs.b, vnames=cliargs.vnames)
    if cliargs.which == 'dat2png':
        run_dat2png(slice_files_pattern=cliargs.b, v=cliargs.v, clabel=cliargs.clabel, levels=cliargs.levels)
    if cliargs.which == 'png2mp4':
        run_png2mp4(slice_files_pattern=cliargs.b, mp4=cliargs.mp4)
    if cliargs.which == 'mat2mp4':
        run_mat2dat(slice_files_pattern=cliargs.b, vnames=cliargs.vnames)
        run_dat2png(slice_files_pattern=cliargs.b, v=cliargs.v, clabel=cliargs.clabel, levels=cliargs.levels)
        run_png2mp4(slice_files_pattern=cliargs.b, mp4=cliargs.mp4)
    sys.exit(0)
