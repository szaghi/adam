#ifndef IO_VTK_H
#define IO_VTK_H

#include "io_vtk.h"
#include "field.h"

void write_vtr(float* f, float *f_2, float *f_3, float *f_4, float *f_5, float *f_6, 
               int nx, int ny, int nz, float *xg, float *yg, float *zg, char* filename);

void write_vtm(int n_blocks, char* filename_all, char* filenames);

void save_vtr(std::unordered_map<std::string, Field> Field_map, int iter, std::vector<int> & np, std::vector<int> & ng);

#endif
