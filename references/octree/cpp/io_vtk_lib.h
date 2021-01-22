#ifndef IO_VTK_LIB_H
#define IO_VTK_LIB_H

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <string>
#include <set>
#include <exception>
#include <iostream>
#include "b64.h"

std::string base64_encode(unsigned char const* , unsigned int len);
std::string base64_decode(std::string const& s);
char * b64_encode (const unsigned char *src, size_t len);

void write_vtr(std::vector<float*> fs, int nx, int ny, int nz, float *xg, float *yg, float *zg, std::string filename, std::vector<std::string> var_names);

void write_vtm(int n_blocks, std::string & filename_all, std::vector<std::string> &filenames);
#endif
