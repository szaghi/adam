#ifndef BC_H
#define BC_H

#include <string>
#include <unordered_map>
#include "geometry.h"
#include "field.h"

void fix_edge_and_corners_block(float *** f, int nx, int ny, int nz, int ngx, int ngy, int ngz);

void find_neighbour(int forest[3], std::string key, std::string direction, std::unordered_map<std::string,Geometry> &geo_map, 
        std::string &key_n_1, std::string &key_n_2, std::string &key_n_3, std::string &key_n_4,
        std::string &n_type, int & coord_x_rel, int & coord_y_rel, int & coord_z_rel);

void update_bc_maps(std::unordered_map<std::string,Geometry> &geo_map, 
        std::unordered_map<std::string,Field> &field_map, int forest[3]);

void update_bc(std::unordered_map<std::string,Geometry> &geo_map, 
        std::unordered_map<std::string,Field> &field_map);

void update_bc_block(std::string face, std::string n_type, float ***f, float ***f_ref, 
                     float ***f_n_1, float ***f_n_2, float ***f_n_3, float ***f_n_4,
                     int nx, int ny, int nz, int ngx, int ngy, int ngz,
                     int coord_x_rel, int coord_y_rel, int coord_z_rel);

#endif
