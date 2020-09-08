#ifndef FIELD_H
#define FIELD_H

#include <unordered_map>
#include "geometry.h"
#include "utils.h"

typedef struct {
    std::array<float,3> center;
    std::array<float,3> normal;
    float area;
} Triangle;

typedef struct {
	float ***tem;
	float ***rhs_tem;
	float ***prhs_tem;
	float ***work_ref;
    std::array<float***, 3> temp; //={NULL};
	int refinement_needed;
	Geometry geo;
    std::string key;
    std::unordered_map<std::string,float ***> dist_map;
    float ***dist;
    std::unordered_map<std::string,std::vector<Triangle>> triangles_map;
    std::string key_n_1[6], key_n_2[6], key_n_3[6], key_n_4[6], n_type[6];
    int coord_x_rel[6], coord_y_rel[6], coord_z_rel[6];
} Field;

void set_simple_immerse(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map, Params & params);

void prhs_compute(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map, Params & params);

void rhs_compute(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map, Params & params);

void linear_compute(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map, Params & params);

void evolve_advanced_immerse(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map);

void correct_advanced_immerse(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map);

void bc_phys(Field &f, std::string face);

void update_bc_block_field(std::string & face, std::string& n_type, Field & f, 
    Field & f_n_1, Field & f_n_2, Field & f_n_3, Field & f_n_4, 
    Geometry & geo, int coord_x_rel, int coord_y_rel, int coord_z_rel);

void derefine_and_merge_field(std::unordered_map<std::string,Field> & field_map, Geometry & geo_first_child, 
    std::string & new_key_str, std::vector<std::string> & key_str_list);

void refine_and_split_field(std::unordered_map<std::string,Field> & field_map, Geometry & geo_parent, 
    std::string & key_str, std::vector<std::string> & key_str_list);

void allocate_field(Field & field, std::vector<int> & np, std::vector<int> & ng);

void deallocate_field(Field & field, std::vector<int> & np, std::vector<int> & ng);

void init_field(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map);

void evolve(Field &f, std::string & immersion_type);

void mark_refinement_needed(Field &f);

std::vector<float> compute_force_field(Triangle & tri, Field & f, Geometry & geo);

void compute_totalforces(int it, std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map, Params & params);

float interp(std::array<float,8> & v, float dx1, float dx2, float dy1, float dy2, float dz1, float dz2);

void fix_edge_and_corners(std::unordered_map<std::string,Field> &field_map);

#endif
