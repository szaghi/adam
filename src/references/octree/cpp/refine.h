#ifndef REFINE_H
#define REFINE_H

#include <unordered_map>
#include "geometry.h"
#include "field.h"
#include "immersed.h"

// SET HERE THE MAXIMUM NUMBER OF LEVELS
#define max_levels 6

void derefine_blocks(std::vector<std::string> key_str_list, std::unordered_map<std::string,Geometry> &geo_map, 
        std::unordered_map<std::string,Field> &field_map,
        std::unordered_map<std::string,Polytree> & immersed_map);

void derefine_and_merge(Geometry & geo, float***f, float***f_ref_glob,
        float *** f_ref_1, float *** f_ref_2, float *** f_ref_3, float *** f_ref_4, 
        float *** f_ref_5, float *** f_ref_6, float *** f_ref_7, float *** f_ref_8);

void derefine(float *** f, float *** f_deref, Geometry & geo);

float derefine_point(float *** f);

void refine_block(std::string key_str, std::unordered_map<std::string,Geometry> &geo_map, 
        std::unordered_map<std::string,Field> &field_map,
        std::unordered_map<std::string,Polytree> & immersed_map);

void refine_and_split(Geometry geo, float***f, float***f_ref_glob,
        float *** f_ref_1, float *** f_ref_2, float *** f_ref_3, float *** f_ref_4, 
        float *** f_ref_5, float *** f_ref_6, float *** f_ref_7, float *** f_ref_8);

void refine(float ***f, float ***f_ref, Geometry geo);

float*** refine_points(float f[3][3][3]);

std::vector<std::vector<std::string>> sanitize_refinement_needed(std::unordered_map<std::string,Geometry> &geo_map);

void update_refine(std::unordered_map<std::string,Geometry> &geo_map, 
        std::unordered_map<std::string,Field> &field_map,
        std::unordered_map<std::string,Polytree> &immersed_map, std::string const & mode, int forest[3]);

#endif
