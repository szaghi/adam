#include <iostream>
#include <vector>
#include "geometry.h"
#include "utils.h"
#include <cmath>

void init_geo(std::unordered_map<std::string, Geometry> & geo_map, std::vector<int> & np, std::vector<int> & ng, 
    float ls[3], float le[3], int forest[3]) {
   
    float ls_f[3], le_f[3];
    float delta;
    for(int i_f=0;i_f<forest[0];i_f++) {
        for(int j_f=0;j_f<forest[1];j_f++) {
            for(int k_f=0;k_f<forest[2];k_f++) {
                std::array<int,7> key_arr = {0,0,0,0,i_f,j_f,k_f};
                std::string key; KeyArrayToString(key_arr, key);
                Geometry geo;
                geo.key    = key;
                geo.coords = key_arr;
                
                ls_f[0] = ls[0] + (le[0]-ls[0])*i_f;
                le_f[0] = ls_f[0] + (le[0]-ls[0]);
                delta = (le[0]-ls[0])/np[0];
                geo.intervals[0] = {0, np[0]-1, np[0], ng[0]};
                geo.lengths[0]   = {ls_f[0], le_f[0],  le[0]-ls[0], delta};
                
                ls_f[1] = ls[1] + (le[1]-ls[1])*j_f;
                le_f[1] = ls_f[1] + (le[1]-ls[1]);
                delta = (le[1]-ls[1])/np[1];
                geo.intervals[1] = {0, np[1]-1, np[1], ng[1]};
                geo.lengths[1]   = {ls_f[1], le_f[1],  le[1]-ls[1], delta};
                
                ls_f[2] = ls[2] + (le[2]-ls[2])*k_f;
                le_f[2] = ls_f[2] + (le[2]-ls[2]);
                delta = (le[2]-ls[2])/np[2];
                geo.intervals[2] = {0, np[2]-1, np[2], ng[2]};
                geo.lengths[2]   = {ls_f[2], le_f[2],  le[2]-ls[2], delta};

                geo_map[key] = geo;
            }
        }
    }
}

void copy_geo(Geometry &geo_source, Geometry &geo_target) {
    geo_target.coords = geo_source.coords;
    geo_target.key = geo_source.key;
    for(int i=0;i<3;i++) {
        geo_target.intervals[i] = geo_source.intervals[i];
        geo_target.lengths[i]   = geo_source.lengths[i];
    }
}

//float mesh_law(float & distance) {
//    float min_delta = 0.05;
//    if(distance < 0.F) return min_delta;
//    float delta = min_delta + distance/3.;
//    return delta;
//}

float mesh_law(float & distance) {
    //1.5 Milioni
    //RIMETTEREfloat min_delta = 0.03;
    //RADOfloat min_delta = 0.06;
    float min_delta = 0.1;
    //3Milioni punti 
    //float min_delta = 0.02;
    //4Milioni punti - max_ref=7
    //float min_delta = 0.01;
    if(distance < 5.*min_delta) return min_delta;
    //RIMETTEREif(distance < 20.*min_delta) return min_delta;
    //    float abs_distance = fabsf(distance);
    float delta = min_delta + 10.*distance;
    return delta;
}
