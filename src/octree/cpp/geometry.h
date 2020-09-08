#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <array>
#include <vector>
#include <string>
#include <unordered_map>
typedef std::array<int,7> Cs7d;
typedef std::array<std::array<int,4>,3> Cin3d; 
typedef std::array<std::array<float,4>,3> Cfl3d; 

typedef struct{
    Cs7d  coords;
    Cin3d intervals; // start,end,size,size_ghosts
    Cfl3d lengths;   // start,end,size,delta
    std::string key;
	int refinement_needed;
}Geometry;

void init_geo(std::unordered_map<std::string, Geometry> & geo_map, std::vector<int> & np, 
    std::vector<int> & ng, float ls[3], float le[3], int forest[3]);

void copy_geo(Geometry &geo_source, Geometry &geo_target);

float mesh_law(float & distance);

#endif
