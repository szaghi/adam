#include <iostream>
#include <unordered_map>
#include <string>
#include <cstdlib>
#include "geometry.h"
#include "field.h"
#include "bc.h"
#include "refine.h"


void fix_edge_and_corners_block(float *** f, int nx, int ny, int nz, int ngx, int ngy, int ngz) {
    int i,j,k;

    // x edges
    j = ngy-1  ;  k = ngz-1;
    for(i = ngx;i<ngx+nx;i++)
        f[i][j][k] = f[i][j+1][k] + f[i][j][k+1] - f[i][j+1][k+1] ;
    j = ngy+ny  ;  k = ngz-1;
    for(i = ngx;i<ngx+nx;i++)
        f[i][j][k] = f[i][j-1][k] + f[i][j][k+1] - f[i][j-1][k+1] ;
    j = ngy-1  ;  k = ngz+nz;
    for(i = ngx;i<ngx+nx;i++)
        f[i][j][k] = f[i][j+1][k] + f[i][j][k-1] - f[i][j+1][k-1] ;
    j = ngy+ny  ;  k = ngz+nz;
    for(i = ngx;i<ngx+nx;i++)
        f[i][j][k] = f[i][j-1][k] + f[i][j][k-1] - f[i][j-1][k-1] ;

    // y edges
    i = ngx-1  ;  k = ngz-1;
    for(j = ngy;j<ngy+ny;j++) 
        f[i][j][k] = f[i+1][j][k] + f[i][j][k+1] - f[i+1][j][k+1] ;
    i = ngx+nx  ;  k = ngz-1;
    for(j = ngy;j<ngy+ny;j++) 
        f[i][j][k] = f[i-1][j][k] + f[i][j][k+1] - f[i-1][j][k+1] ;
    i = ngx-1  ;  k = ngz+nz;
    for(j = ngy;j<ngy+ny;j++) 
        f[i][j][k] = f[i+1][j][k] + f[i][j][k-1] - f[i+1][j][k-1] ;
    i = ngx+nx  ;  k = ngz+nz;
    for(j = ngy;j<ngy+ny;j++) 
        f[i][j][k] = f[i-1][j][k] + f[i][j][k-1] - f[i-1][j][k-1] ;

    // z edges
    j = ngy-1  ;  i = ngx-1;
    for(k = ngz;k<ngz+nz;k++)
        f[i][j][k] = f[i][j+1][k] + f[i+1][j][k] - f[i+1][j+1][k] ;
    j = ngy+ny  ;  i = ngx-1;
    for(k = ngz;k<ngz+nz;k++)
        f[i][j][k] = f[i][j-1][k] + f[i+1][j][k] - f[i+1][j-1][k] ;
    j = ngy-1  ;  i = ngx+nx;
    for(k = ngz;k<ngz+nz;k++)
        f[i][j][k] = f[i][j+1][k] + f[i-1][j][k] - f[i-1][j+1][k] ;
    j = ngy+ny  ;  i = ngx+nx;
    for(k = ngz;k<ngz+nz;k++)
        f[i][j][k] = f[i][j-1][k] + f[i-1][j][k] - f[i-1][j-1][k] ;

    // corners (again with planar interpolation, plane arbitrarily chosen)
    i = ngx-1       ;  j = ngy-1       ;  k = ngz-1;
    f[i][j][k] = f[i+1][j][k] + f[i][j+1][k] - f[i+1][j+1][k] ;
    i = ngx-1       ;  j = ngy-1       ;  k = ngz+nz;
    f[i][j][k] = f[i+1][j][k] + f[i][j+1][k] - f[i+1][j+1][k] ;
    i = ngx-1       ;  j = ngy+ny  ;  k = ngz-1;
    f[i][j][k] = f[i+1][j][k] + f[i][j][k+1] - f[i+1][j][k+1] ;
    i = ngx+nx  ;  j = ngy-1       ;  k = ngz-1;
    f[i][j][k] = f[i][j+1][k] + f[i][j][k+1] - f[i][j+1][k+1] ;
    i = ngx-1       ;  j = ngy+ny  ;  k = ngz+nz;
    f[i][j][k] = f[i][j-1][k] + f[i][j][k-1] - f[i][j-1][k-1] ;
    i = ngx+nx  ;  j = ngy+ny  ;  k = ngz+nz;
    f[i][j][k] = f[i][j-1][k] + f[i][j][k-1] - f[i][j-1][k-1] ;
    i = ngx+nx  ;  j = ngy-1       ;  k = ngz+nz;
    f[i][j][k] = f[i][j+1][k] + f[i][j][k-1] - f[i][j+1][k-1] ;
    i = ngx+nx  ;  j = ngy+ny  ;  k = ngz-1;
    f[i][j][k] = f[i][j-1][k] + f[i][j][k+1] - f[i][j-1][k+1] ;

}

void find_neighbour(int forest[3], std::string key, std::string direction, std::unordered_map<std::string,Geometry> &geo_map, 
        std::string &key_n_1, std::string &key_n_2, std::string &key_n_3, std::string &key_n_4,
        std::string &n_type, int & coord_x_rel, int & coord_y_rel, int & coord_z_rel) {

    std::array<int,7> key_arr; KeyStringToArray(key, key_arr);
    int level = key_arr[0]; int coord_x = key_arr[1]; int coord_y = key_arr[2]; int coord_z = key_arr[3];
    int forest_x = key_arr[4]; int forest_y = key_arr[5]; int forest_z = key_arr[6]; 
    int size_xyz;
    if(level < max_levels) { // cannot shift with negative values
        size_xyz = (2 << (max_levels-level-1)); 
    } else {
        size_xyz = 1;
    }
    int start_x = coord_x; int end_x = coord_x + size_xyz - 1;
    int start_y = coord_y; int end_y = coord_y + size_xyz - 1;
    int start_z = coord_z; int end_z = coord_z + size_xyz - 1;
    std::array<int,7> key_arr_n_1, key_arr_n_2, key_arr_n_3, key_arr_n_4;;
    bool found;
    int level_n, cx_n, cy_n, cz_n = coord_z;

    if(direction == "left") {
        if(start_x == 0) {
            if(forest_x == 0) {
                n_type =  "phys";
            } else {
                // 1st option: same level
                level_n = level; cx_n = (2 << (max_levels-1)) - size_xyz; cy_n = coord_y; cz_n = coord_z;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x-1, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                auto t = geo_map.find(key_n_1);
                found = false;
                // 1st option: neighbour is same refined
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND SAME REFINED NEIGHBOUR MIN-X" << std::endl;
                    found = true; n_type = "n_same_refined";
                    key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
                }
                // 2nd option: neighbour is more refined, there are 4 of them
                if(!found) {
                    level_n = level+1; cx_n = (2 << (max_levels-1)) - size_xyz/2; cy_n = coord_y; cz_n = coord_z;
                    key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x-1, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                    t = geo_map.find(key_n_1);
                    if (t != geo_map.end()) { // key_n_1 found
                        //std::cout << "FOUND MORE REFINED NEIGHBOUR MIN-X" << std::endl;
                        key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x-1, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                        key_arr_n_2 = {level_n, cx_n, cy_n+size_xyz/2, cz_n, forest_x-1, forest_y, forest_z}; KeyArrayToString(key_arr_n_2, key_n_2);
                        key_arr_n_3 = {level_n, cx_n, cy_n, cz_n+size_xyz/2, forest_x-1, forest_y, forest_z}; KeyArrayToString(key_arr_n_3, key_n_3);
                        key_arr_n_4 = {level_n, cx_n, cy_n+size_xyz/2, cz_n+size_xyz/2, forest_x-1, forest_y, forest_z}; KeyArrayToString(key_arr_n_4, key_n_4);
                        found = true; n_type =  "n_more_refined";
                    }
                }
                // 3rd option: neighbour is less refined
                if(!found) {
                    // getting the position of myself wrt the neighbouring coarser interface
                    coord_y_rel = coord_y%(2*size_xyz); coord_z_rel = coord_z%(2*size_xyz);
                    level_n = level-1; cx_n = (2 << (max_levels-1)) - 2*size_xyz; cy_n = coord_y-coord_y_rel; cz_n = coord_z-coord_z_rel;
                    coord_y_rel = coord_y_rel/size_xyz ; coord_z_rel = coord_z_rel/size_xyz;
                    key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x-1, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                    t = geo_map.find(key_n_1);
                    if (t != geo_map.end()) { // key_n_1 found
                        //std::cout << "FOUND LESS REFINED NEIGHBOUR MIN-X" << std::endl;
                        found = true; n_type =  "n_less_refined";
                        key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
                    }
                }
                if(!found) {
                    std::cout << "CANNOT FIND FOREST_X NEIGHBOUR MIN-X! Stopping execution..." << std::endl; exit(EXIT_FAILURE);
                }
            }
        } else {
            // 1st option: same level
            level_n = level; cx_n = coord_x - size_xyz; cy_n = coord_y; cz_n = coord_z;
            key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
            auto t = geo_map.find(key_n_1);
            found = false;
            // 1st option: neighbour is same refined
            if (t != geo_map.end()) { // key_n_1 found
                //std::cout << "FOUND SAME REFINED NEIGHBOUR MIN-X" << std::endl;
                found = true; n_type = "n_same_refined";
                key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
            }
            // 2nd option: neighbour is more refined, there are 4 of them
            if(!found) {
                level_n = level+1; cx_n = coord_x - size_xyz/2; cy_n = coord_y; cz_n = coord_z;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                t = geo_map.find(key_n_1);
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND MORE REFINED NEIGHBOUR MIN-X" << std::endl;
                    key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                    key_arr_n_2 = {level_n, cx_n, cy_n+size_xyz/2, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_2, key_n_2);
                    key_arr_n_3 = {level_n, cx_n, cy_n, cz_n+size_xyz/2, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_3, key_n_3);
                    key_arr_n_4 = {level_n, cx_n, cy_n+size_xyz/2, cz_n+size_xyz/2, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_4, key_n_4);
                    found = true; n_type =  "n_more_refined";
                }
            }
            // 3rd option: neighbour is less refined
            if(!found) {
                // getting the position of myself wrt the neighbouring coarser interface
                coord_y_rel = coord_y%(2*size_xyz); coord_z_rel = coord_z%(2*size_xyz);
                level_n = level-1; cx_n = coord_x-2*size_xyz; cy_n = coord_y-coord_y_rel; cz_n = coord_z-coord_z_rel;
                coord_y_rel = coord_y_rel/size_xyz ; coord_z_rel = coord_z_rel/size_xyz;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                t = geo_map.find(key_n_1);
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND LESS REFINED NEIGHBOUR MIN-X" << std::endl;
                    found = true; n_type =  "n_less_refined";
                    key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
                }
            }
            if(!found) {
                std::cout << "CANNOT FIND NEIGHBOUR MIN-X! Stopping execution..." << std::endl; exit(EXIT_FAILURE);
            }
        }
    }

    if(direction == "right") {
        if(end_x+1 ==  (2 << (max_levels-1))) {
            if(forest_x == (forest[0]-1)) {  // TODO GENERALIZE RIMETTERE
                n_type =  "phys";
            } else {
                // 1st option: same level
                level_n = level; cx_n = 0; cy_n = coord_y; cz_n = coord_z;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x+1, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                auto t = geo_map.find(key_n_1);
                found = false;
                // 1st option: neighbour is same refined
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND SAME REFINED NEIGHBOUR MAX-X" << std::endl;
                    found = true; n_type = "n_same_refined";
                    key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
                }
                // 2nd option: neighbour is more refined, there are 4 of them
                if(!found) {
                    level_n = level+1; cx_n = 0; cy_n = coord_y; cz_n = coord_z;
                    key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x+1, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                    t = geo_map.find(key_n_1);
                    if (t != geo_map.end()) { // key_n_1 found
                        //std::cout << "FOUND MORE REFINED NEIGHBOUR MAX-X" << std::endl;
                        key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x+1, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                        key_arr_n_2 = {level_n, cx_n, cy_n+size_xyz/2, cz_n, forest_x+1, forest_y, forest_z}; KeyArrayToString(key_arr_n_2, key_n_2);
                        key_arr_n_3 = {level_n, cx_n, cy_n, cz_n+size_xyz/2, forest_x+1, forest_y, forest_z}; KeyArrayToString(key_arr_n_3, key_n_3);
                        key_arr_n_4 = {level_n, cx_n, cy_n+size_xyz/2, cz_n+size_xyz/2, forest_x+1, forest_y, forest_z}; KeyArrayToString(key_arr_n_4, key_n_4);
                        found = true; n_type =  "n_more_refined";
                    }
                }
                // 3rd option: neighbour is less refined
                if(!found) {
                    // getting the position of myself wrt the neighbouring coarser interface
                    coord_y_rel = coord_y%(2*size_xyz); coord_z_rel = coord_z%(2*size_xyz);
                    level_n = level-1; cx_n = 0; cy_n = coord_y-coord_y_rel; cz_n = coord_z-coord_z_rel;
                    coord_y_rel = coord_y_rel/size_xyz ; coord_z_rel = coord_z_rel/size_xyz;
                    key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x+1, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                    //std::cout << "MAX-X: SEARCHING LESS REFINED FOR KEY: " << key << " AS: " << key_n_1 << std::endl;
                    t = geo_map.find(key_n_1);
                    if (t != geo_map.end()) { // key_n_1 found
                        //std::cout << "FOUND LESS REFINED NEIGHBOUR MAX-X" << std::endl;
                        found = true; n_type =  "n_less_refined";
                        key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
                    }
                }
                if(!found) {
                    std::cout << "CANNOT FIND FOREST NEIGHBOUR MAX-X! " << key << " Stopping execution..." << std::endl; exit(EXIT_FAILURE);
                }
            }
        } else {
            // 1st option: same level
            level_n = level; cx_n = coord_x + size_xyz; cy_n = coord_y; cz_n = coord_z;
            key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
            auto t = geo_map.find(key_n_1);
            found = false;
            // 1st option: neighbour is same refined
            if (t != geo_map.end()) { // key_n_1 found
                //std::cout << "FOUND SAME REFINED NEIGHBOUR MAX-X" << std::endl;
                found = true; n_type = "n_same_refined";
                key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
            }
            // 2nd option: neighbour is more refined, there are 4 of them
            if(!found) {
                level_n = level+1; cx_n = coord_x + size_xyz; cy_n = coord_y; cz_n = coord_z;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                t = geo_map.find(key_n_1);
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND MORE REFINED NEIGHBOUR MAX-X" << std::endl;
                    key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                    key_arr_n_2 = {level_n, cx_n, cy_n+size_xyz/2, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_2, key_n_2);
                    key_arr_n_3 = {level_n, cx_n, cy_n, cz_n+size_xyz/2, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_3, key_n_3);
                    key_arr_n_4 = {level_n, cx_n, cy_n+size_xyz/2, cz_n+size_xyz/2, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_4, key_n_4);
                    found = true; n_type =  "n_more_refined";
                }
            }
            // 3rd option: neighbour is less refined
            if(!found) {
                // getting the position of myself wrt the neighbouring coarser interface
                coord_y_rel = coord_y%(2*size_xyz); coord_z_rel = coord_z%(2*size_xyz);
                level_n = level-1; cx_n = coord_x+size_xyz; cy_n = coord_y-coord_y_rel; cz_n = coord_z-coord_z_rel;
                coord_y_rel = coord_y_rel/size_xyz ; coord_z_rel = coord_z_rel/size_xyz;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                //std::cout << "MAX-X: SEARCHING LESS REFINED FOR KEY: " << key << " AS: " << key_n_1 << std::endl;
                t = geo_map.find(key_n_1);
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND LESS REFINED NEIGHBOUR MAX-X" << std::endl;
                    found = true; n_type =  "n_less_refined";
                    key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
                }
            }
            if(!found) {
                std::cout << "CANNOT FIND NEIGHBOUR MAX-X! " << key << " Stopping execution..." << std::endl; exit(EXIT_FAILURE);
            }
        }
    }

    if(direction == "bottom") {
        if(start_y == 0) {
            n_type =  "phys";
        } else {
            // 1st option: same level
            level_n = level; cx_n = coord_x; cy_n = coord_y-size_xyz; cz_n = coord_z;
            key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
            auto t = geo_map.find(key_n_1);
            found = false;
            // 1st option: neighbour is same refined
            if (t != geo_map.end()) { // key_n_1 found
                //std::cout << "FOUND SAME REFINED NEIGHBOUR MIN-Y" << std::endl;
                found = true; n_type = "n_same_refined";
                key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
            }
            // 2nd option: neighbour is more refined, there are 4 of them
            if(!found) {
                level_n = level+1; cx_n = coord_x; cy_n = coord_y - size_xyz/2; cz_n = coord_z;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                t = geo_map.find(key_n_1);
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND MORE REFINED NEIGHBOUR MIN-Y" << std::endl;
                    key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                    key_arr_n_2 = {level_n, cx_n+size_xyz/2, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_2, key_n_2);
                    key_arr_n_3 = {level_n, cx_n, cy_n, cz_n+size_xyz/2, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_3, key_n_3);
                    key_arr_n_4 = {level_n, cx_n+size_xyz/2, cy_n, cz_n+size_xyz/2, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_4, key_n_4);
                    found = true; n_type =  "n_more_refined";
                }
            }
            // 3rd option: neighbour is less refined
            if(!found) {
                // getting the position of myself wrt the neighbouring coarser interface
                coord_x_rel = coord_x%(2*size_xyz); coord_z_rel = coord_z%(2*size_xyz);
                level_n = level-1; cx_n = coord_x-coord_x_rel; cy_n = coord_y-2*size_xyz; cz_n = coord_z-coord_z_rel;
                coord_x_rel = coord_x_rel/size_xyz ; coord_z_rel = coord_z_rel/size_xyz;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                //std::cout << "MIN-Y: SEARCHING LESS REFINED FOR KEY: " << key << " AS: " << key_n_1 << std::endl;
                t = geo_map.find(key_n_1);
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND LESS REFINED NEIGHBOUR MIN-Y" << std::endl;
                    found = true; n_type =  "n_less_refined";
                    key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
                }
            }
            if(!found) {
                std::cout << "CANNOT FIND NEIGHBOUR MIN-Y! Stopping execution..." << std::endl; exit(EXIT_FAILURE);
            }
        }
    }

    if(direction == "top") {
        if(end_y+1 ==  (2 << (max_levels-1))) {
            n_type =  "phys";
        } else {
            // 1st option: same level
            level_n = level; cx_n = coord_x; cy_n = coord_y+size_xyz; cz_n = coord_z;
            key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
            auto t = geo_map.find(key_n_1);
            found = false;
            // 1st option: neighbour is same refined
            if (t != geo_map.end()) { // key_n_1 found
                //std::cout << "FOUND SAME REFINED NEIGHBOUR MAX-Y" << std::endl;
                found = true; n_type = "n_same_refined";
                key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
            }
            // 2nd option: neighbour is more refined, there are 4 of them
            if(!found) {
                level_n = level+1; cx_n = coord_x; cy_n = coord_y + size_xyz; cz_n = coord_z;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                t = geo_map.find(key_n_1);
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND MORE REFINED NEIGHBOUR MAX-Y" << std::endl;
                    key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                    key_arr_n_2 = {level_n, cx_n+size_xyz/2, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_2, key_n_2);
                    key_arr_n_3 = {level_n, cx_n, cy_n, cz_n+size_xyz/2, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_3, key_n_3);
                    key_arr_n_4 = {level_n, cx_n+size_xyz/2, cy_n, cz_n+size_xyz/2, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_4, key_n_4);
                    found = true; n_type =  "n_more_refined";
                }
            }
            // 3rd option: neighbour is less refined
            if(!found) {
                // getting the position of myself wrt the neighbouring coarser interface
                coord_x_rel = coord_x%(2*size_xyz); coord_z_rel = coord_z%(2*size_xyz);
                level_n = level-1; cx_n = coord_x-coord_x_rel; cy_n = coord_y+size_xyz; cz_n = coord_z-coord_z_rel;
                coord_x_rel = coord_x_rel/size_xyz ; coord_z_rel = coord_z_rel/size_xyz;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                //std::cout << "MAX-Y: SEARCHING LESS REFINED FOR KEY: " << key << " AS: " << key_n_1 << std::endl;
                t = geo_map.find(key_n_1);
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND LESS REFINED NEIGHBOUR MAX-Y" << std::endl;
                    found = true; n_type =  "n_less_refined";
                    key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
                }
            }
            if(!found) {
                std::cout << "CANNOT FIND NEIGHBOUR MAX-Y! " << key << " Stopping execution..." << std::endl; exit(EXIT_FAILURE);
            }
        }
    }

    if(direction == "back") {
        if(start_z == 0) {
            n_type =  "phys";
        } else {
            // 1st option: same level
            level_n = level; cx_n = coord_x; cy_n = coord_y; cz_n = coord_z-size_xyz;
            key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
            //std::cout << "BACK " << key << "SEARCHING SAME LEVEL KEY: " << key_n_1 << std::endl;
            auto t = geo_map.find(key_n_1);
            found = false;
            // 1st option: neighbour is same refined
            if (t != geo_map.end()) { // key_n_1 found
                //std::cout << "FOUND SAME REFINED NEIGHBOUR MIN-Z" << std::endl;
                found = true; n_type = "n_same_refined";
                key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
            }
            // 2nd option: neighbour is more refined, there are 4 of them
            if(!found) {
                level_n = level+1; cx_n = coord_x; cy_n = coord_y; cz_n = coord_z - size_xyz/2;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                //std::cout << "BACK " << key << "SEARCHING MORE LEVEL KEY: " << key_n_1 << std::endl;
                t = geo_map.find(key_n_1);
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND MORE REFINED NEIGHBOUR MIN-Z" << std::endl;
                    key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                    key_arr_n_2 = {level_n, cx_n+size_xyz/2, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_2, key_n_2);
                    key_arr_n_3 = {level_n, cx_n, cy_n+size_xyz/2, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_3, key_n_3);
                    key_arr_n_4 = {level_n, cx_n+size_xyz/2, cy_n+size_xyz/2, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_4, key_n_4);
                    found = true; n_type =  "n_more_refined";
                }
            }
            // 3rd option: neighbour is less refined
            if(!found) {
                // getting the position of myself wrt the neighbouring coarser interface
                coord_x_rel = coord_x%(2*size_xyz); coord_y_rel = coord_y%(2*size_xyz);
                level_n = level-1; cx_n = coord_x-coord_x_rel; cy_n = coord_y-coord_y_rel; cz_n = coord_z-2*size_xyz;
                coord_x_rel = coord_x_rel/size_xyz ; coord_y_rel = coord_y_rel/size_xyz;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                //std::cout << "BACK " << key << "SEARCHING LESS LEVEL KEY: " << key_n_1 << std::endl;
                t = geo_map.find(key_n_1);
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND LESS REFINED NEIGHBOUR MIN-Z" << std::endl;
                    found = true; n_type =  "n_less_refined";
                    key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
                }
            }
            if(!found) {
                std::cout << "CANNOT FIND NEIGHBOUR MIN-Z! " << key << " Stopping execution..." << std::endl; exit(EXIT_FAILURE);
            }
        }
    }

    if(direction == "front") {
        if(end_z+1 ==  (2 << (max_levels-1))) {
            n_type =  "phys";
        } else {
            // 1st option: same level
            level_n = level; cx_n = coord_x; cy_n = coord_y; cz_n = coord_z+size_xyz;
            key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
            auto t = geo_map.find(key_n_1);
            found = false;
            // 1st option: neighbour is same refined
            if (t != geo_map.end()) { // key_n found
                //std::cout << "FOUND SAME REFINED NEIGHBOUR MAX-Z" << std::endl;
                found = true; n_type = "n_same_refined";
                key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
            }
            // 2nd option: neighbour is more refined, there are 4 of them
            if(!found) {
                level_n = level+1; cx_n = coord_x; cy_n = coord_y; cz_n = coord_z + size_xyz;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                t = geo_map.find(key_n_1);
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND MORE REFINED NEIGHBOUR MAX-Z" << std::endl;
                    key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                    key_arr_n_2 = {level_n, cx_n+size_xyz/2, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_2, key_n_2);
                    key_arr_n_3 = {level_n, cx_n, cy_n+size_xyz/2, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_3, key_n_3);
                    key_arr_n_4 = {level_n, cx_n+size_xyz/2, cy_n+size_xyz/2, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_4, key_n_4);
                    found = true; n_type =  "n_more_refined";
                }
            }
            // 3rd option: neighbour is less refined
            if(!found) {
                // getting the position of myself wrt the neighbouring coarser interface
                coord_x_rel = coord_x%(2*size_xyz); coord_y_rel = coord_y%(2*size_xyz);
                level_n = level-1; cx_n = coord_x-coord_x_rel; cy_n = coord_y-coord_y_rel; cz_n = coord_z+size_xyz;
                coord_x_rel = coord_x_rel/size_xyz ; coord_y_rel = coord_y_rel/size_xyz;
                key_arr_n_1 = {level_n, cx_n, cy_n, cz_n, forest_x, forest_y, forest_z}; KeyArrayToString(key_arr_n_1, key_n_1);
                t = geo_map.find(key_n_1);
                if (t != geo_map.end()) { // key_n_1 found
                    //std::cout << "FOUND LESS REFINED NEIGHBOUR MAX-Z" << std::endl;
                    found = true; n_type =  "n_less_refined";
                    key_n_2 = key_n_3 = key_n_4 = key_n_1 ; // for next convenience
                }
            }
            if(!found) {
                std::cout << "CANNOT FIND NEIGHBOUR MAX-Z! Stopping execution..." << std::endl; exit(EXIT_FAILURE);
            }
        }
    }
}

void update_bc_block(std::string face, std::string n_type, float ***f, float ***f_ref, 
                     float ***f_n_1, float ***f_n_2, float ***f_n_3, float ***f_n_4,
                     int nx, int ny, int nz, int ngx, int ngy, int ngz,
                     int coord_x_rel, int coord_y_rel, int coord_z_rel) {

    int i, j, k, i_ref, j_ref, k_ref;

    if(face == "left") {
        if(n_type == "n_same_refined") {
            //std::cout << "MIN-X SAME REFINED" << std::endl;
            for(i=0;i<ngx;i++) 
                for(j=ngy;j<ngy+ny;j++) 
                    for(k=ngz;k<ngz+nz;k++) 
                        f[i][j][k] = f_n_1[nx+i][j][k];
        } else if (n_type == "n_more_refined") {
            //std::cout << "MIN-X MORE REFINED" << std::endl;
            for(i=0;i<2*ngx;i++) // need double ngx layer to de-refine and get ngx layer
                for(j=ngy;j<ngy+ny;j++) 
                    for(k=ngz;k<ngz+nz;k++) {
                        f_ref[i][j][k]       = f_n_1[i+nx-ngx][j][k];
                        f_ref[i][j+ny][k]    = f_n_2[i+nx-ngx][j][k];
                        f_ref[i][j][k+nz]    = f_n_3[i+nx-ngx][j][k];
                        f_ref[i][j+ny][k+nz] = f_n_4[i+nx-ngx][j][k];
                    }
            for(i=0;i<ngx;i++) // de-refine
                for(j=ngy;j<ngy+ny;j++) 
                    for(k=ngz;k<ngz+nz;k++) {
                        //BUGi_ref = 2*i; j_ref = 2*j; k_ref = 2*k;
                        i_ref = 2*i; j_ref = 2*(j-ngy)+ngy; k_ref = 2*(k-ngz)+ngz;
                        f[i][j][k] = 0.125*(
                           f_ref[i_ref][j_ref][k_ref]     + f_ref[i_ref+1][j_ref][k_ref] + 
                           f_ref[i_ref][j_ref+1][k_ref]   + f_ref[i_ref+1][j_ref+1][k_ref] + 
                           f_ref[i_ref][j_ref][k_ref+1]   + f_ref[i_ref+1][j_ref][k_ref+1] + 
                           f_ref[i_ref][j_ref+1][k_ref+1] + f_ref[i_ref+1][j_ref+1][k_ref+1]);
                    }
        } else if (n_type == "n_less_refined") {
            //std::cout << "MIN-X LESS REFINED" << std::endl;
            for(i=ngx/2+nx;i<ngx+nx;i++) // loop over neighbour points
                for(j=ngy+coord_y_rel*ny/2;j<ngy+(coord_y_rel+1)*ny/2;j++) 
                    for(k=ngz+coord_z_rel*nz/2;k<ngz+(coord_z_rel+1)*nz/2;k++) {
                        i_ref = 2*(i-ngx/2-nx); // DA VERIFICARE TUTTI E 3
                        j_ref = 2*(j-ngy-coord_y_rel*ny/2)+ngy;
                        k_ref = 2*(k-ngz-coord_z_rel*nz/2)+ngz;
                        f[i_ref][j_ref][k_ref]       = f_n_1[i][j][k];
                        f[i_ref+1][j_ref][k_ref]     = f_n_1[i][j][k];
                        f[i_ref][j_ref+1][k_ref]     = f_n_1[i][j][k];
                        f[i_ref+1][j_ref+1][k_ref]   = f_n_1[i][j][k];
                        f[i_ref][j_ref][k_ref+1]     = f_n_1[i][j][k];
                        f[i_ref+1][j_ref][k_ref+1]   = f_n_1[i][j][k];
                        f[i_ref][j_ref+1][k_ref+1]   = f_n_1[i][j][k];
                        f[i_ref+1][j_ref+1][k_ref+1] = f_n_1[i][j][k];
                    }
        }
    }

    if(face == "right") {
        if(n_type == "n_same_refined") {
            //std::cout << "MAX-X SAME REFINED" << std::endl;
            for(int i=0;i<ngx;i++) 
                for(int j=ngy;j<ngy+ny;j++) 
                    for(int k=ngz;k<ngz+nz;k++) 
                        f[i+nx+ngx][j][k] = f_n_1[ngx+i][j][k];
        } else if (n_type == "n_more_refined") {
            //std::cout << "MAX-X MORE REFINED" << std::endl;
            for(i=0;i<2*ngx;i++) // need double ngx layer to de-refine and get ngx layer
                for(j=ngy;j<ngy+ny;j++) 
                    for(k=ngz;k<ngz+nz;k++) {
                        f_ref[i][j][k]       = f_n_1[i+ngx][j][k];
                        f_ref[i][j+ny][k]    = f_n_2[i+ngx][j][k];
                        f_ref[i][j][k+nz]    = f_n_3[i+ngx][j][k];
                        f_ref[i][j+ny][k+nz] = f_n_4[i+ngx][j][k];
                    }
            for(i=0;i<ngx;i++) // de-refine
                for(j=ngy;j<ngy+ny;j++) 
                    for(k=ngz;k<ngz+nz;k++) {
                        i_ref = 2*i; j_ref = 2*(j-ngy)+ngy; k_ref = 2*(k-ngz)+ngz;
                        f[i+nx+ngx][j][k] = 0.125*(
                           f_ref[i_ref][j_ref][k_ref]     + f_ref[i_ref+1][j_ref][k_ref] + 
                           f_ref[i_ref][j_ref+1][k_ref]   + f_ref[i_ref+1][j_ref+1][k_ref] + 
                           f_ref[i_ref][j_ref][k_ref+1]   + f_ref[i_ref+1][j_ref][k_ref+1] + 
                           f_ref[i_ref][j_ref+1][k_ref+1] + f_ref[i_ref+1][j_ref+1][k_ref+1]);
                    }
        } else if (n_type == "n_less_refined") {
            //std::cout << "MAX-X LESS REFINED" << std::endl;
            for(i=ngx;i<ngx+ngx/2;i++) // loop over neighbour points
                for(j=ngy+coord_y_rel*ny/2;j<ngy+(coord_y_rel+1)*ny/2;j++) 
                    for(k=ngz+coord_z_rel*nz/2;k<ngz+(coord_z_rel+1)*nz/2;k++) {
                        i_ref = 2*(i-ngx)+nx+ngx; // DA VERIFICARE TUTTI E 3
                        j_ref = 2*(j-ngy-coord_y_rel*ny/2)+ngy;
                        k_ref = 2*(k-ngz-coord_z_rel*nz/2)+ngz;
                        f[i_ref][j_ref][k_ref]       = f_n_1[i][j][k];
                        f[i_ref+1][j_ref][k_ref]     = f_n_1[i][j][k];
                        f[i_ref][j_ref+1][k_ref]     = f_n_1[i][j][k];
                        f[i_ref+1][j_ref+1][k_ref]   = f_n_1[i][j][k];
                        f[i_ref][j_ref][k_ref+1]     = f_n_1[i][j][k];
                        f[i_ref+1][j_ref][k_ref+1]   = f_n_1[i][j][k];
                        f[i_ref][j_ref+1][k_ref+1]   = f_n_1[i][j][k];
                        f[i_ref+1][j_ref+1][k_ref+1] = f_n_1[i][j][k];
                    }
        }
    }

    if(face == "bottom") {
        if(n_type == "n_same_refined") {
            //std::cout << "MIN-Y SAME REFINED" << std::endl;
            for(int i=ngx;i<ngx+nx;i++) 
                for(int j=0;j<ngy;j++) 
                    for(int k=ngz;k<ngz+nz;k++) 
                        f[i][j][k] = f_n_1[i][ny+j][k];
        } else if (n_type == "n_more_refined") {
            //std::cout << "MIN-Y MORE REFINED" << std::endl;
            for(i=ngx;i<ngx+nx;i++) // need double ngx layer to de-refine and get ngx layer
                for(j=0;j<2*ngy;j++) 
                    for(k=ngz;k<ngz+nz;k++) {
                        f_ref[i][j][k]       = f_n_1[i][j+ny-ngy][k];
                        f_ref[i+nx][j][k]    = f_n_2[i][j+ny-ngy][k];
                        f_ref[i][j][k+nz]    = f_n_3[i][j+ny-ngy][k];
                        f_ref[i+nx][j][k+nz] = f_n_4[i][j+ny-ngy][k];
                    }
            for(i=ngx;i<ngx+nx;i++) // de-refine
                for(j=0;j<ngy;j++) 
                    for(k=ngz;k<ngz+nz;k++) {
                        //BUGi_ref = 2*i; j_ref = 2*j; k_ref = 2*k;
                        i_ref = 2*(i-ngx)+ngx; j_ref = 2*j; k_ref = 2*(k-ngz)+ngz;
                        f[i][j][k] = 0.125*(
                           f_ref[i_ref][j_ref][k_ref]     + f_ref[i_ref+1][j_ref][k_ref] + 
                           f_ref[i_ref][j_ref+1][k_ref]   + f_ref[i_ref+1][j_ref+1][k_ref] + 
                           f_ref[i_ref][j_ref][k_ref+1]   + f_ref[i_ref+1][j_ref][k_ref+1] + 
                           f_ref[i_ref][j_ref+1][k_ref+1] + f_ref[i_ref+1][j_ref+1][k_ref+1]);
                    }
        } else if (n_type == "n_less_refined") { // TO UPDATE
            //std::cout << "MIN-Y LESS REFINED" << std::endl;
            for(j=ngy/2+ny;j<ngy+ny;j++) // loop over neighbour points
                for(i=ngx+coord_x_rel*nx/2;i<ngx+(coord_x_rel+1)*nx/2;i++) 
                    for(k=ngz+coord_z_rel*nz/2;k<ngz+(coord_z_rel+1)*nz/2;k++) {
                        i_ref = 2*(i-ngx-coord_x_rel*nx/2)+ngx;
                        j_ref = 2*(j-ngy/2-ny); // DA VERIFICARE TUTTI E 3
                        k_ref = 2*(k-ngz-coord_z_rel*nz/2)+ngz;
                        f[i_ref][j_ref][k_ref]       = f_n_1[i][j][k];
                        f[i_ref+1][j_ref][k_ref]     = f_n_1[i][j][k];
                        f[i_ref][j_ref+1][k_ref]     = f_n_1[i][j][k];
                        f[i_ref+1][j_ref+1][k_ref]   = f_n_1[i][j][k];
                        f[i_ref][j_ref][k_ref+1]     = f_n_1[i][j][k];
                        f[i_ref+1][j_ref][k_ref+1]   = f_n_1[i][j][k];
                        f[i_ref][j_ref+1][k_ref+1]   = f_n_1[i][j][k];
                        f[i_ref+1][j_ref+1][k_ref+1] = f_n_1[i][j][k];
                    }
        }
    }

    if(face == "top") {
        if(n_type == "n_same_refined") {
            //std::cout << "MAX-Y SAME REFINED" << std::endl;
            for(int i=ngx;i<ngx+nx;i++) 
                for(int j=0;j<ngy;j++) 
                    for(int k=ngz;k<ngz+nz;k++) 
                        f[i][j+ny+ngy][k] = f_n_1[i][ngy+j][k];
        } else if (n_type == "n_more_refined") {
            //std::cout << "MAX-Y MORE REFINED" << std::endl;
            for(j=0;j<2*ngy;j++) // need double ngx layer to de-refine and get ngx layer
                for(i=ngx;i<ngx+nx;i++) 
                    for(k=ngz;k<ngz+nz;k++) {
                        f_ref[i][j][k]       = f_n_1[i][j+ngy][k];
                        f_ref[i+nx][j][k]    = f_n_2[i][j+ngy][k];
                        f_ref[i][j][k+nz]    = f_n_3[i][j+ngy][k];
                        f_ref[i+nx][j][k+nz] = f_n_4[i][j+ngy][k];
                    }
            for(j=0;j<ngy;j++) // de-refine
                for(i=ngx;i<ngx+nx;i++) 
                    for(k=ngz;k<ngz+nz;k++) {
                        i_ref = 2*(i-ngx)+ngx; j_ref = 2*j; k_ref = 2*(k-ngz)+ngz;
                        f[i][j+ny+ngy][k] = 0.125*(
                           f_ref[i_ref][j_ref][k_ref]     + f_ref[i_ref+1][j_ref][k_ref] + 
                           f_ref[i_ref][j_ref+1][k_ref]   + f_ref[i_ref+1][j_ref+1][k_ref] + 
                           f_ref[i_ref][j_ref][k_ref+1]   + f_ref[i_ref+1][j_ref][k_ref+1] + 
                           f_ref[i_ref][j_ref+1][k_ref+1] + f_ref[i_ref+1][j_ref+1][k_ref+1]);
                    }
        } else if (n_type == "n_less_refined") {
            //std::cout << "MAX-Y LESS REFINED" << std::endl;
            for(j=ngy;j<ngy+ngy/2;j++) // loop over neighbour points
                for(i=ngx+coord_x_rel*nx/2;i<ngx+(coord_x_rel+1)*nx/2;i++) 
                    for(k=ngz+coord_z_rel*nz/2;k<ngz+(coord_z_rel+1)*nz/2;k++) {
                        j_ref = 2*(j-ngy)+ny+ngy; // DA VERIFICARE TUTTI E 3
                        i_ref = 2*(i-ngx-coord_x_rel*nx/2)+ngx;
                        k_ref = 2*(k-ngz-coord_z_rel*nz/2)+ngz;

            //for(j=ngy/2+ny;j<ngy+ny;j++) // loop over neighbour points
            //    for(i=ngx+coord_x_rel*nx/2;i<ngx+(coord_x_rel+1)*nx/2;i++) 
            //        for(k=ngz+coord_z_rel*nz/2;k<ngz+(coord_z_rel+1)*nz/2;k++) {
            //            i_ref = 2*(i-ngx-coord_x_rel*nx/2)+ngx;
            //            j_ref = 2*(j-ngy/2-ny); // DA VERIFICARE TUTTI E 3
            //            k_ref = 2*(k-ngz-coord_z_rel*nz/2)+ngz;

                        f[i_ref][j_ref][k_ref]       = f_n_1[i][j][k];
                        f[i_ref+1][j_ref][k_ref]     = f_n_1[i][j][k];
                        f[i_ref][j_ref+1][k_ref]     = f_n_1[i][j][k];
                        f[i_ref+1][j_ref+1][k_ref]   = f_n_1[i][j][k];
                        f[i_ref][j_ref][k_ref+1]     = f_n_1[i][j][k];
                        f[i_ref+1][j_ref][k_ref+1]   = f_n_1[i][j][k];
                        f[i_ref][j_ref+1][k_ref+1]   = f_n_1[i][j][k];
                        f[i_ref+1][j_ref+1][k_ref+1] = f_n_1[i][j][k];
                    }
        }
    }

    if(face == "back") {
        if(n_type == "n_same_refined") {
            //std::cout << "MIN-Z SAME REFINED" << std::endl;
            for(int i=ngx;i<ngx+nx;i++) 
                for(int j=ngy;j<ngy+ny;j++) 
                    for(int k=0;k<ngz;k++) 
                        f[i][j][k] = f_n_1[i][j][nz+k];
        } else if (n_type == "n_more_refined") {
            //std::cout << "MIN-Z MORE REFINED" << std::endl;
            for(k=0;k<2*ngz;k++) // need double ngx layer to de-refine and get ngx layer
                for(j=ngy;j<ngy+ny;j++) 
                    for(i=ngx;i<ngx+nx;i++) {
                        f_ref[i][j][k]       = f_n_1[i][j][k+nz-ngz];
                        f_ref[i+nx][j][k]    = f_n_2[i][j][k+nz-ngz];
                        f_ref[i][j+ny][k]    = f_n_3[i][j][k+nz-ngz];
                        f_ref[i+nx][j+ny][k] = f_n_4[i][j][k+nz-ngz];
                    }
            for(k=0;k<ngz;k++) // de-refine
                for(j=ngy;j<ngy+ny;j++) 
                    for(i=ngx;i<ngx+nx;i++) {
                        i_ref = 2*(i-ngx)+ngx; j_ref = 2*(j-ngy)+ngy; k_ref = 2*k;
                        f[i][j][k] = 0.125*(
                           f_ref[i_ref][j_ref][k_ref]     + f_ref[i_ref+1][j_ref][k_ref] + 
                           f_ref[i_ref][j_ref+1][k_ref]   + f_ref[i_ref+1][j_ref+1][k_ref] + 
                           f_ref[i_ref][j_ref][k_ref+1]   + f_ref[i_ref+1][j_ref][k_ref+1] + 
                           f_ref[i_ref][j_ref+1][k_ref+1] + f_ref[i_ref+1][j_ref+1][k_ref+1]);
                    }
        } else if (n_type == "n_less_refined") {  // TO UPDATE
            //std::cout << "MIN-Z LESS REFINED" << std::endl;
            for(k=ngz/2+nz;k<ngz+nz;k++) // loop over neighbour points
                for(j=ngy+coord_y_rel*ny/2;j<ngy+(coord_y_rel+1)*ny/2;j++) 
                    for(i=ngx+coord_x_rel*nx/2;i<ngx+(coord_x_rel+1)*nx/2;i++) {
                        k_ref = 2*(k-ngz/2-nz); // DA VERIFICARE TUTTI E 3
                        j_ref = 2*(j-ngy-coord_y_rel*ny/2)+ngy;
                        i_ref = 2*(i-ngx-coord_x_rel*nx/2)+ngx;
                        f[i_ref][j_ref][k_ref]       = f_n_1[i][j][k];
                        f[i_ref+1][j_ref][k_ref]     = f_n_1[i][j][k];
                        f[i_ref][j_ref+1][k_ref]     = f_n_1[i][j][k];
                        f[i_ref+1][j_ref+1][k_ref]   = f_n_1[i][j][k];
                        f[i_ref][j_ref][k_ref+1]     = f_n_1[i][j][k];
                        f[i_ref+1][j_ref][k_ref+1]   = f_n_1[i][j][k];
                        f[i_ref][j_ref+1][k_ref+1]   = f_n_1[i][j][k];
                        f[i_ref+1][j_ref+1][k_ref+1] = f_n_1[i][j][k];
                    }
        }
    }

    if(face == "front") {
        if(n_type == "n_same_refined") {
            //std::cout << "MAX-Z SAME REFINED" << std::endl;
            for(int i=ngx;i<ngx+nx;i++) 
                for(int j=ngy;j<ngy+ny;j++) 
                    for(int k=0;k<ngz;k++) 
                        f[i][j][k+nz+ngz] = f_n_1[i][j][ngz+k];
        } else if (n_type == "n_more_refined") {
            //std::cout << "MAX-Z MORE REFINED" << std::endl;
            for(k=0;k<2*ngz;k++) // need double ngx layer to de-refine and get ngx layer
                for(j=ngy;j<ngy+ny;j++) 
                    for(i=ngx;i<ngx+nx;i++) {
                        f_ref[i][j][k]       = f_n_1[i][j][k+ngz];
                        f_ref[i+nx][j][k]    = f_n_2[i][j][k+ngz];
                        f_ref[i][j+ny][k]    = f_n_3[i][j][k+ngz];
                        f_ref[i+nx][j+ny][k] = f_n_4[i][j][k+ngz];
                    }
            for(k=0;k<ngz;k++) // de-refine
                for(j=ngy;j<ngy+ny;j++) 
                    for(i=ngx;i<ngx+nx;i++) {
                        i_ref = 2*(i-ngx)+ngx; j_ref = 2*(j-ngy)+ngy; k_ref = 2*k;
                        f[i][j][k+nz+ngz] = 0.125*(
                           f_ref[i_ref][j_ref][k_ref]     + f_ref[i_ref+1][j_ref][k_ref] + 
                           f_ref[i_ref][j_ref+1][k_ref]   + f_ref[i_ref+1][j_ref+1][k_ref] + 
                           f_ref[i_ref][j_ref][k_ref+1]   + f_ref[i_ref+1][j_ref][k_ref+1] + 
                           f_ref[i_ref][j_ref+1][k_ref+1] + f_ref[i_ref+1][j_ref+1][k_ref+1]);
                    }
        } else if (n_type == "n_less_refined") {
            //std::cout << "MAX-Z LESS REFINED" << std::endl;
            for(k=ngz;k<ngz+ngz/2;k++) // loop over neighbour points
                for(j=ngy+coord_y_rel*ny/2;j<ngy+(coord_y_rel+1)*ny/2;j++) 
                    for(i=ngx+coord_x_rel*nx/2;i<ngx+(coord_x_rel+1)*nx/2;i++) {
                        k_ref = 2*(k-ngz)+nz+ngz; // DA VERIFICARE TUTTI E 3
                        j_ref = 2*(j-ngy-coord_y_rel*ny/2)+ngy;
                        i_ref = 2*(i-ngx-coord_x_rel*nx/2)+ngx;
                        f[i_ref][j_ref][k_ref]       = f_n_1[i][j][k];
                        f[i_ref+1][j_ref][k_ref]     = f_n_1[i][j][k];
                        f[i_ref][j_ref+1][k_ref]     = f_n_1[i][j][k];
                        f[i_ref+1][j_ref+1][k_ref]   = f_n_1[i][j][k];
                        f[i_ref][j_ref][k_ref+1]     = f_n_1[i][j][k];
                        f[i_ref+1][j_ref][k_ref+1]   = f_n_1[i][j][k];
                        f[i_ref][j_ref+1][k_ref+1]   = f_n_1[i][j][k];
                        f[i_ref+1][j_ref+1][k_ref+1] = f_n_1[i][j][k];
                    }
        }
    }

}

void update_bc_maps(std::unordered_map<std::string,Geometry> &geo_map, 
        std::unordered_map<std::string,Field> &field_map, int forest[3]) {

    std::string key_n_1, key_n_2, key_n_3, key_n_4, n_type;
    int coord_x_rel, coord_y_rel, coord_z_rel;
    std::vector<std::string> faces = {"left","right","bottom","top","back","front"};

    for(auto &n: field_map) {
        std::string key = n.first;
        Field &f = n.second;
        Geometry &geo = geo_map[key];
 
        int iface=0;
        for(auto &face: faces) {
            find_neighbour(forest, key, face, geo_map, key_n_1, key_n_2, key_n_3, key_n_4, n_type, coord_x_rel, coord_y_rel, coord_z_rel);
            f.key_n_1[iface] = key_n_1; f.key_n_2[iface] = key_n_2; f.key_n_3[iface] = key_n_3; f.key_n_4[iface] = key_n_4;
            f.n_type[iface]  = n_type;
            f.coord_x_rel[iface] = coord_x_rel; f.coord_y_rel[iface] = coord_y_rel; f.coord_z_rel[iface] = coord_z_rel;
            iface++;
        }
    }

}

void update_bc(std::unordered_map<std::string,Geometry> &geo_map, 
        std::unordered_map<std::string,Field> &field_map) {

    std::string key_n_1, key_n_2, key_n_3, key_n_4, n_type;
    int coord_x_rel, coord_y_rel, coord_z_rel;
    std::vector<std::string> faces = {"left","right","bottom","top","back","front"};

    for(auto &n: field_map) {
        std::string key = n.first;
        Field &f = n.second;
        Geometry &geo = geo_map[key];
 
        int iface=0;
        for(auto &face: faces) {
            key_n_1 = f.key_n_1[iface]; key_n_2 = f.key_n_2[iface]; key_n_3 = f.key_n_3[iface]; key_n_4 = f.key_n_4[iface]; 
            n_type = f.n_type[iface];
            coord_x_rel = f.coord_x_rel[iface]; coord_y_rel = f.coord_y_rel[iface]; coord_z_rel = f.coord_z_rel[iface];
            //std::cout << "face: " << face << "n_type: " << n_type << std::endl;
            if(n_type == "phys") {
                //std::cout << "UPDATE_BC: " << "face: " << face << " - n_type: " << n_type << std::endl;
                bc_phys(f, face);
            } else {
                Field &f_n_1 = field_map[key_n_1]; Field &f_n_2 = field_map[key_n_2]; 
                Field &f_n_3 = field_map[key_n_3]; Field &f_n_4 = field_map[key_n_4]; 
                update_bc_block_field(face, n_type, f, f_n_1, f_n_2, f_n_3, f_n_4, geo, coord_x_rel, coord_y_rel, coord_z_rel);
                //update_bc_block(face, n_type, f.u, f.work_ref, f_n_1.u, f_n_2.u, f_n_3.u, f_n_4.u,
                //    nx, ny, nz, ngx, ngy, ngz, coord_x_rel, coord_y_rel, coord_z_rel);
            }
            iface++;
        }
    }

}
