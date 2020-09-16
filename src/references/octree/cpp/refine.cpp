// ----------------------------------------------------------------------------------
// Lagrangian 3D interpolator from fine to coarse and from coarse to fine grids
// Francesco Salvadore, 2018
// ----------------------------------------------------------------------------------
// Lagrangian interpolation direction-by-direction
// f(x,1,1) = f[1,1,1]*(x-x2)/(x1-x2)*(x-x3)/(x1-x3)+
//            f[2,1,1]*(x-x1)/(x2-x1)*(x-x3)/(x2-x3)+
//            f[3,1,1]*(x-x1)/(x3-x1)*(x-x2)/(x3-x2);
//
// Derefine (1,2,3) => 3/2
// f_3_2_1_1 = 3./8.*f[1,1,1]+3./4.*f[2,1,1]-1./8.*f[3,1,1]
//
// Refine (1,2,3) => (7/4,9/4)
// f_7_4_1_1 = 5./32.*f[1,1,1]+15./16.*f[2,1,1]-3./32.*f[3,1,1]
// f_9_4_1_1 = -3./32.*f[1,1,1]+15./16.*f[2,1,1]+5./32.*f[3,1,1]
//
// From fine to coarse:
// i=1:n/2 ==> i_ref = 2*i-1 interpolate from 
//     [i_ref:i_ref+2]x[j_ref:j_ref+2]x[k_ref:k_ref+2] to 1 value i
//
// From coarse to fine:
// i=1:n/2 ==> i_ref = 2*i-1 interpolate from [i-1:i+1, j-1:j+1, k-1:k+1] to i
//     8 values [i_ref:i_ref+1,j_ref:j_ref+1,k_ref:k_ref+1]
// ----------------------------------------------------------------------------------

#include <iostream>
#include <array>
#include <unordered_map>
#include <vector>
#include <set>
#include "geometry.h"
#include "field.h"
#include "refine.h"
#include "utils.h"
#include "bc.h"
#include "immersed.h"

void derefine_blocks(std::vector<std::string> key_str_list, std::unordered_map<std::string,Geometry> &geo_map, 
        std::unordered_map<std::string,Field> &field_map,
        std::unordered_map<std::string,Polytree> & immersed_map) {

    std::array<std::array<int,7>,8> keys;
    int ik=0;
    for(auto &key_str: key_str_list) {
        KeyStringToArray(key_str, keys[ik]); ik++;
    }
    Geometry geo_first_child = geo_map[key_str_list[0]];
    Field field_child = field_map[key_str_list[0]];

    // Geo
    std::array<int,7> new_key;
    int level = new_key[0] = keys[0][0] - 1;
    Geometry s_temp;
    new_key[1] = keys[0][1]; new_key[2] = keys[0][2]; new_key[3] = keys[0][3]; new_key[4] = keys[0][4]; new_key[5] = keys[0][5]; new_key[6] = keys[0][6];
    s_temp.coords = {new_key[0],new_key[1],new_key[2],new_key[3],new_key[4],new_key[5],new_key[6]};
    std::string new_key_str; KeyArrayToString(new_key, new_key_str);
    for(int c=0 ;c<3;c++){
        int npoints  = geo_first_child.intervals[c][2];
        int nghosts  = geo_first_child.intervals[c][3];
        s_temp.intervals[c][0] = 0;
        s_temp.intervals[c][1] = npoints-1;
        s_temp.intervals[c][2] = npoints;
        s_temp.intervals[c][3] = nghosts;
    }
    float le_x = geo_first_child.lengths[0][2] * 2.;
    s_temp.lengths[0][0] = geo_first_child.lengths[0][0];
    s_temp.lengths[0][1] = s_temp.lengths[0][0] + le_x;
    s_temp.lengths[0][2] = le_x;
    s_temp.lengths[0][3] = geo_first_child.lengths[0][3] * 2.;
    float le_y = geo_first_child.lengths[1][2] * 2.;
    s_temp.lengths[1][0] = geo_first_child.lengths[1][0];
    s_temp.lengths[1][1] = s_temp.lengths[1][0] + le_y;
    s_temp.lengths[1][2] = le_y;
    s_temp.lengths[1][3] = geo_first_child.lengths[1][3] * 2.;
    float le_z = geo_first_child.lengths[2][2] * 2.;
    s_temp.lengths[2][0] = geo_first_child.lengths[2][0];
    s_temp.lengths[2][1] = s_temp.lengths[2][0] + le_z;
    s_temp.lengths[2][2] = le_z;
    s_temp.lengths[2][3] = geo_first_child.lengths[2][3] * 2.;
    s_temp.key = new_key_str;
    geo_map[new_key_str] = s_temp;

    // Field
    Field f_temp;
    copy_geo(s_temp, f_temp.geo);
    f_temp.refinement_needed = 0;
    std::vector<int> np = {s_temp.intervals[0][2], s_temp.intervals[1][2], s_temp.intervals[2][2]};
    std::vector<int> ng = {s_temp.intervals[0][3], s_temp.intervals[1][3], s_temp.intervals[2][3]};
    allocate_field(f_temp, np, ng);
    f_temp.key = new_key_str;
    allocate_distances_block(immersed_map, f_temp);
    compute_distances_block(immersed_map, f_temp);
    field_map[new_key_str] = f_temp;

    // Interpolate data from parent field
    //derefine_and_merge(geo_first_child, field_map[new_key_str].u, field_map[new_key_str].work_ref,
    //        field_map[key_str_list[0]].u, field_map[key_str_list[1]].u, 
    //        field_map[key_str_list[2]].u, field_map[key_str_list[3]].u, 
    //        field_map[key_str_list[4]].u, field_map[key_str_list[5]].u, 
    //        field_map[key_str_list[6]].u, field_map[key_str_list[7]].u);

    derefine_and_merge_field(field_map, geo_first_child, new_key_str, key_str_list);

    // Deallocate children fields and delete their keys from maps
    for(auto &key_str: key_str_list) {
        std::cout << "derefine removing key: " << key_str << std::endl;
        Field &field = field_map[key_str];
        deallocate_field(field, np, ng);
        deallocate_distances_block(field);
        geo_map.erase(key_str);
        field_map.erase(key_str);
    }
}

// Takes 8 blocks and return the derefined block
void derefine_and_merge(Geometry & geo, float***f, float***f_ref_glob,
        float *** f_ref_1, float *** f_ref_2, float *** f_ref_3, float *** f_ref_4, 
        float *** f_ref_5, float *** f_ref_6, float *** f_ref_7, float *** f_ref_8) {

    int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
    int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];

    for(int i=0; i<nx+ngx; i++) for(int j=0; j<ny+ngy; j++) for(int k=0; k<nz+ngz; k++) {
        f_ref_glob[i][j][k]                      = f_ref_1[i][j][k];
        f_ref_glob[i+nx+ngx][j][k]               = f_ref_2[i+ngx][j][k];
        f_ref_glob[i][j+ny+ngy][k]               = f_ref_3[i][j+ngy][k];
        f_ref_glob[i+nx+ngx][j+ny+ngy][k]        = f_ref_4[i+ngx][j+ngy][k];
        f_ref_glob[i][j][k+nz+ngz]               = f_ref_5[i][j][k+ngz];
        f_ref_glob[i+nx+ngx][j][k+nz+ngz]        = f_ref_6[i+ngx][j][k+ngz];
        f_ref_glob[i][j+ny+ngy][k+nz+ngz]        = f_ref_7[i][j+ngy][k+ngz];
        f_ref_glob[i+nx+ngx][j+ny+ngy][k+nz+ngz] = f_ref_8[i+ngx][j+ngy][k+ngz];
    }

    derefine(f_ref_glob, f, geo);

}

void derefine(float *** f, float *** f_deref, Geometry & geo) {

    int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
    int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];
    int i_ref, j_ref, k_ref;

    float *** f_temp = allocate3d({3,3,3});
    // loop over coarse grid
    for(int i=ngx; i<nx+ngx; i++)
        for(int j=ngy; j<ny+ngy; j++)
            for(int k=ngz; k<nz+ngz; k++) {
// get the corresponding fine grid point "below"
                i_ref = (i-ngx)*2+ngx; j_ref = (j-ngy)*2+ngy; k_ref = (k-ngz)*2+ngz; 
// derefine from 3*3*3 fine points (i_ref:i_ref+2 could be changed to i_ref-1:i_ref+1)
                for(int ii=0;ii<3;ii++)
                    for(int jj=0;jj<3;jj++)
                        for(int kk=0;kk<3;kk++) {
                            f_temp[ii][jj][kk] = f[i_ref+ii][j_ref+jj][k_ref+kk];
                        }
                f_deref[i][j][k] = derefine_point(f_temp);
            }
    deallocate3d(&f_temp, {3,3,3});

}

float derefine_point(float *** f) {

    //RIMETTEREfloat f_3_2_1_1     = 3./8.*f[0][0][0]+3./4.*f[1][0][0]-1./8.*f[2][0][0];
    //RIMETTEREfloat f_3_2_2_1     = 3./8.*f[0][1][0]+3./4.*f[1][1][0]-1./8.*f[2][1][0];
    //RIMETTEREfloat f_3_2_3_1     = 3./8.*f[0][2][0]+3./4.*f[1][2][0]-1./8.*f[2][2][0];
    //RIMETTEREfloat f_3_2_1_2     = 3./8.*f[0][0][1]+3./4.*f[1][0][1]-1./8.*f[2][0][1];
    //RIMETTEREfloat f_3_2_2_2     = 3./8.*f[0][1][1]+3./4.*f[1][1][1]-1./8.*f[2][1][1];
    //RIMETTEREfloat f_3_2_3_2     = 3./8.*f[0][2][1]+3./4.*f[1][2][1]-1./8.*f[2][2][1];
    //RIMETTEREfloat f_3_2_1_3     = 3./8.*f[0][0][2]+3./4.*f[1][0][2]-1./8.*f[2][0][2];
    //RIMETTEREfloat f_3_2_2_3     = 3./8.*f[0][1][2]+3./4.*f[1][1][2]-1./8.*f[2][1][2];
    //RIMETTEREfloat f_3_2_3_3     = 3./8.*f[0][2][2]+3./4.*f[1][2][2]-1./8.*f[2][2][2];

    //RIMETTEREfloat f_3_2_3_2_1   = 3./8.*f_3_2_1_1+3./4.*f_3_2_2_1-1./8.*f_3_2_3_1;
    //RIMETTEREfloat f_3_2_3_2_2   = 3./8.*f_3_2_1_2+3./4.*f_3_2_2_2-1./8.*f_3_2_3_2;
    //RIMETTEREfloat f_3_2_3_2_3   = 3./8.*f_3_2_1_3+3./4.*f_3_2_2_3-1./8.*f_3_2_3_3;

    //RIMETTEREfloat f_3_2_3_2_3_2 = 3./8.*f_3_2_3_2_1+3./4.*f_3_2_3_2_2-1./8.*f_3_2_3_2_3;

    //RIMETTEREreturn f_3_2_3_2_3_2;
    
    float f_interp = 0.125*(f[0][0][0]+f[1][0][0]+f[0][1][0]+f[1][1][0]+
                            f[0][0][1]+f[1][0][1]+f[0][1][1]+f[1][1][1]);
    return f_interp;
}

void refine_block(std::string key_str, std::unordered_map<std::string,Geometry> &geo_map, 
        std::unordered_map<std::string,Field> &field_map,
        std::unordered_map<std::string,Polytree> & immersed_map) {

    std::array<int,7> key;
    KeyStringToArray(key_str, key);
    Geometry geo_parent = geo_map[key_str];
    Field field_parent = field_map[key_str];
    std::vector<std::string> key_str_list;

    // Create 2x2x2 descendants from parent blocks (geo and field)
    for(int k=0 ;k<2;k++)
        for(int j=0 ;j<2;j++)
            for(int i=0 ;i<2;i++)
            {
                // Geo
                std::array<int,7> new_key;
                int level = new_key[0] = key[0] + 1;
                // warning! 2 << -1 is not accepted, and gives wrong results with intel compiler
                if(level < max_levels) {  
                    new_key[1] = key[1] + i * (2 << (max_levels-level-1));
                    new_key[2] = key[2] + j * (2 << (max_levels-level-1));
                    new_key[3] = key[3] + k * (2 << (max_levels-level-1));
                } else {
                    new_key[1] = key[1] + i;
                    new_key[2] = key[2] + j;
                    new_key[3] = key[3] + k;
                }
                new_key[4] = key[4];
                new_key[5] = key[5];
                new_key[6] = key[6];
                Geometry s_temp;
                s_temp.coords = {new_key[0],new_key[1],new_key[2],new_key[3],new_key[4],new_key[5],new_key[6]};
                std::string new_key_str;
                KeyArrayToString(new_key, new_key_str);

                key_str_list.push_back(new_key_str);

                s_temp.key = new_key_str;
                int levcoords[3] = {i,j,k};
                int start[3];
                int end[3];
                int size[3];
                for(int c=0 ;c<3;c++){
                    int npoints  = geo_parent.intervals[c][2];
                    int nghosts  = geo_parent.intervals[c][3];

                    s_temp.intervals[c][0] = 0;
                    s_temp.intervals[c][1] = npoints-1;
                    s_temp.intervals[c][2] = npoints;
                    s_temp.intervals[c][3] = nghosts;
                }
                float le_x = geo_parent.lengths[0][2] / 2.;
                s_temp.lengths[0][0] = geo_parent.lengths[0][0] + i*le_x;
                s_temp.lengths[0][1] = s_temp.lengths[0][0] + le_x;
                s_temp.lengths[0][2] = le_x;
                s_temp.lengths[0][3] = geo_parent.lengths[0][3] / 2.;
                float le_y = geo_parent.lengths[1][2] / 2.;
                s_temp.lengths[1][0] = geo_parent.lengths[1][0] + j*le_y;
                s_temp.lengths[1][1] = s_temp.lengths[1][0] + le_y;
                s_temp.lengths[1][2] = le_y;
                s_temp.lengths[1][3] = geo_parent.lengths[1][3] / 2.;
                float le_z = geo_parent.lengths[2][2] / 2.;
                s_temp.lengths[2][0] = geo_parent.lengths[2][0] + k*le_z;
                s_temp.lengths[2][1] = s_temp.lengths[2][0] + le_z;
                s_temp.lengths[2][2] = le_z;
                s_temp.lengths[2][3] = geo_parent.lengths[2][3] / 2.;

                s_temp.key = new_key_str;
                geo_map[new_key_str] = s_temp;

                std::cout << "adding refined key: " << new_key_str << std::endl;

                // Field
                Field f_temp;
                copy_geo(s_temp, f_temp.geo);
                f_temp.refinement_needed = 0;
                std::vector<int> np = {s_temp.intervals[0][2], s_temp.intervals[1][2], s_temp.intervals[2][2]};
                std::vector<int> ng = {s_temp.intervals[0][3], s_temp.intervals[1][3], s_temp.intervals[2][3]};
                allocate_field(f_temp, np, ng);
                f_temp.key = new_key_str;
                allocate_distances_block(immersed_map, f_temp);
                compute_distances_block(immersed_map, f_temp);
                field_map[new_key_str] = f_temp;
            }

    // Interpolate data from parent field
    refine_and_split_field(field_map, geo_parent, key_str, key_str_list);
    //refine_and_split(geo_parent, field_parent.u, field_parent.work_ref,
    //        field_map[key_str_list[0]].u, field_map[key_str_list[1]].u, 
    //        field_map[key_str_list[2]].u, field_map[key_str_list[3]].u, 
    //        field_map[key_str_list[4]].u, field_map[key_str_list[5]].u, 
    //        field_map[key_str_list[6]].u, field_map[key_str_list[7]].u);

    // Deallocate parent fields and delete parent key from maps
    std::vector<int> np = {geo_parent.intervals[0][2], geo_parent.intervals[1][2], geo_parent.intervals[2][2]};
    std::vector<int> ng = {geo_parent.intervals[0][3], geo_parent.intervals[1][3], geo_parent.intervals[2][3]};
    deallocate_field(field_parent, np, ng);
    deallocate_distances_block(field_parent);
    geo_map.erase(key_str);
    field_map.erase(key_str);
}

// Takes one block and returns 8 blocks result of refinement and split
void refine_and_split(Geometry geo, float***f, float***f_ref_glob,
        float *** f_ref_1, float *** f_ref_2, float *** f_ref_3, float *** f_ref_4, 
        float *** f_ref_5, float *** f_ref_6, float *** f_ref_7, float *** f_ref_8) {

    int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
    int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];

    refine(f, f_ref_glob, geo);

    for(int i=0; i<nx+2*ngx; i++)
        for(int j=0; j<ny+2*ngy; j++)
            for(int k=0; k<nz+2*ngz; k++) {
                f_ref_1[i][j][k] = f_ref_glob[i][j][k];
                f_ref_2[i][j][k] = f_ref_glob[i+nx][j][k];
                f_ref_3[i][j][k] = f_ref_glob[i][j+ny][k];
                f_ref_4[i][j][k] = f_ref_glob[i+nx][j+ny][k];
                f_ref_5[i][j][k] = f_ref_glob[i][j][k+nz];
                f_ref_6[i][j][k] = f_ref_glob[i+nx][j][k+nz];
                f_ref_7[i][j][k] = f_ref_glob[i][j+ny][k+nz];
                f_ref_8[i][j][k] = f_ref_glob[i+nx][j+ny][k+nz];
            }
}

void refine(float ***f, float ***f_ref, Geometry geo) {

    int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
    int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];
    float f_temp[3][3][3]; //, f_temp_2[2][2][2]; // ;
    float *** f_temp_2 = allocate3d({2,2,2});
    // loop over coarse grid
    for(int ia=0 ;ia<nx;ia++)
        for(int ja=0 ;ja<ny;ja++)
            for(int ka=0 ;ka<nz;ka++) {
                int i=ia+ngx;      int j=ja+ngy;      int k=ka+ngz;
                // get the corresponding fine grid point "below"
                int i_ref=2*i-ngx; int j_ref=2*j-ngy; int k_ref=2*k-ngz; 
                // refine from 3*3*3 coarse points to 2*2*2 fine points surrounding the coarse central point
                // f_ref[i_ref:i_ref+1][j_ref:j_ref+1][k_ref:k_ref+1] = refine_points(f[i-1:i+1][j-1:j+1][k-1:k+1])
                for(int ii=0;ii<3;ii++)
                    for(int jj=0;jj<3;jj++)
                        for(int kk=0;kk<3;kk++)
                            f_temp[ii][jj][kk] = f[i-1+ii][j-1+jj][k-1+kk];

                f_temp_2 = refine_points(f_temp);

                for(int ii=0;ii<2;ii++)
                    for(int jj=0;jj<2;jj++)
                        for(int kk=0;kk<2;kk++)
                            f_ref[i_ref+ii][j_ref+jj][k_ref+kk] = f_temp_2[ii][jj][kk];
            }
    deallocate3d(&f_temp_2, {2,2,2});
}

float*** refine_points(float f[3][3][3]) {

    float *** ref_points = allocate3d({2,2,2});
    for(int i=0 ;i<2;i++)
        for(int j=0 ;j<2;j++)
            for(int k=0 ;k<2;k++)
                ref_points[i][j][k] = 0.F;

    for(int ii=0 ;ii<2;ii++)
        for(int jj=0 ;jj<2;jj++)
            for(int kk=0 ;kk<2;kk++) {
                //RIMETTEREfloat coeff_i_1 = 5./32.*(1-ii)+(-3./32.)*ii;  // might be optimized
                //RIMETTEREfloat coeff_i_2 = 15./16.;
                //RIMETTEREfloat coeff_i_3 = 1. - coeff_i_1 - coeff_i_2;
                //RIMETTEREfloat coeff_j_1 = 5./32.*(1-jj)+(-3./32.)*jj;  // might be optimized
                //RIMETTEREfloat coeff_j_2 = 15./16.;
                //RIMETTEREfloat coeff_j_3 = 1. - coeff_j_1 - coeff_j_2;
                //RIMETTEREfloat coeff_k_1 = 5./32.*(1-kk)+(-3./32.)*kk;  // might be optimized
                //RIMETTEREfloat coeff_k_2 = 15./16.;
                //RIMETTEREfloat coeff_k_3 = 1. - coeff_k_1 - coeff_k_2;

                //RIMETTEREfloat f_p_1_1 = coeff_i_1*f[0][0][0]+coeff_i_2*f[1][0][0]+coeff_i_3*f[2][0][0];
                //RIMETTEREfloat f_p_2_1 = coeff_i_1*f[0][1][0]+coeff_i_2*f[1][1][0]+coeff_i_3*f[2][1][0];
                //RIMETTEREfloat f_p_3_1 = coeff_i_1*f[0][2][0]+coeff_i_2*f[1][2][0]+coeff_i_3*f[2][2][0];
                //RIMETTEREfloat f_p_1_2 = coeff_i_1*f[0][0][1]+coeff_i_2*f[1][0][1]+coeff_i_3*f[2][0][1];
                //RIMETTEREfloat f_p_2_2 = coeff_i_1*f[0][1][1]+coeff_i_2*f[1][1][1]+coeff_i_3*f[2][1][1];
                //RIMETTEREfloat f_p_3_2 = coeff_i_1*f[0][2][1]+coeff_i_2*f[1][2][1]+coeff_i_3*f[2][2][1];
                //RIMETTEREfloat f_p_1_3 = coeff_i_1*f[0][0][2]+coeff_i_2*f[1][0][2]+coeff_i_3*f[2][0][2];
                //RIMETTEREfloat f_p_2_3 = coeff_i_1*f[0][1][2]+coeff_i_2*f[1][1][2]+coeff_i_3*f[2][1][2];
                //RIMETTEREfloat f_p_3_3 = coeff_i_1*f[0][2][2]+coeff_i_2*f[1][2][2]+coeff_i_3*f[2][2][2];

                //RIMETTEREfloat f_p_p_1 = coeff_j_1*f_p_1_1+coeff_j_2*f_p_2_1+coeff_j_3*f_p_3_1;
                //RIMETTEREfloat f_p_p_2 = coeff_j_1*f_p_1_2+coeff_j_2*f_p_2_2+coeff_j_3*f_p_3_2;
                //RIMETTEREfloat f_p_p_3 = coeff_j_1*f_p_1_3+coeff_j_2*f_p_2_3+coeff_j_3*f_p_3_3;

                //RIMETTEREfloat f_p_p_p = coeff_k_1*f_p_p_1+coeff_k_2*f_p_p_2+coeff_k_3*f_p_p_3;

                //RIMETTEREref_points[ii][jj][kk] = f_p_p_p;

                ref_points[ii][jj][kk] = f[1][1][1];
            }

    return (float***) ref_points;
}

std::vector<std::vector<std::string>> sanitize_refinement_needed(std::unordered_map<std::string,Geometry> &geo_map, int forest[3]) {

    std::vector<std::vector<std::string>> derefine_octuplets;

    const int max_sanitize_iterations = 10;
    bool sanitize_completed = false;
    std::string key, key_s, key_n_1, key_n_2, key_n_3, key_n_4, n_type;
    std::array<int,7> key_arr, key_arr_s, key_arr_n_1, key_arr_n_2, key_arr_n_3, key_arr_n_4;
    int coord_temp_a, coord_temp_b, coord_temp_c;
    int level, size_xyz, coord_x, coord_y, coord_z, coord_x_s, coord_y_s, coord_z_s;
    int coord_x_rel, coord_y_rel, coord_z_rel;
    Geometry geo, geo_n_1, geo_n_2, geo_n_3, geo_n_4;
    int new_level, new_level_n_1, new_level_n_2, new_level_n_3, new_level_n_4;
    std::vector<std::string> faces = {"left","right","bottom","top","back","front"};

    // First avoid new level being > max_levels or < 0
    for(auto &n: geo_map) {
        key = n.first; KeyStringToArray(key, key_arr);
        geo = n.second;
        new_level = key_arr[0] + geo.refinement_needed;
        if(new_level > max_levels || new_level < 0) {
            std::cout << "preventing excessive refinement: " << key << " " << geo.refinement_needed << std::endl;
            geo_map[key].refinement_needed = 0;
        }
    }

    // Then iterate checking that:
    // (a) all blocks to de-refine have all the siblings to de-refine too
    // (b) neighbouring blocks have level difference by 1 or 0 or -1
    for(int i_san_iter=0;i_san_iter<max_sanitize_iterations;i_san_iter++) {
        std::cout << "sanitizing iteration: " << i_san_iter << std::endl;
        sanitize_completed = true;

        std::vector<std::string> derefine_analyzed_keys; // replace with set instead of vector
        derefine_octuplets = {};
        for(auto &n: geo_map) {
            key = n.first; KeyStringToArray(key, key_arr);
            geo = n.second;

            // (a)
            std::vector<std::string> sibling_keys;
            if( (geo.refinement_needed == -1) 
                && (std::find(derefine_analyzed_keys.begin(), derefine_analyzed_keys.end(), key) == derefine_analyzed_keys.end()) )
                {
                std::cout << "possible derefine for key= " << key << std::endl;
                level = key_arr[0]; coord_x = key_arr[1]; coord_y = key_arr[2]; coord_z = key_arr[3];
                if(level < max_levels) { // cannot shift with negative values
                    size_xyz = (2 << (max_levels-level-1)); 
                } else {
                    size_xyz = 1;
                }
                coord_x_rel = coord_x%(2*size_xyz); coord_y_rel = coord_y%(2*size_xyz); coord_z_rel = coord_z%(2*size_xyz);
                if(coord_x_rel == 0) { coord_x_s= coord_x+size_xyz; } else { coord_x_s= coord_x-size_xyz; }
                if(coord_y_rel == 0) { coord_y_s= coord_y+size_xyz; } else { coord_y_s= coord_y-size_xyz; }
                if(coord_z_rel == 0) { coord_z_s= coord_z+size_xyz; } else { coord_z_s= coord_z-size_xyz; }
                key_arr_s = {level, std::min(coord_x,coord_x_s), std::min(coord_y,coord_y_s), std::min(coord_z,coord_z_s)};
                KeyArrayToString(key_arr_s, key_s); sibling_keys.push_back(key_s); 
                key_arr_s = {level, std::max(coord_x,coord_x_s), std::min(coord_y,coord_y_s), std::min(coord_z,coord_z_s)};
                KeyArrayToString(key_arr_s, key_s); sibling_keys.push_back(key_s); 
                key_arr_s = {level, std::min(coord_x,coord_x_s), std::max(coord_y,coord_y_s), std::min(coord_z,coord_z_s)};
                KeyArrayToString(key_arr_s, key_s); sibling_keys.push_back(key_s); 
                key_arr_s = {level, std::max(coord_x,coord_x_s), std::max(coord_y,coord_y_s), std::min(coord_z,coord_z_s)};
                KeyArrayToString(key_arr_s, key_s); sibling_keys.push_back(key_s); 
                key_arr_s = {level, std::min(coord_x,coord_x_s), std::min(coord_y,coord_y_s), std::max(coord_z,coord_z_s)};
                KeyArrayToString(key_arr_s, key_s); sibling_keys.push_back(key_s); 
                key_arr_s = {level, std::max(coord_x,coord_x_s), std::min(coord_y,coord_y_s), std::max(coord_z,coord_z_s)};
                KeyArrayToString(key_arr_s, key_s); sibling_keys.push_back(key_s); 
                key_arr_s = {level, std::min(coord_x,coord_x_s), std::max(coord_y,coord_y_s), std::max(coord_z,coord_z_s)};
                KeyArrayToString(key_arr_s, key_s); sibling_keys.push_back(key_s); 
                key_arr_s = {level, std::max(coord_x,coord_x_s), std::max(coord_y,coord_y_s), std::max(coord_z,coord_z_s)};
                KeyArrayToString(key_arr_s, key_s); sibling_keys.push_back(key_s); 
                bool cancel_derefine = false;
                for(auto &sk: sibling_keys) {
                    auto t = geo_map.find(sk);
                    if (t == geo_map.end()) { // sk not found
                        cancel_derefine = true;
                        break;
                    } else if (geo_map[sk].refinement_needed != -1) {
                        cancel_derefine = true;
                        break;
                    }
                }
                if(cancel_derefine) {
                    for(auto &sk: sibling_keys) {
                        auto t = geo_map.find(sk);
                        if (t != geo_map.end()) { // sk found
                            if(geo_map[sk].refinement_needed == -1) {
                                std::cout << "canceling derefine for key= " << key << std::endl;
                                sanitize_completed = false;
                                geo_map[sk].refinement_needed = 0;
                            }
                        }
                    }
                } else {  // derefine active (up to now)
                    derefine_octuplets.push_back(sibling_keys);
                    for(auto sk: sibling_keys) {
                        std::cout << "confirming derefine for key= " << sk << std::endl;
                        derefine_analyzed_keys.push_back(sk);
                    }
                }
            }
        }

        for(auto &n: geo_map) {
            key = n.first; KeyStringToArray(key, key_arr);
            geo = n.second;
            // (b)
            new_level = key_arr[0] + geo.refinement_needed;
            for(auto &face: faces) {
                //std::cout << "fast check start: " << std::endl;
                find_neighbour(forest, key, face, geo_map, key_n_1, key_n_2, key_n_3, key_n_4, n_type, coord_temp_a, coord_temp_b, coord_temp_c);
                //std::cout << "fast check end: " << std::endl;
                if(n_type != "phys") {
                    if(n_type == "n_more_refined") {
                        geo_n_1 = geo_map[key_n_1]; KeyStringToArray(key_n_1, key_arr_n_1);
                        new_level_n_1 = key_arr_n_1[0] + geo_n_1.refinement_needed;
                        geo_n_2 = geo_map[key_n_2]; KeyStringToArray(key_n_2, key_arr_n_2);
                        new_level_n_2 = key_arr_n_2[0] + geo_n_2.refinement_needed;
                        geo_n_3 = geo_map[key_n_3]; KeyStringToArray(key_n_3, key_arr_n_3);
                        new_level_n_3 = key_arr_n_3[0] + geo_n_3.refinement_needed;
                        geo_n_4 = geo_map[key_n_4]; KeyStringToArray(key_n_4, key_arr_n_4);
                        new_level_n_4 = key_arr_n_4[0] + geo_n_4.refinement_needed;
                        if((new_level_n_1 > new_level+1) || (new_level_n_2 > new_level+1) || 
                                (new_level_n_3 > new_level+1) || (new_level_n_4 > new_level+1)) {
                            std::cout << "sanitizing key: " << key << " due to key (or siblings): " << key_n_1 << " " << key_n_2 << " " << key_n_3 << " " << key_n_4 << std::endl;
                            sanitize_completed = false;
                            geo_map[key].refinement_needed = std::min(geo_map[key].refinement_needed+1,1);
                            //geo_map[key].refinement_needed += 1;
                        }
                    } else {
                        geo_n_1 = geo_map[key_n_1]; KeyStringToArray(key_n_1, key_arr_n_1);
                        new_level_n_1 = key_arr_n_1[0] + geo_n_1.refinement_needed;
                        if(new_level_n_1 > new_level+1) {
                            std::cout << "sanitizing key: " << key << " due to key: " << key_n_1 << std::endl;
                            sanitize_completed = false;
                            geo_map[key].refinement_needed = std::min(geo_map[key].refinement_needed+1,1);
                            //geo_map[key].refinement_needed += 1;
                        }
                    }
                }
            }
            if(geo_map[key].refinement_needed > 1) {
                std::cout << "CANNOT REFINE TWICE IN A ROW. SOMETHING WENT TERRIBLY WRONG. EXITING!" << std::endl; exit(EXIT_FAILURE);
            }
            if(key_arr[0] + geo.refinement_needed > max_levels) {
                std::cout << "CANNOT REFINE MORE. SOMETHING WENT TERRIBLY WRONG" << std::endl; exit(EXIT_FAILURE);
            }
        }
        if(sanitize_completed) {
            std::cout << "Sanitize completed at i_san_iter: " << i_san_iter << std::endl;
            break;
        }
    }

    return derefine_octuplets;

}

void update_refine(std::unordered_map<std::string,Geometry> &geo_map, 
        std::unordered_map<std::string,Field> &field_map,
        std::unordered_map<std::string,Polytree> &immersed_map, std::string const & mode, int forest[3]) {

    for(auto &n: field_map) {
        if(mode=="mark") {
            mark_refinement_needed(n.second);
            std::cout << "key: " << n.first << " - ref_needed 0/1: " << n.second.refinement_needed << std::endl;
        } else if(mode=="force") {
            n.second.refinement_needed = 1;
        }
    }

    for(auto &n: field_map) {
        std::string k = n.first;
        Field & f = n.second;
        geo_map[k].refinement_needed = f.refinement_needed;
    }

    std::vector<std::vector<std::string>> derefine_octuplets = sanitize_refinement_needed(geo_map, forest);

    std::vector<std::string> keys;
    int ii=-1;
    for(auto &n: field_map) {
        ii++;
        if(geo_map[n.first].refinement_needed == 1) {
            //if(ii%4 == 0) {
            keys.push_back(n.first);
            std::cout << "Adding key to refine: " << n.first << std::endl;
        }
    }
    for(auto &k:keys)
        refine_block(k, geo_map, field_map, immersed_map);

     //fix_edge_and_corners(field_map); // METTERLO ANCHE ALTROVE FORSE

    for(auto &o: derefine_octuplets) {
        std::cout << "derefine octuplet: " << o[0] << " " << o[1] << " " << o[2] << " " << o[3] << 
                     " " << o[4] << " " << o[5] << " " << o[6] << " " << o[7] << std::endl;
        derefine_blocks(o, geo_map, field_map, immersed_map);
    }
}
