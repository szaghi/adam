#include <limits>
#include <iostream>
#include <cmath>        // std::abs
#include "geometry.h"
#include "field.h"
//NONRIMETTERESPERO#include "field_cuda.h"
#include "refine.h"
#include "bc.h"
#include "der_inter.h"

void bc_phys(Field &f, std::string face){
    Geometry geo = f.geo;
    int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
    int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];
    if(face == "left") {
#pragma omp parallel for
        for(int i=0;i<ngx;i++) 
            for(int j=0;j<2*ngy+ny;j++) 
#pragma omp simd
                for(int k=0;k<2*ngz+nz;k++) {
                    f.tem[i][j][k] = 1.;
                }
    } else if(face == "right") {
#pragma omp parallel for
        for(int i=ngx+nx;i<2*ngx+nx;i++) 
            for(int j=0;j<2*ngy+ny;j++) 
#pragma omp simd
                for(int k=0;k<2*ngz+nz;k++) {
                    f.tem[i][j][k] = f.tem[ngx+nx-1][j][k]; 
                }
    } else if(face == "bottom") {
#pragma omp parallel for
        for(int i=0;i<2*ngx+nx;i++) 
            for(int j=0;j<ngy;j++) 
#pragma omp simd
                for(int k=0;k<2*ngz+nz;k++) {
                    f.tem[i][j][k] = f.tem[i][ngy][k];
                }
    } else if(face == "top") {
#pragma omp parallel for
        for(int i=0;i<2*ngx+nx;i++) 
            for(int j=ngy+ny;j<2*ngy+ny;j++) 
#pragma omp simd
                for(int k=0;k<2*ngz+nz;k++) {
                    f.tem[i][j][k] = f.tem[i][ngy+ny-1][k];
                }
    } else if(face == "back") {
#pragma omp parallel for
        for(int i=0;i<2*ngx+nx;i++) 
            for(int j=0;j<2*ngy+ny;j++) 
#pragma omp simd
                for(int k=0;k<ngz;k++) {
                    f.tem[i][j][k] = f.tem[i][j][ngz]; 
                }
    } else if(face == "front") {
#pragma omp parallel for
        for(int i=0;i<2*ngx+nx;i++) 
            for(int j=0;j<2*ngy+ny;j++) 
#pragma omp simd
                for(int k=ngz+nz;k<2*ngz+nz;k++) {
                    f.tem[i][j][k] = f.tem[i][j][ngz+nz-1];
                }
    }
}

void fix_edge_and_corners(std::unordered_map<std::string,Field> &field_map) {
    for(auto &n: field_map) {
        std::string key = n.first;
        Field &field = n.second;

        Geometry geo = field.geo;
        int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
        int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];

        fix_edge_and_corners_block(field.tem, nx, ny, nz, ngx, ngy, ngz);
    }
}

void update_bc_block_field(std::string & face, std::string& n_type, Field & f, 
        Field & f_n_1, Field & f_n_2, Field & f_n_3, Field & f_n_4, 
        Geometry & geo, int coord_x_rel, int coord_y_rel, int coord_z_rel) {

    int nx   = geo.intervals[0][2]; int ny   = geo.intervals[1][2]; int nz   = geo.intervals[2][2];
    int ngx  = geo.intervals[0][3]; int ngy  = geo.intervals[1][3]; int ngz  = geo.intervals[2][3];

    update_bc_block(face, n_type, f.tem, f.work_ref, f_n_1.tem, f_n_2.tem, f_n_3.tem, f_n_4.tem,
            nx, ny, nz, ngx, ngy, ngz, coord_x_rel, coord_y_rel, coord_z_rel);

}

void derefine_and_merge_field(std::unordered_map<std::string,Field> & field_map, Geometry & geo_first_child, 
        std::string & new_key_str, std::vector<std::string> & key_str_list) {

    derefine_and_merge(geo_first_child, field_map[new_key_str].tem, field_map[new_key_str].work_ref,
            field_map[key_str_list[0]].tem, field_map[key_str_list[1]].tem, 
            field_map[key_str_list[2]].tem, field_map[key_str_list[3]].tem, 
            field_map[key_str_list[4]].tem, field_map[key_str_list[5]].tem, 
            field_map[key_str_list[6]].tem, field_map[key_str_list[7]].tem);

}

void refine_and_split_field(std::unordered_map<std::string,Field> & field_map, Geometry & geo_parent, 
        std::string & key_str, std::vector<std::string> & key_str_list) {

    refine_and_split(geo_parent, field_map[key_str].tem, field_map[key_str].work_ref,
            field_map[key_str_list[0]].tem, field_map[key_str_list[1]].tem, 
            field_map[key_str_list[2]].tem, field_map[key_str_list[3]].tem, 
            field_map[key_str_list[4]].tem, field_map[key_str_list[5]].tem, 
            field_map[key_str_list[6]].tem, field_map[key_str_list[7]].tem);
}

void allocate_field(Field & field, std::vector<int> & np, std::vector<int> & ng) {
    std::vector<int> alloc_size = (np+ng); alloc_size = alloc_size + ng;
    std::vector<int> alloc_size_ref = np+np; alloc_size_ref = alloc_size_ref+ng; alloc_size_ref = alloc_size_ref+ng;
    field.tem      = allocate3d(alloc_size); field.rhs_tem  = allocate3d(alloc_size); field.prhs_tem = allocate3d(alloc_size);

    for(auto &t: field.temp) {
        t = allocate3d(alloc_size); 
    }

    field.work_ref = allocate3d(alloc_size_ref);

    int nx  = np[0] ; int ny  = np[1] ; int nz  = np[2];
    int ngx = ng[0] ; int ngy = ng[1] ; int ngz = ng[2];
    for(int i=0;i<2*ngx+nx;i++) 
        for(int j=0;j<2*ngy+ny;j++) 
            for(int k=0;k<2*ngz+nz;k++) {
                field.tem[i][j][k]      = 0.; field.rhs_tem[i][j][k]  = 0.; field.prhs_tem[i][j][k] = 0.;
            }
//#ifdef USE_CUDA
//    allocate_field_dev(field, alloc_size, alloc_size_ref);
//#endif
}

void deallocate_field(Field & field, std::vector<int> & np, std::vector<int> & ng) {
    std::vector<int> alloc_size = (np+ng); alloc_size = alloc_size + ng;
    std::vector<int> alloc_size_ref = np+np; alloc_size_ref = alloc_size_ref+ng; alloc_size_ref = alloc_size_ref+ng;
    deallocate3d(&(field.tem), alloc_size); deallocate3d(&(field.rhs_tem), alloc_size); deallocate3d(&(field.prhs_tem), alloc_size);

    for(auto &t: field.temp) {
        deallocate3d(&t, alloc_size);
    }

    deallocate3d(&(field.work_ref), alloc_size_ref);
//#ifdef USE_CUDA
//    deallocate_field_dev(field, alloc_size, alloc_size_ref);
//#endif
}

void init_field(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map) {

    for(auto &n: geo_map) {
        Field field;
        Geometry & geo = n.second;
        int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
        int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];
        std::vector<int> np = {nx, ny, nz};
        std::vector<int> ng = {ngx, ngy, ngz};
        allocate_field(field, np, ng);
        field.refinement_needed = 0;
        for(int i=0;i<2*ngx+nx;i++) 
            for(int j=0;j<2*ngy+ny;j++) 
                for(int k=0;k<2*ngz+nz;k++) {
                    field.tem[i][j][k] = 0.;

                    field.rhs_tem[i][j][k] = 0.F; field.prhs_tem[i][j][k] = 0.F;
                }
        copy_geo(geo, field.geo);
        field_map[n.first] = field;
    }

}

void evolve_advanced_immerse(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map) {

    for(auto &n: field_map) {
        Field & f= n.second;
        Geometry & geo = geo_map[n.first];
        evolve_advanced_immerse_block(f.tem, f.temp[0], f.dist, geo);
    }
}

void correct_advanced_immerse(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map) {

    for(auto &n: field_map) {
        Field & f= n.second;
        Geometry & geo = geo_map[n.first];
        correct_advanced_immerse_block(1., f.tem, f.temp[0], f.dist, geo);
    }
}

void prhs_compute(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map, Params & params) {

    float qnrk=(params.dt)*(params.brk[params.n]);
    //std::cout << "prhs_compute qnrk= " << qnrk << std::endl;
    for(auto &n: field_map) {
        Field & f= n.second;
        Geometry & geo = geo_map[n.first];
        int nx   = geo.intervals[0][2]; int ny   = geo.intervals[1][2]; int nz   = geo.intervals[2][2];
        int ngx  = geo.intervals[0][3]; int ngy  = geo.intervals[1][3]; int ngz  = geo.intervals[2][3];
#pragma omp parallel for
        for(int i=ngx;i<ngx+nx;i++) 
            for(int j=ngy;j<ngy+ny;j++) 
#pragma omp simd
                for(int k=ngz;k<ngz+nz;k++) {
                    if(f.dist[i][j][k] > 0.) {
                        f.prhs_tem[i][j][k] = f.tem[i][j][k] + qnrk*f.rhs_tem[i][j][k];
                    }
                }
    }
}

void linear_compute(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map, Params & params) {

//#ifdef USE_CUDA
//    linear_compute_dev(geo_map, field_map, params);
//#else
    float qnrk=(params.dt)*(params.ark[params.n]);
    //std::cout << "linear_compute qnrk= " << qnrk << std::endl;
    for(auto &n: field_map) {
        Field & f= n.second;
        Geometry & geo = geo_map[n.first];
        int nx   = geo.intervals[0][2]; int ny   = geo.intervals[1][2]; int nz   = geo.intervals[2][2];
        int ngx  = geo.intervals[0][3]; int ngy  = geo.intervals[1][3]; int ngz  = geo.intervals[2][3];
#pragma omp parallel for
        for(int i=ngx;i<ngx+nx;i++) 
            for(int j=ngy;j<ngy+ny;j++) 
#pragma omp simd
                for(int k=ngz;k<ngz+nz;k++) {
                    if(f.dist[i][j][k] > 0.) {
                        f.tem[i][j][k] = f.prhs_tem[i][j][k] + qnrk*f.rhs_tem[i][j][k];
                    }
                }
    }
//#endif

}

void set_simple_immerse(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map, Params & params) {

    float dt = params.dt;
    std::string & immersion_type = params.immersion_type;
    for(auto &n: field_map) {
        Field & f= n.second;
        Geometry & geo = geo_map[n.first];
        int nx   = geo.intervals[0][2]; int ny   = geo.intervals[1][2]; int nz   = geo.intervals[2][2];
        int ngx  = geo.intervals[0][3]; int ngy  = geo.intervals[1][3]; int ngz  = geo.intervals[2][3];
        float dx = geo.lengths[0][3]  ; float dy = geo.lengths[1][3]; float dz = geo.lengths[2][3];

#pragma omp parallel for
        for(int i=0;i<2*ngx+nx;i++) 
            for(int j=0;j<2*ngy+ny;j++) 
#pragma omp simd
                for(int k=0;k<2*ngz+nz;k++) {
                    if(f.dist[i][j][k] < 0.) {
                        f.tem[i][j][k] = 0.F;
                    }
                }
    }
}

void rhs_compute(std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map, Params & params) {

    float dt = params.dt;

    std::string & immersion_type = params.immersion_type;
    for(auto &n: field_map) {
        Field & f= n.second;
        Geometry & geo = geo_map[n.first];
        int nx   = geo.intervals[0][2]; int ny   = geo.intervals[1][2]; int nz   = geo.intervals[2][2];
        int ngx  = geo.intervals[0][3]; int ngy  = geo.intervals[1][3]; int ngz  = geo.intervals[2][3];
        float dx = geo.lengths[0][3]  ; float dy = geo.lengths[1][3];   float dz = geo.lengths[2][3];
        int nt[3] = {nx, ny, nz}; int ng[3] = {ngx, ngy, ngz}; float delta[3] = {dx, dy, dz};
        //std::cout << "dx, dy, dz: " << dx << " " << dy << " " << dz << std::endl;

        std::vector<int> npp = {nx, ny, nz};
        std::vector<int> ngg = {ngx,ngy,ngz};
        std::vector<int> alloc_size = (npp+ngg); alloc_size = alloc_size + ngg;

        float ***d2t_dx2  = f.temp[0];  float ***d2t_dy2  = f.temp[1];  float ***d2t_dz2  = f.temp[2];  

        der2(f.tem, d2t_dx2, 1, delta, nt, ng);
        der2(f.tem, d2t_dy2, 2, delta, nt, ng);
        der2(f.tem, d2t_dz2, 3, delta, nt, ng);

        //check_nan(f.tem, alloc_size);

#pragma omp parallel for collapse(2)
        for(int i=ngx;i<ngx+nx;i++) 
            for(int j=ngy;j<ngy+ny;j++) 
#pragma omp simd
                for(int k=ngz;k<ngz+nz;k++) {
                    if(f.dist[i][j][k] >= 0.) {
                        f.rhs_tem[i][j][k] = d2t_dx2[i][j][k] + d2t_dy2[i][j][k] + d2t_dz2[i][j][k];
                    }
                }

        check_nan(f.rhs_tem, alloc_size);
    }

}

void mark_refinement_needed(Field &f){

    std::string ref_type = "distance"; //gradient";
    Geometry &geo = f.geo;

    if(ref_type == "gradient") {
//        float max_grad_allowed=0.02F;
//        float max_grad=0.F;
//        float cur_grad;
//        int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
//        int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];
//        for(int i=ngx;i<ngx+nx;i++) 
//            for(int j=ngy;j<ngy+ny;j++) 
//                for(int k=ngz;k<ngz+nz;k++) 
//                {
//                    cur_grad = (f.u[i+1][j][k]-f.u[i-1][j][k])*(f.u[i+1][j][k]-f.u[i-1][j][k]) + 
//                        (f.u[i][j+1][k]-f.u[i][j-1][k])*(f.u[i][j+1][k]-f.u[i][j-1][k]) + 
//                        (f.u[i][j][k+1]-f.u[i][j][k-1])*(f.u[i][j][k+1]-f.u[i][j][k-1]) ;
//                    if(cur_grad > max_grad) max_grad = cur_grad;
//                    //            std::cout << "cur_grad: " << cur_grad << std::endl;
//                }
//        std::cout << "max_grad: " << max_grad << std::endl;
//        if(max_grad > max_grad_allowed) {
//            f.refinement_needed = +1;
//            std::cout << "marking refinement needed: " << std::endl;
//        } else {
//            f.refinement_needed = 0;
//        }
//
    } else if(ref_type == "distance") {
        int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
        int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];
        float dx = geo.lengths[0][3] ; float dy = geo.lengths[1][3];  float dz = geo.lengths[2][3];
        float sx = geo.lengths[0][0] ; float sy = geo.lengths[1][0];  float sz = geo.lengths[2][0];
        float ex = geo.lengths[0][1] ; float ey = geo.lengths[1][1];  float ez = geo.lengths[2][1];
        float mx = 0.5*(sx+ex);        float my = 0.5*(sy+ey);        float mz = 0.5*(sz+ez);
        float min_distance = std::numeric_limits<float>::max();
#pragma omp parallel for
        for(int i=ngx;i<ngx+nx;i++) 
            for(int j=ngy;j<ngy+ny;j++) 
                for(int k=ngz;k<ngz+nz;k++) {
                    min_distance = std::min(f.dist[i][j][k], min_distance);
                }
        float max_dx_allowed = mesh_law(min_distance);

        //RIMETTEREfloat miny = std::min(std::abs(sy), std::abs(ey));
        //RIMETTEREfloat minz = std::min(std::abs(sz), std::abs(ez));
        //RIMETTEREif(miny*miny+minz*minz < 1.) {
        //RIMETTERE    max_dx_allowed = std::min(0.1F, max_dx_allowed);
        //RIMETTERE    //RIMETTEREmax_dx_allowed = std::min(0.05F, max_dx_allowed);
        //RIMETTERE}
        std::cout << "dx, max_dx_allowed, min_distance: " << dx << " ; " << max_dx_allowed << " " << min_distance << std::endl;
        if(dx > max_dx_allowed) {
            f.refinement_needed = +1;
            std::cout << "marking refinement needed: " << std::endl;
        } else if(2.2*dx < max_dx_allowed) {
            f.refinement_needed = -1;
            std::cout << "DISABLED marking derefinement needed: " << std::endl;
        } else {
            f.refinement_needed = 0;
        }
    }

}

std::vector<float> compute_force_field(Triangle & tri, Field & f, Geometry & geo) {
    // TODO
}
