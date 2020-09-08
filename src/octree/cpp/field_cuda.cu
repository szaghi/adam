#include "cuda.h"
#include "geometry.h"
#include "field.h"
#include "field_cuda.h"
#include "stdio.h"
#include <iostream>
#include <cstdlib>

void allocate_dev(std::unordered_map<std::string, Field> & field_map, Fields_cuda & fields_cuda, 
    std::vector<int> & np, std::vector<int> & ng, float ls[3], float le[3]) {

    int device_num = 0;
    float safety_mem_factor = 0.9;

    int nx = np[0];
    int ny = np[1];
    int nz = np[2];
    int ngx = ng[0];
    int ngy = ng[1];
    int ngz = ng[2];
    float dx = (le[0]-ls[0])/np[0];
    float dy = (le[1]-ls[1])/np[1];
    float dz = (le[2]-ls[2])/np[2];
    fields_cuda.nx  = nx;
    fields_cuda.ny  = ny;
    fields_cuda.nz  = nz;
    fields_cuda.ngx = ngx;
    fields_cuda.ngy = ngy;
    fields_cuda.ngz = ngz;

    int f1_size = (nx+2*ngx)*(ny+2*ngy)*(nz+2*ngz);
    int f1_size_ref = (2*nx+2*ngx)*(2*ny+2*ngy)*(2*nz+2*ngz);
    fields_cuda.f1_size = f1_size;
    fields_cuda.f1_size_ref = f1_size_ref;
    // 8 is the number of variables with size f1_size
    size_t block_size_f1 = (f1_size)*sizeof(float); 
    size_t block_size_f1_ref = (f1_size_ref)*sizeof(float); 
    size_t block_size_allfields = 7*block_size_f1+1*block_size_f1_ref+4*sizeof(float);

    struct cudaDeviceProp prop;
    cudaGetDeviceProperties (&prop, device_num);
    std::cout << "Total global memory: " << prop.totalGlobalMem << std::endl;
    size_t totalGlobalMem = prop.totalGlobalMem;

    int max_blocks = totalGlobalMem * safety_mem_factor / block_size_allfields;
    fields_cuda.max_blocks = max_blocks;
    std::cout << "Max number of blocks on GPU: " << max_blocks << std::endl;

    cudaMalloc((void **)&(fields_cuda.dx_dev), max_blocks*sizeof(float));
    cudaMalloc((void **)&(fields_cuda.dy_dev), max_blocks*sizeof(float));
    cudaMalloc((void **)&(fields_cuda.dz_dev), max_blocks*sizeof(float));
    cudaMalloc((void **)&(fields_cuda.level_dev), max_blocks*sizeof(int));

    size_t allblocks_size_f1 = block_size_f1*max_blocks;
    size_t allblocks_size_f1_ref = block_size_f1_ref*max_blocks;

    cudaMalloc((void **)&(fields_cuda.tem_dev), allblocks_size_f1);
    cudaMemset(fields_cuda.tem_dev, 0., allblocks_size_f1);

    cudaMalloc((void **)&(fields_cuda.rhs_tem_dev), allblocks_size_f1);
    cudaMemset(fields_cuda.rhs_tem_dev, 0., allblocks_size_f1);

    cudaMalloc((void **)&(fields_cuda.prhs_tem_dev), allblocks_size_f1);
    cudaMemset(fields_cuda.prhs_tem_dev, 0., allblocks_size_f1);

    cudaMalloc((void **)&(fields_cuda.dist_dev), allblocks_size_f1);
    cudaMemset(fields_cuda.dist_dev, 0., allblocks_size_f1);

    for(auto &t_dev: fields_cuda.temp_dev) {
        cudaMalloc((void **)&(t_dev), allblocks_size_f1);
        cudaMemset(t_dev, 0., allblocks_size_f1);
    }

    cudaMalloc((void **)&(fields_cuda.work_ref_dev), allblocks_size_f1_ref);
    cudaMemset(fields_cuda.work_ref_dev, 0., allblocks_size_f1_ref);

    size_t n_size = 6*max_blocks*sizeof(int);
    fields_cuda.i_n_1           = (int*)malloc(n_size);
    fields_cuda.i_n_2           = (int*)malloc(n_size);
    fields_cuda.i_n_3           = (int*)malloc(n_size);
    fields_cuda.i_n_4           = (int*)malloc(n_size);
    fields_cuda.n_type          = (int*)malloc(n_size);
    fields_cuda.coord_x_rel     = (int*)malloc(n_size);
    fields_cuda.coord_y_rel     = (int*)malloc(n_size);
    fields_cuda.coord_z_rel     = (int*)malloc(n_size);
    cudaMalloc((void **)&(fields_cuda.i_n_1_dev), n_size);
    cudaMalloc((void **)&(fields_cuda.i_n_2_dev), n_size);
    cudaMalloc((void **)&(fields_cuda.i_n_3_dev), n_size);
    cudaMalloc((void **)&(fields_cuda.i_n_4_dev), n_size);
    cudaMalloc((void **)&(fields_cuda.n_type_dev), n_size);
    cudaMalloc((void **)&(fields_cuda.coord_x_rel_dev), n_size);
    cudaMalloc((void **)&(fields_cuda.coord_y_rel_dev), n_size);
    cudaMalloc((void **)&(fields_cuda.coord_z_rel_dev), n_size);

    // start cuda streams
    for (int i = 0; i < 9; i++)
        cudaStreamCreate(&(fields_cuda.streams[i]));

#ifdef NVIDIA_GPU
    fields_cuda.cb_x    = 32;
    fields_cuda.cb_y    = 32;
    fields_cuda.cb_z    = 1;
    fields_cuda.cb_bc_x = 32;
    fields_cuda.cb_bc_y = 16;
    fields_cuda.cb_bc_z = 1;
#endif
#ifdef AMD_GPU
    fields_cuda.cb_x    = 64;
    fields_cuda.cb_y    = 4;
    fields_cuda.cb_z    = 4;
    fields_cuda.cb_bc_x = 64;
    fields_cuda.cb_bc_y = 16;
    fields_cuda.cb_bc_z = 1;
#endif

}

void update_cpu_gpu(std::unordered_map<std::string, Field> & field_map, Fields_cuda & fields_cuda, 
    std::unordered_map<std::string, Geometry> & geo_map, std::string dir) {

    float *T_tem_dev      = fields_cuda.tem_dev;
    float *T_rhs_tem_dev  = fields_cuda.rhs_tem_dev;
    float *T_prhs_tem_dev = fields_cuda.prhs_tem_dev;
    float *T_dist_dev     = fields_cuda.dist_dev;
    float *T_dx_dev       = fields_cuda.dx_dev;
    float *T_dy_dev       = fields_cuda.dy_dev;
    float *T_dz_dev       = fields_cuda.dz_dev;
    int *T_level_dev      = fields_cuda.level_dev;

    int n_blocks = 0;
    int ib=0;
    fields_cuda.index_map.clear();
    for(auto &n: field_map) {
        n_blocks++;
        fields_cuda.index_map[n.first] = ib;
        ib++;
    }
    if(n_blocks > fields_cuda.max_blocks) {
        std::cout << "n_blocks > max_blocks: " << n_blocks << " > " << fields_cuda.max_blocks << std::endl;
        exit(EXIT_FAILURE);
    }

    int sf1     = fields_cuda.f1_size;
    size_t sf1b = sf1*sizeof(float);
    size_t sf   = sizeof(float);

    ib = 0;
    std::string ts;
    for(auto &n: field_map) {
        std::array<int,7> key_arr; KeyStringToArray(n.first, key_arr);
        int level = key_arr[0];
        Field & field = n.second;
        if(dir == "cpu_to_gpu") {
            // --------------------------------------------------------------------------------------------------
            //  TODO - PADDING PER FAR INIZIARE OGNI BLOCCO IN MULTIPLO DI 128B (CUDA) E ??? (AMD) PER AVERE COALESCENCE
            // --------------------------------------------------------------------------------------------------
            cudaMemcpy(T_tem_dev,      &(field.tem[0][0][0]),      sf1b, cudaMemcpyHostToDevice);
            cudaMemcpy(T_rhs_tem_dev,  &(field.rhs_tem[0][0][0]),  sf1b, cudaMemcpyHostToDevice);
            cudaMemcpy(T_prhs_tem_dev, &(field.prhs_tem[0][0][0]), sf1b, cudaMemcpyHostToDevice);
            cudaMemcpy(T_dist_dev,     &(field.dist[0][0][0]),     sf1b, cudaMemcpyHostToDevice);
            cudaMemcpy(T_dx_dev,       &(field.geo.lengths[0][3]), sf,   cudaMemcpyHostToDevice);
            cudaMemcpy(T_dy_dev,       &(field.geo.lengths[1][3]), sf,   cudaMemcpyHostToDevice);
            cudaMemcpy(T_dz_dev,       &(field.geo.lengths[2][3]), sf,   cudaMemcpyHostToDevice);
            cudaMemcpy(T_level_dev,    &level,                     sf,   cudaMemcpyHostToDevice);
            for(int iface=0;iface<6;iface++) {
                ts = field.key_n_1[iface]; fields_cuda.i_n_1[ib*6+iface] = fields_cuda.index_map[ts];
                ts = field.key_n_2[iface]; fields_cuda.i_n_2[ib*6+iface] = fields_cuda.index_map[ts];
                ts = field.key_n_3[iface]; fields_cuda.i_n_3[ib*6+iface] = fields_cuda.index_map[ts];
                ts = field.key_n_4[iface]; fields_cuda.i_n_4[ib*6+iface] = fields_cuda.index_map[ts];
                fields_cuda.coord_x_rel[ib*6+iface] = field.coord_x_rel[iface];
                fields_cuda.coord_y_rel[ib*6+iface] = field.coord_y_rel[iface];
                fields_cuda.coord_z_rel[ib*6+iface] = field.coord_z_rel[iface];
                ts = field.n_type[iface];
                if(ts == "phys") {
                    fields_cuda.n_type[ib*6+iface] = 0;
                } else if(ts == "n_same_refined") {
                    fields_cuda.n_type[ib*6+iface] = 1;
                } else if(ts == "n_more_refined") {
                    fields_cuda.n_type[ib*6+iface] = 2;
                } else if(ts == "n_less_refined") {
                    fields_cuda.n_type[ib*6+iface] = 3;
                }
            }
        } else if(dir == "gpu_to_cpu") {
            cudaMemcpy(&(field.tem[0][0][0]),      T_tem_dev,      sf1b, cudaMemcpyDeviceToHost);
            cudaMemcpy(&(field.rhs_tem[0][0][0]),  T_rhs_tem_dev,  sf1b, cudaMemcpyDeviceToHost);
            cudaMemcpy(&(field.prhs_tem[0][0][0]), T_prhs_tem_dev, sf1b, cudaMemcpyDeviceToHost);
            cudaMemcpy(&(field.dist[0][0][0]),     T_dist_dev,     sf1b, cudaMemcpyDeviceToHost);
        }
        T_tem_dev      += sf1;
        T_rhs_tem_dev  += sf1;
        T_prhs_tem_dev += sf1;
        T_dist_dev     += sf1;
        T_dx_dev       += 1;
        T_dy_dev       += 1;
        T_dz_dev       += 1;
        T_level_dev    += 1;
        ib++;
    }
    if(dir == "cpu_to_gpu") {
        fields_cuda.n_blocks = n_blocks;
        size_t sf1c = 6*n_blocks*sizeof(int);
        cudaMemcpy(fields_cuda.i_n_1_dev, fields_cuda.i_n_1, sf1c, cudaMemcpyHostToDevice);
        cudaMemcpy(fields_cuda.i_n_2_dev, fields_cuda.i_n_2, sf1c, cudaMemcpyHostToDevice);
        cudaMemcpy(fields_cuda.i_n_3_dev, fields_cuda.i_n_3, sf1c, cudaMemcpyHostToDevice);
        cudaMemcpy(fields_cuda.i_n_4_dev, fields_cuda.i_n_4, sf1c, cudaMemcpyHostToDevice);
        cudaMemcpy(fields_cuda.n_type_dev, fields_cuda.n_type, sf1c, cudaMemcpyHostToDevice);
        cudaMemcpy(fields_cuda.coord_x_rel_dev, fields_cuda.coord_x_rel, sf1c, cudaMemcpyHostToDevice);
        cudaMemcpy(fields_cuda.coord_y_rel_dev, fields_cuda.coord_y_rel, sf1c, cudaMemcpyHostToDevice);
        cudaMemcpy(fields_cuda.coord_z_rel_dev, fields_cuda.coord_z_rel, sf1c, cudaMemcpyHostToDevice);
    }
}

__global__ void prhs_compute_dev_kernel(float *tem_dev, float *prhs_tem_dev, float *rhs_tem_dev, float qnrk, float * dist_dev, 
                                        int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size) {
    int k  = threadIdx.x+blockIdx.x*blockDim.x + ngz; // coalescence
    int j  = threadIdx.y+blockIdx.y*blockDim.y + ngy;
    int ib = threadIdx.z+blockIdx.z*blockDim.z;
    size_t ind=0;
    int i=0;
    //printf("GGG FROM CUDA: %d %d %d %d %d %d %d %d %d\n:",threadIdx.x, blockIdx.x, blockDim.x, threadIdx.y, blockIdx.y, blockDim.y, j, k, ib);
    if(j < (ngy+ny) && k < (ngz+nz) && ib < n_blocks) {
        //for(ib=0;ib<n_blocks;ib++) {
        for(i=ngx;i<ngx+nx;i++) {
            ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
            if(dist_dev[ind] > 0.) {
                prhs_tem_dev[ind] = tem_dev[ind] + qnrk * rhs_tem_dev[ind];
            }
        }
        //}
    }
}

void prhs_compute_dev(Fields_cuda & fields_cuda, Params & params) {

    float qnrk=(params.dt)*(params.brk[params.n]);
    //std::cout << "linear_compute qnrk= " << qnrk << std::endl;

    int nx  = fields_cuda.nx;  int ny  = fields_cuda.ny;  int nz  = fields_cuda.nz;
    int ngx = fields_cuda.ngx; int ngy = fields_cuda.ngy; int ngz = fields_cuda.ngz;
    int n_blocks = fields_cuda.n_blocks; int f1_size  = fields_cuda.f1_size;

    //dim3 cuda_blocks = {32, 16, 1};
    // Number of threads per cuda_block must be less or equal than 1024!!!
    //dim3 cuda_blocks = {32, 4, 8};
    dim3 cuda_blocks = {std::min(fields_cuda.cb_x, nz), std::min(fields_cuda.cb_y, ny), std::min(fields_cuda.cb_z, n_blocks)};
    unsigned int grid_x = (nz+cuda_blocks.x-1)/cuda_blocks.x;
    unsigned int grid_y = (ny+cuda_blocks.y-1)/cuda_blocks.y;
    unsigned int grid_z = (n_blocks+cuda_blocks.z-1)/cuda_blocks.z;
    //dim3 cuda_grid  = {grid_x, grid_y, 1};
    dim3 cuda_grid  = {grid_x, grid_y, grid_z};

    //std::cout << "CUDA grid: " << cuda_grid.x << " " << cuda_grid.y << " " << cuda_grid.z << " " << 
    //    cuda_blocks.x << " " << cuda_blocks.y << " " << cuda_blocks.z << std::endl;

    prhs_compute_dev_kernel<<<cuda_grid, cuda_blocks>>>(fields_cuda.tem_dev, fields_cuda.prhs_tem_dev, 
        fields_cuda.rhs_tem_dev, qnrk, fields_cuda.dist_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size);
    MY_CUDA_CHECK( cudaDeviceSynchronize() );

}

__global__ void linear_compute_dev_kernel(float *tem_dev, float *prhs_tem_dev, float *rhs_tem_dev, float qnrk, float * dist_dev, 
                                          int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size) {
    int k = threadIdx.x+blockIdx.x*blockDim.x + ngz; // coalescence
    int j = threadIdx.y+blockIdx.y*blockDim.y + ngy;
    int ib = threadIdx.z+blockIdx.z*blockDim.z;

    //KILLSCACHEf[i][j][k][ib]
    //KILLSCACHEint ib = threadIdx.x+blockIdx.x*blockDim.x; // coalescence
    //KILLSCACHEint k  = threadIdx.y+blockIdx.y*blockDim.y + ngz;
    //KILLSCACHEint j  = threadIdx.z+blockIdx.z*blockDim.z + ngy;

    size_t ind=0;
    int i=0;
    //printf("FROM CUDA: %d %d %d %d %d %d %d %d %d %d\n:",threadIdx.x, blockIdx.x, blockDim.x, threadIdx.y, blockIdx.y, blockDim.y, i, j, k, ind);
    if(j < (ngy+ny) && k < (ngz+nz) && ib < n_blocks) {
        //for(ib=0;ib<n_blocks;ib++) {
            for(i=ngx;i<ngx+nx;i++) {
                ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                if(dist_dev[ind] > 0.) {
                    tem_dev[ind] = prhs_tem_dev[ind] + qnrk * rhs_tem_dev[ind];
                }
            }
        //}
    }
}

void linear_compute_dev(Fields_cuda & fields_cuda, Params & params) {

    float qnrk=(params.dt)*(params.ark[params.n]);
    //std::cout << "linear_compute qnrk= " << qnrk << std::endl;

    int nx  = fields_cuda.nx;  int ny  = fields_cuda.ny;  int nz  = fields_cuda.nz;
    int ngx = fields_cuda.ngx; int ngy = fields_cuda.ngy; int ngz = fields_cuda.ngz;
    int n_blocks = fields_cuda.n_blocks; int f1_size  = fields_cuda.f1_size;

    //dim3 cuda_blocks = {32, 4, 8};
    dim3 cuda_blocks = {std::min(fields_cuda.cb_x, nz), std::min(fields_cuda.cb_y, ny), std::min(fields_cuda.cb_z, n_blocks)};
    unsigned int grid_x = (nz+cuda_blocks.x-1)/cuda_blocks.x;
    unsigned int grid_y = (ny+cuda_blocks.y-1)/cuda_blocks.y;
    unsigned int grid_z = (n_blocks+cuda_blocks.z-1)/cuda_blocks.z;
    dim3 cuda_grid  = {grid_x, grid_y, grid_z};
    //std::cout << "CUDA grid: " << cuda_grid.x << " " << cuda_grid.y << " " << cuda_grid.z << " " << 
    //    cuda_blocks.x << " " << cuda_blocks.y << " " << cuda_blocks.z << std::endl;

    linear_compute_dev_kernel<<<cuda_grid, cuda_blocks>>>(fields_cuda.tem_dev, fields_cuda.prhs_tem_dev, 
        fields_cuda.rhs_tem_dev, qnrk, fields_cuda.dist_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size);
    cudaDeviceSynchronize();

}

__global__ void rhs_compute_dev_kernel(float *tem_dev, float *rhs_tem_dev, float *dist_dev, 
        float *dx_dev, float *dy_dev, float *dz_dev, 
        int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size) {

    int k = threadIdx.x+blockIdx.x*blockDim.x + ngz; // coalescence
    int j = threadIdx.y+blockIdx.y*blockDim.y + ngy;
    int ib = threadIdx.z+blockIdx.z*blockDim.z;
    size_t ind=0; size_t ind_px; size_t ind_mx; size_t ind_py; size_t ind_my; size_t ind_pz; size_t ind_mz;
    int i=0;
    float d2t_dx2; float d2t_dy2; float d2t_dz2;
    float inv_dx2; float inv_dy2; float inv_dz2;
    //printf("FROM CUDA: %d %d %d %d %d %d %d %d %d %d\n:",threadIdx.x, blockIdx.x, blockDim.x, threadIdx.y, blockIdx.y, blockDim.y, i, j, k, ind);
    if(j < (ngy+ny) && k < (ngz+nz) && ib < n_blocks) {
        //for(ib=0;ib<n_blocks;ib++) {
            inv_dx2 = 1./(dx_dev[ib]*dx_dev[ib]); 
            inv_dy2 = 1./(dy_dev[ib]*dy_dev[ib]); 
            inv_dz2 = 1./(dz_dev[ib]*dz_dev[ib]); 
            //float dxl = dx
            for(i=ngx;i<ngx+nx;i++) {
                ind    = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                ind_px = k+j*(nz+2*ngz)+(i+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                ind_mx = k+j*(nz+2*ngz)+(i-1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                ind_py = k+(j+1)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                ind_my = k+(j-1)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                ind_pz = k+1+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                ind_mz = k-1+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                if(dist_dev[ind] > 0.) {
                    d2t_dx2 = inv_dx2 * (tem_dev[ind_px]+tem_dev[ind_mx]-2.*tem_dev[ind]);
                    d2t_dy2 = inv_dy2 * (tem_dev[ind_py]+tem_dev[ind_my]-2.*tem_dev[ind]);
                    d2t_dz2 = inv_dz2 * (tem_dev[ind_pz]+tem_dev[ind_mz]-2.*tem_dev[ind]);
                    rhs_tem_dev[ind] = d2t_dx2 + d2t_dy2 + d2t_dz2;
                }
            }
        //}
    }

}

void rhs_compute_dev(Fields_cuda & fields_cuda, Params & params) {

    int nx  = fields_cuda.nx;  int ny  = fields_cuda.ny;  int nz  = fields_cuda.nz;
    int ngx = fields_cuda.ngx; int ngy = fields_cuda.ngy; int ngz = fields_cuda.ngz;
    int n_blocks = fields_cuda.n_blocks; int f1_size  = fields_cuda.f1_size;

    //dim3 cuda_blocks = {32, 4, 8};
    dim3 cuda_blocks = {std::min(fields_cuda.cb_x, nz), std::min(fields_cuda.cb_y, ny), std::min(fields_cuda.cb_z, n_blocks)};
    unsigned int grid_x = (nz+cuda_blocks.x-1)/cuda_blocks.x;
    unsigned int grid_y = (ny+cuda_blocks.y-1)/cuda_blocks.y;
    unsigned int grid_z = (n_blocks+cuda_blocks.z-1)/cuda_blocks.z;
    dim3 cuda_grid  = {grid_x, grid_y, grid_z};
    //std::cout << "CUDA grid: " << cuda_grid.x << " " << cuda_grid.y << " " << cuda_grid.z << " " << 
    //    cuda_blocks.x << " " << cuda_blocks.y << " " << cuda_blocks.z << std::endl;

    rhs_compute_dev_kernel<<<cuda_grid, cuda_blocks>>>(fields_cuda.tem_dev, fields_cuda.rhs_tem_dev, 
        fields_cuda.dist_dev, fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev,
        nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size);
    cudaDeviceSynchronize();

}

//void evolve_advanced_immerse_dev(std::unordered_map<std::string,Geometry> & geo_map,
//        std::unordered_map<std::string,Field> & field_map) {

__global__ void evolve_advanced_immerse_dev_kernel_p1(float *f_dev, float *f_dev_new, float *dist_dev, 
        float *dx_dev, float *dy_dev, float *dz_dev, 
        int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size) {

    int k = threadIdx.x+blockIdx.x*blockDim.x + ngz; // coalescence
    int j = threadIdx.y+blockIdx.y*blockDim.y + ngy;
    int ib = threadIdx.z+blockIdx.z*blockDim.z;
    size_t ind=0; size_t ind_px; size_t ind_mx; size_t ind_py; size_t ind_my; size_t ind_pz; size_t ind_mz;
    int i=0;
    float dist_x; float dist_y; float dist_z; float mod_dist; float f_x; float f_y; float f_z;
    float inv_dx2; float inv_dy2; float inv_dz2;
    //printf("FROM CUDA: %d %d %d %d %d %d %d %d %d %d\n:",threadIdx.x, blockIdx.x, blockDim.x, threadIdx.y, blockIdx.y, blockDim.y, i, j, k, ind);
    if(j < (ngy+ny) && k < (ngz+nz) && ib < n_blocks) {
        //for(ib=0;ib<n_blocks;ib++) {
            for(i=ngx;i<ngx+nx;i++) {
                ind    = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                if(dist_dev[ind] < 0.) {
                    ind_px = k+j*(nz+2*ngz)+(i+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_mx = k+j*(nz+2*ngz)+(i-1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_py = k+(j+1)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_my = k+(j-1)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_pz = k+1+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_mz = k-1+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    dist_x = dist_dev[ind_px]-dist_dev[ind_mx];
                    dist_y = dist_dev[ind_py]-dist_dev[ind_my];
                    dist_z = dist_dev[ind_pz]-dist_dev[ind_mz];
                    //mod_dist = sqrt(dist_x*dist_x+dist_y*dist_y+dist_z*dist_z);
                    mod_dist = std::abs(dist_x)+std::abs(dist_y)+std::abs(dist_z);
                    dist_x = dist_x/mod_dist;
                    dist_y = dist_y/mod_dist;
                    dist_z = dist_z/mod_dist;
                    if(dist_x > 0) {
                        f_x = f_dev[ind_px]-f_dev[ind];
                    } else {
                        f_x = f_dev[ind]-f_dev[ind_mx];
                    }
                    if(dist_y > 0) {
                        f_y = f_dev[ind_py]-f_dev[ind];
                    } else {
                        f_y = f_dev[ind]-f_dev[ind_my];
                    }
                    if(dist_z > 0) {
                        f_z = f_dev[ind_pz]-f_dev[ind];
                    } else {
                        f_z = f_dev[ind]-f_dev[ind_mz];
                    }
                    f_dev_new[ind] = f_dev[ind] + 0.9*(dist_x*f_x+dist_y*f_y+dist_z*f_z);
                }
            }
        //}
    }

}

// --------------------------------------------------
// TODO fare la prossima come switch puntatori
// --------------------------------------------------
__global__ void evolve_advanced_immerse_dev_kernel_p2(float *f_dev, float *f_dev_new, float *dist_dev, 
        float *dx_dev, float *dy_dev, float *dz_dev, 
        int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size) {

    int k = threadIdx.x+blockIdx.x*blockDim.x + ngz; // coalescence
    int j = threadIdx.y+blockIdx.y*blockDim.y + ngy;
    int ib = threadIdx.z+blockIdx.z*blockDim.z;
    size_t ind=0; size_t ind_px; size_t ind_mx; size_t ind_py; size_t ind_my; size_t ind_pz; size_t ind_mz;
    int i=0;
    float dist_x; float dist_y; float dist_z; float mod_dist; float f_x; float f_y; float f_z;
    float inv_dx2; float inv_dy2; float inv_dz2;
    //printf("FROM CUDA: %d %d %d %d %d %d %d %d %d %d\n:",threadIdx.x, blockIdx.x, blockDim.x, threadIdx.y, blockIdx.y, blockDim.y, i, j, k, ind);
    if(j < (ngy+ny) && k < (ngz+nz) && ib < n_blocks) {
        for(i=ngx;i<ngx+nx;i++) {
            ind    = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
            if(dist_dev[ind] < 0.) {
                f_dev[ind] = f_dev_new[ind];
            }
        }
    }

}

void evolve_advanced_immerse_dev(Fields_cuda & fields_cuda, Params & params) {

    int nx  = fields_cuda.nx;  int ny  = fields_cuda.ny;  int nz  = fields_cuda.nz;
    int ngx = fields_cuda.ngx; int ngy = fields_cuda.ngy; int ngz = fields_cuda.ngz;
    int n_blocks = fields_cuda.n_blocks; int f1_size  = fields_cuda.f1_size;

    //dim3 cuda_blocks = {32, 4, 8};
    dim3 cuda_blocks = {std::min(fields_cuda.cb_x, nz), std::min(fields_cuda.cb_y, ny), std::min(fields_cuda.cb_z, n_blocks)};
    unsigned int grid_x = (nz+cuda_blocks.x-1)/cuda_blocks.x;
    unsigned int grid_y = (ny+cuda_blocks.y-1)/cuda_blocks.y;
    unsigned int grid_z = (n_blocks+cuda_blocks.z-1)/cuda_blocks.z;
    dim3 cuda_grid  = {grid_x, grid_y, grid_z};
    //std::cout << "CUDA grid: " << cuda_grid.x << " " << cuda_grid.y << " " << cuda_grid.z << " " << 
    //    cuda_blocks.x << " " << cuda_blocks.y << " " << cuda_blocks.z << std::endl;
    //cudaStream_t streams[num_streams];
    //float *data[num_streams];
    //for (int i = 0; i < num_streams; i++) {
    //    cudaStreamCreate(&streams[i]);

    evolve_advanced_immerse_dev_kernel_p1<<<cuda_grid, cuda_blocks>>>(fields_cuda.tem_dev,
        fields_cuda.temp_dev[0], fields_cuda.dist_dev, fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev,
        nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size);

    evolve_advanced_immerse_dev_kernel_p2<<<cuda_grid, cuda_blocks>>>(fields_cuda.tem_dev,
        fields_cuda.temp_dev[0], fields_cuda.dist_dev, fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev,
        nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size);

    cudaDeviceSynchronize();

}

__global__ void correct_advanced_immerse_dev_kernel(float imposed_val, float *f_dev, float *f_dev_new, float *dist_dev, 
        float *dx_dev, float *dy_dev, float *dz_dev, 
        int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size) {

    int k = threadIdx.x+blockIdx.x*blockDim.x + ngz; // coalescence
    int j = threadIdx.y+blockIdx.y*blockDim.y + ngy;
    int ib = threadIdx.z+blockIdx.z*blockDim.z;
    size_t ind=0; size_t ind_px; size_t ind_mx; size_t ind_py; size_t ind_my; size_t ind_pz; size_t ind_mz;
    int i=0;
    float dist_x; float dist_y; float dist_z; float mod_dist; float f_x; float f_y; float f_z;
    float inv_dx2; float inv_dy2; float inv_dz2;
    //printf("FROM CUDA: %d %d %d %d %d %d %d %d %d %d\n:",threadIdx.x, blockIdx.x, blockDim.x, threadIdx.y, blockIdx.y, blockDim.y, i, j, k, ind);
    if(j < (ngy+ny) && k < (ngz+nz) && ib < n_blocks) {
        //for(ib=0;ib<n_blocks;ib++) {
            for(i=ngx;i<ngx+nx;i++) {
                ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                if(dist_dev[ind] < 0.) {
                    //COLPADILUCAf_dev_new[ind] = 2.*imposed_val-f_dev[ind];
                    f_dev[ind] = 2.*imposed_val-f_dev[ind];
                }
            }
       // }
        //for(ib=0;ib<n_blocks;ib++) {
        //COLPADILUCA    for(i=ngx;i<ngx+nx;i++) {
        //COLPADILUCA        ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
        //COLPADILUCA        if(dist_dev[ind] < 0.) {
        //COLPADILUCA            f_dev[ind] = f_dev_new[ind];
        //COLPADILUCA        }
        //COLPADILUCA    }
       // }
    }

}

void correct_advanced_immerse_dev(Fields_cuda & fields_cuda, Params & params) {

    int nx  = fields_cuda.nx;  int ny  = fields_cuda.ny;  int nz  = fields_cuda.nz;
    int ngx = fields_cuda.ngx; int ngy = fields_cuda.ngy; int ngz = fields_cuda.ngz;
    int n_blocks = fields_cuda.n_blocks; int f1_size  = fields_cuda.f1_size;

    //dim3 cuda_blocks = {32, 4, 8};
    dim3 cuda_blocks = {std::min(fields_cuda.cb_x, nz), std::min(fields_cuda.cb_y, ny), std::min(fields_cuda.cb_z, n_blocks)};
    unsigned int grid_x = (nz+cuda_blocks.x-1)/cuda_blocks.x;
    unsigned int grid_y = (ny+cuda_blocks.y-1)/cuda_blocks.y;
    unsigned int grid_z = (n_blocks+cuda_blocks.z-1)/cuda_blocks.z;
    dim3 cuda_grid  = {grid_x, grid_y, grid_z};
    //std::cout << "CUDA grid: " << cuda_grid.x << " " << cuda_grid.y << " " << cuda_grid.z << " " << 
    //    cuda_blocks.x << " " << cuda_blocks.y << " " << cuda_blocks.z << std::endl;

    float imposed_val = 1.F;
    correct_advanced_immerse_dev_kernel<<<cuda_grid, cuda_blocks>>>(imposed_val, fields_cuda.tem_dev,
        fields_cuda.temp_dev[0], fields_cuda.dist_dev, fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev,
        nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size);
    cudaDeviceSynchronize();

}

// TODO LEVARE GHOST NEI THREADS
__global__ void update_bc_phys_same_lr_dev_kernel(float *f_dev, float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size) {

    int k = threadIdx.x+blockIdx.x*blockDim.x + ngz; // coalescence
    int j = threadIdx.y+blockIdx.y*blockDim.y + ngy;
    int ib, ib_n_1, i;
    size_t ind, ind_p;
    if(j < (ngy+ny) && k < (ngz+nz)) {
        for(ib=0; ib < n_blocks; ib++) {
            if(n_type_dev[6*ib] == 0) { // left phys
                //printf("FROM CUDA IMPOSING LEFT PHYSICS FOR BLOCK: %d\n:",ib);
                for(i=0;i<ngx;i++) {
                    ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind] = 1.;
                }
            }
            if(n_type_dev[6*ib] == 1) { // left n_same_refined
                ib_n_1 = i_n_1_dev[6*ib];
                for(i=0;i<ngx;i++) {
                    ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_p = k+j*(nz+2*ngz)+(nx+i)*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    f_dev[ind] = f_dev[ind_p];
                }
            }
            if(n_type_dev[6*ib+1] == 0) { // right phys
                for(i=ngx+nx;i<2*ngx+nx;i++) {
                    ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_p = k+j*(nz+2*ngz)+(ngx+nx-1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind] = f_dev[ind_p];
                }
            }
            if(n_type_dev[6*ib+1] == 1) { // right n_same_refined
                ib_n_1 = i_n_1_dev[6*ib+1];
                for(i=0;i<ngx;i++) {
                    ind = k+j*(nz+2*ngz)+(i+nx+ngx)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_p = k+j*(nz+2*ngz)+(ngx+i)*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    f_dev[ind] = f_dev[ind_p];
                }
            }
        }
    }
}

__global__ void update_bc_phys_same_bt_dev_kernel(float *f_dev, float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size) {

    int k = threadIdx.x+blockIdx.x*blockDim.x + ngz; // coalescence
    int i = threadIdx.y+blockIdx.y*blockDim.y + ngx;
    int ib, ib_n_1, j;
    size_t ind, ind_p;
    if(i < (ngx+nx) && k < (ngz+nz)) {
        for(ib=0; ib < n_blocks; ib++) {
            if(n_type_dev[6*ib+2] == 0) { // bottom phys
                //printf("FROM CUDA IMPOSING LEFT PHYSICS FOR BLOCK: %d\n:",ib);
                for(j=0;j<ngy;j++) {
                    ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_p = k+ngy*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind] = f_dev[ind_p];
                }
            }
            if(n_type_dev[6*ib+2] == 1) { // bottom n_same_refined
                ib_n_1 = i_n_1_dev[6*ib+2];
                for(j=0;j<ngy;j++) {
                    ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_p = k+(ny+j)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    f_dev[ind] = f_dev[ind_p];
                }
            }
            if(n_type_dev[6*ib+3] == 0) { // top phys
                for(j=ngy+ny;j<2*ngy+ny;j++) {
                    ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_p = k+(ngy+ny-1)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind] = f_dev[ind_p];
                }
            }
            if(n_type_dev[6*ib+3] == 1) { // top n_same_refined
                ib_n_1 = i_n_1_dev[6*ib+3];
                for(j=0;j<ngy;j++) {
                    ind = k+(j+ny+ngy)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_p = k+(j+ngy)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    f_dev[ind] = f_dev[ind_p];
                }
            }
        }
    }
}

__global__ void update_bc_phys_same_bf_dev_kernel_OPT(float *f_dev, float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size) {

    int k = threadIdx.x; // coalescence
    int j = threadIdx.y+blockIdx.y*blockDim.y + ngy;
    int ib = blockIdx.z;

    int ib_n_1, i;
    size_t ind, ind_p;
    if(j < (ngy+ny) && ib < n_blocks) {
        if(n_type_dev[6*ib+4] == 0) { // back phys
            //printf("FROM CUDA IMPOSING LEFT PHYSICS FOR BLOCK: %d\n:",ib);
            for(i=ngx;i<ngx+nx;i++) {
                ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                ind_p = ngz+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                f_dev[ind] = f_dev[ind_p];
            }
        }
        if(n_type_dev[6*ib+4] == 1) { // back n_same_refined
            ib_n_1 = i_n_1_dev[6*ib+4];
            for(i=ngx;i<ngx+nx;i++) {
                ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                ind_p = (nz+k)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                f_dev[ind] = f_dev[ind_p];
            }
        }
        if(n_type_dev[6*ib+5] == 0) { // front phys
            k += ngz+nz;
            for(i=ngx;i<ngx+nx;i++) {
                ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                ind_p = (ngz+nz-1)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size; 
                f_dev[ind] = f_dev[ind_p];
            }
        }
        if(n_type_dev[6*ib+5] == 1) { // front n_same_refined
            ib_n_1 = i_n_1_dev[6*ib+5];
            for(i=ngx;i<ngx+nx;i++) {
                ind = (k+nz+ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                ind_p = (k+ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                f_dev[ind] = f_dev[ind_p];
            }
        }
    }
}
__global__ void update_bc_phys_same_bf_dev_kernel(float *f_dev, float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size) {

    int j = threadIdx.x+blockIdx.x*blockDim.x + ngy; // NOT coalescence
    int i = threadIdx.y+blockIdx.y*blockDim.y + ngx;

    int ib, ib_n_1, k;
    size_t ind, ind_p;
    if(j < (ngy+ny) && i < (ngx+nx)) {
        for(ib=0; ib < n_blocks; ib++) {
            if(n_type_dev[6*ib+4] == 0) { // back phys
                //printf("FROM CUDA IMPOSING LEFT PHYSICS FOR BLOCK: %d\n:",ib);
                for(k=0;k<ngz;k++) {
                    ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_p = ngz+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind] = f_dev[ind_p];
                }
            }
            if(n_type_dev[6*ib+4] == 1) { // back n_same_refined
                ib_n_1 = i_n_1_dev[6*ib+4];
                for(k=0;k<ngz;k++) {
                    ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_p = (nz+k)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    f_dev[ind] = f_dev[ind_p];
                }
            }
            if(n_type_dev[6*ib+5] == 0) { // front phys
                for(k=ngz+nz;k<2*ngz+nz;k++) {
                    ind = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_p = (ngz+nz-1)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size; 
                    f_dev[ind] = f_dev[ind_p];
                }
            }
            if(n_type_dev[6*ib+5] == 1) { // front n_same_refined
                ib_n_1 = i_n_1_dev[6*ib+5];
                for(k=0;k<ngz;k++) {
                    ind = (k+nz+ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    ind_p = (k+ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    f_dev[ind] = f_dev[ind_p];
                }
            }
        }
    }
}


__global__ void update_bc_less_lr_dev_kernel(float *f_dev, float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size,
        int *coord_x_rel, int *coord_y_rel, int *coord_z_rel) {

    int k_base = threadIdx.x+blockIdx.x*blockDim.x; // coalescence
    int j_base = threadIdx.y+blockIdx.y*blockDim.y;
    int i, j, k, ib, ib_n_1, i_ref, j_ref, k_ref;
    size_t ind=0, ind_ref;
    if(j_base < (ny/2) && k_base < (nz/2)) {
        for(ib=0; ib < n_blocks; ib++) {
            if(n_type_dev[6*ib] == 3) { // left less
                k = k_base + ngz+coord_z_rel[6*ib]*nz/2;
                j = j_base + ngy+coord_y_rel[6*ib]*ny/2;
                ib_n_1 = i_n_1_dev[6*ib];
                //printf("FROM CUDA IMPOSING LEFT PHYSICS FOR BLOCK: %d\n:",ib);
                for(i=ngx/2+nx;i<ngx+nx;i++) {
                    i_ref   = 2*(i-ngx/2-nx);
                    j_ref   = 2*(j-ngy-coord_y_rel[6*ib]*ny/2)+ngy;
                    k_ref   = 2*(k-ngz-coord_z_rel[6*ib]*nz/2)+ngz;
                    ind     = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    ind_ref = k_ref+j_ref*(nz+2*ngz)+i_ref*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+j_ref*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+(j_ref+1)*(nz+2*ngz)+(i_ref)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+(j_ref+1)*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+j_ref*(nz+2*ngz)+i_ref*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+j_ref*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+(j_ref+1)*(nz+2*ngz)+(i_ref)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+(j_ref+1)*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                }
            }
            if(n_type_dev[6*ib+1] == 3) { // right less
                k = k_base + ngz+coord_z_rel[6*ib+1]*nz/2;
                j = j_base + ngy+coord_y_rel[6*ib+1]*ny/2;
                ib_n_1 = i_n_1_dev[6*ib+1];
                for(i=ngx;i<ngx+ngx/2;i++) { 
                    i_ref = 2*(i-ngx)+nx+ngx; 
                    j_ref = 2*(j-ngy-coord_y_rel[6*ib+1]*ny/2)+ngy;
                    k_ref = 2*(k-ngz-coord_z_rel[6*ib+1]*nz/2)+ngz;
                    ind     = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    ind_ref = k_ref+j_ref*(nz+2*ngz)+i_ref*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+j_ref*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+(j_ref+1)*(nz+2*ngz)+(i_ref)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+(j_ref+1)*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+j_ref*(nz+2*ngz)+i_ref*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+j_ref*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+(j_ref+1)*(nz+2*ngz)+(i_ref)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+(j_ref+1)*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                }
            }
        }
    }
}

__global__ void update_bc_less_bt_dev_kernel(float *f_dev, float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size,
        int *coord_x_rel, int *coord_y_rel, int *coord_z_rel) {

    int k_base = threadIdx.x+blockIdx.x*blockDim.x; // coalescence
    int i_base = threadIdx.y+blockIdx.y*blockDim.y;
    int i, j, k, ib, ib_n_1, i_ref, j_ref, k_ref;
    size_t ind=0, ind_ref;
    if(k_base < (nz/2) && i_base < (nx/2)) {
        for(ib=0; ib < n_blocks; ib++) {
            if(n_type_dev[6*ib+2] == 3) { // bottom less
                k = k_base + ngz+coord_z_rel[6*ib+2]*nz/2;
                i = i_base + ngx+coord_x_rel[6*ib+2]*nx/2;
                ib_n_1 = i_n_1_dev[6*ib+2];
                //printf("FROM CUDA IMPOSING LEFT PHYSICS FOR BLOCK: %d\n:",ib);
                for(j=ngy/2+ny;j<ngy+ny;j++) {
                    j_ref   = 2*(j-ngy/2-ny);
                    k_ref   = 2*(k-ngz-coord_z_rel[6*ib+2]*nz/2)+ngz;
                    i_ref   = 2*(i-ngx-coord_x_rel[6*ib+2]*nx/2)+ngx;
                    ind     = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    ind_ref = k_ref+j_ref*(nz+2*ngz)+i_ref*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+j_ref*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+(j_ref+1)*(nz+2*ngz)+(i_ref)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+(j_ref+1)*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+j_ref*(nz+2*ngz)+i_ref*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+j_ref*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+(j_ref+1)*(nz+2*ngz)+(i_ref)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+(j_ref+1)*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                }
            }
            if(n_type_dev[6*ib+3] == 3) { // top less
                k = k_base + ngz+coord_z_rel[6*ib+3]*nz/2;
                i = i_base + ngx+coord_x_rel[6*ib+3]*nx/2;
                ib_n_1 = i_n_1_dev[6*ib+3];
                for(j=ngy;j<ngy+ngy/2;j++) { 
                    j_ref = 2*(j-ngy)+ny+ngy; 
                    k_ref = 2*(k-ngz-coord_z_rel[6*ib+3]*nz/2)+ngz;
                    i_ref = 2*(i-ngx-coord_x_rel[6*ib+3]*nx/2)+ngx;
                    ind     = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    ind_ref = k_ref+j_ref*(nz+2*ngz)+i_ref*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+j_ref*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+(j_ref+1)*(nz+2*ngz)+(i_ref)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+(j_ref+1)*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+j_ref*(nz+2*ngz)+i_ref*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+j_ref*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+(j_ref+1)*(nz+2*ngz)+(i_ref)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+(j_ref+1)*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                }
            }
        }
    }
}
__global__ void update_bc_less_bf_dev_kernel(float *f_dev, float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size,
        int *coord_x_rel, int *coord_y_rel, int *coord_z_rel) {

    int j_base = threadIdx.x+blockIdx.x*blockDim.x; // NOT coalescence
    int i_base = threadIdx.y+blockIdx.y*blockDim.y;
    int i, j, k, ib, ib_n_1, i_ref, j_ref, k_ref;
    size_t ind=0, ind_ref;
    if(j_base < (ny/2) && i_base < (nx/2)) {
        for(ib=0; ib < n_blocks; ib++) {
            if(n_type_dev[6*ib+4] == 3) { // back less
                j = j_base + ngy+coord_y_rel[6*ib+4]*ny/2;
                i = i_base + ngx+coord_x_rel[6*ib+4]*nx/2;
                ib_n_1 = i_n_1_dev[6*ib+4];
                for(k=ngz/2+nz;k<ngz+nz;k++) {
                    k_ref   = 2*(k-ngz/2-nz);
                    j_ref   = 2*(j-ngy-coord_y_rel[6*ib+4]*ny/2)+ngy;
                    i_ref   = 2*(i-ngx-coord_x_rel[6*ib+4]*nx/2)+ngx;
                    ind     = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    ind_ref = k_ref+j_ref*(nz+2*ngz)+i_ref*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+j_ref*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+(j_ref+1)*(nz+2*ngz)+(i_ref)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+(j_ref+1)*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+j_ref*(nz+2*ngz)+i_ref*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+j_ref*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+(j_ref+1)*(nz+2*ngz)+(i_ref)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+(j_ref+1)*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                }
            }
            if(n_type_dev[6*ib+5] == 3) { // front less
                j = j_base + ngy+coord_y_rel[6*ib+5]*ny/2;
                i = i_base + ngx+coord_x_rel[6*ib+5]*nx/2;
                ib_n_1 = i_n_1_dev[6*ib+5];
                for(k=ngz;k<ngz+ngz/2;k++) { 
                    k_ref = 2*(k-ngz)+nz+ngz; 
                    j_ref = 2*(j-ngy-coord_y_rel[6*ib+5]*ny/2)+ngy;
                    i_ref = 2*(i-ngx-coord_x_rel[6*ib+5]*nx/2)+ngx;
                    ind     = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    ind_ref = k_ref+j_ref*(nz+2*ngz)+i_ref*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+j_ref*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+(j_ref+1)*(nz+2*ngz)+(i_ref)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = k_ref+(j_ref+1)*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+j_ref*(nz+2*ngz)+i_ref*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+j_ref*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+(j_ref+1)*(nz+2*ngz)+(i_ref)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                    ind_ref = (k_ref+1)+(j_ref+1)*(nz+2*ngz)+(i_ref+1)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    f_dev[ind_ref] = f_dev[ind];
                }
            }
        }
    }
}

__global__ void update_bc_more_p1_lr_dev_kernel(float *f_dev, float *f_ref_dev,
        float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *i_n_2_dev, int *i_n_3_dev, int *i_n_4_dev, 
        int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size,
        int f1_size_ref, int *coord_x_rel, int *coord_y_rel, int *coord_z_rel) {

    int k = threadIdx.x+blockIdx.x*blockDim.x+ngz; // coalescence
    int j = threadIdx.y+blockIdx.y*blockDim.y+ngy;
    int i, ib, ib_n_1, ib_n_2, ib_n_3, ib_n_4;
    size_t ind, ind_ref;
    if(j < (ngy+ny) && k < (ngz+nz)) {
        for(ib=0; ib < n_blocks; ib++) {
            if(n_type_dev[6*ib] == 2) { // left more
                ib_n_1 = i_n_1_dev[6*ib];
                ib_n_2 = i_n_2_dev[6*ib];
                ib_n_3 = i_n_3_dev[6*ib];
                ib_n_4 = i_n_4_dev[6*ib];
                for(i=0;i<2*ngx;i++) {
                    ind     = k+j*(nz+2*ngz)+(i+nx-ngx)*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    ind_ref = k+j*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = k+j*(nz+2*ngz)+(i+nx-ngx)*(nz+2*ngz)*(ny+2*ngy)+ib_n_2*f1_size;
                    ind_ref = k+(j+ny)*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = k+j*(nz+2*ngz)+(i+nx-ngx)*(nz+2*ngz)*(ny+2*ngy)+ib_n_3*f1_size;
                    ind_ref = (k+nz)+j*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = k+j*(nz+2*ngz)+(i+nx-ngx)*(nz+2*ngz)*(ny+2*ngy)+ib_n_4*f1_size;
                    ind_ref = (k+nz)+(j+ny)*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                }
            }
            if(n_type_dev[6*ib+1] == 2) { // right more
                ib_n_1 = i_n_1_dev[6*ib+1];
                ib_n_2 = i_n_2_dev[6*ib+1];
                ib_n_3 = i_n_3_dev[6*ib+1];
                ib_n_4 = i_n_4_dev[6*ib+1];
                for(i=0;i<2*ngx;i++) {
                    ind     = k+j*(nz+2*ngz)+(i+ngx)*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    ind_ref = k+j*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = k+j*(nz+2*ngz)+(i+ngx)*(nz+2*ngz)*(ny+2*ngy)+ib_n_2*f1_size;
                    ind_ref = k+(j+ny)*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = k+j*(nz+2*ngz)+(i+ngx)*(nz+2*ngz)*(ny+2*ngy)+ib_n_3*f1_size;
                    ind_ref = (k+nz)+j*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = k+j*(nz+2*ngz)+(i+ngx)*(nz+2*ngz)*(ny+2*ngy)+ib_n_4*f1_size;
                    ind_ref = (k+nz)+(j+ny)*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                }
            }
        }
    }
}

__global__ void update_bc_more_p1_bt_dev_kernel(float *f_dev, float *f_ref_dev,
        float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *i_n_2_dev, int *i_n_3_dev, int *i_n_4_dev, 
        int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size,
        int f1_size_ref, int *coord_x_rel, int *coord_y_rel, int *coord_z_rel) {

    int k = threadIdx.x+blockIdx.x*blockDim.x+ngz; // coalescence
    int i = threadIdx.y+blockIdx.y*blockDim.y+ngx;
    int j, ib, ib_n_1, ib_n_2, ib_n_3, ib_n_4;
    size_t ind, ind_ref;
    if(i < (ngx+nx) && k < (ngz+nz)) {
        for(ib=0; ib < n_blocks; ib++) {
            if(n_type_dev[6*ib+2] == 2) { // bottom more
                ib_n_1 = i_n_1_dev[6*ib+2];
                ib_n_2 = i_n_2_dev[6*ib+2];
                ib_n_3 = i_n_3_dev[6*ib+2];
                ib_n_4 = i_n_4_dev[6*ib+2];
                for(j=0;j<2*ngy;j++) {
                    ind     = k+(j+ny-ngy)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    ind_ref = k+j*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = k+(j+ny-ngy)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_2*f1_size;
                    ind_ref = k+j*(2*nz+2*ngz)+(i+nx)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = k+(j+ny-ngy)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_3*f1_size;
                    ind_ref = (k+nz)+j*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = k+(j+ny-ngy)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_4*f1_size;
                    ind_ref = (k+nz)+j*(2*nz+2*ngz)+(i+nx)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                }
            }
            if(n_type_dev[6*ib+3] == 2) { // top more
                ib_n_1 = i_n_1_dev[6*ib+3];
                ib_n_2 = i_n_2_dev[6*ib+3];
                ib_n_3 = i_n_3_dev[6*ib+3];
                ib_n_4 = i_n_4_dev[6*ib+3];
                for(j=0;j<2*ngy;j++) {
                    ind     = k+(j+ngy)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    ind_ref = k+j*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = k+(j+ngy)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_2*f1_size;
                    ind_ref = k+j*(2*nz+2*ngz)+(i+nx)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = k+(j+ngy)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_3*f1_size;
                    ind_ref = (k+nz)+j*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = k+(j+ngy)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_4*f1_size;
                    ind_ref = (k+nz)+j*(2*nz+2*ngz)+(i+nx)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                }
            }
        }
    }
}

__global__ void update_bc_more_p1_bf_dev_kernel(float *f_dev, float *f_ref_dev,
        float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *i_n_2_dev, int *i_n_3_dev, int *i_n_4_dev, 
        int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size,
        int f1_size_ref, int *coord_x_rel, int *coord_y_rel, int *coord_z_rel) {

    int j = threadIdx.x+blockIdx.x*blockDim.x+ngy; // NOT coalescence
    int i = threadIdx.y+blockIdx.y*blockDim.y+ngx;
    int k, ib, ib_n_1, ib_n_2, ib_n_3, ib_n_4;
    size_t ind, ind_ref;
    if(j < (ngy+ny) && i < (ngx+nx)) {
        for(ib=0; ib < n_blocks; ib++) {
            if(n_type_dev[6*ib+4] == 2) { // back more
                ib_n_1 = i_n_1_dev[6*ib+4];
                ib_n_2 = i_n_2_dev[6*ib+4];
                ib_n_3 = i_n_3_dev[6*ib+4];
                ib_n_4 = i_n_4_dev[6*ib+4];
                for(k=0;k<2*ngz;k++) {
                    ind     = (k+nz-ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    ind_ref = k+j*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = (k+nz-ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_2*f1_size;
                    ind_ref = k+j*(2*nz+2*ngz)+(i+nx)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = (k+nz-ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_3*f1_size;
                    ind_ref = k+(j+ny)*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = (k+nz-ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_4*f1_size;
                    ind_ref = k+(j+ny)*(2*nz+2*ngz)+(i+nx)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                }
            }
            if(n_type_dev[6*ib+5] == 2) { // front more
                ib_n_1 = i_n_1_dev[6*ib+5];
                ib_n_2 = i_n_2_dev[6*ib+5];
                ib_n_3 = i_n_3_dev[6*ib+5];
                ib_n_4 = i_n_4_dev[6*ib+5];
                for(k=0;k<2*ngz;k++) {
                    ind     = (k+ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_1*f1_size;
                    ind_ref = k+j*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = (k+ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_2*f1_size;
                    ind_ref = k+j*(2*nz+2*ngz)+(i+nx)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = (k+ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_3*f1_size;
                    ind_ref = k+(j+ny)*(2*nz+2*ngz)+i*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                    ind     = (k+ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib_n_4*f1_size;
                    ind_ref = k+(j+ny)*(2*nz+2*ngz)+(i+nx)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_ref_dev[ind_ref] = f_dev[ind];
                }
            }
        }
    }
}

__global__ void update_bc_more_p2_lr_dev_kernel(float *f_dev, float *f_ref_dev,
        float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *i_n_2_dev, int *i_n_3_dev, int *i_n_4_dev, 
        int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size, int f1_size_ref, 
        int *coord_x_rel, int *coord_y_rel, int *coord_z_rel) {

    int k = threadIdx.x+blockIdx.x*blockDim.x+ngz; // coalescence
    int j = threadIdx.y+blockIdx.y*blockDim.y+ngy;
    int i, ib, ib_n_1, ib_n_2, ib_n_3, ib_n_4;
    int i_ref, j_ref, k_ref; 
    int ind_ref_1, ind_ref_2, ind_ref_3, ind_ref_4, ind_ref_5, ind_ref_6, ind_ref_7, ind_ref_8;
    size_t ind, ind_ref;
    if(j < (ngy+ny) && k < (ngz+nz)) {
        for(ib=0; ib < n_blocks; ib++) {
            if(n_type_dev[6*ib] == 2) { // left more
                for(i=0;i<ngx;i++) {
                    ind   = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    i_ref = 2*i; j_ref = 2*(j-ngy)+ngy; k_ref = 2*(k-ngz)+ngz;
                    ind_ref_1 = k_ref+j_ref*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_2 = k_ref+j_ref*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_3 = k_ref+(j_ref+1)*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_4 = k_ref+(j_ref+1)*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_5 = k_ref+1+j_ref*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_6 = k_ref+1+j_ref*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_7 = k_ref+1+(j_ref+1)*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_8 = k_ref+1+(j_ref+1)*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_dev[ind] = 0.125*(f_ref_dev[ind_ref_1]+f_ref_dev[ind_ref_2]+f_ref_dev[ind_ref_3]+f_ref_dev[ind_ref_4]+
                                      f_ref_dev[ind_ref_5]+f_ref_dev[ind_ref_6]+f_ref_dev[ind_ref_7]+f_ref_dev[ind_ref_8]);
                }
            }
            if(n_type_dev[6*ib+1] == 2) { // right more
                for(i=0;i<ngx;i++) {
                    ind   = k+j*(nz+2*ngz)+(i+nx+ngx)*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    i_ref = 2*i; j_ref = 2*(j-ngy)+ngy; k_ref = 2*(k-ngz)+ngz;
                    ind_ref_1 = k_ref+j_ref*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_2 = k_ref+j_ref*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_3 = k_ref+(j_ref+1)*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_4 = k_ref+(j_ref+1)*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_5 = k_ref+1+j_ref*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_6 = k_ref+1+j_ref*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_7 = k_ref+1+(j_ref+1)*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_8 = k_ref+1+(j_ref+1)*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_dev[ind] = 0.125*(f_ref_dev[ind_ref_1]+f_ref_dev[ind_ref_2]+f_ref_dev[ind_ref_3]+f_ref_dev[ind_ref_4]+
                                      f_ref_dev[ind_ref_5]+f_ref_dev[ind_ref_6]+f_ref_dev[ind_ref_7]+f_ref_dev[ind_ref_8]);
                }
            }
        }
    }
}

__global__ void update_bc_more_p2_bt_dev_kernel(float *f_dev, float *f_ref_dev,
        float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *i_n_2_dev, int *i_n_3_dev, int *i_n_4_dev, 
        int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size, int f1_size_ref, 
        int *coord_x_rel, int *coord_y_rel, int *coord_z_rel) {

    int k = threadIdx.x+blockIdx.x*blockDim.x+ngz; // coalescence
    int i = threadIdx.y+blockIdx.y*blockDim.y+ngx;
    int j, ib, ib_n_1, ib_n_2, ib_n_3, ib_n_4;
    int i_ref, j_ref, k_ref; 
    int ind_ref_1, ind_ref_2, ind_ref_3, ind_ref_4, ind_ref_5, ind_ref_6, ind_ref_7, ind_ref_8;
    size_t ind, ind_ref;
    if(i < (ngx+nx) && k < (ngz+nz)) {
        for(ib=0; ib < n_blocks; ib++) {
            if(n_type_dev[6*ib+2] == 2) { // bottom more
                for(j=0;j<ngy;j++) {
                    ind   = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    i_ref = 2*(i-ngx)+ngx; j_ref = 2*j; k_ref = 2*(k-ngz)+ngz;
                    ind_ref_1 = k_ref+j_ref*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_2 = k_ref+j_ref*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_3 = k_ref+(j_ref+1)*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_4 = k_ref+(j_ref+1)*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_5 = k_ref+1+j_ref*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_6 = k_ref+1+j_ref*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_7 = k_ref+1+(j_ref+1)*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_8 = k_ref+1+(j_ref+1)*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_dev[ind] = 0.125*(f_ref_dev[ind_ref_1]+f_ref_dev[ind_ref_2]+f_ref_dev[ind_ref_3]+f_ref_dev[ind_ref_4]+
                                      f_ref_dev[ind_ref_5]+f_ref_dev[ind_ref_6]+f_ref_dev[ind_ref_7]+f_ref_dev[ind_ref_8]);
                }
            }
            if(n_type_dev[6*ib+3] == 2) { // top more
                for(j=0;j<ngy;j++) {
                    ind   = k+(j+ny+ngy)*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    i_ref = 2*(i-ngx)+ngx; j_ref = 2*j; k_ref = 2*(k-ngz)+ngz;
                    ind_ref_1 = k_ref+j_ref*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_2 = k_ref+j_ref*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_3 = k_ref+(j_ref+1)*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_4 = k_ref+(j_ref+1)*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_5 = k_ref+1+j_ref*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_6 = k_ref+1+j_ref*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_7 = k_ref+1+(j_ref+1)*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_8 = k_ref+1+(j_ref+1)*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_dev[ind] = 0.125*(f_ref_dev[ind_ref_1]+f_ref_dev[ind_ref_2]+f_ref_dev[ind_ref_3]+f_ref_dev[ind_ref_4]+
                                      f_ref_dev[ind_ref_5]+f_ref_dev[ind_ref_6]+f_ref_dev[ind_ref_7]+f_ref_dev[ind_ref_8]);
                }
            }
        }
    }
}

__global__ void update_bc_more_p2_bf_dev_kernel(float *f_dev, float *f_ref_dev,
        float *dx_dev, float *dy_dev, float *dz_dev, 
        int *i_n_1_dev, int *i_n_2_dev, int *i_n_3_dev, int *i_n_4_dev, 
        int *n_type_dev, int nx, int ny, int nz, int ngx, int ngy, int ngz, int n_blocks, int f1_size, int f1_size_ref,
        int *coord_x_rel, int *coord_y_rel, int *coord_z_rel) {

    int j = threadIdx.x+blockIdx.x*blockDim.x+ngy; // NOT coalescence
    int i = threadIdx.y+blockIdx.y*blockDim.y+ngx;
    int k, ib, ib_n_1, ib_n_2, ib_n_3, ib_n_4;
    int i_ref, j_ref, k_ref; 
    int ind_ref_1, ind_ref_2, ind_ref_3, ind_ref_4, ind_ref_5, ind_ref_6, ind_ref_7, ind_ref_8;
    size_t ind, ind_ref;
    if(j < (ngy+ny) && i < (ngx+nx)) {
        for(ib=0; ib < n_blocks; ib++) {
            if(n_type_dev[6*ib+4] == 2) { // back more
                for(k=0;k<ngz;k++) {
                    ind   = k+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    i_ref = 2*(i-ngx)+ngx; j_ref = 2*(j-ngy)+ngy; k_ref = 2*k;
                    ind_ref_1 = k_ref+j_ref*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_2 = k_ref+j_ref*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_3 = k_ref+(j_ref+1)*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_4 = k_ref+(j_ref+1)*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_5 = k_ref+1+j_ref*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_6 = k_ref+1+j_ref*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_7 = k_ref+1+(j_ref+1)*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_8 = k_ref+1+(j_ref+1)*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_dev[ind] = 0.125*(f_ref_dev[ind_ref_1]+f_ref_dev[ind_ref_2]+f_ref_dev[ind_ref_3]+f_ref_dev[ind_ref_4]+
                                      f_ref_dev[ind_ref_5]+f_ref_dev[ind_ref_6]+f_ref_dev[ind_ref_7]+f_ref_dev[ind_ref_8]);
                }
            }
            if(n_type_dev[6*ib+5] == 2) { // front more
                for(k=0;k<ngz;k++) {
                    ind   = (k+nz+ngz)+j*(nz+2*ngz)+i*(nz+2*ngz)*(ny+2*ngy)+ib*f1_size;
                    i_ref = 2*(i-ngx)+ngx; j_ref = 2*(j-ngy)+ngy; k_ref = 2*k;
                    ind_ref_1 = k_ref+j_ref*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_2 = k_ref+j_ref*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_3 = k_ref+(j_ref+1)*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_4 = k_ref+(j_ref+1)*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_5 = k_ref+1+j_ref*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_6 = k_ref+1+j_ref*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_7 = k_ref+1+(j_ref+1)*(2*nz+2*ngz)+i_ref*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    ind_ref_8 = k_ref+1+(j_ref+1)*(2*nz+2*ngz)+(i_ref+1)*(2*nz+2*ngz)*(2*ny+2*ngy)+ib*f1_size_ref;
                    f_dev[ind] = 0.125*(f_ref_dev[ind_ref_1]+f_ref_dev[ind_ref_2]+f_ref_dev[ind_ref_3]+f_ref_dev[ind_ref_4]+
                                      f_ref_dev[ind_ref_5]+f_ref_dev[ind_ref_6]+f_ref_dev[ind_ref_7]+f_ref_dev[ind_ref_8]);
                }
            }
        }
    }
}

void update_bc_dev(std::unordered_map<std::string,Geometry> &geo_map, 
        std::unordered_map<std::string,Field> &field_map,
        Fields_cuda & fields_cuda) {

    int nx  = fields_cuda.nx;  int ny  = fields_cuda.ny;  int nz  = fields_cuda.nz;
    int ngx = fields_cuda.ngx; int ngy = fields_cuda.ngy; int ngz = fields_cuda.ngz;
    int n_blocks = fields_cuda.n_blocks; 
    int f1_size  = fields_cuda.f1_size; int f1_size_ref = fields_cuda.f1_size_ref;

    dim3 cuda_blocks, cuda_grid;
    unsigned int grid_x, grid_y, grid_z;

    // phys and same
    //cuda_blocks = {32, 16, 1};
    cuda_blocks = {std::min(fields_cuda.cb_bc_x, nz), std::min(fields_cuda.cb_bc_y, ny), fields_cuda.cb_bc_z};
    grid_x = (nz+cuda_blocks.x-1)/cuda_blocks.x;
    grid_y = (ny+cuda_blocks.y-1)/cuda_blocks.y;
    cuda_grid  = {grid_x, grid_y, 1};
    update_bc_phys_same_lr_dev_kernel<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[0]>>>(fields_cuda.tem_dev, 
        fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, fields_cuda.i_n_1_dev, 
        fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size);

    //cuda_blocks = {32, 16, 1};
    cuda_blocks = {std::min(fields_cuda.cb_bc_x, nz), std::min(fields_cuda.cb_bc_y, nx), fields_cuda.cb_bc_z};
    grid_x = (nz+cuda_blocks.x-1)/cuda_blocks.x;
    grid_y = (nx+cuda_blocks.y-1)/cuda_blocks.y;
    cuda_grid  = {grid_x, grid_y, 1};
    update_bc_phys_same_bt_dev_kernel<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[1]>>>(fields_cuda.tem_dev, 
        fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, fields_cuda.i_n_1_dev, 
        fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size);

    // ------------------------------------------------------------------------------------------
    // scegliere uno tra i due successivi
    // ------------------------------------------------------------------------------------------
    //cuda_blocks = {32, 16, 1};
    //RIMETTEREcuda_blocks = {std::min(fields_cuda.cb_bc_x, ny), std::min(fields_cuda.cb_bc_y, nx), fields_cuda.cb_bc_z};
    //RIMETTEREgrid_x = (ny+cuda_blocks.x-1)/cuda_blocks.x;
    //RIMETTEREgrid_y = (nx+cuda_blocks.y-1)/cuda_blocks.y;
    //RIMETTEREcuda_grid  = {grid_x, grid_y, 1};
    //RIMETTEREupdate_bc_phys_same_bf_dev_kernel<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[2]>>>(fields_cuda.tem_dev, 
    //RIMETTERE    fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, fields_cuda.i_n_1_dev, 
    //RIMETTERE    fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size);

    cuda_blocks = {ngz, std::min(fields_cuda.cb_bc_y, ny), 1};
    grid_x = 1;
    grid_y = (ny+cuda_blocks.y-1)/cuda_blocks.y;
    grid_z = (n_blocks+cuda_blocks.z-1)/cuda_blocks.z;
    cuda_grid  = {grid_x, grid_y, grid_z};
    update_bc_phys_same_bf_dev_kernel_OPT<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[2]>>>(fields_cuda.tem_dev, 
        fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, fields_cuda.i_n_1_dev, 
        fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size);
    // ------------------------------------------------------------------------------------------

    // less
    //cuda_blocks = {32, 16, 1};
    cuda_blocks = {std::min(fields_cuda.cb_bc_x, nz), std::min(fields_cuda.cb_bc_y, ny), fields_cuda.cb_bc_z};
    grid_x = (nz/2+cuda_blocks.x-1)/cuda_blocks.x;
    grid_y = (ny/2+cuda_blocks.y-1)/cuda_blocks.y;
    cuda_grid  = {grid_x, grid_y, 1};
    update_bc_less_lr_dev_kernel<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[3]>>>(fields_cuda.tem_dev, 
        fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, fields_cuda.i_n_1_dev, 
        fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size,
        fields_cuda.coord_x_rel_dev, fields_cuda.coord_y_rel_dev, fields_cuda.coord_z_rel_dev);

    //cuda_blocks = {32, 16, 1};
    cuda_blocks = {std::min(fields_cuda.cb_bc_x, nz), std::min(fields_cuda.cb_bc_y, nx), fields_cuda.cb_bc_z};
    grid_x = (nz/2+cuda_blocks.x-1)/cuda_blocks.x;
    grid_y = (nx/2+cuda_blocks.y-1)/cuda_blocks.y;
    cuda_grid  = {grid_x, grid_y, 1};
    update_bc_less_bt_dev_kernel<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[4]>>>(fields_cuda.tem_dev, 
        fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, fields_cuda.i_n_1_dev, 
        fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size,
        fields_cuda.coord_x_rel_dev, fields_cuda.coord_y_rel_dev, fields_cuda.coord_z_rel_dev);

    //cuda_blocks = {32, 16, 1};
    cuda_blocks = {std::min(fields_cuda.cb_bc_x, ny), std::min(fields_cuda.cb_bc_y, nx), fields_cuda.cb_bc_z};
    grid_x = (ny/2+cuda_blocks.x-1)/cuda_blocks.x;
    grid_y = (nx/2+cuda_blocks.y-1)/cuda_blocks.y;
    cuda_grid  = {grid_x, grid_y, 1};
    update_bc_less_bf_dev_kernel<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[5]>>>(fields_cuda.tem_dev, 
        fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, fields_cuda.i_n_1_dev, 
        fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size,
        fields_cuda.coord_x_rel_dev, fields_cuda.coord_y_rel_dev, fields_cuda.coord_z_rel_dev);

    // more part 1
    //cuda_blocks = {32, 16, 1};
    cuda_blocks = {std::min(fields_cuda.cb_bc_x, nz), std::min(fields_cuda.cb_bc_y, ny), fields_cuda.cb_bc_z};
    grid_x = (nz+cuda_blocks.x-1)/cuda_blocks.x;
    grid_y = (ny+cuda_blocks.y-1)/cuda_blocks.y;
    cuda_grid  = {grid_x, grid_y, 1};
    update_bc_more_p1_lr_dev_kernel<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[6]>>>(fields_cuda.tem_dev, 
        fields_cuda.work_ref_dev, fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, 
        fields_cuda.i_n_1_dev, fields_cuda.i_n_2_dev, fields_cuda.i_n_3_dev, fields_cuda.i_n_4_dev, 
        fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size, f1_size_ref, 
        fields_cuda.coord_x_rel_dev, fields_cuda.coord_y_rel_dev, fields_cuda.coord_z_rel_dev);

    //cuda_blocks = {32, 16, 1};
    cuda_blocks = {std::min(fields_cuda.cb_bc_x, nz), std::min(fields_cuda.cb_bc_y, nx), fields_cuda.cb_bc_z};
    grid_x = (nz+cuda_blocks.x-1)/cuda_blocks.x;
    grid_y = (nx+cuda_blocks.y-1)/cuda_blocks.y;
    cuda_grid  = {grid_x, grid_y, 1};
    update_bc_more_p1_bt_dev_kernel<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[7]>>>(fields_cuda.tem_dev, 
        fields_cuda.work_ref_dev, fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, 
        fields_cuda.i_n_1_dev, fields_cuda.i_n_2_dev, fields_cuda.i_n_3_dev, fields_cuda.i_n_4_dev, 
        fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size, f1_size_ref, 
        fields_cuda.coord_x_rel_dev, fields_cuda.coord_y_rel_dev, fields_cuda.coord_z_rel_dev);

    //cuda_blocks = {32, 16, 1};
    cuda_blocks = {std::min(fields_cuda.cb_bc_x, ny), std::min(fields_cuda.cb_bc_y, nx), fields_cuda.cb_bc_z};
    grid_x = (ny+cuda_blocks.x-1)/cuda_blocks.x;
    grid_y = (nx+cuda_blocks.y-1)/cuda_blocks.y;
    cuda_grid  = {grid_x, grid_y, 1};
    update_bc_more_p1_bf_dev_kernel<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[8]>>>(fields_cuda.tem_dev, 
        fields_cuda.work_ref_dev, fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, 
        fields_cuda.i_n_1_dev, fields_cuda.i_n_2_dev, fields_cuda.i_n_3_dev, fields_cuda.i_n_4_dev, 
        fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size, f1_size_ref,
        fields_cuda.coord_x_rel_dev, fields_cuda.coord_y_rel_dev, fields_cuda.coord_z_rel_dev);

    //cudaDeviceSynchronize();

    // more part 2
    //cuda_blocks = {32, 16, 1};
    cuda_blocks = {std::min(fields_cuda.cb_bc_x, nz), std::min(fields_cuda.cb_bc_y, ny), fields_cuda.cb_bc_z};
    grid_x = (nz+cuda_blocks.x-1)/cuda_blocks.x;
    grid_y = (ny+cuda_blocks.y-1)/cuda_blocks.y;
    cuda_grid  = {grid_x, grid_y, 1};
    update_bc_more_p2_lr_dev_kernel<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[6]>>>(fields_cuda.tem_dev, 
        fields_cuda.work_ref_dev, fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, 
        fields_cuda.i_n_1_dev, fields_cuda.i_n_2_dev, fields_cuda.i_n_3_dev, fields_cuda.i_n_4_dev, 
        fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size, f1_size_ref, 
        fields_cuda.coord_x_rel_dev, fields_cuda.coord_y_rel_dev, fields_cuda.coord_z_rel_dev);

    //cuda_blocks = {32, 16, 1};
    cuda_blocks = {std::min(fields_cuda.cb_bc_x, nz), std::min(fields_cuda.cb_bc_y, nx), fields_cuda.cb_bc_z};
    grid_x = (nz+cuda_blocks.x-1)/cuda_blocks.x;
    grid_y = (nx+cuda_blocks.y-1)/cuda_blocks.y;
    cuda_grid  = {grid_x, grid_y, 1};
    update_bc_more_p2_bt_dev_kernel<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[7]>>>(fields_cuda.tem_dev, 
        fields_cuda.work_ref_dev, fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, 
        fields_cuda.i_n_1_dev, fields_cuda.i_n_2_dev, fields_cuda.i_n_3_dev, fields_cuda.i_n_4_dev, 
        fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size, f1_size_ref, 
        fields_cuda.coord_x_rel_dev, fields_cuda.coord_y_rel_dev, fields_cuda.coord_z_rel_dev);

    //cuda_blocks = {32, 16, 1};
    cuda_blocks = {std::min(fields_cuda.cb_bc_x, ny), std::min(fields_cuda.cb_bc_y, nx), fields_cuda.cb_bc_z};
    grid_x = (ny+cuda_blocks.x-1)/cuda_blocks.x;
    grid_y = (nx+cuda_blocks.y-1)/cuda_blocks.y;
    cuda_grid  = {grid_x, grid_y, 1};
    update_bc_more_p2_bf_dev_kernel<<<cuda_grid, cuda_blocks, 0, fields_cuda.streams[8]>>>(fields_cuda.tem_dev, 
        fields_cuda.work_ref_dev, fields_cuda.dx_dev, fields_cuda.dy_dev, fields_cuda.dz_dev, 
        fields_cuda.i_n_1_dev, fields_cuda.i_n_2_dev, fields_cuda.i_n_3_dev, fields_cuda.i_n_4_dev, 
        fields_cuda.n_type_dev, nx, ny, nz, ngx, ngy, ngz, n_blocks, f1_size, f1_size_ref, 
        fields_cuda.coord_x_rel_dev, fields_cuda.coord_y_rel_dev, fields_cuda.coord_z_rel_dev);

    cudaDeviceSynchronize();

    //cudaDeviceReset();

}
