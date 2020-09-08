#ifdef USE_CUDA

#ifdef NVIDIA_GPU
#include "cuda.h"
#include "cuda_runtime.h"
#endif

#ifdef AMD_GPU
#include "hip/hip_runtime.h"
#endif

#include "field.h"
#include <vector>
#include <unordered_map>

#define MY_CUDA_CHECK( call) {                                    \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }

#define MY_CHECK_ERROR(errorMessage) {                                    \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    }

typedef struct {
    int n_blocks;
    int max_blocks;
    int nx, ny, nz;
    int ngx, ngy, ngz;
    int *level_dev;
    float *dx_dev, *dy_dev, *dz_dev;
    int f1_size; // size of device arrays
    int f1_size_ref; // size of device arrays ref
	float *tem_dev;
	float *rhs_tem_dev;
	float *prhs_tem_dev;
	float *work_ref_dev;
	float *dist_dev;
    std::array<float*, 3> temp_dev; //={NULL};
    int *i_n_1, *i_n_2, *i_n_3, *i_n_4, *n_type;
    int *coord_x_rel, *coord_y_rel, *coord_z_rel;
    std::unordered_map<std::string, int> index_map;
    int *i_n_1_dev, *i_n_2_dev, *i_n_3_dev, *i_n_4_dev, *n_type_dev;
    int *coord_x_rel_dev, *coord_y_rel_dev, *coord_z_rel_dev;
#ifdef NVIDIA_GPU
    cudaStream_t streams[9];
#endif
#ifdef AMD_GPU
    hipStream_t streams[9];
#endif
    int cb_x, cb_y, cb_z;
    int cb_bc_x, cb_bc_y, cb_bc_z;
} Fields_cuda;

void allocate_dev(std::unordered_map<std::string, Field> & field_map, Fields_cuda & fields_cuda, 
    std::vector<int> & np, std::vector<int> & ng, float ls[3], float le[3]);

void update_cpu_gpu(std::unordered_map<std::string, Field> & field_map, Fields_cuda & fields_cuda, 
    std::unordered_map<std::string, Geometry> & geo_map, std::string dir);

void prhs_compute_dev(Fields_cuda & fields_cuda, Params & params);

void linear_compute_dev(Fields_cuda & fields_cuda, Params & params);

void rhs_compute_dev(Fields_cuda & fields_cuda, Params & params);

void evolve_advanced_immerse_dev(Fields_cuda & fields_cuda, Params & params);

void correct_advanced_immerse_dev(Fields_cuda & fields_cuda, Params & params);

void update_bc_dev(std::unordered_map<std::string,Geometry> &geo_map, 
        std::unordered_map<std::string,Field> &field_map,
        Fields_cuda & fields_cuda);
#endif
