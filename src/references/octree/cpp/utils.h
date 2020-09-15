#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <array>
#include <vector>

typedef struct {
	float t;
    float dt;
    int n;
    int n_step;
    std::array<float,4> ark;
    std::array<float,4> brk;
    std::string immersion_type;
} Params;

std::vector<int> operator+(std::vector<int> &a, std::vector<int> &b);

void KeyStringToArray(const std::string& key_str, std::array<int,7> &key_arr);

void KeyArrayToString(std::array<int,7> &key_arr, std::string& key_str);

void allocate3d_old(float **** v,std::vector<int> dims);

float *** allocate3d(std::vector<int> dims);

void deallocate3d(float ****v, std::vector<int> dims);

void check_nan(float ***v, std::vector<int> dims);

//#define IND(i,j,k) (i)*np[1]*np[2]+(j)*np[2]+(k)

//#define allocate_cube(type, var, x, y, z)                            \
//    do {                                                             \
//        (var) = malloc(sizeof(type **) * (x));                       \
//        if ((var) == NULL) {                                         \
//            perror("Error allocating memory.\n");                    \
//        }                                                            \
//        int i, j;                                                    \
//        for (i = 0; i < x; i++) {                                    \
//            (var)[i] = malloc(sizeof(type*) * (y));                  \
//            if ((var)[i] == NULL) {                                  \
//                perror("Error allocating memory.\n");                \
//            }                                                        \
//            for (j = 0; j < y; j++) {                                \
//                (var)[i][j] = malloc(sizeof(type) * (z));            \
//                if ((var)[i][j] == NULL) {                           \
//                    perror("Error allocating memory.\n");            \
//                }                                                    \
//            }                                                        \
//        }                                                            \
//    } while(0)
//
#endif
