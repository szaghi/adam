#include <string>
#include <array>
#include <vector>
#include <iostream>
#include <iterator>
#include <sstream>
#include <cmath>
#include <cstdlib>

#include "utils.h"

#ifdef DEBALLOC
#define DEBUG_MEMORY
#include "debAlloc.h"
#endif

std::vector<int> operator+(std::vector<int> &a, std::vector<int> &b)
{
  std::vector<int> c;
  c.resize(a.size());
  for(int i = 0 ; i< a.size();i++)
  {
    c[i] = a[i] + b[i];
  }
  return c;
}

void KeyStringToArray(const std::string& key_str, std::array<int,7> &key_arr)
{
	std::vector<std::string> k;
	std::istringstream iss(key_str);
	std::copy(std::istream_iterator<std::string>(iss),
			std::istream_iterator<std::string>(),
			std::back_inserter(k));
	for(int i=0;i<7;i++)
	{
		//std::cout << k[i] << std::endl;
		key_arr[i] = std::stoi(k[i]);
	}
}

// key_arr cannot be a reference. why? yes it can now...
void KeyArrayToString(std::array<int,7> &key_arr, std::string &key_str)
{
	std::stringstream ks;
	ks << key_arr[0] << " " << key_arr[1] << " " << key_arr[2] << " " << key_arr[3]
       << " " << key_arr[4] << " " << key_arr[5] << " " << key_arr[6];
	key_str = ks.str();
}

void allocate3d_old(float **** v, std::vector<int> dims){
	size_t size= dims[0]*dims[1]*dims[2];
	float *data =(float*) malloc(size*sizeof(float)); 
	*v = (float ***)malloc(dims[0]*sizeof(float ***));
	for(int i=0;i<dims[0];i++) 
	{
		(*v)[i] = (float **) malloc(dims[1]*sizeof(float**));
		for(int j=0;j<dims[1];j++) 
		{
			(*v)[i][j] = data+dims[1]*dims[2]*i + dims[2]*j; 
		} 
	}
}

float*** allocate3d(std::vector<int> dims){
    float ***v;
	size_t size= dims[0]*dims[1]*dims[2];
	float *data =(float*) malloc(size*sizeof(float)); 
    if(data == NULL) {
        std::cout << "malloc data failed" << std::endl; exit(EXIT_FAILURE);
    }
	v = (float ***)malloc(dims[0]*sizeof(float ***));
    if(v == NULL) {
        std::cout << "malloc v failed" << std::endl; exit(EXIT_FAILURE);
    }
	for(int i=0;i<dims[0];i++) 
	{
		v[i] = (float **) malloc(dims[1]*sizeof(float**));
        if(v[i] == NULL) {
            std::cout << "malloc v[i] failed" << std::endl; exit(EXIT_FAILURE);
        }
		for(int j=0;j<dims[1];j++) 
		{
			v[i][j] = data+dims[1]*dims[2]*i + dims[2]*j; 
		} 
	}
    return v;
}

void deallocate3d(float ****v, std::vector<int> dims){
    // free data pointer
    free(&((*v)[0][0][0]));
    // free pointer to lines
	for(int i=0;i<dims[0];i++) 
	{
        //std::cout << "freeing: " << i << std::endl;
		free((*v)[i]);
	}
    // free pointer to planes
    free(*v);
}

void check_nan(float ***v, std::vector<int> dims){
    for(int i=0;i<dims[0];i++) 
        for(int j=0;j<dims[1];j++) 
            for(int k=0;k<dims[2];k++) {
                if(std::isnan(v[i][j][k])) {
                    std::cout << "NaN found on i,j,k: " << i << " " << j << " " << k << " " << std::endl;
                    exit(EXIT_FAILURE);
                }
            }
}
