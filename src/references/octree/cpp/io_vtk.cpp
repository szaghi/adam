#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <vector>
#include "io_vtk.h"
#include "field.h"
#include "geometry.h"
#include "io_vtk_lib.h"

void save_vtr(std::unordered_map<std::string, Field> Field_map, int iter, std::vector<int> & np, std::vector<int> & ng) {

    std::vector<std::string> var_names = {"rho", "rhu", "rhv", "rhw", "rhe", "dist"};

    float *xg, *yg, *zg;

    std::string filename_list = "";
    std::vector<std::string> filenames;
    std::vector<float *> fs;
    fs.resize(6);
    for(auto &f: fs) { 
        f  = (float*)malloc(np[0]*np[1]*np[2]* sizeof(float));
    }
    int nx  = np[0]; int ny  = np[1]; int nz  = np[2]; 
    int ngx = ng[0]; int ngy = ng[1]; int ngz = ng[2];
    xg = (float*)malloc((nx+1) * sizeof(float));
    yg = (float*)malloc((ny+1) * sizeof(float));
    zg = (float*)malloc((nz+1) * sizeof(float));

    int i_field=-1;
    for(auto &n: Field_map) {
        i_field += 1;
        Field &field = n.second; // reference to field block

        Geometry geo = field.geo;
    
        float dx  = geo.lengths[0][3]  ; float dy  = geo.lengths[1][3]  ; float dz  = geo.lengths[2][3]  ; 
        float sx  = geo.lengths[0][0]  ; float sy  = geo.lengths[1][0]  ; float sz  = geo.lengths[2][0]  ; 
        std::cout << "save_vtr n* g*: => " << geo.key << " - " << nx << ny << nz << ngx << ngy << ngz << std::endl;
        std::cout << "save_vtr d* s*: => " << geo.key << " - " << dx << dy << dz << sx  << sy  << sz  << std::endl;
    
        // init grid
        for(int i=0;i<nx+1;i++) {
            xg[i] = sx + i *dx;
            //std::cout << "xg[i]: " << xg[i] << std::endl;
        }
        for(int j=0;j<ny+1;j++)
            yg[j] = sy + j *dy;
        for(int k=0;k<nz+1;k++)
            zg[k] = sz + k *dz;
    
        // print file
        std::string key_underscore = geo.key;
        std::replace(key_underscore.begin(), key_underscore.end(), ' ', '_'); // replace all ' ' to '_'
        std::string filename;

        filename = "output_"+key_underscore+"_b"+std::to_string(i_field)+"_i"+std::to_string(iter)+".vtr";
        filenames.push_back(filename);
        const char* filename_c = filename.c_str();
        filename_list += filename;

        // init field
        for(int i=0;i<nx;i++)
        for(int j=0;j<ny;j++)
        for(int k=0;k<nz;k++) {
            // Transpose data to match Fortran style
            //f[k*nx*ny+j*nx+i]   = field.rho[i+ngx][j+ngy][k+ngz];
            //f_2[k*nx*ny+j*nx+i] = field.rhu[i+ngx][j+ngy][k+ngz];
            //f_3[k*nx*ny+j*nx+i] = field.rhv[i+ngx][j+ngy][k+ngz];
            //f_4[k*nx*ny+j*nx+i] = field.rhw[i+ngx][j+ngy][k+ngz];
            //f_5[k*nx*ny+j*nx+i] = field.rhe[i+ngx][j+ngy][k+ngz];
            //f_6[k*nx*ny+j*nx+i] = field.dist[i+ngx][j+ngy][k+ngz];
            fs[0][k*nx*ny+j*nx+i] = field.tem[i+ngx][j+ngy][k+ngz];
            fs[1][k*nx*ny+j*nx+i] = field.tem[i+ngx][j+ngy][k+ngz];
            fs[2][k*nx*ny+j*nx+i] = field.tem[i+ngx][j+ngy][k+ngz];
            fs[3][k*nx*ny+j*nx+i] = field.tem[i+ngx][j+ngy][k+ngz];
            fs[4][k*nx*ny+j*nx+i] = field.tem[i+ngx][j+ngy][k+ngz];
            fs[5][k*nx*ny+j*nx+i] = field.dist[i+ngx][j+ngy][k+ngz];
        }

    //    std::cout << "RR: " << i_field << "tem[1+ngx][4+ngy][2+ngz]: " << field.tem[1+ngx][4+ngy][2+ngz] << std::endl;
    
        write_vtr(fs, nx, ny, nz, xg, yg, zg, filename, var_names);
    }

    std::string filename_all = "output_all_i"+std::to_string(iter)+".vtm";
    const char* filename_all_c = filename_all.c_str();
    const char* filename_list_c = filename_list.c_str();
    std::cout << "OUTPUT VTM" << filename_all_c << " " << filename_list_c << std::endl;
    write_vtm(i_field + 1, filename_all, filenames);

    for(auto &f: fs) { 
        free(f);
    }
    free(xg); free(yg); free(zg);
}
