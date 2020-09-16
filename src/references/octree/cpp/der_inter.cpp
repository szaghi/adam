#include "der_inter.h"

void der1(float *** v_in, float *** v_out, int dir, float delta[3], int nt[3], int ng[3]) {
    int nx = nt[0]  ; int ny = nt[1]  ; int nz = nt[2] ;
    int ngx = ng[0] ; int ngy = ng[1] ; int ngz = ng[2];
    if(dir == 1) {
        float invdelta = 1./(2.*delta[0]);
#pragma omp parallel for collapse(2)
        for(int i=ngx;i<ngx+nx;i++) 
            for(int j=ngy;j<ngy+ny;j++) 
                for(int k=ngz;k<ngz+nz;k++)
                    v_out[i][j][k] = invdelta * (v_in[i+1][j][k] - v_in[i-1][j][k]);
    } else if(dir == 2) {
        float invdelta = 1./(2.*delta[1]);
#pragma omp parallel for collapse(2)
        for(int i=ngx;i<ngx+nx;i++) 
            for(int j=ngy;j<ngy+ny;j++) 
                for(int k=ngz;k<ngz+nz;k++)
                    v_out[i][j][k] = invdelta * (v_in[i][j+1][k] - v_in[i][j-1][k]);
    } else if(dir == 3) {
        float invdelta = 1./(2.*delta[2]);
#pragma omp parallel for collapse(2)
        for(int i=ngx;i<ngx+nx;i++) 
            for(int j=ngy;j<ngy+ny;j++) 
                for(int k=ngz;k<ngz+nz;k++)
                    v_out[i][j][k] = invdelta * (v_in[i][j][k+1] - v_in[i][j][k-1]);
    }
}

void der2(float *** v_in, float *** v_out, int dir, float *delta, int *nt, int *ng) {
    int nx = nt[0]  ; int ny = nt[1]  ; int nz = nt[2] ;
    int ngx = ng[0] ; int ngy = ng[1] ; int ngz = ng[2];
    if(dir == 1) {
        float invdelta2 = 1./(delta[0]*delta[0]);
#pragma omp parallel for collapse(2)
        for(int i=ngx;i<ngx+nx;i++) 
            for(int j=ngy;j<ngy+ny;j++) 
                for(int k=ngz;k<ngz+nz;k++)
                    v_out[i][j][k] = invdelta2 * (v_in[i+1][j][k] + v_in[i-1][j][k] - 2.*v_in[i][j][k]);
    } else if(dir == 2) {
        float invdelta2 = 1./(delta[1]*delta[1]);
#pragma omp parallel for collapse(2)
        for(int i=ngx;i<ngx+nx;i++) 
            for(int j=ngy;j<ngy+ny;j++) 
                for(int k=ngz;k<ngz+nz;k++)
                    v_out[i][j][k] = invdelta2 * (v_in[i][j+1][k] + v_in[i][j-1][k] - 2.*v_in[i][j][k]);
    } else if(dir == 3) {
        float invdelta2 = 1./(delta[2]*delta[2]);
#pragma omp parallel for collapse(2)
        for(int i=ngx;i<ngx+nx;i++) 
            for(int j=ngy;j<ngy+ny;j++) 
                for(int k=ngz;k<ngz+nz;k++)
                    v_out[i][j][k] = invdelta2 * (v_in[i][j][k+1] + v_in[i][j][k-1] - 2.*v_in[i][j][k]);
    }
}

void extrapolate(float *** v_inout, float delta[3], int nt[3], int ng[3]) {
    int nx = nt[0]  ; int ny = nt[1]  ; int nz = nt[2] ;
    int ngx = ng[0] ; int ngy = ng[1] ; int ngz = ng[2];
    int i, j, k;
#pragma omp parallel for private(k,i) collapse(1) 
    for(j=ngy;j<ngy+ny;j++) 
        for(k=ngz;k<ngz+nz;k++) 
            for(i=0;i<ngx;i++) {
                v_inout[i][j][k] = v_inout[ngx][j][k];
                v_inout[nx+ngx+i][j][k] = v_inout[nx+ngx-1][j][k];
            }
#pragma omp parallel for private(k,j) collapse(1)
    for(i=ngx;i<ngx+nx;i++) 
        for(k=ngz;k<ngz+nz;k++) 
            for(j=0;j<ngy;j++) {
                v_inout[i][j][k] = v_inout[i][ngy][k];
                v_inout[i][ny+ngy+j][k] = v_inout[i][ny+ngy-1][k];
            }
#pragma omp parallel for private(j,k) collapse(1)
    for(i=ngx;i<ngx+nx;i++) 
        for(j=ngy;j<ngy+ny;j++) 
            for(k=0;k<ngz;k++) {
                v_inout[i][j][k] = v_inout[i][j][ngz];
                v_inout[i][j][nz+ngz+k] = v_inout[i][j][nz+ngz-1];
            }
}
