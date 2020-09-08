#include <unordered_map>
#include <string>
#include <vector>
#include <math.h>
#include <cmath>        // std::abs
#include <limits>
#include <algorithm>    // std::min
#include "immersed.h"

void init_immersed(std::unordered_map<std::string,Polytree> & immersed_map, 
        std::vector<std::string> filenames) {
    for(auto &f: filenames) {
        Polytree ptree;
        polyhedron_from_file(ptree, f.c_str());
        immersed_map[f] = ptree; 
    }
    get_bb(immersed_map);
}

void get_bb(std::unordered_map<std::string,Polytree> & immersed_map) {

    double bmin_x, bmin_y, bmin_z, bmax_x, bmax_y, bmax_z;
    for(auto &i: immersed_map) {
        Polytree & p = i.second;
        polyhedron_bbox(p, &bmin_x, &bmin_y, &bmin_z, &bmax_x, &bmax_y, &bmax_z);
        p.boundingbox[0][0] = bmin_x;
        p.boundingbox[0][1] = bmax_x;
        p.boundingbox[1][0] = bmin_y;
        p.boundingbox[1][1] = bmax_y;
        p.boundingbox[2][0] = bmin_z;
        p.boundingbox[2][1] = bmax_z;
        std::cout << "bb ==> " << bmin_x << " " << bmin_y << " " << bmin_z << " " << bmax_x << " " << bmax_y << " " << bmax_z << std::endl;
    }

}

bool check_bb_interact(Geometry &geo, Polytree &ptree) {
    //TODO implement real bounding box interaction check
    bool interact;
    float gsx = geo.lengths[0][0]       ; float gex = geo.lengths[0][1];
    float gsy = geo.lengths[1][0]       ; float gey = geo.lengths[1][1];
    float gsz = geo.lengths[2][0]       ; float gez = geo.lengths[2][1];
    float psx = ptree.boundingbox[0][0] ; float pex = ptree.boundingbox[0][1];
    float psy = ptree.boundingbox[1][0] ; float pey = ptree.boundingbox[1][1];
    float psz = ptree.boundingbox[2][0] ; float pez = ptree.boundingbox[2][1];
    interact = true;
    return interact;
}

void allocate_distances(std::unordered_map<std::string,Polytree> & immersed_map,
        std::unordered_map<std::string,Field> & field_map) {
    for(auto &ifi: field_map) {
        Field & field = ifi.second;
        allocate_distances_block(immersed_map, field);
    }
}

void deallocate_distances_block(Field & field) {
    Geometry &geo = field.geo;
    int a_alloc_x = geo.intervals[0][2] + 2*geo.intervals[0][3];
    int a_alloc_y = geo.intervals[1][2] + 2*geo.intervals[1][3];
    int a_alloc_z = geo.intervals[2][2] + 2*geo.intervals[2][3];
    std::vector<int> np = {a_alloc_x, a_alloc_y, a_alloc_z};
    bool atleastoneimmersedbody = false;
    for(auto &iim: field.dist_map) {
        deallocate3d(&(iim.second), np);
        atleastoneimmersedbody = true;
    }
    if(atleastoneimmersedbody) deallocate3d(&(field.dist), np);

    field.triangles_map.clear(); //USELESS AUTOMATICALLY CALLED BY BLOCK DESTRUCTOR
}

void allocate_distances_block(std::unordered_map<std::string,Polytree> & immersed_map, Field & field) {
    Geometry &geo = field.geo;
    int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
    int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];
    float dx = geo.lengths[0][3]  ; float dy = geo.lengths[1][3]; float dz = geo.lengths[2][3];
    int a_alloc_x = geo.intervals[0][2] + 2*geo.intervals[0][3];
    int a_alloc_y = geo.intervals[1][2] + 2*geo.intervals[1][3];
    int a_alloc_z = geo.intervals[2][2] + 2*geo.intervals[2][3];
    std::vector<int> np = {a_alloc_x, a_alloc_y, a_alloc_z};
    bool atleastoneimmersedbody = false;
    for(auto &iim: immersed_map) {
        const std::string & imm_key = iim.first;
        Polytree & ptree = immersed_map[imm_key];
        bool boundingboxinteract = check_bb_interact(geo, ptree);
        if(boundingboxinteract) {
            field.dist_map[imm_key] = allocate3d(np);
            field.triangles_map[imm_key] = {};
            atleastoneimmersedbody = true;
        }
    }
    if(atleastoneimmersedbody) {
        field.dist = allocate3d(np);
    }
}

void compute_distances(std::unordered_map<std::string,Polytree> & immersed_map,
        std::unordered_map<std::string,Field> & field_map) {

    for(auto &ifi: field_map) {
        Field & field = ifi.second;
        std::cout << "block computing distances" << std::endl;
        compute_distances_block(immersed_map, field);
    }
}

void compute_distances_block(std::unordered_map<std::string,Polytree> & immersed_map, Field & field) {
    double near_x, near_y, near_z, query_x, query_y, query_z, distance, min_distance;
    Geometry &geo = field.geo;
    int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
    int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];
    float dx = geo.lengths[0][3]  ; float dy = geo.lengths[1][3]; float dz = geo.lengths[2][3];
    float sx  = geo.lengths[0][0]  ; float sy  = geo.lengths[1][0]  ; float sz  = geo.lengths[2][0]  ; 
    float ex  = geo.lengths[0][1]  ; float ey  = geo.lengths[1][1]  ; float ez  = geo.lengths[2][1]  ; 

    for(int i=0;i<2*ngx+nx;i++) {
        std::cout << "compute_distances_block first loop i=" << i << std::endl;
        query_x = (i-ngx)*dx+0.5*dx + sx;
        for(int j=0;j<2*ngy+ny;j++) {
            query_y = (j-ngy)*dy+0.5*dy + sy;
            for(int k=0;k<2*ngz+nz;k++) {
                query_z = (k-ngz)*dz+0.5*dz + sz;
                min_distance = std::numeric_limits<float>::max();
                for(auto &iim: field.dist_map) {
                    const std::string & imm_key = iim.first;
                    Polytree & ptree = immersed_map[imm_key];
                    polyhedron_closest(ptree, query_x, query_y, query_z, &near_x, &near_y, &near_z);
                    distance = sqrt((query_x-near_x)*(query_x-near_x)+
                            (query_y-near_y)*(query_y-near_y)+
                            (query_z-near_z)*(query_z-near_z));
                    if(polyhedron_inside(ptree, query_x, query_y, query_z)) distance = -distance;
                    min_distance = std::min(distance, min_distance);
                    iim.second[i][j][k] = distance;
                }
                field.dist[i][j][k] = min_distance;
            }
        }
    }

    std::cout << "compute_distances_block first part concluded" << std::endl;

    for(auto &iim: field.triangles_map) {
        const std::string & imm_key = iim.first;
        Polytree & ptree = immersed_map[imm_key];
        for(Facet_iterator f = ptree.poly->facets_begin(); f != ptree.poly->facets_end(); ++f) { 
            std::cout << "compute_distances_block second loop" << std::endl;
            Halfedge_facet_circulator j = f->facet_begin();
            //std::cout << CGAL::circulator_size(j) << std::endl;
            Point p1 = j->vertex()->point(); j++;
            Point p2 = j->vertex()->point(); j++;
            Point p3 = j->vertex()->point();
            float x1 = p1.x(); float x2 = p2.x(); float x3 = p3.x();
            float y1 = p1.y(); float y2 = p2.y(); float y3 = p3.y();
            float z1 = p1.z(); float z2 = p2.z(); float z3 = p3.z();
            std::array<float,3> center = {(float)((x1+x2+x3)/3.),
                                          (float)((y1+y2+y3)/3.),
                                          (float)((z1+z2+z3)/3.)};
            if(center[0] >= sx && center[0] < ex && 
               center[1] >= sy && center[1] < ey && 
               center[2] >= sz && center[2] < ez) {
                float ux = x2-x1;  float uy = y2-y1;  float uz = z2-z1;
                float vx = x3-x1;  float vy = y3-y1;  float vz = z3-z1;
                float nx = uy*vz-uz*vy; float ny = uz*vx-ux*vz; float nz = ux*vy-uy*vx;
                float nmod = sqrt(nx*nx+ny*ny+nz*nz);
                nx = nx/nmod; ny = ny/nmod; nz = nz/nmod;
                std::array<float,3> normal = {nx,ny,nz};
                //float area = 0.5*sqrt((x2*y3-x3*y2)*(x2*y3-x3*y2)+
                //                      (x3*y1-x1*y3)*(x3*y1-x1*y3)+
                //                      (x1*y2-x2*y1)*(x1*y2-x2*y1));
                float edge12 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
                float edge13 = sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+(z3-z1)*(z3-z1));
                float edge23 = sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2));
                float semip = 0.5*(edge12+edge13+edge23);
                float area = sqrt(semip*(semip-edge12)*(semip-edge13)*(semip-edge23));
                Triangle triangle;
                triangle.area = area;
                triangle.normal = normal;
                triangle.center = center;

                iim.second.push_back(triangle);
            }
            //std::cout << "faccia: " << f.size_of_vertexes() << std::endl;
        } 
    } 

    //for(auto &iim: field.triangles_map) {
    //    std::cout << "n triangles body=" << iim.first << " is=" << iim.second.size() << std::endl;
    //}
    //exit(EXIT_FAILURE);


}

void move_immersed(std::string type, std::array<float,3> & amount, 
    std::unordered_map<std::string,Polytree> &immersed_map, const std::string imm_key) {
    Polytree & ptree = immersed_map[imm_key];
    if(type == "translate") {
        polyhedron_move(ptree, amount);
    }
}

void polyhedron_move(Polytree &ptree, std::array<float,3> & amount) {
    const Vector_3 transvec(amount[0], amount[1], amount[2]);
    Aff_transformation_3 transl(CGAL::TRANSLATION, transvec);
    std::transform(ptree.poly->points_begin(), ptree.poly->points_end(), ptree.poly->points_begin(), transl);
    delete ptree.tree;
    Tree *tree = new Tree(faces(*(ptree.poly)).first, faces(*(ptree.poly)).second, *(ptree.poly));
    tree->accelerate_distance_queries();
    ptree.tree = tree;
}

void polyhedron_from_file (Polytree  &ptree, const char *fname){

    Polyhedron *P = new Polyhedron;

    std::ifstream input(fname, std::ios::in | std::ios::binary);

    //-----------------------------------------------------------
    // START - STL read
    //-----------------------------------------------------------
    //CGAL::read_STL( input,
    //    points,
    //    faces,
    //    true);
    //std::cout.precision(17);
    //std::cout << "OFF\n" << points.size() << " " << faces.size()  << " 0" << std::endl;
    //for(std::size_t i=0; i < points.size(); i++){
    //    std::cout << points[i][0] << " " << points[i][1] << " " << points[i][2]<< std::endl;
    //}
    //for(std::size_t i=0; i < faces.size(); i++){
    //    std::cout << "3 " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << std::endl;
    //}
    //-----------------------------------------------------------
    // END - STL read
    //-----------------------------------------------------------

    std::ifstream stream(fname); 
    if(!stream) { 
        std::cerr << "Cannot open file!" << std::endl; 
    } 
    stream >> *P; 
    if(!stream) { 
        std::cerr << "test.off is not a polyhedron" << std::endl; 
    } 

    Tree *tree = new Tree(faces(*P).first, faces(*P).second, *P);
    tree->accelerate_distance_queries();

    std::cout << " facets: "    << P->size_of_facets()    << std::endl;
    std::cout << " halfedges: " << P->size_of_halfedges() << std::endl;
    std::cout << " vertices: "  << P->size_of_vertices()  << std::endl;

    ptree.poly = P;
    ptree.tree = tree;

}

void polyhedron_closest (Polytree &ptree, double query_x, double query_y, double query_z, 
        double *near_x, double *near_y, double *near_z){

    Point query_point(query_x,query_y,query_z);

    Point closest = ptree.tree->closest_point(query_point);

    *near_x = closest.x();
    *near_y = closest.y();
    *near_z = closest.z();
}

bool polyhedron_inside(Polytree &ptree, double query_x, double query_y, double query_z) {

    Point_inside inside_tester(*(ptree.tree));
    Point query_point = Point(query_x,query_y,query_z);

    // Determine the side and return true if inside!
    int ret=0;
    bool is_inside = (inside_tester(query_point) == CGAL::ON_BOUNDED_SIDE);
    if(is_inside) 
        ret = 1;
    return ret;
}

void polyhedron_bbox(Polytree &ptree, 
        double *bmin_x, double *bmin_y, double *bmin_z,
        double *bmax_x, double *bmax_y, double *bmax_z){
    Bbox_3 bbox = ptree.tree->bbox();
    *bmin_x = bbox.xmin();
    *bmin_y = bbox.ymin();
    *bmin_z = bbox.zmin();
    *bmax_x = bbox.xmax();
    *bmax_y = bbox.ymax();
    *bmax_z = bbox.zmax();
}

void polyhedron_finalize(Polytree &ptree){
    delete ptree.tree; ptree.tree = NULL;
    delete ptree.poly; ptree.poly = NULL;
}

void correct_advanced_immerse_block(float val, float ***f, float ***f_new, float ***dist, Geometry & geo) {
    int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
    int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];
    float dx = geo.lengths[0][3]  ; float dy = geo.lengths[1][3]; float dz = geo.lengths[2][3];
#pragma omp parallel for
    for(int i=0;i<2*ngx+nx;i++) 
        for(int j=0;j<2*ngy+ny;j++) 
#pragma omp simd
            for(int k=0;k<2*ngz+nz;k++) {
                if(dist[i][j][k] < 0.) {
                    f_new[i][j][k] = 2.*val-f[i][j][k];
                }
            }
#pragma omp parallel for
    for(int i=0;i<2*ngx+nx;i++) 
        for(int j=0;j<2*ngy+ny;j++) 
#pragma omp simd
            for(int k=0;k<2*ngz+nz;k++) {
                if(dist[i][j][k] < 0.) {
                    f[i][j][k] = f_new[i][j][k];
                }
            }
}

void evolve_advanced_immerse_block(float ***f, float ***f_new, float ***dist, Geometry & geo) {
    int nx  = geo.intervals[0][2]; int ny  = geo.intervals[1][2]; int nz  = geo.intervals[2][2];
    int ngx = geo.intervals[0][3]; int ngy = geo.intervals[1][3]; int ngz = geo.intervals[2][3];
    float dx = geo.lengths[0][3]  ; float dy = geo.lengths[1][3]; float dz = geo.lengths[2][3];
    float dist_x, dist_y, dist_z, f_x, f_y, f_z, mod_dist;

#pragma omp parallel for
    for(int i=ngx;i<ngx+nx;i++) 
        for(int j=ngy;j<ngy+ny;j++) 
#pragma omp simd
            for(int k=ngz;k<ngz+nz;k++) {
                if(dist[i][j][k] < 0.) {
                    dist_x = dist[i+1][j][k]-dist[i-1][j][k];
                    dist_y = dist[i][j+1][k]-dist[i][j-1][k];
                    dist_z = dist[i][j][k+1]-dist[i][j][k-1];
                    //mod_dist = sqrt(dist_x*dist_x+dist_y*dist_y+dist_z*dist_z);
                    mod_dist = std::abs(dist_x)+std::abs(dist_y)+std::abs(dist_z);
                    dist_x = dist_x/mod_dist;
                    dist_y = dist_y/mod_dist;
                    dist_z = dist_z/mod_dist;
                    if(dist_x > 0) {
                        f_x = f[i+1][j][k]-f[i][j][k];
                    } else {
                        f_x = f[i][j][k]-f[i-1][j][k];
                    }
                    if(dist_y > 0) {
                        f_y = f[i][j+1][k]-f[i][j][k];
                    } else {
                        f_y = f[i][j][k]-f[i][j-1][k];
                    }
                    if(dist_z > 0) {
                        f_z = f[i][j][k+1]-f[i][j][k];
                    } else {
                        f_z = f[i][j][k]-f[i][j][k-1];
                    }
                    f_new[i][j][k] = f[i][j][k] + 0.9*(dist_x*f_x+dist_y*f_y+dist_z*f_z);
                    //f_new[i][j][k] = f[i][j][k] + 0.3*(dist_x*f_x+dist_y*f_y+dist_z*f_z);
                    //f_new[i][j][k] = f[i][j][k] + 0.002*(dist_x*f_x+dist_y*f_y+dist_z*f_z);
                }
            }

#pragma omp parallel for
    for(int i=ngx;i<ngx+nx;i++) 
        for(int j=ngy;j<ngy+ny;j++) 
#pragma omp simd
            for(int k=ngz;k<ngz+nz;k++) {
                if(dist[i][j][k] < 0.) {
                    f[i][j][k] = f_new[i][j][k];
                }
    }
}

void compute_forces(int it, std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map, 
        std::unordered_map<std::string,Polytree> &immersed_map, Params & params) {

    for(auto &iim: immersed_map) {
        iim.second.forces = {0.,0.,0.,0.,0.,0.};
    }

    for(auto &n: field_map) {
        Field & f = n.second;
        Geometry & geo = geo_map[n.first];
        for(auto &tm: f.triangles_map) {
            const std::string & body = tm.first;
            std::vector<float> & imm_forces = immersed_map[body].forces;
            for(auto &tri: tm.second) {
                std::vector<float> forces = compute_force_field(tri, f, geo);
                //TESTstd::cout << "triangle forces=(" << 
                //TEST    forces[0] << ";" << forces[1] << ";" << forces[2] << ";" << 
                //TEST    forces[3] << ";" << forces[4] << ";" << forces[5] << ")" << std::endl;
                int iforce=0;
                for(auto &force: forces) {
                    imm_forces[iforce] += force;
                    iforce++;
                }
            }
        }
    }

    std::ofstream forces_file;
    forces_file.open("forces.dat", std::ios_base::app);
    for(auto &iim: immersed_map) {
        std::vector<float> & ff = iim.second.forces;
        //forces_file << "body=" << iim.first << " - forces=(" << 
        forces_file << it << " " << ff[0] << " " << ff[1] << " " << ff[2] << " " << 
            ff[3] << " " << ff[4] << " " << ff[5] << std::endl;
    }
    forces_file.close();
    
}
