#ifndef IMMERSED_H
#define IMMERSED_H

#include <unordered_map>
#include <fstream>
#include <iostream>
#include <vector>
#include <array>
#include <string>

// the next is special: it contains implementation and therefore has to be included once
// #include <CGAL/IO/STL_reader.h>
// the previous is special: it contains implementation and therefore has to be included once

//#include <CGAL/array.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h> // not found in cineca opt installation

#include "geometry.h"
#include "field.h"

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef CGAL::Side_of_triangle_mesh<Polyhedron, K> Point_inside;
typedef CGAL::Bbox_3 Bbox_3;

// see https://stackoverflow.com/questions/30039937/syntax-error-for-cgal-affine-transformation
typedef K::Vector_3 Vector_3;
typedef CGAL::Aff_transformation_3<K> Aff_transformation_3;

// for iterating over triangles
typedef typename Polyhedron::Facet_iterator Facet_iterator; 
typedef Polyhedron::Halfedge_around_facet_circulator Halfedge_facet_circulator;

typedef struct {
    Polyhedron *poly; 
    Tree *tree;
    std::array<std::array<int,2>,3> boundingbox;
    std::vector<float> forces;
} Polytree;

void init_immersed(std::unordered_map<std::string,Polytree> & immersed_map, 
    std::vector<std::string> filenames); 

void get_bb(std::unordered_map<std::string,Polytree> & immersed_map);

bool check_bb_interact(Geometry &geo, Polytree &ptree);

void allocate_distances(std::unordered_map<std::string,Polytree> & immersed_map,
    std::unordered_map<std::string,Field> & field_map);

void allocate_distances_block(std::unordered_map<std::string,Polytree> & immersed_map, Field & field);

void deallocate_distances_block(Field & field);

void compute_distances(std::unordered_map<std::string,Polytree> & immersed_map,
    std::unordered_map<std::string,Field> & field_map);

void compute_distances_block(std::unordered_map<std::string,Polytree> & immersed_map, Field & field);

void move_immersed(std::string type, std::array<float,3> & amount, 
    std::unordered_map<std::string,Polytree> &immersed_map, const std::string imm_key);

void polyhedron_move(Polytree &ptree, std::array<float,3> & amount);

void polyhedron_from_file (Polytree  &ptree, const char *fname);

void polyhedron_closest (Polytree &ptree, double query_x, double query_y, double query_z, 
        double *near_x, double *near_y, double *near_z);

bool polyhedron_inside(Polytree &ptree, double query_x, double query_y, double query_z);

void polyhedron_bbox(Polytree &ptree, double *bmin_x, double *bmin_y, double *bmin_z, double *bmax_x, double *bmax_y, double *bmax_z);

void polyhedron_finalize(Polytree &ptree);

void evolve_advanced_immerse_block(float ***f, float ***f_new, float ***dist, Geometry & geo);

void correct_advanced_immerse_block(float val, float ***f, float ***f_new, float ***dist, Geometry & geo);

void compute_forces(int it, std::unordered_map<std::string,Geometry> & geo_map,
        std::unordered_map<std::string,Field> & field_map, 
        std::unordered_map<std::string,Polytree> &immersed_map, Params & params);

#endif
