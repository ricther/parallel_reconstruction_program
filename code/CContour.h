#pragma once
//athour xu liu 07/10/13
//the class of contour
#include <vector>
#include <fstream>
#include <string>
#include "CMedialMap.h"
#include "CMap.h"
#include "initial.h"
class CPoint;
class CLayer;


class CContour
{

 public:
  CContour(const float ID,CLayer*);
  CLayer* m_layer;
  std::string filename;
  float LayerID;
  float contourID;
  int PointNum;//the number of points in the layer
  std::vector<CPoint*> vec_points;// this is the primitive points;
  std::vector<CPoint*> vec_points_orderd;// the point after check the sequnce
  std::vector<CPoint*> vec_Points_Origin;// the points after normalize;
  std::vector<CPoint*> vec_new_points;//change after bspline_update
  std::vector<CPoint*> vec_Points_Inter;
  std::vector<CPoint*> vec_Points_Vicinity;
  std::vector<CPoint*> vec_new_points_vicinity;//change after bspline_update
  std::vector<CPoint*> vec_Points_project;//points project from adjacent layer's medial_axis;
  void operator=(CContour &temp);
  bool read_single_layer_without_z(std::fstream&);
  bool read_contour_with_z(std::fstream& fin);
  CMap* m_Map;
  CMedialMap* m_medial_map;
  CMap* m_temp_map;
  double** lattice_x,**lattice_y;//as NXB NYB will be initialize in CRegister
  double** new_lattice_x,**new_lattice_y;//as NNXB NNYB the next control point coordinates position
  void normalize(CPoint); 
  void InitMap();
  std::vector<double> vec_intensity_old;
  std::vector<double> vec_intensity_new;
  CPoint* center_point;
  CPoint* moment_one_point;
  double max_x,max_y,min_x,min_y;
  double length;
  void reset();
  int use_as_higher_contour_count,use_as_lower_contour_count;
  void smooth();
  void calculate_medial_map(double**);
  void swap_map_medialmap();
  void swap_medialmap_map();
private:
  void check_edge(double,double);
  void get_center();
  void calculate_one_moment();
  void check_point_seq();
  float sum_x,sum_y;
  double last_x,last_y;
  CPoint* get_the_nearest_unused_point(CPoint* point);
};
