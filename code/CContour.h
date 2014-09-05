#pragma once
//athour xu liu 07/10/13
//the class of contour
#include <vector>
#include <fstream>
#include <string>
#include "CMedialMap.h"
#include "CMap.h"
#include "initial.h"
#include <map>

class CVirtualContour;
class CPoint;
class CLayer;
class CContourArrangement;

class CContour
{

  friend class CVirtualContour;
  friend class CContourArrangement;
 public:
  static float static_max_x;
  static float static_min_x;
  static float static_max_y;
  static float static_min_y;
  CContour(const float ID,CLayer*);
  int contour_index;
  CLayer* m_layer;
  std::string filename;
  float LayerID;
  float contourID;
  int PointNum;//the number of points in the layer
  std::vector<CPoint*> vec_points;// this is the primitive points;
  std::vector<CPoint*> vec_points_orderd;// the point after check the sequnce
  std::vector<CPoint*> vec_Points_Origin;// the points after normalize;
  std::vector<CPoint*> vec_Points_Inter; //just used in CMap. can clear after calculate Distancemap
  std::vector<CVirtualContour*> vec_virtual_contour;
  void operator=(CContour &temp);
  bool read_single_layer_without_z(std::fstream&);
  bool read_contour_with_z(std::fstream& fin);

  //  CMedialMap* m_medial_map;
  //  CMap* m_temp_map;
  void normalize(CPoint); 
  void InitMap();
  CPoint* center_point;
  CPoint* moment_one_point;
  float max_x,max_y,min_x,min_y, original_max_x,original_max_y,original_min_x,original_min_y;
  float length,original_length;         /**< the length is the contour length after scale, the original_length is the original length from the data */
  void reset();
  //  int use_as_higher_contour_count,use_as_lower_contour_count;
  void smooth();
  /* void calculate_medial_map(float**); */
  /* void swap_map_medialmap(); */
  /* void swap_medialmap_map(); */
  std::multimap<float,CContour*> map_neighbour; /**< record the contour intersection with this contour, which in the neighbour layer ,the map_neighbour must in the CContour class, can't put in the CVirtualContour class, because in the CVirtualContour class only record the information for one pair contour, can't record the neighbour*/
  int use_counter;
private:
  void check_edge(float,float);
  void check_edge_before_scale(float,float);
  void get_center();
  void calculate_one_moment();
  void check_point_seq();
  float sum_x,sum_y;
  float last_x,last_y,original_last_x,original_last_y;
  CPoint* get_the_nearest_unused_point(CPoint* point);
  static int static_contour_index;             /**< the unique id for contour */
  CMap* m_Map;
};
