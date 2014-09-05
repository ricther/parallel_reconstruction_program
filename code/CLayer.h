#pragma once
//athour xu liu 07/10/13
//the class of layer
#include <vector>
#include <fstream>
#include <string>
#include <map>
#include "initial.h"
class CMap;
class CPoint;
class CContour;

class CLayer
{
 public:
  CLayer(const float ID);
  std::string filename;
  float LayerID;
  int contour_Num;
  
  void operator=(CLayer &temp);
  bool read_layer_multi_contour_with_z(std::fstream& fin);
  void setup(CPoint);




  CPoint* center_point;
  CPoint* moment_one_point;
  float max_x,max_y,min_x,min_y, original_max_x,original_max_y,original_min_x,original_min_y;

  void reset();
  std::map<int,int> map_contourID;//contour_Num,contourID
  std::map<int,CContour*> map_contour;//contourID,CContour*
  float total_length;
private:
  void check_edge(float,float);
  void check_edge_before_scale(float tempx,float tempy);
  void get_center();
  void check_length(float);
  void calculate_one_moment();
  float sum_x,sum_y;

};
