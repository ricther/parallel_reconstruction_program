#ifndef _CCONTOURARRANGEMENT_H_
#define _CCONTOURARRANGEMENT_H_
#include "CLayer.h"
#include "CContour.h"
#include <vector>
#include "initial.h"
#include "CVirtualContour.h"
class CShape;
class contour_couple
{
public:
  contour_couple();
  CVirtualContour* vcontour1;
  CVirtualContour* vcontour2;
  bool b_intersection;
  float distance;
  bool b_couple;
  float boundary_coverage_rate;
  bool operator == (const contour_couple& right);

};

class CContourArrangement
{
public:
CContourArrangement(CShape* shape):m_shape(shape){};

  CShape* m_shape;
  std::vector<contour_couple*> vec_contour_couple;
  void setup();
  void setup_use_omp();
  float check_intersection(CContour*,CContour*);
  float cal_distance(CContour*,CContour*);
  float check_boundary_box_coverage(CContour* u, CContour*d);
  float check_boundary_box_distance(float num1,float num2,float num3, float num4);
  void get_couple(CLayer* first_layer, CLayer* second_layer,  std::vector<contour_couple*>& vec_couple );
  float check_intersection_use_level_set(CContour* u, CContour* d);
  float check_coverage_level_set(CContour* u, CContour*d);
  float check_intersection_level_set(CContour* u, CContour*d);
  void get_union_couple(  std::vector<contour_couple*>&,  std::vector<contour_couple*>&);
};

#endif /* _CCONTOURARRANGEMENT_H_ */
