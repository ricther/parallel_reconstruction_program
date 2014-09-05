/**
 * @file   CVirtualContour.h
 * @author xuliu <xl@animation-Precision-T7610>
 * @date   Fri Aug  8 08:05:41 2014
 * 
 * @brief  used for the free form algorithm, it's the virtual contour , instead of the CContour. Sinc *e one contour may used serval times in performing the FFD. but one CContour can't store serval data *for the serval FFD. So use virtual contour to instead of the CContour.
 * 
 * 
 */

#ifndef _CVIRTUALCONTOUR_H_
#define _CVIRTUALCONTOUR_H_

#include "CContour.h"
#include "Memory.h"

class CPoint;
class CVirtualContour
{
public:
  CVirtualContour(CContour* father_contour)
  {
    m_contour=father_contour;
    LayerID=m_contour->LayerID;
    m_contour->vec_virtual_contour.push_back(this);
    lattice_x=NULL;lattice_y=NULL;new_lattice_x=NULL;new_lattice_y=NULL;
    XB=NULL;YB=NULL;dXB=NULL;dYB=NULL;
    lamda=0.02;
    medial_axis=NULL;
    medial_axis_count=0;
    m_medial_map=NULL;
    m_Map=NULL;
    m_temp_map=NULL;
    m_Map=m_contour->m_Map;
    vec_medial_points.clear();
  }
  ~CVirtualContour()
  {
     free_2D_float_array(XB); //can't release XB, scince the XB may used in the incremental FFD.
     free_2D_float_array(YB);
  }
  void release_array()
  {
    free_2D_float_array(lattice_x);
    free_2D_float_array(lattice_y);
    free_2D_float_array(new_lattice_x);
    free_2D_float_array(new_lattice_y);
    //    free_2D_float_array(XB); //can't release XB, scince the XB may used in the incremental FFD.
    //    free_2D_float_array(YB);
    free_2D_float_array(dXB);
    free_2D_float_array(dYB);
    vec_Points_Vicinity.clear();
    vec_Points_Vicinity.resize(0);
    vec_new_points_vicinity.clear();
    vec_new_points_vicinity.resize(0);
    vec_intensity_old.clear();
    vec_intensity_old.resize(0);
    vec_intensity_new.clear();
    vec_intensity_new.resize(0);
  }
  CContour* m_contour;
  float LayerID;
  std::vector<CPoint*> vec_new_points;  /**< the final points positon, record the CPoints position after FFD ;used in FFD and Display*/
  std::vector<CPoint*> vec_Points_Vicinity; /**< the narrow band for ffd. */
  std::vector<CPoint*> vec_new_points_vicinity; /**< the narrow band for ffd, after update */
  std::vector<CPoint*> vec_Points_project;/**<points project from adjacent layer's medial_axis*/
  std::vector<float> vec_intensity_old; /**< temp intensity for FFD */
  std::vector<float> vec_intensity_new;/**< temp intensity for FFD */

  float** lattice_x,**lattice_y;/** the pre-step lattice position ,as NXB NYB will be initialize in CRegister;only used in FFD*/
  float** new_lattice_x,**new_lattice_y;/** next-step lattice position ,as NNXB NNYB the next control point coordinates position;only used in FFD*/
  float** XB;/**the original lattice position of every-FFD, when used in incremental FFD, XB will differed for every FFD, each times the XB record the last time's lattice's position */
  float** YB;/**like XB*/
  float** dXB;                          /**< record the value of  delta lattice, use the dXB to change the lattice, when parallel , one contour may processed in many thread, so transfer dXB from CRegistration to here, every virtual contour has a copy.   */
  float** dYB;                          /**< like dXB */
  float lamda;                          /**< control the step in gradient descent!, each contour has his own step, so must put lamda in virtual contour */


  float** medial_axis;
  int medial_axis_count;
  void calculate_medial_axis(float layerID);
  std::vector<CPoint*> vec_medial_points;
  void calculate_medial_map(float**);
  void swap_map_medialmap();
  void swap_medialmap_map();
  CMedialMap* m_medial_map;
  CMap* m_Map;
  CMap* m_temp_map;
};


#endif /* _CVIRTUALCONTOUR_H_ */
