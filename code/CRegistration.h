#pragma once
#include <vector>
#include "CContour.h"
#include "initial.h"

//function of freeform code

class CPoint;
class CShape;
class CLayer;
class CCorrespond;
class CVirtualContour;

class ParamRecord
{
  public:
  ParamRecord()
  {
    lattice_x=0;
    lattice_y=0;
    errorE=0;
  }
  float **lattice_x;
  float **lattice_y;
  float errorE;
  void operator=(const ParamRecord &);
  void initial(float errorE,CVirtualContour* higher);
};

class CRegistration
{
  public:
  CRegistration(CShape*);
  CShape *SourceShape;
  void freeform_res1(CVirtualContour*,CVirtualContour*);//first is lower_layer equal target, second is higher_layer equal source
  int kNumberOfIteration;
  void Register();
  void Register_use_openmp();
  private:
  void reset(CVirtualContour*);
  void shape_vicinity(CVirtualContour*,int);//first is higher_layer as source,second is the narrow band width
  int narrow_band;
  int vicinity_points_num;
  //compute the intensity of the sample points in distance map
  void compute_intensity(std::vector<CPoint*>&,float**&,std::vector<float>&);
  //compute the intensity of single point in distance map
  float compute_intensity_by_point(float**& dist,float tx, float ty);
  //the grill in distance map
  void init_lattice(CVirtualContour*);
  
  //bspline_update
  void bspline_update(CVirtualContour*,int mode,std::vector<CPoint*>&,std::vector<CPoint*>&);// mode = 0,use the normal lattice, mode = 1 ,use the new lattice.p
  //compute energy
  float energy_func_square_diff(CVirtualContour*);
  //drive the eqution
  void gradientdescent_smoother(CVirtualContour* lower,CVirtualContour*higher,float errorE);
  //calculate the lattice
  void calculate_dlattice_by_point(float**&,float**&,CVirtualContour*,CVirtualContour*);
  
  bool update_lattice(float**&,float**&m,CVirtualContour*,CVirtualContour*,float&);
  ParamRecord current_params;
  float cubic_spline(float u,int o);
  void get_correspondence(CCorrespond*,CVirtualContour*higher,CVirtualContour*lower);
  int get_closest_point(std::vector<CPoint*>&,CPoint);
  int get_closest_point(std::vector<CPoint*>&vec_points,std::vector<CPoint*>&medial_points,CPoint,bool&);
  //fill the hole on the boundary
  void fill_the_hole(CCorrespond* corres,int& count,int first,int size,int last_index,int index, CVirtualContour *lower,CPoint* point2,int gap,bool isforeward);
  void fill_the_hole(CCorrespond* corres,int& count,int first,int size,int last_index,int index, CPoint* point1, CVirtualContour *higher,int gap,bool isforeward);
  int get_relate_index(int first, int size, int old_index);
  void get_gap(int first, int size, int index, int last_index,int &real_gap,bool& isforeward );

  float* xvb;//the lattice's x value
  float* yvb;//the lattice's y value


  float kappa;//weight factor controlling the smoothness term

  const int NumberRows,NumberCols;
  inline float max(float* val,int size)
  {
    float max = val[0];
    for (int i = 0; i < size; ++i)
    {
      if (val[i]>max)
      {
        max=val[i];
      }
    }
    return max;
  }
  inline float min(float*val,int size)
  {
    float min =val[0];
    for (int i = 0; i < size; ++i)
    {
      if (val[i]<min)
      {
        min=val[i];
      }
    }
    return min;
  }
  inline float dot(float*v1,float*v2)
  {
    return v1[0]*v2[0]+v1[1]*v2[1];
  }
  inline float spline_deriv(float u, int o)
  {
    
	float b;
	switch (o)
	{
	case 0:
		b = (float) 1.0/ (float) 6.0*3.0*(1.0-u)*(1.0-u)*(-1);
		break;
	case 1:
		b = (float) 1.0/(float) 6.0*(3.0*3.0*u*u-6.0*2.0*u);
		break;
	case 2:
		b = (float) 1.0/(float) 6.0*(-3.0*3.0*u*u+3.0*2.0*u+3.0);
		break;
	case 3:
		b = (float)1.0/(float) 6.0*3.0*u*u;
		break;
	}
	return b;
  }
  //find the nearest contour except used
  CVirtualContour* find_nearest_contour(CVirtualContour*,CLayer*,int);
  bool check_contour_distance(CVirtualContour* lower_contour,CVirtualContour* higher_contour);
  void regist_lower_long_higher_short(CVirtualContour*,CVirtualContour*);

};






