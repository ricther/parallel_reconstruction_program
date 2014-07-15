#pragma once
#include <vector>
#include "CContour.h"
#include "initial.h"
//function of freeform code

class CPoint;
class CShape;
class CLayer;
class CCorrespond;

class ParamRecord
{
  public:
  ParamRecord()
  {
    lattice_x=0;
    lattice_y=0;
    errorE=0;
  }
  double **lattice_x;
  double **lattice_y;
  double errorE;
  void operator=(const ParamRecord &);
  void initial(double errorE,CContour* higher);
};

class CRegistration
{
  public:
  CRegistration(CShape*);
  CShape *SourceShape;
  void freeform_res1(CContour*,CContour*);//first is lower_layer equal target, second is higher_layer equal source
  int kNumberOfIteration;
  void Register();
  void Register_use_openmp();
  private:
  void reset();
  void shape_vicinity(CContour*,int);//first is higher_layer as source,second is the narrow band width
  int narrow_band;
  int vicinity_points_num;
  //compute the intensity of the sample points in distance map
  void compute_intensity(std::vector<CPoint*>&,double**&,std::vector<double>&);
  //compute the intensity of single point in distance map
  double compute_intensity_by_point(double**& dist,double tx, double ty);
  //the grill in distance map
  void init_lattice(CContour*);
  
  double* xvb;//the lattice's x value
  double* yvb;//the lattice's y value
  //bspline_update
  void bspline_update(CContour*,int mode,std::vector<CPoint*>&,std::vector<CPoint*>&);// mode = 0,use the normal lattice, mode = 1 ,use the new lattice.p
  //compute energy
  double energy_func_square_diff(CContour*);
  //drive the eqution
  void gradientdescent_smoother(CContour* lower,CContour*higher,double errorE);
  //calculate the lattice
  void calculate_dlattice_by_point(double**&,double**&,CContour*,CContour*);
  
  bool update_lattice(double**&,double**&m,CContour*,CContour*,double&);
  ParamRecord current_params;
  double cubic_spline(double u,int o);
  void get_correspondence(CCorrespond*,CContour*higher,CContour*lower);
  int get_closest_point(std::vector<CPoint*>&,CPoint);
  int get_closest_point(std::vector<CPoint*>&vec_points,std::vector<CPoint*>&medial_points,CPoint,bool&);
  //fill the hole on the boundary
  void fill_the_hole(CCorrespond* corres,int& count,int first,int size,int last_index,int index, CContour *lower,CPoint* point2,int gap,bool isforeward);
  void fill_the_hole(CCorrespond* corres,int& count,int first,int size,int last_index,int index, CPoint* point1, CContour *higher,int gap,bool isforeward);
  int get_relate_index(int first, int size, int old_index);
  void get_gap(int first, int size, int index, int last_index,int &real_gap,bool& isforeward );
  double** XB;// original control point coordinates;
  double** YB;
  double** dXB;
  double** dYB;

  float lamda;// step size
  float kappa;//weight factor controlling the smoothness term

  const int NumberRows,NumberCols;
  inline double max(double* val,int size)
  {
    double max = val[0];
    for (int i = 0; i < size; ++i)
    {
      if (val[i]>max)
      {
        max=val[i];
      }
    }
    return max;
  }
  inline double min(double*val,int size)
  {
    double min =val[0];
    for (int i = 0; i < size; ++i)
    {
      if (val[i]<min)
      {
        min=val[i];
      }
    }
    return min;
  }
  inline double dot(double*v1,double*v2)
  {
    return v1[0]*v2[0]+v1[1]*v2[1];
  }
  inline double spline_deriv(double u, int o)
  {
    
	double b;
	switch (o)
	{
	case 0:
		b = (double) 1.0/ (double) 6.0*3.0*(1.0-u)*(1.0-u)*(-1);
		break;
	case 1:
		b = (double) 1.0/(double) 6.0*(3.0*3.0*u*u-6.0*2.0*u);
		break;
	case 2:
		b = (double) 1.0/(double) 6.0*(-3.0*3.0*u*u+3.0*2.0*u+3.0);
		break;
	case 3:
		b = (double)1.0/(double) 6.0*3.0*u*u;
		break;
	}
	return b;
  }
  //find the nearest contour except used
  CContour* find_nearest_contour(CContour*,CLayer*,int);
  bool check_contour_distance(CContour* lower_contour,CContour* higher_contour);
  void regist_lower_long_higher_short(CContour*,CContour*);

};






