#pragma once

#include <string>
#include <vtkSmartPointer.h>
#include <vtkContourWidget.h>
#include <vtkOrientedGlyphContourRepresentation.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkRenderer.h>
#include <vtkImagePlaneWidget.h>
#include <map>
#include "CSkeletonDisplay.h"
#include "CContourDisplay.h"
#include "CContourSetDisplay.h"
#include "COriginContour.h"
#include "CNormalContour.h"
#include "vtkTextActor.h"
#include "CCorespondDisplay.h"
#include "CTransformedContour.h"
#include "CTriangleDisplay.h"
#include "CContourSurface.h"
class CLayer;
class CPoint;
class CShape;

extern  const int MatrixRes;
extern  const int NumRows;
extern  const int NumCols;
class CShapeDisplay
{
public:
CShapeDisplay(CShape *,vtkSmartPointer<vtkRenderWindow> temp_window,vtkSmartPointer<vtkRenderWindowInteractor> temp_interactor);


private:
  CShape* m_shape;
  CSkeletonDisplay* m_skeleton;
  CContourDisplay* m_skeleton_contour;
  CCorespondDisplay* m_correspond;
  CTriangleDisplay* m_triangle;
  CContourSurface* m_contour_surface;
  
  void draw_medial_axis(vtkSmartPointer<vtkRenderer>);
  void draw_lattice(vtkSmartPointer<vtkRenderer>);
  void draw_lattice_content(vtkSmartPointer<vtkRenderer>,CContour*,int);
      
  CContourSetDisplay<COriginContour>* m_origin_contours;
  CContourSetDisplay<CNormalContour>* m_normal_contours;
  CContourSetDisplay<CTransformedContour>* m_transformed_contours;
  
  vtkSmartPointer<vtkRenderWindow> m_render_window;
  vtkSmartPointer<vtkRenderWindowInteractor> m_interactor;
  std::map<int,vtkSmartPointer<vtkContourWidget> > map_Contour;
  int count_contour;
public:
  void Draw();
  void draw_origin_points_widget(vtkSmartPointer<vtkRenderWindowInteractor>,vtkSmartPointer<vtkRenderer>);
  void draw_transformed_points_widget(vtkSmartPointer<vtkRenderWindowInteractor>,vtkSmartPointer<vtkRenderer>);
  void draw_compared_pionts_widget(vtkSmartPointer<vtkRenderWindowInteractor>,vtkSmartPointer<vtkRenderer>);

  void draw_origin_points(vtkSmartPointer<vtkRenderer>);
  void draw_transformed_points(vtkSmartPointer<vtkRenderer>);
  void draw_compared_points(vtkSmartPointer<vtkRenderer>);
  void draw_correspond( vtkSmartPointer<vtkRenderer>);
  void draw_normal_points(vtkSmartPointer<vtkRenderer>);
  void draw_skeleton(vtkSmartPointer<vtkRenderer>,vtkSmartPointer<vtkRenderWindowInteractor>);
  void draw_triangle(vtkSmartPointer<vtkRenderer>,vtkSmartPointer<vtkRenderWindowInteractor>);

  float current_LayerID;
  int current_lineID;
  int current_correspondID;
  void key_down();
  void key_up();
  void key_minus();
  void key_equal();
  void key_r();
  void key_f();
  void key_bracketleft();
  void key_bracketright();
  void key_change_viewport(std::string);
  void key_o();
  void key_p();

  
private:
    void init_contour_represetation_widget(vtkSmartPointer<vtkRenderWindowInteractor> interactor,vtkSmartPointer<vtkRenderer>,std::vector<CPoint*>,int scale);
  void init_contour_represetation(vtkSmartPointer<vtkRenderer>,std::vector<CPoint*>,int scale,double level);
  vtkSmartPointer<vtkActor> initialActor(std::vector<CPoint*>::iterator,std::vector<CPoint*>::iterator,int,double);
  

 
  
  std::map<float,CLayer*>::iterator itr,etr;
  void draw_text(vtkSmartPointer<vtkRenderer>,std::string);
  void draw_layerID(vtkSmartPointer<vtkRenderer> renderer ,float layerID);
  void draw_contour_for_skeleton(int layerID,vtkSmartPointer<vtkRenderer>);
  void update_layer_text(float layerID);
  void update_viewport(std::string key);
  vtkSmartPointer<vtkTextActor> m_layer_text_actor;
  vtkSmartPointer<vtkActor> m_contour_actor;
  vtkSmartPointer<vtkPolyData> m_contour_PolyData;
  
  vtkSmartPointer<vtkImagePlaneWidget> draw_plane(vtkSmartPointer<vtkRenderer>,vtkSmartPointer<vtkRenderWindowInteractor>);
  double line_space;
  void update_contour_for_skeleton(int layerID);
  bool b_draw_contour_for_skeleton;

  float min_z,max_z;
};









