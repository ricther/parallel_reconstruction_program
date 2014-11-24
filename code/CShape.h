#pragma once
//athour xu liu 01/10/13
//the class of layer
#include <map>
#include"CCorespond.h"
#include<vector>
#include<string>

#include<vtkSmartPointer.h>
#include<vtkRenderWindow.h>
#include<vtkRenderWindowInteractor.h>
#include"CShapeDisplay.h"
#include"CContourArrangement.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "initial.h"
class CLayer;
class CRegistration;
class CSkeleton;

extern const float layer_interval_scale;
class CShape
{
  friend class CShapeDisplay;
public:
 CShape();
 std::vector<CCorrespond*> vec_Cor;
 std:: map<float,CLayer*> map_Layer;

 void Registration();
 void initial(std::string);//if size of vec is biger than one the section will read multi files and one file contain one //layer else the section will read wohle shape from one file
 void Setup();
 void Setup_use_openmp();
 void initial_display(vtkSmartPointer<vtkRenderWindow>,vtkSmartPointer<vtkRenderWindowInteractor>);
 CShapeDisplay* m_display;
 CSkeleton *m_skeleton;
 float get_next_layer(float now_layerID,int direction);//direction 0 down 1 up
 std::vector<float> vec_layerID;
 CContourArrangement* contour_arrangement;


 std::map<int,int> map_CPointsIndex_vtkIndex;
 vtkSmartPointer<vtkPoints> m_total_points; /**< the m_total_points for the shape */
 int get_vtk_points_index(){return ++vtk_points_counter;};
 
 void insert_map_CPointsIndex_vtkIndex(CPoint* point,int vtkIndex)
 {
   m_total_points->InsertPoint(static_cast<vtkIdType>(vtkIndex),point->x,point->y,point->z);
   map_CPointsIndex_vtkIndex.insert(std::make_pair(point->index,vtkIndex));
 }

 void initial_vtk_points();
private:
 CRegistration * m_registration;
 std::string m_filename;
 void read();
 void check_edge(float ,float ,float);
 void check_edge_before_scale(float,float);
 void get_center();
 void calculate_one_moment();
 CPoint center_point;
 CPoint moment_one_point;
 float max_x,min_x,max_y,min_y,max_z,min_z,sum_x,sum_y;
 float original_max_x,original_max_y,original_min_x,original_min_y;
 void check_point_scale();
 int vtk_points_counter;
 void insert_vtk_points(std::vector<CPoint*>& vec_points);
 void write_vtk_points();
 int total_points_num;
};
















