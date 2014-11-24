#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkContourWidget.h>
#include <vtkOrientedGlyphContourRepresentation.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCommand.h>
#include <vtkDebugLeaks.h>
#include <vtkCamera.h>
#include <vtkPlane.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkMath.h>
#include <vtkWidgetEvent.h>
#include <vtkWidgetEventTranslator.h>
#include <vtkPolyDataMapper.h>
#include <vtkLightCollection.h>
#include "CShape.h"
#include "fstream"
#include "string"
#include <unistd.h>
#include "CPoint.h"
#include "CLayer.h"
#include <vtkLight.h>
#include "initial.h"
#include "CDataReader.h"
using namespace std;

static CDataReader m_data_reader("../config.xml");
const int MatrixRes=m_data_reader.get_value_int("MatrixRes"); /**< the resolution for the lattice */
const int NumRows = m_data_reader.get_value_int ("NumRows"); /**< the numbers of the rows for the distancemap */
const int NumCols = m_data_reader.get_value_int ("NumCols"); /**< the numbers of the cols for the distancemap */
const int interval_of_point = m_data_reader.get_value_int ("interval_of_point"); /**< if the interval of point is n, there will be n-1 points not read */
const float layer_interval_scale = m_data_reader.get_value_float ("layer_interval_scale"); /**< the z value will multip the scale value when read */
const int interval_of_line_layer = m_data_reader.get_value_int ("interval_of_line_layer");/**< the interval of the lines within a layer use for correspondence display*/
const int distance_medialaxis_contour = m_data_reader.get_value_int ("distance_medialaxis_contour"); /**< the threshold for the medialaxis */
const int iteration_number = m_data_reader.get_value_int ("iteration_number"); /**< the iteration number of the gradient descent */
const bool use_medial_axis = m_data_reader.get_value_bool ("use_medial_axis"); /**< whether use the medial axis*/
const bool use_normalize_for_points = m_data_reader.get_value_bool ("use_normalize_for_points");/**<not use now, must use ,  becuse maybe the positon of the points was not in the first section for the cooradinate.*/
const bool show_vertex_on_edge = m_data_reader.get_value_bool ("show_vertex_on_edge"); /**< whether show the vertex on the edge */
const bool use_length_energy = m_data_reader.get_value_bool ("use_length_energy"); /**< whether use the length to decide the deformation sequnce */
const float contour_distance_threshold = m_data_reader.get_value_float ("contour_distance_threshold"); /**< wheter the contours can be a couple depend the their distance smaller than the threshold */
const float point_scale = m_data_reader.get_value_float ("point_scale"); /**< the scale of the point x,y value when read */
const float contour_length_threshold = m_data_reader.get_value_float ("contour_length_threshold");/**<if use this threshold the, the raw data for one  must be continous.like for contourId 1. the data must be in together, can't mixed with otherdata.  */
const float contour_pointnum_threshold = m_data_reader.get_value_float ("contour_pointnum_threshold");/**< the contour's length small than the threshold will be cascaded*/
const float point_distance_scale = m_data_reader.get_value_float ("point_distance_scale"); /**< avoid the disorder of points, the default value is 0. */
const float boundary_coverage_rate_threshold = m_data_reader.get_value_float ("boundary_coverage_rate_threshold");/**< this is a intersection rate threshold. if adjacent contour intersection bigger than the threshold will be arranged.*/
const bool use_open_mp = m_data_reader.get_value_bool ("use_open_mp"); /**< whether use open_mp */
const int thread_num = m_data_reader.get_value_int ("thread_num"); /**< the threads number when using openmp */
//@}
const int gap_threshold = m_data_reader.get_value_float("gap_threshold"); /**< the gap of the hole, the threshold decide whether the hole to fill, the gap is metric by points num */
const bool use_contour_smooth= m_data_reader.get_value_bool("use_contour_smooth"); /**< whether use the smooth in setup the ccontour */
const bool write_shape_triangle_to_data= m_data_reader.get_value_bool("write_shape_triangle_to_data");
const bool write_shape_points_to_data= m_data_reader.get_value_bool("write_shape_points_to_data");

const bool debug_contour_arrangement= m_data_reader.get_value_bool("debug_contour_arrangement",2);
const bool debug_registration= m_data_reader.get_value_bool("debug_registration",2);
const bool show_contour_surface= m_data_reader.get_value_bool("show_contour_surface");


CShape m_source;

void init_data();
void vtk_setup();
void registration();
void vtk_draw_view1(vtkSmartPointer<vtkRenderWindow>,vtkSmartPointer<vtkRenderWindowInteractor>);
void vtk_draw_view2(vtkSmartPointer<vtkRenderWindow>,vtkSmartPointer<vtkRenderWindowInteractor>);
void vtk_draw_view3(vtkSmartPointer<vtkRenderWindow>,vtkSmartPointer<vtkRenderWindowInteractor>);
void vtk_draw_view4(vtkSmartPointer<vtkRenderWindow>,vtkSmartPointer<vtkRenderWindowInteractor>);
void vtk_draw_view5(vtkSmartPointer<vtkRenderWindow>,vtkSmartPointer<vtkRenderWindowInteractor>);
int main( int argc, char *argv[] )
{
  // Create the RenderWindow, Renderer and both Actors
  //

  init_data();

  registration();
  
  vtk_setup();
  
  return EXIT_SUCCESS;
}

void registration()
{
  if (use_open_mp)
  {
    //    m_source.Setup_use_openmp();//can't use the parallel version, cause CPoint::get_index()will introduce the conflic! when invoke the contour->setup();
    m_source.Setup();
  }
  else
  {
    m_source.Setup();
  }
  
  m_source.Registration();
  m_source.initial_vtk_points();
}

void init_data()
{
  char currentPath[200];
  getcwd(currentPath, sizeof(currentPath));
  m_data_reader.read_file(&m_source);
}

void vtk_setup()
{
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  m_source.initial_display(renderWindow,interactor);
  interactor->SetRenderWindow(renderWindow);  
  //  renderWindow->SetFullScreen(true);
  renderWindow->SetSize(600,600);
  int N=5;
  for (int i = 0; i < N; ++i)
  {
    switch(i+1)
    {
      case 1:
        vtk_draw_view1(renderWindow,interactor);
        break;
      case 2:
        vtk_draw_view2(renderWindow,interactor);
        break;
      case 3:
        vtk_draw_view3(renderWindow,interactor);
        break;
      case 4:
        vtk_draw_view4(renderWindow,interactor);
        break;
      case 5:
        vtk_draw_view5(renderWindow,interactor);
        break;
    }
  }
  renderWindow->Render();
  interactor->Initialize();
  interactor->Start();
}

//cordinate for normal
double xmins[5] = {0,.25,0.5,0,.5};
double xmaxs[5] = {0.25,0.50,1,0.5,1};
double ymins[5] = {0,0,0,.5,.5};
double ymaxs[5]= {0.5,0.5,0.5,1,1};

 // double xmins[5] = {0,0,0,0,0};
 // double xmaxs[5] = {1,0,0,0,0};
 // double ymins[5] = {0,0,0,0,0};
 // double ymaxs[5]= {1,0,0,0,0};

//double backgroundcolor={0.1,0.2,0.4};
double backgroundcolor[3]={1,1,1};
double contour_line_color[3]={0.3,0.3,0.3};
double contour_surface_color[3]={0.063,0.48,0.254};
double triangle_surface_color[3]={0.027,0.27,0.12};
double correspondence_color[3]={0,0,0};
int linewidth=5;


vtkSmartPointer<vtkRenderer> Renderers[5];
void vtk_draw_view1(vtkSmartPointer<vtkRenderWindow> renderWindow,vtkSmartPointer<vtkRenderWindowInteractor> interactor)
{
  vtkSmartPointer<vtkRenderer> renderer= vtkSmartPointer<vtkRenderer>::New();
  Renderers[0]=renderer;
  //  renderer->RemoveLight( renderer->GetLights()->GetNextItem());
  renderWindow->AddRenderer(renderer);
  
  renderer->SetBackground(backgroundcolor);

  renderer->SetViewport(xmins[0],ymins[0],xmaxs[0],ymaxs[0]);
  
    m_source.m_display->draw_origin_points(renderer);
  //  m_source.m_display->draw_normal_points(renderer);

  renderer->GetActiveCamera()->SetParallelProjection(1);
  
  renderer->ResetCamera();

  // renderWindow->Render();
}

void vtk_draw_view2(vtkSmartPointer<vtkRenderWindow> renderWindow,vtkSmartPointer<vtkRenderWindowInteractor> interactor)
{
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderWindow->AddRenderer(renderer);
  Renderers[1]=renderer;
  //renderer->RemoveLight( renderer->GetLights()->GetNextItem());
  vtkSmartPointer<vtkLight> light = vtkSmartPointer<vtkLight>::New();
  light->SetLightTypeToSceneLight();
  light->SetPosition(100, 100, 100);
  light->SetFocalPoint(50,50,0); 
  light->SetColor(1,1,1);
  light->SetSpecularColor(1,1,1);
  light->SetPositional(true); // required for vtkLightActor below
  light->SetSwitch(true);

  renderer->UpdateLightsGeometryToFollowCamera();

  // light->SetConeAngle(10);
  // light->SetFocalPoint(lightFocalPoint[0], lightFocalPoint[1], lightFocalPoint[2]);
  // light->SetDiffuseColor(1,0,0);
  // light->SetAmbientColor(0,1,0);
  // light->SetSpecularColor(0,0,1);
  
  renderer->SetBackground(backgroundcolor);

  renderer->SetViewport(xmins[1],ymins[1],xmaxs[1],ymaxs[1]);
  
  //  m_source.m_display->draw_normal_points(renderer);
  m_source.m_display->draw_triangle(renderer,interactor);

  renderer->GetActiveCamera()->SetParallelProjection(1);
    
  renderer->ResetCamera();

  renderWindow->Render();
  renderer->AddLight(light);
}

void vtk_draw_view3(vtkSmartPointer<vtkRenderWindow> renderWindow,vtkSmartPointer<vtkRenderWindowInteractor> interactor)
{
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderWindow->AddRenderer(renderer);
  Renderers[2]=renderer;
  renderer->SetBackground(backgroundcolor);

  renderer->SetViewport(xmins[2],ymins[2],xmaxs[2],ymaxs[2]);
  //if transformed vector size =0 ; will render the origin points
  m_source.m_display->draw_transformed_points(renderer);

  renderer->GetActiveCamera()->SetParallelProjection(1);
  
  renderer->ResetCamera();

  // renderWindow->Render();
}
#include "vtkSphereSource.h"
void vtk_draw_view4(vtkSmartPointer<vtkRenderWindow> renderWindow,vtkSmartPointer<vtkRenderWindowInteractor> interactor)
{
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderWindow->AddRenderer(renderer);
  Renderers[3]=renderer;
  renderer->SetBackground(backgroundcolor);

  renderer->SetViewport(xmins[3],ymins[3],xmaxs[3],ymaxs[3]);
  
  m_source.m_display->draw_skeleton(renderer,interactor);

  renderer->GetActiveCamera()->SetParallelProjection(1);
  
  renderer->ResetCamera();

  // renderWindow->Render();
}

void vtk_draw_view5(vtkSmartPointer<vtkRenderWindow> renderWindow,vtkSmartPointer<vtkRenderWindowInteractor> interactor)
{
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderWindow->AddRenderer(renderer);
  Renderers[4]=renderer;
  renderer->SetBackground(backgroundcolor);

  renderer->SetViewport(xmins[4],ymins[4],xmaxs[4],ymaxs[4]);
  
  m_source.m_display->draw_correspond(renderer);

  renderer->GetActiveCamera()->SetParallelProjection(1);
  
  renderer->ResetCamera();

  // renderWindow->Render();
      // Create a sphere
}
