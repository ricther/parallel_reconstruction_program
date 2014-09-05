#include "CRegistration.h"
#include "CShape.h"
#include "string.h"
#include "Memory.h"
#include "CLayer.h"
#include "CPoint.h"
#include "CMap.h"
#include "CCorespond.h"
#include "CFileDebug.h"
#include "initial.h"
#include "CContourArrangement.h"
#include <omp.h>
using namespace std;
void ParamRecord::operator=(const ParamRecord & right)
{
  int N =MatrixRes;
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      this->lattice_x[i][j] = right.lattice_x[i][j];
      this->lattice_y[i][j] = right.lattice_y[i][j];
    }
  }
  this->errorE = right.errorE;
}

void ParamRecord::initial(float errorE,CVirtualContour* higher)
{
  this->lattice_x = make_2D_float_array(MatrixRes,MatrixRes);//initial the new lattice
  this->lattice_y = make_2D_float_array(MatrixRes,MatrixRes);
  for (int j = 0; j < MatrixRes; ++j)
  {
    for (int k = 0; k < MatrixRes; ++k)
    {
      this->lattice_x[j][k]=higher->lattice_x[j][k];
      this->lattice_y[j][k]=higher->lattice_y[j][k];
    }
  }
  this->errorE = errorE;
}

//source is the higher layer, target is lower layer.
//call source & target is consistency with the registration convention

CRegistration::CRegistration(CShape* source):narrow_band(10),NumberRows(NumRows),NumberCols(NumCols)
{
  SourceShape=source;
  
  kappa = 5.0;

  xvb=new float[MatrixRes];
  yvb=new float[MatrixRes];
  
  int interval = NumberRows/MatrixRes;
  for (int i = 0; i < MatrixRes; ++i)
  {
    xvb[i]= -interval+i*interval -1;
    yvb[i]= xvb[i];
  }
}

void CRegistration::reset(CVirtualContour* higher)
{
  for (int i = 0; i < MatrixRes; ++i)
  {
    for (int j = 0; j < MatrixRes; ++j)
    {
      higher->XB[i][j]=0;higher->YB[i][j]=0;
      higher->dXB[i][j]=0;higher->dYB[i][j]=0;
    }
  }
  higher->lamda = 0.02;
}

void CRegistration::freeform_res1(CVirtualContour* lower,CVirtualContour*higher)
{
  float *xv;
  float *yv;
  int count;
  init_lattice(higher);
  shape_vicinity(higher,narrow_band);
  ////////////////
  // string filestr;
  //filestr=higher->filename+"intensity_data_1";
  // CFileDebug m_file(filestr);
  //filestr=higher->filename+"intensity_data_2";
  // CFileDebug m_file2("synthetic3_higher_1_dismap");
  ///////////////
  //first time calculate the intensity before bespline as though the points do not changed
  //calculate the intensity in self distancesmap
  compute_intensity(higher->vec_Points_Vicinity,higher->m_Map->get_distance_map(),higher->vec_intensity_old);
  // m_file.Output_intensity_data(higher->vec_Points_Vicinity,higher->vec_intensity_old);
  // m_file2.Output(higher->m_Map->DistancsMap);
  //  init_lattice(lower);

  bspline_update(higher,0,higher->vec_Points_Vicinity,higher->vec_new_points_vicinity);
  //second time calculate the intensity after bespline update as though the points was changed
  //calculate the intensity in lower distancemap
  compute_intensity(higher->vec_new_points_vicinity,lower->m_Map->get_distance_map(),higher->vec_intensity_new);
  //m_file2.Output_intensity_data(higher->vec_new_points_vicinity,higher->vec_intensity_new);
  float errorE=energy_func_square_diff(higher);
  gradientdescent_smoother(lower,higher,errorE);
}
#include "iostream"
void CRegistration::Register()
{
  float start_time=omp_get_wtime();
  int src_layer_num=SourceShape->map_Layer.size();
  SourceShape->contour_arrangement->setup();
  //  contour_arrangement

  std::vector<contour_couple*>::iterator itr=SourceShape->contour_arrangement->vec_contour_couple.begin(),etr=SourceShape->contour_arrangement->vec_contour_couple.end();
  int size=SourceShape->contour_arrangement->vec_contour_couple.size();
  int counter=0;
  for (;itr!=etr;++itr)
  {
    counter++;
    std::cout<<"processing couple:"<<counter<<"\t total couple num:"<<size<<"\t complete:"<<counter*1.0/(size*1.0)<<"\n";
    regist_lower_long_higher_short((*itr)->vcontour1,(*itr)->vcontour2);
  }
  
  // std::map<float,CLayer*>::iterator itr=SourceShape->map_Layer.begin(),etr=SourceShape->map_Layer.end(),nitr;
  // int count=1;
  // contour_couple* temp_couple;
  
  // for (;itr!=etr; ++itr)
  // {

  //   nitr=itr;
  //   ++nitr;
  //   cout<<"processing layer:"<<itr->first<<"-"<<nitr->first<<"; count:"<<count<<"\ttotal layer num:"<<src_layer_num<<"\t"<<"complete:"<<count*1.0/(src_layer_num*1.0)<<"\n";
  //   count++;
  //   if (nitr==etr)
  //   {
  //     break;
  //   }
  //   lower_layer=itr->second;
  //   higher_layer=nitr->second;
  //   lower_layer->reset();
  //   higher_layer->reset();


  //   int couple_num=contour_arrangement.vec_contour_couple.size();
  //   if (couple_num==0)
  //   {
  //     continue;
  //   }
    
  //   for (int i = 0; i < couple_num; ++i)
  //   {
  //    temp_couple=contour_arrangement.vec_contour_couple[i];
  //    regist_lower_long_higher_short(temp_couple->contour1,temp_couple->contour2);
  //   }
  // }
  float end_time=omp_get_wtime();
  std::cout<<"***************************\n"<<"registration cost time:"<<start_time-end_time<<endl<<"***************************\n";

}

void CRegistration::Register_use_openmp()
{
  float start_time=omp_get_wtime();
  int src_layer_num=SourceShape->map_Layer.size();
  SourceShape->contour_arrangement->setup_use_omp();

  std::vector<contour_couple*>::iterator itr=SourceShape->contour_arrangement->vec_contour_couple.begin(),etr=SourceShape->contour_arrangement->vec_contour_couple.end();
  int size=SourceShape->contour_arrangement->vec_contour_couple.size();
  int counter=0;
  
#pragma omp parallel for num_threads(thread_num)
  for (int i=0;i<size;++i)
  {
    counter++;
    std::cout<<"processing couple:"<<counter<<"\t total couple num:"<<size<<"\t complete:"<<counter*1.0/(size*1.0)<<"\n";
    regist_lower_long_higher_short(SourceShape->contour_arrangement->vec_contour_couple[i]->vcontour1,SourceShape->contour_arrangement->vec_contour_couple[i]->vcontour2);
  }
  


  
//   std::map<float,CLayer*>::iterator itr=SourceShape->map_Layer.begin(),etr=SourceShape->map_Layer.end(),nitr;
//   std::vector<float> vec_layerID;
//   int count=1;
//   contour_couple* temp_couple;

//   for (;itr!=etr ; ++itr)
//   {
//     vec_layerID.push_back(itr->first);
//   }

//   int size=vec_layerID.size();

//   std::cout<<"use openmp to acc"<<std::endl;
//   for (int ii = 0; ii < 2; ++ii)
//   {
// #pragma omp parallel for num_threads(thread_num)
//     for (int ly=ii;ly<size-1-2;ly+=2)
//     {
//     CLayer *lower_layer=NULL;
//     CLayer *higher_layer=NULL;
//     CContour* lower_contour=NULL;
//     CContour* higher_contour=NULL;


//     CContourArrangement contour_arrangement(SourceShape);
//     contour_arrangement.setup();

//     lower_layer=SourceShape->map_Layer[vec_layerID[ly]];
//     higher_layer=SourceShape->map_Layer[vec_layerID[ly+1]];
//     #pragma omp critical
//     {
//     cout<<"layer_itration:"<<ly<<";layer_id"<<lower_layer->LayerID<<endl;
//     cout<<"processing layer:"<<lower_layer->LayerID<<"-"<<higher_layer->LayerID<<"; count:"<<count<<"\ttotal layer num:"<<src_layer_num<<"\t"<<"complete:"<<count*1.0/(src_layer_num*1.0)<<"\n";
//     }
//     count++;

//     lower_layer->reset();
//     higher_layer->reset();

//     int couple_num=contour_arrangement.vec_contour_couple.size();
//     if (couple_num==0)
//     {
//       continue;
//     }
    
//     for (int i = 0; i < couple_num; ++i)
//     {
//      temp_couple=contour_arrangement.vec_contour_couple[i];
//      regist_lower_long_higher_short(temp_couple->contour1,temp_couple->contour2);
//     }
//   }
//   }
  float end_time=omp_get_wtime();
  std::cout<<"***************************\n"<<"registration cost time with use openmp:"<<start_time-end_time<<endl<<"***************************\n";
}





void CRegistration:: regist_lower_long_higher_short(CVirtualContour* lower,CVirtualContour* higher)
{
  if (lower->LayerID>higher->LayerID)
  {
    CVirtualContour* temp = lower;
    lower=higher;
    higher=temp;
  }
  CVirtualContour* longer,*shorter;
  if (use_length_energy)
  {
    if (lower->m_contour->length>=higher->m_contour->length)
    {
      longer=lower;shorter=higher;
    }
    else
    {
      longer=higher;shorter=lower;
    }
  }
  else
  {
    longer=lower;shorter=higher;
  }

  if (use_medial_axis)
  {
    longer->calculate_medial_axis(shorter->LayerID);
    if(longer->medial_axis_count>0)
    {
      longer->calculate_medial_map(longer->medial_axis);
      longer->swap_map_medialmap();
      freeform_res1(longer,shorter);
      longer->swap_medialmap_map();
    }
    else
    {
      freeform_res1(longer,shorter);
    }
  }
  else
  {
    freeform_res1(longer,shorter);
  }
  CCorrespond* cor = new CCorrespond();
  cor->CorrespondLayer1=lower->LayerID;
  cor->CorrespondLayer2=higher->LayerID;    
  get_correspondence(cor,shorter,longer);
  
  if(debug_registration)
  cout<<lower->LayerID<<"-"<<higher->LayerID<<"have correspondence points:"<<cor->map_cor.size()<<endl;
  
  SourceShape->vec_Cor.push_back(cor);
  //only the higher contour need release the array, the lower don't alloc the memory
  if(shorter->lattice_x!=NULL)
  {
     shorter->release_array();
  }
  else
  {
     longer->release_array();
  }

}

CVirtualContour* CRegistration::find_nearest_contour(CVirtualContour* lower_contour,CLayer* higher_layer,int diffecence)
{
  assert(false);
  // float dis=9999999;
  // map<int,CVirtualContour*>::iterator itr,etr;
  CVirtualContour* temp_contour=NULL;
  // itr=higher_layer->map_contour.begin();etr=higher_layer->map_contour.end();
  // for (;itr!=etr;++itr)
  // {
  //   //    CPoint* point1= lower_contour->center_point;
  //   //    CPoint* point2= (itr->second)->center_point;
  //   CPoint* point1= lower_contour->moment_one_point;
  //   CPoint* point2= (itr->second)->moment_one_point;
  //   if (itr->second->use_as_higher_contour_count<=diffecence)
  //   {
  //     float tempdis=(point1->x-point2->x)*(point1->x-point2->x)+(point1->y-point2->y)*(point1->y-point2->y);
  //     if (tempdis<dis)
  //     {
  //       dis=tempdis;
  //       temp_contour=(itr->second);
  //     }
  //   }
  // }
  // temp_contour->use_as_higher_contour_count++;
  return temp_contour;
}
#include "math.h"
bool CRegistration::check_contour_distance(CVirtualContour* lower_contour,CVirtualContour* higher_contour)
{
  if(contour_distance_threshold<=0)
  {
    return true;
  }
  // if have medial points, then the two layer have intersections, so needn't to calculate the distance
  if (lower_contour->vec_medial_points.size()>0||higher_contour->vec_medial_points.size()>0)
  {
    return true;
  }
  CPoint* point1= lower_contour->m_contour->moment_one_point;
  CPoint* point2= higher_contour->m_contour->moment_one_point;

  float tempdis=(point1->x-point2->x)*(point1->x-point2->x)+(point1->y-point2->y)*(point1->y-point2->y);
  tempdis= sqrt(tempdis);
  std::cout<<point1->z<<"-"<<point2->z<<":"<<"contour_distance:"<<tempdis<<"\n";
  if (contour_distance_threshold>tempdis)
  {
    return true;
  }
  else
  {
    return false;
  }
}


void CRegistration::shape_vicinity(CVirtualContour* source,int band)
{
  int size = source->m_contour->vec_Points_Origin.size();
  if (source->vec_Points_Vicinity.size()!=0)
  {
    return;
  }
  source->vec_Points_Vicinity.clear();
  source->vec_Points_Vicinity.reserve(size*band);
  int k=0;
  for (int i = 0; i < size; ++i)
  {
    int cx=(int)source->m_contour->vec_Points_Origin[i]->x+1;
    int cy=(int)source->m_contour->vec_Points_Origin[i]->y+1;
    //TODO    int cz
    for (int j = 0; j < (int)(2.0*band/3.0); ++j)
    {
      CPoint* temp=new CPoint();
      temp->x=cx-band+j*3-1;
      temp->y=cy-band+j*3-1;
      //TODO temp->z
      source->vec_Points_Vicinity.push_back(temp);
      ++k;
    }
  }
  vicinity_points_num=k;
}

void CRegistration::compute_intensity(std::vector<CPoint*>&vec_points,float**&map,std::vector<float>&vec_intensity)
{
  vec_intensity.clear();
  int size = vec_points.size();
  for (int i = 0; i < size; ++i)
  {//TODO here x,y ,a ,b exchange position . consistency with xiaolei's code but , maybe not right ...have change back , exchange is not right/// for test now exchange ,because already change the distmap x,y 
    float px,py;
    px=vec_points[i]->y;
    py=vec_points[i]->x;
    int x = (int)px;
    int y= (int)py;
    float a = (float)(px-x);
    float b = (float)(py-y);
    if ((x>=0)&&(y>=0)&&(x<(NumberRows-1))&&(y<(NumberCols-1)))
    {
      float e=(float)((1.0-a)*(1.0-b)*map[x][y]+(1.0-a)*b*map[x][y+1]+a*b*map[x+1][y+1]+a*(1.0-b)*map[x+1][y]);
      vec_intensity.push_back(e);
    }
    else
    {
      if(x<0)x=0;
      if (y<0)
      {
        y=0;
      }
      if (x>=(NumberRows-2))
      {
        x=NumberRows-1;
      }
      if (y>=(NumberCols-2))
      {
        y=NumberCols-1;
      }
      vec_intensity.push_back(map[x][y]);
    }
  }
}


void CRegistration::init_lattice(CVirtualContour* contour)
{
  if (contour->XB==NULL)
  {
    contour->XB=make_2D_float_array(MatrixRes,MatrixRes);
    contour->YB=make_2D_float_array(MatrixRes,MatrixRes);
    //    reset(contour);
  }
  contour->dXB=make_2D_float_array(MatrixRes,MatrixRes);
  contour->dYB=make_2D_float_array(MatrixRes,MatrixRes);

  for (int i = 0; i < MatrixRes; ++i)
  {
    for (int j = 0;j < MatrixRes; ++j)
    {
      contour->XB[j][i]=xvb[i];
      contour->YB[j][i]=yvb[j];
    }
    
  }

  contour->lattice_x=make_2D_float_array(MatrixRes,MatrixRes);
  contour->lattice_y=make_2D_float_array(MatrixRes,MatrixRes);
  contour->new_lattice_x=make_2D_float_array(MatrixRes,MatrixRes);
  contour->new_lattice_y=make_2D_float_array(MatrixRes,MatrixRes);
  
  for (int i = 0; i < MatrixRes; ++i)
  {
    for (int j = 0; j < MatrixRes; ++j)
    {
      contour->lattice_x[i][j]= contour->XB[i][j]+contour->dXB[i][j];
      contour->lattice_y[i][j]= contour->YB[i][j]+contour->dYB[i][j];
    }
  }
}
#include <vector>
//mode=0 use old_lattice mode=1 use new_lattice
void CRegistration::bspline_update(CVirtualContour* contour,int mode,std::vector<CPoint*>&vec_points,std::vector<CPoint*>&vec_new_points)
{
  int res = 1;
  int nx = MatrixRes;
  int ny = MatrixRes;
  int size = vec_points.size();
  
  vec_new_points.clear();
  vec_new_points.reserve(size);
  for (int i = 0; i < size; ++i)
  {
    float x = vec_points[i]->x;
    float y = vec_points[i]->y;

    float x_index = (x- min(xvb,nx))/(max(xvb,nx)-min(xvb,nx))*(nx-1);
    float y_index = (y- min(yvb,ny))/(max(yvb,ny)-min(yvb,ny))*(ny-1);

    int i_index=(int)(x_index)-1;
    int j_index=(int)(y_index)-1;
    
    float u_index=(x_index)-(int)(x_index);
    float v_index=(y_index)-(int)(y_index);

    if (i_index<0||i_index>(nx-3-1)||j_index<0||j_index>(ny-3-1))
    {
      float xb = x; float nxb= x;
      float yb = y; float nyb=y;
      CPoint* temp =new CPoint();
      temp->x=nxb;temp->y=nyb;temp->z=vec_points[i]->z;
      vec_new_points.push_back(temp);
    }
    else
    {
      float nxb=0.0,nyb=0.0;
      for (int m = 0; m < 4; ++m)
       {
        for (int n = 0; n < 4; ++n)
        {
          if(mode == 0)
          {
            nxb=nxb+cubic_spline(u_index,m)*cubic_spline(v_index,n)*contour->lattice_x[j_index+n][i_index+m];
            nyb=nyb+cubic_spline(u_index,m)*cubic_spline(v_index,n)*contour->lattice_y[j_index+n][i_index+m];
          }
          else if (mode == 1)
          {
            nxb=nxb+cubic_spline(u_index,m)*cubic_spline(v_index,n)*contour->new_lattice_x[j_index+n][i_index+m];
            nyb=nyb+cubic_spline(u_index,m)*cubic_spline(v_index,n)*contour->new_lattice_y[j_index+n][i_index+m];
          }
        }
      }
      CPoint* temp = new CPoint();
      temp->x=nxb;
      temp->y=nyb;
      temp->z=vec_points[i]->z;
      vec_new_points.push_back(temp);
    }
  }
}

float CRegistration::energy_func_square_diff(CVirtualContour* contour)
{
  int N=contour->vec_intensity_old.size();
  float temp=0.0;
  for (int i = 0; i < N; ++i)
  {
    float tempp=( contour->vec_intensity_old[i]-contour->vec_intensity_new[i]);
    temp+=tempp*tempp;
  }
  return temp;
}

void CRegistration::gradientdescent_smoother(CVirtualContour* lower,CVirtualContour* higher,float errorE)
{
  kNumberOfIteration=iteration_number;
  int magicnumber=MatrixRes;

  current_params.lattice_x=make_2D_float_array(MatrixRes,MatrixRes);
  current_params.lattice_y=make_2D_float_array(MatrixRes,MatrixRes);
  ParamRecord record[kNumberOfIteration+1];
  //  Use Levenberg-Marquardt method to implement the minimization
  //	float lamda = 0.05;  // initialize lamda
  //	float kappa = 2.0;   // the weight parameter between image data and smoothness term
  /////////////////////////////////////////////////////////////////
  float** Dcpx = make_2D_float_array(MatrixRes,MatrixRes);
  float** Dcpy = make_2D_float_array(MatrixRes,MatrixRes);


  ////////////////////////////////////////////////////////////////
  bool stop=0;//this stop works on the iteration
  int it=0;
  for (it = 0; it <= kNumberOfIteration; ++it)
  {
    //  after each iteration, save a copy
    //  save record_hand_freeForm_RES1 record;
    //  save record_hand_freeForm_vicinity record;
    // Set up the parameters (x/y coordinates of the control points)
    memset(Dcpx[0],0,sizeof(float)*MatrixRes*MatrixRes);
    memset(Dcpy[0],0,sizeof(float)*MatrixRes*MatrixRes);
    record[it].initial(errorE,higher);
    calculate_dlattice_by_point(Dcpx,Dcpy,higher,lower);
    //////////////////////////////////////////////////////
    // fstream fileout("Dcpx_Dcpy",ios_base::out|ios_base::app);
    // for (int i = 0; i < MatrixRes; ++i)
    // {
    //   for (int j = 0; j < MatrixRes; ++j)
    //   {
    //     fileout<<"Dcpx["<<i<<"]["<<j<<"]:"<<Dcpx[i][j]<<"\t"<<"Dcpy["<<i<<"]["<<j<<"]:"<<Dcpy[i][j]<<"\r\n";
    //   }
    // }
    // fileout<<"****************************************************\n";
    // fileout<<it<<"end\n";
    // fileout<<"****************************************************\n";
    // fileout.close();
    //////////////////////////////////////////////////////
    //TODO here ingnore  the constrains segment from 709 to 872 in the xiaolei's code (FreeForm.cpp)
    stop = update_lattice(Dcpx,Dcpy,higher,lower,errorE);
    if(stop)
    {
      break;
    }
    //end of for (int it=0...)
  }
  if (it<kNumberOfIteration)
  {
    current_params = record[it];
  }
  else
  {
    current_params = record[kNumberOfIteration];
  }
  std::cout<<"higher_ID:"<<higher->LayerID<<"(contourID:"<<higher->m_contour->contourID<<")"<<"-lower_ID:"<<lower->LayerID<<"(contourID:"<<lower->m_contour->contourID<<")"<<"\tIteration:"<<it<<"\terrorE:"<<errorE<<"\n";
  // for (int xx = 0; xx < MatrixRes; ++xx)
  // {
  //   for (int yy = 0; yy < MatrixRes; ++yy)
  //   {
  //     float x=Dcpx[xx][yy],y=Dcpy[xx][yy];
  //     if (x!=0&&y!=0)
  //     {
  //       std::cout<<"Dcpxy["<<xx<<","<<yy<<"]:"<<x<<","<<y<<";";
  //     }
  //   }
  //   std::cout<<"\n";
  // }
  free_2D_float_array(Dcpx);
  free_2D_float_array(Dcpy);
}

float CRegistration::compute_intensity_by_point(float **& dist, float tx,float ty)
{
  //for text ,exchange x,y .because already exchange the x,y int distmap
  float px=ty;
  float py=tx;
  int x = (int) px;
  int y = (int) py;
  
  float a = px - x;
  float b = py -y;
  int N= NumberRows;
  if (( x >= 0) && (y >= 0) && (x < (N-1)) && (y < (N-1))) 
    return  (float) ((1.0 - a) * (1.0 - b) * dist[x][y] + (1.0 - a) * b * dist[x][y+1] +
                      a * b * dist[x+1][y+1] + a * (1.0 - b) * dist[x+1][y]);
  else 
  {
    if (x < 0) x = 0;
    if (y < 0) y = 0;
    if (x >= (N-2)) x = N-1;
    if (y >= (N-2)) y = N-1;
    return  dist[x][y];
  }
}

float CRegistration::cubic_spline(float u, int o)

	// u: t, the time
	// o: the order
{
	float b;
	switch (o)
	{
	case 0:
		b = (float) (1.0-u)*(1.0-u)*(1.0-u)/ (float) 6.0;
		break;
	case 1: 
		b = (float) (3.0*u*u*u-6*u*u+4)/(float) 6.0;
		break;
	case 2: 
		b = (float) (-3*u*u*u+3*u*u+3*u+1)/(float) 6.0;
		break;
	case 3: 
		b = (float) u*u*u/(float) 6.0;
		break;
	}
	return b;
}

void CRegistration::calculate_dlattice_by_point(float **&Dcpx,float**&Dcpy,CVirtualContour*higher,CVirtualContour*lower)
{
  float fx=0,gtx=0;
  float Dg[2],Dp[2];
  float dux=0,duy=0;
  float dLx[2],dLy[2];
  float dxdLx[2],dydLx[2],dxdLy[2],dydLy[2];
  float dx[2],dy[2];
  int size = higher->vec_Points_Vicinity.size();
  //fstream filefx("fx_gtx_Dg",ios_base::out|ios_base::app);
  for(int k=0;k<size;++k)
  {
    float x = higher->vec_Points_Vicinity[k]->x;
    float y = higher->vec_Points_Vicinity[k]->y;
    float xt= higher->vec_new_points_vicinity[k]->x;
    float yt= higher->vec_new_points_vicinity[k]->y;
    //TODO why y and x exchange position?

    fx = compute_intensity_by_point(higher->m_Map->get_distance_map(),x,y);
    gtx = compute_intensity_by_point(lower->m_Map->get_distance_map(),xt,yt);
    //distance between a transfored feature point and its corresponding target feature point
    //TODO why yt xt exchange their position? don't need exchange
    Dg[0] = compute_intensity_by_point(lower->m_Map->get_gx(),xt,yt);
    Dg[1] = compute_intensity_by_point(lower->m_Map->get_gy(),xt,yt);
    // filefx<<x<<"\t"<<y<<"\t"<<xt<<"\t"<<yt<<"\t"<<fx<<"\t"<<gtx<<"\t"<<Dg[0]<<"\t"<<Dg[1]<<"\r\n";
    //    filefx<<"******************************************************"<<"\r\n";
    //find the index of current pixel
    // find the index of current pixel
    float xIdx = 1+(x-min(xvb, MatrixRes))/(max(xvb,MatrixRes)-min(xvb,MatrixRes))*(MatrixRes-1);  // the index in the control point lattice
    float yIdx = 1+(y-min(yvb,MatrixRes))/(max(yvb,MatrixRes)-min(yvb,MatrixRes))*(MatrixRes-1);
    
    
    int i = (int) (xIdx)-2;
    int j = (int) (yIdx)-2;
    
    float u = xIdx - (int) (xIdx);
    float v = yIdx - (int) (yIdx);
      
    // compute partial u/partial x and partial v/partial x
    dux = 1.0/(max(xvb, MatrixRes)-min(xvb,MatrixRes))*(MatrixRes-1.0);
    duy = 1.0/(max(yvb,MatrixRes)-min(yvb,MatrixRes))*(MatrixRes-1.0);
    
    dLx[0]=(float)0.0;
    dLx[1]=(float)0.0;
    dLy[0]=(float)0.0;
    dLy[1]=(float)0.0;
    
    // cuurent pixel only affects sixteen control points , thus only 2*16 parameters
    if (!(i<0||i>(MatrixRes-3-1)||j<0||j>(MatrixRes-3-1)))
    {
      int m=0;
      for(m=0;m<=3;++m)
      {
        for(int n=0;n<=3;++n)
        {
         //the derivative for control point lattice(j+n,i+m) derived from the image term
          //to do the y x position exchanged!!
          dx[0] = cubic_spline(u,m)*cubic_spline(v,n);
          dx[1] = (float)0.0;
          
          dy[0] = (float)0.0;
          dy[1] = cubic_spline(u,m)*cubic_spline(v,n);
          
          Dcpx[j+n][i+m] = Dcpx[j+n][i+m] - 2*(fx-gtx)*dot(Dg,dx);
          Dcpy[j+n][i+m] = Dcpy[j+n][i+m] - 2*(fx-gtx)*dot(Dg,dy);
          
            // the smoothness term related
          dLx[0] = dLx[0] + spline_deriv(u,m)*dux*cubic_spline(v,n)*higher->dXB[j+n][i+m];
          dLx[1] = dLx[1] + spline_deriv(u,m)*dux*cubic_spline(v, n)*higher->dYB[j+n][i+m];
          dLy[0] = dLy[0] + cubic_spline(u,m)*spline_deriv(v, n)*duy*higher->dXB[j+n][i+m];
          dLy[1] = dLy[1] + cubic_spline(u,m)*spline_deriv(v, n)*duy*higher->dYB[j+n][i+m];
        }
      }
      // the smoothness term increment
      for (m = 0; m <= 3; m++)
      {
        for (int n = 0; n <= 3; n++)
        {
          
          // the derivative for control point NXB(j+n, i+m) derived from the image term
          dxdLx[0] = spline_deriv(u, m)*dux*cubic_spline(v, n);     // derivative with respect to the x coordinate
          dxdLx[1] = (float)0.0;    // derivative with respect to the x coordinate
          dydLx[0] = (float)0.0;
          dydLx[1] = spline_deriv(u, m)*dux*cubic_spline(v, n);
          
          //  	dydLx = [0 spline_deriv(u, m)*dux*cubic_spline(v, n)]';
          
          dxdLy[0] = cubic_spline(u, m)*spline_deriv(v, n)*duy;
          dxdLy[1] = (float) 0.0;     // derivative with respect to the x coordinate
          dydLy[0] = (float) 0.0;
          dydLy[1] = cubic_spline(u, m)*spline_deriv(v, n)*duy;
          
          // update based on the smoothness term
          Dcpx[j+n][i+m] = Dcpx[j+n][i+m] + kappa*((float) 2.0*dot(dLx, dxdLx) + (float) 2.0*dot(dLy, dxdLy));
          Dcpy[j+n][i+m] = Dcpy[j+n][i+m] + kappa*((float) 2.0*dot(dLx, dydLx) + (float) 2.0*dot(dLy, dydLy));
        }
      }
      //end of if
    }
    //end of  for(int k=0;k<size;++k) 
  }
}

bool CRegistration::update_lattice(float**&Dcpx,float**&Dcpy,CVirtualContour*higher,CVirtualContour*lower,float& errorE)
{
  float **ddXB = make_2D_float_array(MatrixRes,MatrixRes);
  float **ddYB = make_2D_float_array(MatrixRes,MatrixRes);

  //because the direction (Dcpx,Dcpy) is the max incremental direction of E, so we use the oppse direction of the (Dcpx,Dcpy) 
  for (int i = 0; i < MatrixRes; ++i)
  {
    for (int j = 0; j < MatrixRes; ++j)
    {
      ddXB[i][j]= -Dcpx[i][j]*higher->lamda;
      ddYB[i][j]= -Dcpy[i][j]*higher->lamda;
    }
  }

  ///////////////////////////////////////////////////
    //     fstream filout_dd("ddxB_ddyB",ios_base::out|ios_base::app);
    // for (int i = 0; i <14; ++i)
    // {
    //   for (int j = 0; j < 14; ++j)
    //   {
    //     filout_dd<<"Dcpx["<<i<<"]["<<j<<"]:"<<ddXB[i][j]<<"\t"<<"Dcpy["<<i<<"]["<<j<<"]:"<<ddYB[i][j]<<"\r\n";
    //   }
    // }
    // filout_dd<<"****************************************************\n";
    // //filout_dd<<it<<"end\n";
    // filout_dd<<"****************************************************\n";
    // filout_dd.close();
    // //////////////////////////
    //////////////////
  bool stop=false;
  int numTries = 0;
  bool boolTry=1; //update at a new point, will try several step sizes
  while(boolTry)
  {
    int N=MatrixRes;
    for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        higher->new_lattice_x[i][j] = higher->lattice_x[i][j] + ddXB[i][j];
        higher->new_lattice_y[i][j] = higher->lattice_y[i][j] + ddYB[i][j];
      }
    }
    // fstream filout_dd("new_lattice",ios_base::out|ios_base::app);
    // for (int i = 0; i <14; ++i)
    // {
    //   for (int j = 0; j < 14; ++j)
    //   {
    //     filout_dd<<"new_lattice_x["<<i<<"]["<<j<<"]:"<< higher->new_lattice_x[i][j]<<"\t"<<"new_lattcie_y["<<i<<"]["<<j<<"]:"<< higher->new_lattice_y[i][j]<<"\r\n";
    //   }
    // }
    // filout_dd<<"****************************************************\n";
    // //filout_dd<<it<<"end\n";
    // filout_dd<<"****************************************************\n";
    // filout_dd.close();
    
    bspline_update(higher,1,higher->vec_Points_Vicinity,higher->vec_new_points_vicinity);//mode =1 ,use the new lattice.
    // fstream filout_p("newpoint",ios_base::out|ios_base::app);
    // int size_p=higher->vec_new_points_vicinity.size();
    // for (int i = 0; i < size_p; ++i)
    // {
    //   filout_p<<i<<":"<<higher->vec_new_points_vicinity[i]->x<<"\t"<<higher->vec_new_points_vicinity[i]->y<<"\r\n";
    // }
    // filout_p<<"****************************************************\n";
    // //filout_dd<<it<<"end\n";
    // filout_p<<"****************************************************\n";
    // filout_p.close();
    compute_intensity(higher->vec_new_points_vicinity,lower->m_Map->get_distance_map(),higher->vec_intensity_new);
    float nextE = energy_func_square_diff(higher);
    if (nextE < errorE)
    {
      for (int i = 0; i < MatrixRes; ++i)
      {
        for (int j = 0; j < MatrixRes; ++j)
        {
          higher->lattice_x[i][j] = higher->new_lattice_x[i][j];
          higher->lattice_y[i][j] = higher->new_lattice_y[i][j];
          higher->dXB[i][j] = higher->dXB[i][j]+ddXB[i][j];
          higher->dYB[i][j] = higher->dYB[i][j]+ddYB[i][j];
        }
      }
      // fstream filout_dxb("dXB",ios_base::out|ios_base::app);

      // for (int i = 0; i < 14; ++i)
      // {
      //   for (int j = 0; j < 14; ++j)
      //   {
      //     filout_dxb<<i<<"\t"<<j<<":"<<higher->lattice_x[i][j]<<"\t"<<higher->lattice_y[i][j]<<"\t"<<dXB[i][j]<<"\t"<<dYB[i][j]<<"\r\n";
      //   }
      // }
      // filout_dxb<<"****************************************************\n";
      // //filout_dd<<it<<"end\n";
      // filout_dxb<<"****************************************************\n";
      // filout_dxb.close();


      
      errorE = nextE;
      // higher->vec_intensity_old.clear();
      //higher->vec_intensity_old = higher->vec_intensity_new;
      //   higher->vec_Points_Vicinity.clear();
      // higher->vec_Points_Vicinity = higher->vec_new_points_vicinity;
      numTries=0;
      boolTry=0;
    }
    else
    {
      numTries = numTries +1;
      higher->lamda = higher->lamda /2;
      if (numTries>10)
      {
        stop = true;// the stop works on the iteratrations
        boolTry=0;
      }
      int N = MatrixRes;
      for (int i = 0; i < N; ++i)
      {
        for (int j = 0; j < N; ++j)
        {
          ddXB[i][j] = -Dcpx[i][j]*higher->lamda;
          ddYB[i][j] = -Dcpy[i][j]*higher->lamda;
        }
      }
    }
  }
  free_2D_float_array(ddXB);
  free_2D_float_array(ddYB);
  return stop;
}

void CRegistration::get_correspondence(CCorrespond* corres,CVirtualContour* shorter,CVirtualContour* longer)
{
  // int N=lower->vec_Points_Origin.size();
  // bspline_update(higher,1,higher->vec_Points_Origin,higher->vec_new_points);
  // //  string filestr;
  // // filestr=higher->filename+"_Original_Points";
  // // CFileDebug m_file(filestr);
  // //m_file.Output_points(higher->vec_Points_Origin);
  // // filestr=higher->filename+"_New_Points";
  // //CFileDebug m_file2(filestr);
  // //m_file2.Output_points(higher->vec_new_points);
  // int index=-1;
  // CPoint temp;
  // for (int i = 0; i < N; ++i)
  // {
  //   temp =*(lower->vec_Points_Origin[i]);
  //   index=get_closest_point(higher->vec_new_points,temp);
  //   Parapoint * new_para_point= new Parapoint();
  //   new_para_point->point1 = temp;
  //   new_para_point->point2 = *(higher->vec_points[index]);
  //   corres->map_cor.insert(std::make_pair(i,new_para_point));
  // }
  bool use_medial_axis_point=false;
  int N=shorter->m_contour->vec_Points_Origin.size();
  int longer_size=longer->m_contour->vec_Points_Origin.size();
  int real_gap_threshold=gap_threshold<longer_size?gap_threshold:longer_size;
  bspline_update(shorter,1,shorter->m_contour->vec_Points_Origin,shorter->vec_new_points);
  //  string filestr;
  // filestr=higher->filename+"_Original_Points";
  // CFileDebug m_file(filestr);
  //m_file.Output_points(higher->vec_Points_Origin);
  // filestr=higher->filename+"_New_Points";
  //CFileDebug m_file2(filestr);
  //m_file2.Output_points(higher->vec_new_points);
  int index=-1;
  int last_index=-1;
  int first_index=-1;
  CPoint temp;
  int offset=0;
  for (int i = 0; i < N; ++i)
  {
    temp =*(shorter->vec_new_points[i]);
    int shorter_layerID=shorter->LayerID;
    if (longer->medial_axis_count>0&&use_medial_axis)
    {
      index=get_closest_point(longer->m_contour->vec_Points_Origin,longer->vec_medial_points,temp, use_medial_axis_point);
      //      std::cout<<"use medial_axis"<<"\n";
    }
    else
    {
      index=get_closest_point(longer->m_contour->vec_Points_Origin,temp);
      if (i==0)
      {
        first_index=index;
      }
    }

    Parapoint * new_para_point= new Parapoint();

    if (corres->CorrespondLayer1==shorter_layerID)
    {
      // new_para_point->point1 = *(shorter->vec_points[i]);
      //new_para_point->point2 = *(longer->vec_points[index]);
       new_para_point->point1 = *(shorter->m_contour->vec_Points_Origin[i]);
       if (use_medial_axis_point&&use_medial_axis)
       {
         CPoint* project_point=new CPoint();
         project_point->index=project_point->get_index();
         int temp_index = project_point->index;
         *project_point=*(longer->vec_medial_points[index]);
         project_point->index=temp_index;
         project_point->z=longer->m_contour->vec_Points_Origin[0]->z;
         new_para_point->point2 = *project_point;
         longer->vec_Points_project.push_back(project_point);
       }
       else
       {
         int gap=0;bool isforeward=true;
         get_gap(first_index,longer_size,index,last_index,gap,isforeward);
         if (gap>1&&gap<real_gap_threshold&&last_index!=-1)
         {
           fill_the_hole(corres,offset,first_index,longer_size,last_index,index,shorter->m_contour->vec_Points_Origin[i],longer,gap,isforeward);
         }
         else if(abs(gap)>1)
         {
           if (debug_registration)
           {
             cout<<"gap:"<<gap<<"\tindex:"<<index<<"\tlast_index:"<<last_index<<"\tfirst_index:"<<first_index<<"\tsize:"<<longer_size<<"\n";
           }
         }
         new_para_point->point2 = *(longer->m_contour->vec_Points_Origin[index]);
       }
    }
    else
    {
      //new_para_point->point2 = *(shorter->vec_points[i]);
      //new_para_point->point1 = *(longer->vec_points[index]);
      new_para_point->point2 = *(shorter->m_contour->vec_Points_Origin[i]);
      if (use_medial_axis_point&&use_medial_axis)
      {
        CPoint* project_point=new CPoint();
        project_point->index=project_point->get_index();
        int temp_index = project_point->index;
        *project_point=*(longer->vec_medial_points[index]);
        project_point->index=temp_index;
        project_point->z=longer->m_contour->vec_Points_Origin[0]->z;
        new_para_point->point1 = *project_point;
        longer->vec_Points_project.push_back(project_point);
      }
      else
      {
        int gap=0;bool isforeward=true;
        get_gap(first_index,longer_size,index,last_index,gap,isforeward);
        if (gap>1&&gap<real_gap_threshold&&last_index!=-1)
        {
          fill_the_hole(corres,offset,first_index,longer_size,last_index,index,longer,shorter->m_contour->vec_Points_Origin[i],gap,isforeward);
        }
        else if(abs(gap)>1)
        {
          if(debug_registration)
          cout<<"gap:"<<gap<<"\tindex:"<<index<<"\tlast_index:"<<last_index<<"\tfirst_index:"<<first_index<<"\tsize:"<<longer_size<<"\n";
        }

        new_para_point->point1 = *(longer->m_contour->vec_Points_Origin[index]);
      }
    }
    last_index=index;
    corres->map_cor.insert(std::make_pair(offset++,new_para_point));
  }
}

int CRegistration::get_relate_index(int first, int size, int old_index)
{
  int relate_index=(old_index-first)%size;
  if (relate_index<0)
  {
    relate_index+=size;
  }
  return relate_index;
}

void CRegistration::get_gap(int first, int size, int index, int last_index,int &real_gap,bool& isforeward )
{
  isforeward=true;
  int gap=index-last_index;
  isforeward=gap>0?true:false;
  
  int reverse_gap=size-abs(gap);
  real_gap=abs(gap);
  if (abs(gap)>reverse_gap)
  {
    real_gap=reverse_gap;
    isforeward=!isforeward;
  }
}

void CRegistration::fill_the_hole(CCorrespond *corres,int& count,int first,int size,int last_index,int index, CPoint* point1, CVirtualContour *higher,int gap,bool isforeward)
{
  if (isforeward)
  {
    for (int i = 1; i < gap; ++i)
    {
      Parapoint * new_para_point= new Parapoint();
      new_para_point->point1=*point1;
      new_para_point->point2=*(higher->m_contour->vec_Points_Origin[(last_index+i)%size]);
      corres->map_cor.insert(std::make_pair(count++,new_para_point));
    }
  }
  else
  {
    for (int i = 1; i < gap; ++i)
    {
      Parapoint * new_para_point= new Parapoint();
      new_para_point->point1=*point1;
      new_para_point->point2=*(higher->m_contour->vec_Points_Origin[(index+i)%size]);
      corres->map_cor.insert(std::make_pair(count++,new_para_point));
    }
  }
}


void CRegistration::fill_the_hole(CCorrespond *corres,int& count,int first,int size,int last_index,int index, CVirtualContour *lower,CPoint* point2,int gap,bool isforeward)
{
  if (isforeward)
  {
    for (int i = 1; i <gap; ++i)
    {
      Parapoint * new_para_point= new Parapoint();
      new_para_point->point2=*point2;
      new_para_point->point1=*(lower->m_contour->vec_Points_Origin[(last_index+i)%size]);
      corres->map_cor.insert(std::make_pair(count++,new_para_point));
    }
  }
  else
  {
    for (int i = 1; i <gap; ++i)
    {
      Parapoint * new_para_point= new Parapoint();
      new_para_point->point2=*point2;
      new_para_point->point1=*(lower->m_contour->vec_Points_Origin[(index+i)%size]);
      corres->map_cor.insert(std::make_pair(count++,new_para_point));
    }
  }
}



int CRegistration::get_closest_point(std::vector<CPoint*>&vec_points,CPoint point)
{
  int N = vec_points.size();
  CPoint temp;
  int result=-1;
  float distancs=999999999;
  for (int i = 0; i < N; ++i)
  {
    temp =*(vec_points[i]);
    float d=(temp.x-point.x)*(temp.x-point.x)+(temp.y-point.y)*(temp.y-point.y);
    d=sqrt(d);
    if (d<distancs)
    {
      distancs=d;
      result=i;
    }
  }
  return result;
}

int CRegistration::get_closest_point(std::vector<CPoint*>&vec_points,std::vector<CPoint*>&medial_points,CPoint point,   bool& use_medial_axis_point)
{
  use_medial_axis_point=false;
  int N = vec_points.size();
  CPoint temp;
  int result=-1;
  float distancs=999999999;
  for (int i = 0; i < N; ++i)
  {
    temp =*(vec_points[i]);
    float d=(temp.x-point.x)*(temp.x-point.x)+(temp.y-point.y)*(temp.y-point.y);
    d=sqrt(d);
    if (d<distancs)
    {
      distancs=d;
      result=i;
    }
  }
  N = medial_points.size();
  for (int i = 0; i < N; ++i)
  {
    temp =*(medial_points[i]);
    float d=(temp.x-point.x)*(temp.x-point.x)+(temp.y-point.y)*(temp.y-point.y);
    d=sqrt(d);
    if (d<distancs)
    {
      use_medial_axis_point=true;
      distancs=d;
      result=i;
    }
  }
  return result;
}
