//////////////////////////////////////////////////////////
//   Implementation for Memory allocation functions    ///
//////////////////////////////////////////////////////////


#include "stdlib.h"
#include "Memory.h"

// Allocate the space for a 2D array of boolean
double *make_1D_double_array(int size)
{
	double *arr;
	if ((arr=(double *)calloc(size,sizeof(double )))==NULL)
    {return(NULL);}
	return(arr);
} 

// Free the space of a 2D array of boolean
void free_1D_double_array(double *arr)
{
	free((char*)arr);
}

// Allocate the space for a 2D array of boolean
float *make_1D_float_array(int size)
{
	float *arr;
	if ((arr=(float *)calloc(size,sizeof(float)))==NULL)
    {return(NULL);}
	return(arr);
} 

// Free the space of a 2D array of boolean
void free_1D_float_array(float *arr)
{
	free((char*)arr);
}

// Allocate the space for a 2D array of boolean
int *make_1D_int_array(int size)
{
	int *arr;
	if ((arr=(int *)calloc(size,sizeof(int)))==NULL)
    {return(NULL);}
	return(arr);
} 

// Free the space of a 2D array of boolean
void free_1D_int_array(int *arr)
{
	free((char*)arr);
}

// Allocate the space for a 2D array of boolean
bool **make_2D_bool_array(int row,int col)
{
	bool *arr_ptr;
	bool **arr;
	int i;
	if ((arr=(bool **)calloc(row,sizeof(bool *)))==NULL)
    {return(NULL);}
	if ((arr_ptr=(bool *)calloc(row*col,sizeof(bool)))==NULL)
    {return(NULL);}
	for(i=0;i<row;i++) arr[i]=&arr_ptr[i*col];
	return(arr);
} 

// Free the space of a 2D array of boolean
void free_2D_bool_array(bool **arr)
{
	free((char*)arr[0]);
	free((char*)arr);
}

// Allocate the space for a 2D array of unsigned char
unsigned char **make_2D_byte_array(int row,int col)
{
	unsigned char *arr_ptr;
	unsigned char **arr;
	int i;
	if ((arr=(unsigned char **)calloc(row,sizeof(unsigned char *)))==NULL)
    {return(NULL);}
	if ((arr_ptr=(unsigned char *)calloc(row*col,sizeof(unsigned char)))==NULL)
    {return(NULL);}
	for(i=0;i<row;i++) arr[i]=&arr_ptr[i*col];
	return(arr);
} 

// Free the space of a 2D array of unsigned char
void free_2D_byte_array(unsigned char **arr)
{
	free((char*)arr[0]);
	free((char*)arr);
}

// Allocate the space for a 2D array of int
int **make_2D_int_array(int row,int col)
{
	int *arr_ptr;
	int **arr;
	int i;
	if ((arr=(int **)calloc(row,sizeof(int *)))==NULL)
    {return(NULL);}
	if ((arr_ptr=(int *)calloc(row*col,sizeof(int)))==NULL)
    {return(NULL);}
	for(i=0;i<row;i++) arr[i]=&arr_ptr[i*col];
	return(arr);
} 

// Free the space of a 2D array of int
void free_2D_int_array(int **arr)
{
	free((char*)arr[0]);
	free((char*)arr);
}

/********************************************/
// Allocate the space for a 2D array of float
float **make_2D_float_array(int row,int col)
{
	float *arr_ptr;
	float **arr;
	int i;
	if ((arr=(float **)calloc(row,sizeof(float *)))==NULL)
    {return(NULL);}
	if ((arr_ptr=(float *)calloc(row*col,sizeof(float)))==NULL)
    {return(NULL);}
	for(i=0;i<row;i++) arr[i]=&arr_ptr[i*col];
	return(arr);
} 

/***************************************/
// Free the space of a 2D array of float
void free_2D_float_array(float **arr)
{
	free((char*)arr[0]);
	free((char*)arr);
}

/*********************************************/
// Allocate the space for a 2D array of double
double **make_2D_double_array(int row,int col)
{
	double *arr_ptr;
	double **arr;
	int i;
	if ((arr=(double **)calloc(row,sizeof(double *)))==NULL)
    {return(NULL);}
	if ((arr_ptr=(double *)calloc(row*col,sizeof(double)))==NULL)
    {return(NULL);}
	for(i=0;i<row;i++) arr[i]=&arr_ptr[i*col];
	return(arr);
} 

/****************************************/
// Free the space of a 2D array of double
void free_2D_double_array(double **arr)
{
	free((char*)arr[0]);
	free((char*)arr);
}


