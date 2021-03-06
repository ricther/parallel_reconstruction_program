#pragma once
// for algorithm initial
extern const int MatrixRes;
extern const int NumRows;
extern const int NumCols;
extern const int interval_of_point;
extern const float layer_interval_scale;
extern const int interval_of_line_layer;// the interval of the lines within one layer
extern const int distance_medialaxis_contour;
extern const int iteration_number;
extern const bool use_medial_axis;
extern const bool use_length_energy;
extern const float contour_distance_threshold;
extern const float point_scale;
extern const float contour_length_threshold;
extern const float contour_pointnum_threshold;
extern const float point_distance_scale;//if the distance bewteen two point bigger than point_distance_scale*avarge_length. the point is disordered
extern const float boundary_coverage_rate_threshold;
extern const bool use_open_mp;
extern const int thread_num;
extern const int gap_threshold;
extern const bool use_contour_smooth;
extern const bool show_contour_surface;
extern const bool write_shape_points_to_data;
extern const bool write_shape_triangle_to_data;
// for debug

extern const bool debug_contour_arrangement;
extern const bool debug_registration;
