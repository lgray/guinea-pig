#ifndef MESH_SEEN
#define MESH_SEEN
#include "rndmCPP.h"

class MESH
{

  double delta_x_,delta_y_,delta_z_;
  double offset_x_,offset_y_;

 public:

 MESH() : delta_x_(0.0),delta_y_(0.0),delta_z_(0.0),offset_x_(0.0),offset_y_(0.0) {;}
  ~MESH() {;}
MESH(double delta_x, double delta_y, double delta_z, double offset_x, double offset_y) 
{
  delta_x_=delta_x;
  delta_y_=delta_y;
  delta_z_=delta_z;
  offset_x_=offset_x;
  offset_y_=offset_y;
}

 inline double get_offset_x() const {return offset_x_;}
 inline double get_offset_y() const {return offset_y_;}
 inline double get_delta_x() const {return delta_x_;}
 inline double get_delta_y() const {return delta_y_;}
 inline double get_delta_z() const {return delta_z_;}



inline void guess_position_in_cell(int cell_x, int cell_y, double min_z, double& x,double& y, double& z, RNDM& rndm_generator) const
{
    z=(min_z + rndm_generator.rndm())*delta_z_;
    x=(cell_x + rndm_generator.rndm() - offset_x_)*delta_x_;
    y=(cell_y + rndm_generator.rndm() - offset_y_)*delta_y_;
}

inline void pair_guess_position_in_cell(int cell_x, int cell_y, double min_z, double& x,double& y, double& z, RNDM& rndm_generator) const
{
    z=(min_z + rndm_generator.rndm_pairs())*delta_z_;
    x=(cell_x + rndm_generator.rndm_pairs() - offset_x_)*delta_x_;
    y=(cell_y + rndm_generator.rndm_pairs() - offset_y_)*delta_y_;
}

};

#endif
