#ifndef BEAMPARAMETERS_SEEN
#define BEAMPARAMETERS_SEEN


#include <cstdio>
#include <sstream>
#include <string>
#include "abstractParticle.h"
#include "parametersCPP.h"
#include "mathematicalEntities.h"

class BEAM_PARAMETERS : public ABSTRACT_IO_CLASS
{

  double ebeam_,n_particles_;
  double gamma_beam_;
  
  // emittances in m.rad
  double em_x_,em_y_;
  
  double sigma_x_,sigma_y_,sigma_z_;
  int dist_x_, dist_z_;
  
  // in mm
  double beta_x_,beta_y_;
  
  // initial polarization direction 
  // if this vector is not normed (the norme must be <= 1)
  // polarization represents a mixed state
  TRIDVECTOR polar_;
  
  double offset_x_,offset_y_,offset_z_;
  double waist_x_,waist_y_;
  // double couple_xy_;
  double phi_angle_;
  double x_angle_,y_angle_;
  double L_0_,L_;
  //  int   what_distribution_;
  int bunches_per_train_;
  double frep_;
  int trav_focus_;
  char extension_[3];

  void set_polar(double compx, double compy, double compz);

  void create_name(std::string& name, std::string param);

  inline std::string param_with_extension(std::string param) const
    {
      return param.append(extension_);
    }
  
  
  bool acc_test(double& emitt,double& beta,double& sigma);
  
 public:
  
  BEAM_PARAMETERS();
  
  virtual std::string output_flow() const;
  
  inline int label() const 
    {
      std::istringstream stream1;
      stream1.str(std::string(&extension_[1]));
      int lab;
      stream1 >> lab;
      return lab;
    }  
  
  void read(const PARAMETERS& param);
  
  void setLabel(char label);
  
  inline double ebeam() const {return ebeam_;}
  inline double gamma_beam() const {return gamma_beam_;}
  
  // betas in mm
  inline double beta_x() const { return beta_x_;}
  inline double beta_y() const { return beta_y_;}
  
  inline double n_particles()  const  {return   n_particles_;}     
  // emittances in m.rad
  inline double em_x() const {return em_x_;};
  inline double em_y() const {return em_y_;};
  
  inline double sigma_x() const { return sigma_x_;}
  inline double sigma_y() const { return sigma_y_;}
  
  inline double sigma_z() const 
    {
      return sigma_z_;
    }
  
  inline int dist_x() const {return dist_x_;}
  inline int dist_z() const {return dist_z_;}
  inline double offset_x() const {return offset_x_;}
  inline double offset_y() const {return offset_y_;}
  inline double offset_z() const {return offset_z_;}
  inline double waist_x() const {return waist_x_;}
  inline double waist_y() const {return waist_y_;}
  inline double phi_angle() const {return phi_angle_;}
  inline double x_angle() const {return x_angle_;}
  inline double y_angle() const {return y_angle_;}
  inline double L_0() const {return L_0_;}
  inline double L() const {return L_;}
  inline int bunches_per_train() const {return bunches_per_train_;}
  inline double frep() const {return frep_;}
  inline int trav_focus() const  {return trav_focus_;}
  inline TRIDVECTOR get_polar() const { return polar_;}
  // inline double get_polar_rate() const {return polar_rate_;}
};

#endif
