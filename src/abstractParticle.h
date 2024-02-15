#ifndef ABSTRACTPARTICLE_SEEN
#define ABSTRACTPARTICLE_SEEN
#include <iostream>
#include <sstream>
#include <string>

#include "typeDefs.h"
#include "abstractIOclass.h"
#include "mathematicalEntities.h"

class ABSTRACT_BEAM
{

 protected:
 public: 
  ABSTRACT_BEAM() {;}
  virtual ~ABSTRACT_BEAM() {;}
  virtual int sizeOfPhotonSlice(int slice) const = 0;
  virtual int number_of_slices() const =0;
  virtual int numberOfParticlesOfSlice(int slice) const = 0;
  
};


class ABSTRACT_PARTICLE
{
  
 protected: 
  
  // unit of xpos_, ypos_ : nm 
  double xpos_, ypos_;
  
  double zpos_;
  
  // unit of velocity : radian (?)
  double vx_,vy_;

  // energy in GeV
  double energy_;

  inline void set(double energyi, double xi, double yi,double zi,double vxi,double vyi)
    {
      energy_ = energyi;
      xpos_= xi;
      ypos_ = yi;
      zpos_ = zi;
      vx_ = vxi;
      vy_ = vyi;
    }
  
  inline void set(const ABSTRACT_PARTICLE& part)
    {
      set (part.energy_, part.xpos_, part.ypos_, part.zpos_, part.vx_, part.vy_);
    }
  
  
  // internal units seems to be nm for x,y,z
  inline void transform_to_internal_units(double facPosition, double facVelocity)
    {
      xpos_ *= facPosition;
      ypos_ *= facPosition;
      zpos_ *= facPosition;
      
      vx_ *= facVelocity;
      vy_ *= facVelocity;
    }
  
  
  ABSTRACT_PARTICLE (double xi, double yi,double zi, double vxi,double vyi,double energyi)
    {
      set(energyi, xi, yi,zi, vxi, vyi);
    }

  
  public : 
    
    ABSTRACT_PARTICLE () 
    {
      set(0.0, 0.0, 0.0, 0.0, 0.0, 0.0); 
    }
  
  virtual ~ABSTRACT_PARTICLE () {;}
  
  inline double z() const {return zpos_;}
  inline void setZ( double zi) { zpos_ = zi;};
  
  virtual const TRIDVECTOR& getSpin() const = 0;
  
  
  inline void razEnergy() {energy_ = 0.0;}
  
  
  inline void velocities(double& vxout, double& vyout) const 
    {
      vxout = vx_;
      vyout = vy_;
    };
  
  inline void XYposition(double& xpos, double& ypos) const
    {
      xpos = xpos_;
      ypos = ypos_;
    };
  
  inline double energy() const {return energy_;};
  
  inline void setEnergy(double en) {energy_ = en;}; 
  
  inline void advancePosition(double distance)
    {
      xpos_ += vx_*distance;
      ypos_ += vy_*distance;
    };
  
  inline void advanceVelocities(double Fx, double Fy, double step)
  {
    vx_ += step*Fx/energy_;
    vy_ += step*Fy/energy_;
  }
  
  inline void advanceVelocities(double Fx, double Fy, double step, double& deltaVx, double& deltaVy)
    {
      deltaVx = step*Fx/energy_;
      deltaVy = step*Fy/energy_;
      vx_ += deltaVx;
      vy_ += deltaVy;
    }

  
  inline void set_z_for_dump(int istep,double dz0, double max_z, int sign_label ) 
    {
      double z,dz;
      double sign = (double)sign_label;
      
      z = zpos_;
      dz = max_z-istep*dz0;
      //  std::cout << " set_z_for_dump istep= " << istep <<  " dz= " << dz << std::endl;
      z += dz;
      zpos_ = sign*(z + dz0);
    }
  
  inline void set_z_velocity_corrected_for_dump(int slice, int istep,int complement, double dz0, double max_z, int sign_label ) 
    {
      double dz;
      //      double vx, vy unused;
      set_z_for_dump(istep, dz0, max_z, sign_label);
      dz=(slice-istep + complement -1)*dz0;
      advancePosition(-dz);
    }
};



class ABSTRACT_BHABHA_PHOTON_SAMPLES
{
 public :

   ABSTRACT_BHABHA_PHOTON_SAMPLES() {;}
 virtual  ~ABSTRACT_BHABHA_PHOTON_SAMPLES() {;}
 virtual  int get_label() const = 0;
 virtual  unsigned int nb_samples() const = 0;
 
 virtual  void get_parameters_for_output(unsigned int number, int& number_bhabha, double& en,double& vx,double& vy, double&vz) const = 0;
 
 virtual  void add_bhabha_photon(int nbhabha, double px, double py, double pz, double en) = 0;
// virtual  void create_bhabha_photon(int nbhabha, double px, double py, double pz, double en) = 0;
 
 
};

class ABSTRACT_BHABHASAMPLES
{
  public : 
    ABSTRACT_BHABHASAMPLES() {;}
  virtual ~ABSTRACT_BHABHASAMPLES() {;}
  virtual unsigned int nb_samples() const = 0;
  
  virtual  void get_parameters_for_output(unsigned int number, unsigned int& evtIdx, double& eCM, double& mother1_en,double&e1,double&vx1,double& vy1, double&vz1, double& mother2_en, double& e2, double& vx2, double&vy2, double&vz2, int& nbphot) const = 0;
  
  virtual  void add_bhabha(unsigned int evtIdx, double px1, double py1, double pz1, double e1, double px2, double py2, double pz2, double e2, int nbphot) = 0;
  
};

class ABSTRACT_CROSS_DATA
{
  
  public :
    
    ABSTRACT_CROSS_DATA() {;}
  virtual ~ABSTRACT_CROSS_DATA() {;}
  virtual void resize(int n, int nval) = 0;
  virtual void add_data(double ener, const double* data) =0;
  
};


class ABSTRACT_LUMI_HEAP : public ABSTRACT_IO_CLASS
{

 public:
  ABSTRACT_LUMI_HEAP() {;}
  virtual ~ABSTRACT_LUMI_HEAP() {;}

  virtual  int nb_pairs() const = 0;
  virtual  void get_parameters_for_output(unsigned int number, double& e1,double& e2,double& x,double& y,double& z) const = 0;
  virtual  void get_parameters_for_output(unsigned int number, double& e1,double& e2,double& x,double& y,double& z, double& vx1,double& vy1,double& vx2,double& vy2, double& sx1, double& sy1, double& sz1, double& sx2, double& sy2, double& sz2,int& t)  const = 0;
  virtual std::string output_flow() const 
    {
      std::ostringstream out;
      out << " ABSTRACT_LUMI_HEAP:: no data for output file in abstract class " << std::endl;
      return out.str();
    }
};

class PARTICLE_INTERFACE : public ABSTRACT_PARTICLE
{
  friend class BEAM_FROM_FILE;
  protected : 
    
    
    TRIDVECTOR polar_;
  
  
  
 public:
  
  PARTICLE_INTERFACE() {;}
  virtual ~PARTICLE_INTERFACE() {;}
  
  inline bool good() const
    {
    bool test = true;
    // if the energy is 0 or negative the reading of the particle on the
    // file has probably failed. This test is made because the program
    // has to divide by energy.
    if ( energy_ < 1.9e-12 ) test = false;
    return test;
  }

 inline void init_from_input_file(double energy, double xpos, double ypos, double zpos, double vx, double vy, double polx, double poly, double polz)
   {
     set(energy, xpos, ypos,zpos,vx,vy);
     polar_.setComponents((double)polx, (double)poly, (double)polz);
   }



 inline void get_parameters(double& energy, double& xpos, double& ypos, double& zpos, double& vx, double& vy) const
 {
   xpos = xpos_;
   ypos = ypos_;
   zpos = zpos_;
   vx = vx_;
   vy = vy_;
   energy = energy_;
 }



 inline double get_helicity() const { return polar_(0);}
 
 virtual inline const TRIDVECTOR& getSpin() const 
   { 
     return polar_;
   }
 
};


#endif
