#ifndef PARTICLES_SEEN
#define PARTICLES_SEEN
#include <iostream>
#include <cmath>
#include "define.h"
#include "typeDefs.h"
#include "abstractIOclass.h"
#include "abstractParticle.h"
#include "physicalTools.h"
#include "mathematicalEntities.h"
#include "tridentCPP.h"

class PARTICLE : public ABSTRACT_PARTICLE, public ABSTRACT_IO_CLASS
{

  protected :

    double ups_;
    TRIDENT trident_;
  

  PARTICLE (double xi, double yi, double zi, double vxi,double vyi, double energyi) : ABSTRACT_PARTICLE (xi, yi, zi,vxi, vyi,energyi)  {;}
  
     
 public:

  PARTICLE() : ABSTRACT_PARTICLE(),ABSTRACT_IO_CLASS(),ups_(0.0) {;}

  // return the quantity dz/radius
  inline double advanceDueToEBfield(TRIVECTOR EBfield, double distance, double scal_step_local)
    {
      double deltaVx, deltaVy;
      // divide by energy in GeV.
      // the total field (electric+magnetic) is equivalent to 
      // 2 X electric field

      advanceVelocities(-2.0*EBfield(0), -2.0*EBfield(1), distance, deltaVx, deltaVy);
      advancePosition(distance*scal_step_local); 
      double dz = distance*1e-9*scal_step_local;
      double dzOnRadius = sqrt(deltaVx*deltaVx+deltaVy*deltaVy);
      double radius_i = dzOnRadius/dz;
      double upsilon=LAMBDA_BAR/(EMASS*EMASS)*energy_*energy_*radius_i;
      ups_ = upsilon;
      return dzOnRadius;
    }

  inline double advanceTridentDueToEBfield(TRIVECTOR EBfield, double distance, double scal_step_local)
    {
      double deltaVx, deltaVy;
      // divide by energy in GeV.
      // the total field (electric+magnetic) is equivalent to 
      // 2 X electric field
      
      advanceVelocities(-2.0*EBfield(0), -2.0*EBfield(1), distance, deltaVx, deltaVy);
      advancePosition(distance*scal_step_local); 
      double dz = distance*1e-9*scal_step_local;
      double dzOnRadius = sqrt(deltaVx*deltaVx+deltaVy*deltaVy);
      double radius_i = dzOnRadius/dz;
      double upsilon=LAMBDA_BAR/(EMASS*EMASS)*energy_*energy_*radius_i;
      ups_ = upsilon;
      return dzOnRadius;
    }

  inline void apply_position_offset(double dx, double dy)
    {
      xpos_ += dx;
      ypos_ += dy;
    }
  
  
  inline void rotate(double cosa, double sina)
    {
      double x = xpos_;
      double y = ypos_;
      xpos_ = cosa*x + sina*y;
      ypos_ = -sina*x + cosa*y;
      double vx = vx_;
      double vy = vy_;
      vx_ = cosa*vx + sina*vy;
      vy_ = -sina*vx + cosa*vy;
    }


  inline void adjustEnergyToMinimum(double emin)
    {
      if (energy_ < emin)
	{
	  std::cout << "PARTICLE:: e_low2: " << energy_ << std::endl;
	  energy_ = emin;
	}
    };


  virtual inline void setParticle(double x,double y,double z,double vx,double vy,double energy)
    {
      ABSTRACT_PARTICLE::set(energy, x, y,z,vx,vy);
    }
  
  virtual inline void setParticle(const PARTICLE_INTERFACE& part)
    {
      ABSTRACT_PARTICLE::set(part);
    }


  virtual inline void setParticle(double x,double y,double z,double vx,double vy,double energy, const TRIDVECTOR* /*dummy*/)
    {
      setParticle(x,y,z,vx,vy,energy);
    }
  inline double getEnergy() const {return energy_;}
  inline void setUps(double upsilon) { ups_ = upsilon;}
  inline double getUps() const { return ups_;}
  
  
  virtual inline const TRIDVECTOR& getSpin() const
    {
      std::cerr << " PARTICLE::get_spin : a PARTICLE has no spin, use PARTICLE_WITH_SPIN " << std::endl;
      exit(0);
/*       const TRIDVECTOR tv(0.0,0.0,0.0); */
/*       return tv; */
    }


  virtual inline void rotateBMT(TRIDVECTOR /*Efield*/, TRIDVECTOR /*Bfield*/, double /*charge_sign*/, double /*dt*/) {;}
  
  inline void beamstrahlung(TRIDVECTOR /*Efield*/, TRIDVECTOR /*Bfield*/, double dzOnRadius,  double emin, std::vector<double>& photonEnergies, RNDM& rndm_generator)
    {
      energy_ = PHYSTOOLS::synrad_no_spin_flip(ups_, energy_, dzOnRadius ,photonEnergies, rndm_generator);
      adjustEnergyToMinimum(emin);
    }
    
  inline void createTridents(double dz,std::vector<double>* electrons,std::vector<double>* positrons,std::vector<double>* virt,RNDM& rndm_generator)
    {
      trident_.createTridents(&energy_,ups_,dz,electrons,positrons,virt,rndm_generator);
    }

  inline void createTridents(double dz,std::vector<double>* electrons,std::vector<double>* positrons,RNDM& rndm_generator)
    {
      trident_.createTridents(&energy_,ups_,dz,electrons,positrons,rndm_generator);
    }
  
  virtual inline void beamstrahlungSokolov(TRIDVECTOR /*Efield*/, TRIDVECTOR /*Bfield*/, double /*dzOnRadius*/,  double /*emin*/, double /*charge_sign*/, std::vector<double>& /*photonEnergies*/, RNDM& /*rndm_generator*/)
    {
      std::cerr << " illegal call of PARTICLE::beamstrahlungSokolov " << std::endl;
      exit(0);
    }
  

  virtual std::string  output_flow() const
    {
      std::cerr << "  output_flow not programmed for class PARTICLE " << std::endl;
      return std::string(" ");
    }

  virtual std::string persistent_flow() const
    {
      std::ostringstream out;
      out <<  energy_ << " " << xpos_*1e-3 << " " << ypos_*1e-3 << " " << zpos_*1e-3 << " " << vx_*1e6 << " " << vy_*1e6;
      return out.str();
    }
  

  virtual  inline void init_from_input_file(double energy, double xpos, double ypos, double zpos, double vx, double vy) 
    {
      double facPosition = 1.0e3;
      double facVelocity = 1.0e-6;
      set(energy, xpos, ypos,zpos,vx,vy);
      // transform to nm and rad
      transform_to_internal_units(facPosition, facVelocity);
    }

};


class PARTICLE_WITH_SPIN : public PARTICLE
{
  protected :
  TRIDVECTOR spin_;

  PARTICLE_WITH_SPIN (double xi, double yi, double zi, double vxi,double vyi, double energyi, TRIDVECTOR sp) : PARTICLE (xi, yi, zi, vxi, vyi, energyi) 
    {
      spin_ = sp;
    }
  
  
 public:
  
  PARTICLE_WITH_SPIN() : PARTICLE() {;}
  
 
  virtual inline const TRIDVECTOR& getSpin() const
    {
      return spin_;
    }
  
  
  virtual inline void setParticle(double x,double y,double z,double vx,double vy,double energy, const TRIDVECTOR* polar)
    {
      ABSTRACT_PARTICLE::set(energy, x, y,z,vx,vy);
      spin_ = *polar;
    }
  
  
  virtual inline void setParticle(const PARTICLE_INTERFACE& part)
    {
      ABSTRACT_PARTICLE::set(part);
      spin_ = part.getSpin();
    }
  
  
  inline void setSpin(TRIDVECTOR sp)
    {
      //      std::cout << " appel setSpin " << std::endl;
      spin_ = sp;
    }


  virtual  inline void init_from_input_file(double /*energy*/, double /*xpos*/, double /*ypos*/, double /*zpos*/, double /*vx*/, double /*vy*/) 
    {
      std::cerr << " reading from input file not programmed for PARTICLE_WITH_SPIN " << std::endl;
      exit(0);
      
    }

  virtual void rotateBMT(TRIDVECTOR Efield, TRIDVECTOR Bfield, double charge_sign, double dt);
  
  
  virtual inline void beamstrahlungSokolov(TRIDVECTOR Efield, TRIDVECTOR Bfield, double dzOnRadius,  double emin, double charge_sign, std::vector<double>& photonEnergies, RNDM& rndm_generator)
    {
      TRIDVECTOR ev1, ev2, ev3; 
      PHYSTOOLS::referenceSpin((double)vx_, (double)vy_, ev1, ev2, ev3, Efield, Bfield, charge_sign);
      energy_ = PHYSTOOLS::synrad_spin_flip(ups_,energy_, ev1, ev2, ev3, spin_,dzOnRadius,photonEnergies, rndm_generator);
      adjustEnergyToMinimum(emin);
    }
  
  
  virtual std::string persistent_flow() const
    {
      std::ostringstream out;
      out <<  PARTICLE::persistent_flow() << " " << spin_(0) << " " << spin_(1) << " " << spin_(2);
      return out.str();
    }
  
  
};

class PHOTON : public ABSTRACT_PARTICLE, public ABSTRACT_IO_CLASS
{
  TRIDVECTOR stokes_;
  double helicity_;
  
  int no_;
  
  
 public:
  
  PHOTON() 
    {
      helicity_=0.;
      no_=0;
    }

  PHOTON(  double xi, double yi,  double zi, double vxi,double vyi, double energyi, double helicityi, int no) : ABSTRACT_PARTICLE (xi, yi, zi, vxi, vyi,energyi)
    {
      helicity_ = helicityi;
      no_ = no;
    }
    
  PHOTON(double energy, const PARTICLE& particle, double helicityi, int no)
    {
      
      double x,y,vx,vy;
      particle.XYposition(x,y);
      particle.velocities(vx,vy);
      xpos_ = x;
      ypos_ = y;
      zpos_ = particle.z();
      vx_ = vx;
      vy_ = vy;
      energy_ = energy;
      helicity_ = helicityi;
      no_ = no;
    }
  
  virtual std::string  output_flow() const
    {
      std::cerr << "  output_flow not programmed for class PHOTON " << std::endl;
      return std::string(" ");
    }

  virtual std::string persistent_flow() const
    {
      std::ostringstream out;
      /*   out <<  energy_ << " " << xpos_*1e-3 << " " << ypos_*1e-3 << " " << zpos_*1e-3 << " " << vx_*1e6 << " " << vy_*1e6 << " " << helicity_ << " " << no_; */
      // from original, to be modified perhaps
      out <<  energy_ << " " <<  vx_ << " " << vy_<< " ";
      
      return out.str();
    }

  virtual  inline void init_from_input_file(double energy, double xpos, double ypos, double zpos, double vx, double vy) 
    {
      set(energy, xpos, ypos,zpos,vx,vy);
      std::cerr << " reading from input file to be checked for photons " << std::endl;
      exit(0);
    } 
  
  virtual  inline void init_from_input_file(double energy, double xpos, double ypos, double zpos, double vx, double vy, double hel) 
    {
      
      PHOTON::init_from_input_file(energy, xpos, ypos, zpos, vx, vy);
      helicity_ = hel;
      //  no_ = 0;
      std::cerr << " reading from input file to be checked for photons " << std::endl;
      exit(0);
    } 
  
  virtual inline const TRIDVECTOR& getSpin() const 
    { 
      std::cerr << " PHOTON::getSpin : not programmed " << std::endl;
      exit(0);
/*       const TRIDVECTOR tv(0.0,0.0,0.0); */
/*       return tv; */
    }
  
  
  inline int no() const {return no_;} 
  
};

inline bool operator < ( const PARTICLE& p1, const PARTICLE& p2)
{
  return (p1.z() < p2.z());
}


class PAIR_PARTICLE :  public ABSTRACT_PARTICLE, public ABSTRACT_IO_CLASS
{
  
  double velocityz_;
  
  int process_;
  int label_;

  int icharge_;
  double e_inv_;
  bool last_rescaling_ok_;

  // MANU, to retrieve the time :
   int slice1_ ;
   int slice2_ ;

 public:
  
  PAIR_PARTICLE() : ABSTRACT_PARTICLE ()
    {
      velocityz_ =0.0;
      e_inv_ = 0.0;
      icharge_ = 0;
      process_ = -99;
      label_ = -99;
      last_rescaling_ok_ = true;

      // -- MANU :
	slice1_ = -1 ;
	slice2_ =-1 ;
      // ----
    }
  
  PAIR_PARTICLE(int label, int index_of_process, double x,double y,double z, double vx,double vy,double vz, double energy) : ABSTRACT_PARTICLE (x,y,z,vx,vy,energy)
    {
      velocityz_ =vz;
      if (energy > 0.0) 
	{
	  icharge_ = 0;
	}
      else 
	{
	  icharge_ = 1;
	  energy_ = -energy;
	}
      e_inv_ = 1.0/energy_;
      process_ = index_of_process;
      label_ = label;
      last_rescaling_ok_ = true;
    }
  virtual ~PAIR_PARTICLE () {;}
  
  virtual const TRIDVECTOR& getSpin() const
  {
    std::cerr << " PARTICLE::get_spin : a PAIR_PARTICLE has no spin, use PARTICLE_WITH_SPIN " << std::endl;
    exit(0);
/*     const TRIDVECTOR tv(0.0,0.0,0.0); */
/*     return tv; */
  }
    
  virtual std::string  output_flow() const
    {
      std::cerr << "  output_flow not programmed for class PAIR_PARTICLE " << std::endl;
      return std::string(" ");
    }
  
  
  virtual std::string persistent_flow() const
  {
    std::ostringstream out;
    double ener;
    if ( icharge_) ener = -energy_;
    else ener = energy_;

    out <<  ener << " " <<  vx_ << " " << vy_ << " " << velocityz_ << " " << xpos_ << " " << ypos_ << " " << zpos_ << " " << process_ << " " <<  slice1_ <<" " << slice2_  ;
    return out.str();
  }

  inline int get_label() const {return label_;}
  inline int get_process() const {return process_;}
  
  inline bool last_rescaling_ok() const {return  last_rescaling_ok_;}
  inline double unsigned_energy() const { return energy_;}
  
  inline double signed_energy() const
    {
      double ener;
      if ( icharge_) ener = -energy_;
      else ener = energy_;
      return ener;
    }

  inline void set_slices12( int i1, int i2 ){ 
	slice1_ = i1 ;
	slice2_ = i2;
	return ;
  }

  inline void get_slices12( int &i1,  int &i2 ) {
	i1 = slice1_ ;
	i2 = slice2_ ;
	return;
  }



  inline double Zvelocity() const 
    {
      return  velocityz_;
    }
  
  inline double velocity_q() const 
    {
      return vx_*vx_ + vy_*vy_ + velocityz_*velocityz_;
    }
  

  
  inline void advancePosition(double step)
    {
      ABSTRACT_PARTICLE::advancePosition(step);
      zpos_ += velocityz_*step;
    }

  void apply_electric_field_on_pair(double step_2, double ex, double ey); 
  
  double  scale_pair(double vold2, double mass);
  
  inline void scale_velocities_sync(double scal, double& vx0, double& vy0, double& vz0) const
    {
#ifdef SCALE_ENERGY
#ifdef PAIR_SYN
  vx0 *= scal;
  vy0 *= scal;
  vz0 *= scal;
#endif
#endif
    }


  void synchrotron_radiation(double step, double mass, double step_fraction, double vx0, double vy0, double vz0, std::vector<double>* photon_e, RNDM& rndm_generator);

  void synchrotron_radiation(double step, double mass, double step_fraction, double vx0, double vy0, double vz0, RNDM& rndm_generator);
 
  double apply_magnetic_field_on_pair(double fac_theta, double step_q, double bx, double by);

  double apply_final_half_step_fields(double step, double mass, double ex,double ey, double bx, double by, double thetamx, RNDM& rndm_generator);

  double apply_initial_half_step_fields(double step, double mass, double ex,double ey, double bx, double by, std::vector<double>* photon_e, RNDM& rndm_generator);

  double apply_initial_half_step_fields(double step, double mass, double ex,double ey, double bx, double by, RNDM& rndm_generator);
 
  double apply_full_step_fields(double step, double mass, double ex,double ey, double bx, double by, std::vector<double>* photon_e, RNDM& rndm_generator); 

  double apply_full_step_fields(double step, double mass, double ex,double ey, double bx, double by, RNDM& rndm_generator); 

  void update_energy(double mass);


};

class  PHOTON_DATA
{
  int i_,n_,j_;
  double scal_,scal_i_,scal2_,scal2_i_,emin_;
  
  
 public:
  
  PHOTON_DATA()
    {
      n_       = 1;
      j_       = 1;
      scal_    = 1.0*n_;
      scal_i_  = 1.0/scal_;
      emin_    = 50.0;
      scal2_   = 1.e3;
      scal2_i_ = 1.0;
    }
  
  
  inline void updateXX() 
    { 
      if (!(--j_)) j_ = n_; 
      i_ = j_; 
    }

  inline void nTOi() 
    {
      i_ = n_; 
    }

  inline double get_scal2() const 
    {
      std::cout << "PHOTON_DATA : scal2 " << std::endl; 
      return scal2_;
    }

};

class PHOTON_COUNTER
{
  double esum_;
  int number_;
  
 public:
  
  PHOTON_COUNTER()
    {
      esum_ = 0.0;
      number_ = 0;
    }
  inline void addPhoton(double energy)
    {
      esum_ += energy;
      number_++;
    };
  

  inline double getSum() const
    {
      return esum_;
    }
  inline int getNumber() const 
    {
      return number_;
    }
};

// Photons from incoherent pairs
class TERTPHOTON : public ABSTRACT_PARTICLE, public ABSTRACT_IO_CLASS
{
  double vz_;
 
 public:
  
  TERTPHOTON() 
    {
      vz_=0.;
    }

  TERTPHOTON(  double xi, double yi,  double zi, double vxi,double vyi, double vzi, double energyi) : ABSTRACT_PARTICLE (xi, yi, zi, vxi, vyi,energyi)
    {
      vz_=vzi;
    }
    
  TERTPHOTON(double energy, const PAIR_PARTICLE& particle)
    {
      double x,y,vx,vy,norm;
      particle.XYposition(x,y);
      particle.velocities(vx,vy);
      xpos_ = x;
      ypos_ = y;
      zpos_ = particle.z();
      vx_ = vx;
      vy_ = vy;
      vz_=particle.Zvelocity();
      norm=sqrt(vx_*vx_+vy_*vy_+vz_*vz_);
      vx_/=norm;
      vy_/=norm;
      vz_/=norm;
      energy_ = energy;
    }
  
  virtual std::string  output_flow() const
    {
      std::cerr << "  output_flow not programmed for class TERTPHOTON " << std::endl;
      return std::string(" ");
    }

  virtual std::string persistent_flow() const
    {
      std::ostringstream out;
      /*   out <<  energy_ << " " << xpos_*1e-3 << " " << ypos_*1e-3 << " " << zpos_*1e-3 << " " << vx_*1e6 << " " << vy_*1e6 << " " << helicity_ << " " << no_; */
      // from original, to be modified perhaps
      out <<  energy_ << " " <<  vx_ << " " << vy_<< " " << vz_ << " "<<  xpos_ << " " << ypos_<< " " << zpos_ << " ";
      
      return out.str();
    }
  
  virtual inline const TRIDVECTOR& getSpin() const 
    { 
      std::cerr << " PHOTON::getSpin : not programmed " << std::endl;
      exit(0);
/*       const TRIDVECTOR tv(0.0,0.0,0.0); */
/*       return tv; */
    }
};

#endif
