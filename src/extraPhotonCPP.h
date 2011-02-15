#ifndef EXTRAPHOTON_SEEN
#define EXTRAPHOTON_SEEN
#include <vector>
#include <list>
#include <iostream>

using namespace std;


class  EXTRA
{
  float energy_,vx_,vy_,q2_, eorg_;
  //  float weight_;

 public:

  EXTRA() {;}

  EXTRA(float energy,float vx,float vy,float q2,float eorg) 
    {
      energy_ = energy;
      vx_     = vx;
      vy_      = vy;
      q2_      = q2;
      eorg_    = eorg;
      //         weight_  = weight;
    }

  ~EXTRA() {;}
  inline void get_parameters(float& energy,float& vx,float& vy,float& q2,float& eorg) const
  {
      energy = energy_;
      vx     = vx_;
      vy      = vy_;
      q2      = q2_;
      eorg    = eorg_;
      //         weight  = weight_;
  }
    inline float energy() const {return energy_;}
    //       inline float weight() const {return weight_;}
    inline float q2() const {return q2_;}
    inline float eorg() const { return eorg_;}
    inline  void velocities(float& vx, float& vy) const
    {
      vx = vx_;
      vy = vy_;
    }
};


/* class EXTRA_PHOTON */

/* { */


/*  vector< list<EXTRA> > cell_; */

/*  public: */

/*   EXTRA_PHOTON() {;}; */
/* ~EXTRA_PHOTON() */
/*   {;}; */


/* inline void clear_extra_photons() */
/* { */
/*   unsigned int j; */
/*       for (j=0; j < cell_.size(); j++)    cell_[j].clear(); */
/* } */
/* inline void init_extra_photons(int number) */
/* { */
/*   cell_ = vector< list<EXTRA> >(number); */
/* } */

/* inline void store_vir_photon(float energy,float vx,float vy, */
/* 		      float q2,float eorg,float weight,int j) */
/* { */
/*  cell_[j].push_front(EXTRA(energy, vx, vy, q2, eorg, weight)); */
/* } */

/*  inline const list<EXTRA>& get_cell(int j) const {return cell_[j];} */
/* }; */
#endif
