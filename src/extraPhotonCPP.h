#ifndef EXTRAPHOTON_SEEN
#define EXTRAPHOTON_SEEN

class  EXTRA
{
  double energy_,vx_,vy_,q2_, eorg_;
  //  double weight_;

 public:

  EXTRA() {;}

  EXTRA(double energy,double vx,double vy,double q2,double eorg) 
    {
      energy_ = energy;
      vx_     = vx;
      vy_      = vy;
      q2_      = q2;
      eorg_    = eorg;
      //         weight_  = weight;
    }

  ~EXTRA() {;}
  inline void get_parameters(double& energy,double& vx,double& vy,double& q2,double& eorg) const
  {
      energy = energy_;
      vx     = vx_;
      vy      = vy_;
      q2      = q2_;
      eorg    = eorg_;
      //         weight  = weight_;
  }
    inline double energy() const {return energy_;}
    //       inline double weight() const {return weight_;}
    inline double q2() const {return q2_;}
    inline double eorg() const { return eorg_;}
    inline  void velocities(double& vx, double& vy) const
    {
      vx = vx_;
      vy = vy_;
    }
};

#endif
