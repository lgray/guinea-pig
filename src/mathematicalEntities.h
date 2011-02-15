#ifndef MATHEMATICALENTITIES_SEEN
#define MATHEMATICALENTITIES_SEEN

#include <iostream>
#include <sstream>
#include "mathematicalTools.h"
//#include "physicalTools.h"
#include "abstractIOclass.h"
//#include "abstractParticle.h"


/* class QUADRIVECTOR */
/* { */
/*   protected : */

/*   float v1_,v2_,v3_; */
/*   float v4_; */


/*  inline void multByFloat(const float& fac) */
/*    { */
/*      v1_ *= fac;  */
/*      v2_ *= fac;  */
/*      v3_ *= fac;  */
/*      v4_ *= fac;  */
/*    } */


/*   public :  */

/* QUADRIVECTOR() {;} */



/* QUADRIVECTOR(float vx,float vy,float vz, float energy)  */
/*     { */
/*       set(vx, vy, vz, energy); */
/*     } */

/* inline void set (float vx,float vy,float vz, float energy) */
/*   { */
/*     v1_     = vx; */
/*     v2_     = vy; */
/*     v3_     = vz; */
/*     v4_ = energy;          */
/*   } */

/* inline void trivector(float& vx,float& vy,float& vz) const */
/*   { */
/*     vx     =  v1_; */
/*     vy     =  v2_; */
/*     vz     =  v3_; */
/*   } */

/* inline void quadrivector(float& vx,float& vy,float& vz, float& energy) const */
/*   { */
/*     vx     =  v1_; */
/*     vy     =  v2_; */
/*     vz     =  v3_; */
/*     energy =  v4_; */
/*   } */


/*  inline float composante4() const {return v4_;} */
/*  inline float composante3() const {return v3_;} */

/*  inline void switchComp4Sign() {v4_ = -v4_;} */

/*  inline QUADRIVECTOR& operator *= (const float& fac) */
/*    { */
/*      multByFloat(fac); */
/*      return *this; */
/*    } */

/* inline void boost_cecile(float e1, float e2, float beta_x, float beta_y) */
/*   { */
/*     PHYSTOOLS::lorent_pair(e1,e2,v4_, v3_); */
/*     PHYSTOOLS::lorent(v4_ ,v1_ ,beta_x); */
/*     PHYSTOOLS::lorent(v4_ ,v2_ ,beta_y); */
/*   } */
 
/*  inline void momentumToVelocity()  */
/*    { */
/*      if (v4_ == 0.0)  */
/*        { */
/* 	 cout << " mathematicalEntities::momentumToVelocity : null energy! " << endl; */
/* 	 return; */
/*        } */
/*      float factor = 1.0/fabs(v4_); */
/*      v1_ *= factor; */
/*      v2_ *= factor; */
/*      v3_ *= factor; */
/*    } */

/* }; */

class TRIVECTOR
{
  float vec_[3];
  
  public :
    
    TRIVECTOR()
    {
      int k;
      for (k=0; k<3; k++) vec_[k] = 0.0;
    }
  
  TRIVECTOR( const TRIVECTOR& triv)
    {
      setComponents( triv.vec_[0],triv.vec_[1], triv.vec_[2]);
    }
  TRIVECTOR(float x,float y,float z)
    {
      vec_[0] = x;
      vec_[1] = y;
      vec_[2] = z;
    }

  inline float& operator() (int index) { return vec_[index]; }
  inline const float& operator() (int index) const { return vec_[index]; }
  
  inline TRIVECTOR& operator = (const TRIVECTOR& triv)
    {
      setComponents( triv.vec_[0],triv.vec_[1], triv.vec_[2]);
      return *this;
    }
  
  inline TRIVECTOR& operator = (float value)
    {
      setComponents( value,value, value);
      return *this;
    }
  
  inline TRIVECTOR& operator *= (const float& factor)
    {
      int k;
      for (k=0; k < 3; k++) vec_[k] *= factor;
      return *this;
    }
  
  inline float getComponent(int index ) const {return vec_[index];}
  
  inline void setComponent(int index, float value )  {vec_[index] = value;}
  
  inline void getComponents(float& x,float& y,float& z) const
    {
      x = vec_[0];
      y = vec_[1];
      z = vec_[2];
    }
  
  inline void setComponents(float x,float y,float z)
    {
      vec_[0] = x;
      vec_[1] = y;
      vec_[2] = z;
    }
  
  
  inline float norm2() const
    {
      return vec_[0]*vec_[0] + vec_[1]*vec_[1] + vec_[2]*vec_[2];
    }
  inline float norm() const
    {
      return sqrt(abs(norm2()));
    }
  
  inline void renormalize()
    {
      int k;
      float normeInv = 1.0/norm();
      for (k=0; k< 3 ; k++) vec_[k] *= normeInv;
    }
  
  inline void opposite()
    {
      int k;
      for (k=0; k< 3 ; k++) vec_[k] = -vec_[k];
    }
  
  inline void print()
    {
      cout << " x comp. = " << vec_[0] << " y comp. = " << vec_[1] << " z comp. = " << vec_[2] << endl;
    }
};

class TRIDVECTOR
{
  double vec_[3];

  public :
    
    TRIDVECTOR() 
    {
      int k;
      for (k=0; k<3; k++) vec_[k] = 0.0;
    }
  
  TRIDVECTOR( const TRIDVECTOR& triv) 
    {
      setComponents( triv.vec_[0],triv.vec_[1], triv.vec_[2]); 
    }
  TRIDVECTOR( const TRIVECTOR& triv) 
    {
      setComponents( (double)triv(0),(double)triv(1), (double)triv(2)); 
    }
  
  
  TRIDVECTOR(double x,double y,double z)
    {
      vec_[0] = x;
      vec_[1] = y;
      vec_[2] = z;
    }

  inline double& operator() (int index) { return vec_[index]; }
  inline const double& operator() (int index) const { return vec_[index]; }
  
  inline TRIDVECTOR& operator = (const TRIDVECTOR& triv) 
    {
      setComponents( triv.vec_[0],triv.vec_[1], triv.vec_[2]); 
      return *this;
    }
  inline TRIDVECTOR& operator = (const TRIVECTOR& triv) 
    {
      setComponents( (double)triv(0),(double)triv(1), (double)triv(2)); 
      return *this;
    }

  inline TRIDVECTOR& operator = (double value) 
    {
      setComponents( value,value, value); 
      return *this;
    }
  
  inline TRIDVECTOR& operator *= (const double& factor)
    {
      int k;
      for (k=0; k < 3; k++) vec_[k] *= factor;
      return *this;
    }
  
  inline double getComponent(int index ) const {return vec_[index];}
  
  inline void setComponent(int index, double value )  {vec_[index] = value;}
  
  inline void getComponents(double& x,double& y,double& z) const
   {
     x = vec_[0];
     y = vec_[1];
     z = vec_[2];
   }

  inline void setComponents(double x, double y, double z) 
    {
      vec_[0] = x;
      vec_[1] = y;
      vec_[2] = z;
    }
  
  inline void clear()
    {
      vec_[0] = 0.0;
      vec_[1] = 0.0;
      vec_[2] = 0.0;
    }

  inline double norm2() const
    {
      return vec_[0]*vec_[0] + vec_[1]*vec_[1] + vec_[2]*vec_[2];
    }
  inline double norm() const
    {
      return sqrt(abs(norm2()));
    }
  
  inline void renormalize()
    {
      int k;
      double normeInv = 1.0/norm();
      for (k=0; k< 3 ; k++) vec_[k] *= normeInv;
    }

  inline void opposite()
    {
      int k;
      for (k=0; k< 3 ; k++) vec_[k] = -vec_[k];
    }
  
  inline void print()
    {
      cout << " x comp. = " << vec_[0] << " y comp. = " << vec_[1] << " z comp. = " << vec_[2] << endl;
    }
};


class named_int_vector : public ABSTRACT_IO_CLASS
{
  vector<int> vecteur_;
  string name_;

 public :

  named_int_vector() {;}

  inline void put_name(string nom) {name_ = nom;}
  inline void clear()
  {
    name_ = string (" ");
    vecteur_.clear();
  }
  inline void  add_component(int comp) 
  {
    vecteur_.push_back(comp);
  }
  ~named_int_vector() {;}

  virtual string output_flow() const 
  {
    ostringstream sortie;
    string entete = string("vector").append(name_);
    sortie << titre(entete);
    sortie << " vector " << name_ << endl;
    int k;
    for (k=0; k < (int)vecteur_.size(); k++)
      {
	sortie << name_ << " k = " << k << " :" << vecteur_[k] << endl;
      }
    return sortie.str();  
  }
  
  
};


#endif
