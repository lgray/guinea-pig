#ifndef RESULTS_SEEN
#define RESULTS_SEEN
#include <iostream>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <algorithm>

#include "LesDefines.h"
#include "physconst.h"
#include "mathconst.h"
//#include "hitCPP.h"
//#include "meshCPP.h"
//#include "pairsCPP.h"
#include "beamCPP.h"
#include "switchesCPP.h"
#include "abstractParticle.h"

using namespace std;


/* Definition of variables to store the results */

class  RESULTS : public ABSTRACT_IO_CLASS
{
  const SWITCHES* switches_;
  double lumi_0;
  double lumi_fine,lumi_ee,lumi_pp,lumi_eg,lumi_ge,lumi_gg,lumi_gg_high,lumi_ee_high;
  double lumi[3][3];
  double lumi_ecm,lumi_ecm2,lumi_ecm3,lumi_ecm4;
  double hadrons_ee[3],hadrons_eg[3],hadrons_ge[3],hadrons_gg[3];
  double minijets;
  //  double temp_1,temp_2,temp_3;
  double eloss_1,eloss_2;
  double ephot_1,ephot_2;
  double c_vx_1,sig_vx_1,c_vy_1,sig_vy_1;
  double c_vx_2,sig_vx_2,c_vy_2,sig_vy_2;
  double c_vx_1_coh,c_vx_2_coh,c_vy_1_coh,c_vy_2_coh;
  double upsmax;

  
  
 public:
   
  RESULTS();
  
  inline void cumulate3array(double source[3], double cible[3])
    {
      int k;
      for (k=0; k<3; k++) cible[k] += source[k];
    }
  
  inline void connect_switches(const SWITCHES* sw)
    {
      switches_ = sw;
    }
  
  inline void get_velocity_quantities1(double& c_vx, double& sig_vx, double& c_vy, double& sig_vy)
    {
      c_vx = c_vx_1;
      sig_vx = sig_vx_1;
      c_vy = c_vy_1;
      sig_vy = sig_vy_1;
    }
  inline void get_velocity_quantities2(double& c_vx, double& sig_vx, double& c_vy, double& sig_vy)
    {
      c_vx = c_vx_2;
      sig_vx = sig_vx_2;
      c_vy = c_vy_2;
      sig_vy = sig_vy_2;
    }

  //inline void raz_results_temp()
  //{
  //  temp_1=0.0;
  //  temp_2=0.0;
  //  temp_3=0.0;
  //}
  
  inline void add_lumi(int firstIndex, int secondIndex, double x)
    {
      lumi[firstIndex][secondIndex]+= x;
    }
  
  inline void add_lumi_ee(double x)
    {
      lumi_ee+= x;
    }
  
  inline void add_lumi_ge(double x)
    {
      lumi_ge += x;
    }
  
  inline void add_lumi_eg(double x)
    {
      lumi_eg += x;
    }
  
  inline void add_lumi_gg(double x)
    {
      lumi_gg += x;
    }
  
  inline void add_lumi_gg_high(double x)
    {
      lumi_gg_high += x;
    }
  
  
  inline void add_lumi_pp(double x)
    {
      lumi_pp+= x;
    }
  
  inline void add_lumi_ee_high(double x)
    {
      lumi_ee_high += x;
    };

  inline void add_lumi_ecm(double x)
    {
      lumi_ecm += x;
    }
  
  inline void add_lumi_ecm2(double x)
    {
      lumi_ecm2 += x;
    }
  
  /*
    inline void add_temp_1(double x)
    {
    temp_1 += x;
    }
    inline void add_temp_2(double x)
    {
    temp_2 += x;
    }
    inline void add_temp_3(double x)
    {
    temp_3 += x;
    }
  */

  
  
  inline void updateUpsmax(float x) { upsmax=max(upsmax,(double)x);};
  
  // cette methode n'est pas terrible, voir qui s'en sert et comment l'eviter
  inline void cumulate_hadrons_ee(double temp[3])
    {
      cumulate3array(temp, hadrons_ee);
      //  hadrons_ee[0] = temp_1; 
      //  hadrons_ee[1] = temp_2; 
      //  hadrons_ee[2] = temp_3; 
    }
  
  // cette methode n'est pas terrible, voir qui s'en sert et comment l'eviter
  inline void cumulate_hadrons_ge(double temp[3])
    {
      cumulate3array(temp, hadrons_ge);
      //  hadrons_ge[0] = temp_1; 
      //  hadrons_ge[1] = temp_2; 
      //  hadrons_ge[2] = temp_3; 
    }
  
  // cette methode n'est pas terrible, voir qui s'en sert et comment l'eviter
  inline void  cumulate_hadrons_eg(double temp[3])
    {
      cumulate3array(temp, hadrons_eg);
      //  hadrons_eg[0] = temp_1; 
      //  hadrons_eg[1] = temp_2; 
      //  hadrons_eg[2] = temp_3; 
    }
  
  // cette methode n'est pas terrible, voir qui s'en sert et comment l'eviter
  inline void  cumulate_hadrons_gg(double temp[3])
    {
      cumulate3array(temp, hadrons_gg);
      //  hadrons_gg[0] = temp_1; 
      //  hadrons_gg[1] = temp_2; 
      //  hadrons_gg[2] = temp_3; 
    }
  
  
  
  inline void add_lumi_fine(double x)
    {
      lumi_fine += x;
    }

  void bpm_signal(const BEAM& beam);
  
  void bpm_signal_coherent(const BEAM& beam) ;

  
  virtual inline string name_of_class() const 
    {
      return string("RESULTS");
    }
  

  void output_flow(ostringstream& out ) const ;
  virtual string  output_flow() const ;

};

class PAIRS_RESULTS : public ABSTRACT_IO_CLASS
{
  int number_;
  double energy_;
  double eproc_[3],nproc_[3];
  double b1_,b2_,n1_,n2_;
  string name_;

  double highptsum_, highpteng_;


  void set();
  void update_contribution(int composante, double ener,double px,double py,double pz, double wgt);
  
  
 public: 
  
 PAIRS_RESULTS()  {set();}
 
 //  PAIRS_RESULTS() : secondaries_pointer_(NULL) {set();}
 //  PAIRS_RESULTS(PAIR_BEAM* secondaries) {secondaries_pointer_ = secondaries; set();}

 inline void set_name(string name){name_=name;}
 inline int number() const {return number_;}
 inline double energy() const {return energy_;}
 
 void store_full_pair(int composante, double e1,double px1,double py1,double pz1,double e2,double px2,double py2,double pz2, double wgt,  bool luckyPair );
 
 inline void storep_(int composante, double e,double wgt, bool lucky)
   {
     update_contribution(composante, e,0.0,0.0,e, wgt);
     eproc_[composante]+= fabs(e) * wgt;
     nproc_[composante]+= wgt;
     if (lucky)
       {
	 number_++;
	 energy_ += fabs(e);
       }
   } /* storep_ */
 
 virtual inline string name_of_class() const 
   {
     return string("PAIRS_RESULTS");
   }
 
 string output_flow() const ;
 
 
};

class COMPT_RESULTS : public ABSTRACT_IO_CLASS
{
  int number_;
  double energy_;
  double eproc_[2],nproc_[2];
  double b1_,b2_,n1_,n2_;


 public:
  
  COMPT_RESULTS();
  
  int store_compt(int composante,double e,double px,double py,double pz,double wgt,RNDM& hasard);

  
  virtual inline string name_of_class() const 
    {
      return string("COMPT_RESULTS");
    }
  
  virtual  inline string output_flow() const 
    {
      ostringstream out;
      out << title(string("Compton results"));
      out << "compt_npart_0 = " << nproc_[0] << " compt_epart_0 = " << eproc_[0] << endl;
      out << "compt_npart_1 = " << nproc_[1] << " compt_epart_1 = " << eproc_[1] << endl;
      out << "compt_n.1 = " << n1_ << " compt_b.1 = " << b1_ << endl;
      out << "compt_n.2 = " << n2_ << " compt_b.2 = " << b2_ << endl;
      return out.str();
    }
};

/*
  class MUON_RESULTS : public ABSTRACT_IO_CLASS
  {
  double eproc[3],nproc[3];
  
  public:
  
  MUON_RESULTS();
  
  inline string output_flow() const 
  {
  ostringstream out;
  out << title(string("Muons results"));
  out << "muon_n.0 = " << nproc[0] << " muon_e.0 = " << eproc[0] << endl;
  out << "muon_n.1 = " << nproc[1] << " muon_e.1 = " << eproc[1] << endl;
  out << "muon_n.2 = " << nproc[2] << " muon_e.2 = " << eproc[2] << endl;
  
  out << "muon_n_sum = " << nproc[0] + nproc[1]+nproc[2] <<  " muon_e_sum = " << eproc[0]+eproc[1] + eproc[2] << endl;
  
  return out.str();
  }
  };
*/

class JET_RESULTS : public ABSTRACT_IO_CLASS
{
  double sigma_[3];
  
  
 public:
  
  JET_RESULTS()
    {
      int k;
      for (k=0; k<3; k++) sigma_[k] = 0.0;
    }
  
  inline void increment_sigma(int index, double x)
    {
     sigma_[index] += x;
    }
  
  virtual inline string name_of_class() const 
    {
      return string("JET_RESULTS");
    }
  
  
  virtual inline string output_flow() const 
    {
      ostringstream out;
      out << title(string("jets results"));
      out << " nb of minijet evts per bx due to the direct process : jets0 = " << sigma_[0] << endl;
      out << "nb of minijet evts per bx due to the once resolved process : jets1 = " << sigma_[1] << endl;
      out << "nb of minijet evts per bx due to the twice resolved process : jets2 = " << sigma_[2] << endl;
      return out.str();
    }
};


/* class  PHOTON_RESULTS : public ABSTRACT_IO_CLASS */
/* { */
/*     double energy_[2]; */
/*     unsigned long int number_[2];  */
/*  public: */
/* PHOTON_RESULTS() */
/*   { */
/*     energy_[0]=0.0; */
/*     energy_[1]=0.0; */
/*     number_[0]=0; */
/*     number_[1]=0; */
/*   }; */
/*  inline void addPhoton(int nbeam, float energy) */
/*       { */
/* 	number_[nbeam-1]++; */
/* 	energy_[nbeam-1]+=energy; */
/*       }; */
/*  inline string output_flow() const  */
/* { */
/*   ostringstream out; */
/*   out << title(string("photons results")); */
/*   out << "e_phot.1 = " << energy_[0]/max(double(1.0),double(number_[0])) << " e_phot.2 = " << energy_[1]/max(double(1.0),double(number_[1])) << endl; */
/*   return out.str(); */
/* } */
/* }; */

class COHERENT_RESULTS : public ABSTRACT_IO_CLASS 
{
  double sum_,sum2_,sumeng_,upsmax_,probmax_,sumreal_,engreal_;
  double total_energy_;
  long ncall_,count_;
  int do_coherent_;

  public :
    
    COHERENT_RESULTS() 
    {
      sum_          = 0.0;
      sum2_         = 0.0;
      sumeng_       = 0.0;
      ncall_        = 0;
      upsmax_       = 0.0;
      probmax_      = 0.0;
      sumreal_      = 0.0;
      engreal_      = 0.0;
      count_        = 0;
      total_energy_ = 0.0;
      do_coherent_ = 0;
    }

  inline void init() { do_coherent_ = 1;}
  
  inline void updateProbmax(PHI_FLOAT tmp)
    {
      probmax_ = max(probmax_,tmp);
    }
  inline void updateSumeng(float ener)
    {
      sumeng_ += ener;
    }
  inline void addCountEnergy(float ener)
    {
      count_++;
      total_energy_ += fabs(ener);
    }
  
  inline void updateSums(float wt, float ener, float ups)
    {
      sum_ += wt;
      sum2_ += wt*wt;
      sumeng_ += wt*ener;
      upsmax_ = max(upsmax_, double(ups));
    }
  
  virtual inline string name_of_class() const 
    {
      return string("COHERENT_RESULTS");
    }

  virtual string output_flow() const ;
  
};

class TRIDENT_RESULTS : public ABSTRACT_IO_CLASS 
{
  //  double sum_,sum2_,sumeng_,upsmax_,probmax_,sumreal_,engreal_;
  double total_energy_;
  //  long ncall_,
  long count_;
  int do_trident_;

  public :
    
    TRIDENT_RESULTS() 
      {;
      //      sum_          = 0.0;
      //      sum2_         = 0.0;
      //      sumeng_       = 0.0;
      //      ncall_        = 0;
      //      upsmax_       = 0.0;
      //      probmax_      = 0.0;
      //      sumreal_      = 0.0;
      //      engreal_      = 0.0;
      count_        = 0;
      total_energy_ = 0.0;
      //      do_coherent_ = 0;
    }

  inline void init() { do_trident_ = 1;}

  inline void addCountEnergy(float ener)
    {
      count_++;
      total_energy_ += fabs(ener);
    }
  
  virtual inline string name_of_class() const 
    {
      return string("TRIDENT_RESULTS");
    }

  virtual string output_flow() const ;
  
};
#endif
