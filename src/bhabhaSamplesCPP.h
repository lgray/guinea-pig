#ifndef BHABHASAMPLES_SEEN
#define BHABHASAMPLES_SEEN


#include "particlesCPP.h"
#include "pairsCPP.h"
#include "mathematicalEntities.h"

#include "abstractParticle.h"
#include "fileInputOutput.h"

#include <iostream>
#include <vector>

using namespace std;

class BHABHA_PHOTON_SAMPLES : public ABSTRACT_BHABHA_PHOTON_SAMPLES
{

  vector<QUADRIVECTOR> bhabha_photons_;
  vector<int> numero_bhabha_;
  int label_;
  int next_;
  public :
    
    BHABHA_PHOTON_SAMPLES() : label_(-1), next_(0) {;}
  
  
  virtual ~BHABHA_PHOTON_SAMPLES() {;}
  
  virtual inline int get_label() const  { return label_;}
  
  
  virtual inline int nb_samples() const 
    {
      if ( next_ == 0 ) return bhabha_photons_.size();
      else return next_;
    }
  
  
  virtual inline void get_parameters_for_output(int numero, int& numero_bhabha, float& en,float& vx,float& vy, float& vz) const
    {
      bhabha_photons_[numero].trivector(vx, vy, vz);
      en = bhabha_photons_[numero].composante4();
      numero_bhabha = numero_bhabha_[numero];
    }
  
  virtual inline void add_bhabha_photon(float px, float py, float pz, float en) 
    {
      QUADRIVECTOR bhab_phot = QUADRIVECTOR(px,py,pz,en);
      numero_bhabha_.push_back(-1);
      bhabha_photons_.push_back(bhab_phot);
    }
  
  
  inline void set_label(int label)  {label_ = label;}
  
  inline void create_bhabha_photon(int numero_bhabha, float px, float py, float pz, float en) 
    {
      QUADRIVECTOR bhab_phot = QUADRIVECTOR(px,py,pz,en);
      bhabha_photons_.push_back(bhab_phot);
      numero_bhabha_.push_back(numero_bhabha);
    }
  
  inline void set_numero_bhabha(int index, int num)
    {
      numero_bhabha_[index] = num;
    }
  
  inline  int load(string bhabhaPhotonFIleIni)
    {
      FILE_IN_OUT filin;
      filin.open_file(bhabhaPhotonFIleIni,"r");
      filin.read_bhabhaPhotonsamples(this);
      filin.close();
      return bhabha_photons_.size();
    }
  
  bool pick_next(float ecmratio, float& en,float& px,float& py,float& pz, int& found)  ; 
  
  inline  void set(int index, float en,float px,float py,float pz, int numero_bhabha)
    {
      bhabha_photons_[index].set(px, py, pz, en);
      numero_bhabha_[index] = numero_bhabha;
    }
  
  inline void save_on_file(string nameOfOutputFile) const 
    {
      FILE_IN_OUT filout;
      filout.open_file(nameOfOutputFile, "w");
      filout.save_bhabhaPhotonSamples(this);
      filout.close();
    }
};

class BHABHASAMPLES : public ABSTRACT_BHABHASAMPLES
{
  
  typedef struct
  {
    
    QUADRIVECTOR p1, p2;
    float mother1,mother2;
    int nbphot;
  } BHABHA_INI;
  
  vector<BHABHA_INI> bhabha_;
  unsigned long next_;
  unsigned long prod_info_;
  
  public :
    
    BHABHASAMPLES() : next_(0), prod_info_(0) {;}
  
  virtual  ~BHABHASAMPLES() {;}
  
  
  virtual inline void get_parameters_for_output(int numero, unsigned long& rank1_index, float& mother1_en,float&e1,float&vx1,float& vy1, float&vz1, unsigned long& rank2_index, float& mother2_en, float& e2, float& vx2, float&vy2, float&vz2, int& nbphot) const
    {
      float px, py, pz;
      
      if (numero >= next_)
	{
	  cerr << " WARNING : BHABHASAMPLES::get_parameters_for_output: numero of bhabha_prod to save is out of range numero = " << numero << " next_= " << next_  << endl;
	  return;
	}
      rank1_index = 2*numero+1;
      rank2_index = rank1_index+1;
      mother1_en = bhabha_[numero].mother1;
      mother2_en = bhabha_[numero].mother2;
      
      e1 = bhabha_[numero].p1.composante4();
      
      bhabha_[numero].p1.trivector(px, py, pz);
      vx1=px/abs(e1);
      vy1=py/abs(e1);
      vz1=pz/abs(e1);
      
      e2=bhabha_[numero].p2.composante4();
      bhabha_[numero].p2.trivector(px, py, pz);
      vx2=px/abs(e2);
      vy2=py/abs(e2);
      vz2=pz/abs(e2); 
      nbphot = bhabha_[numero].nbphot;
 }
 
  virtual inline int nb_samples() const {return next_;}
  

  bool pick_next_bhabha(float e1, float e2, float ecmratio, float& px1,float& py1,float& pz1, float& en1,float& px2,float& py2,float& pz2,float& en2, int& nbphot, int& numero_bhabha); 
  
  
  inline void save_on_file(string nameOfOutputFile) const 
    {
      FILE_IN_OUT filout;
      filout.open_file(nameOfOutputFile, "w");
      filout.save_bhabhasamples(this);
      filout.close();
    }
  
  virtual inline void add_bhabha(float px1, float py1, float pz1, float e1, float px2, float py2, float pz2, float e2, int nbphot)
    {
      BHABHA_INI bhab;
      bhab.p1 = QUADRIVECTOR(px1,py1,pz1,e1);
      bhab.p2 = QUADRIVECTOR(px2,py2,pz2,e2);
      bhab.nbphot = nbphot;
      bhabha_.push_back(bhab);
    }
  
  inline  int load_bhabha(string bhabhaFIleIni)
    {
      FILE_IN_OUT filin;
      filin.open_file(bhabhaFIleIni,"r");
      filin.read_bhabhasamples(this);
      filin.close();
      return bhabha_.size();
    }
  
};

class BHABHA
{
  BHABHASAMPLES bhabhaReserve_; 
  BHABHA_PHOTON_SAMPLES bhabhaPhotonReserve_;
  BHABHA_PHOTON_SAMPLES boostedBhabhaPhotons_;
  //  PAIR_BEAM bhabhas_;
  int nbhabha_ini_;
  int nbhabha_photon_ini_;
  
  
  /********************************************************/
  /*!boost bhabha from CM frame of e+e- (P1P2) to lab frame*/
  /********************************************************/
  void bhabha_rotation(float theta, float phi, float& px1, float& py1, float& pz1)
    {
      float px,py,pz;
      px=px1;
      py=py1;
      pz=pz1;
      px1 = cos(theta)*cos(phi)*px-sin(phi)*py+sin(theta)*cos(phi)*pz;
      py1 = cos(theta)*sin(phi)*px+cos(phi)*py+sin(theta)*sin(phi)*pz;
      pz1 = -sin(theta)*px+cos(theta)*pz;
    }
  
  void lorent_bhabha(float e1,float e2,float pz1,float pz2,float& e,float& pz)
    {
      float beta, eold, gam,pzold;
      
      beta = -(pz1 + pz2) / (e1 + e2);
      gam = 1.0 / sqrt(1.0 - beta * beta);
      eold = e;
      pzold= pz;
      e = gam * (e - beta * pz);
      pz = gam * (pz - beta * eold);
    }
  
  void lorent_bhabha_back(float& e,float& pl,float beta)
    {
      /* +beta because want to transform back */
      float gamma,eold;
      gamma=1.0/sqrt(1.0-beta*beta);
      eold=e;
      e=gamma*(eold + beta * pl);
      pl=gamma*(pl + beta * eold);
    }
  
  void frame_change_part_of_bhabha(float partVx, float partVy, float e, float& px, float& py,float& theta, float& phi)
    {
      float cosphi,sinphi;
      px=e*partVx;
      py=e*partVy;
      
      theta=asin(sqrt(partVx*partVx+partVy*partVy));
      if(abs(theta)<0.00001)
	{
	  phi=0.;
	}
      else
	{
	  cosphi=partVx/sin(theta);
	  sinphi=partVy/sin(theta);
	  if(sinphi>0.) phi=acos(cosphi);
	  if(sinphi<0.) phi=2*PI-acos(cosphi);
	}
    }
  
  void lorent_bhabha_transformation(float e1, float e2, float pz1, float pz2,float beta_x, float beta_y, float theta, float phi, float& pxin,float& pyin,float& pzin,float& ein)
    {
      
      lorent_bhabha(e1,e2,pz1,pz2, ein,pzin);
      lorent_bhabha_back( ein,pxin,beta_x);
      lorent_bhabha_back( ein,pyin,beta_y);
      bhabha_rotation(theta,phi,pxin,pyin,pzin);
    }
  void boost_bhabha(float part1Vx, float part1Vy, float part2Vx, float part2Vy,float e1, float e2, float& px1in,float& py1in,float& pz1in,float& e1in,float& px2in,float& py2in,float& pz2in,float& e2in, int nphot, float ecmratio,  int do_bhabha, int numero_bhabha);
  
 public:
  
  BHABHA() : nbhabha_ini_(0),nbhabha_photon_ini_(0) {;}
  ~BHABHA() {;}
  
  inline void load_samples(int do_bhabhas, string bhabha_samples, string bhabha_photon_samples)
    {
      nbhabha_ini_ = bhabhaReserve_.load_bhabha(bhabha_samples);
      if (do_bhabhas > 1) 
	{
	  nbhabha_photon_ini_ =  bhabhaPhotonReserve_.load(bhabha_photon_samples);
	}
    }
  
  inline void make_bhabha(PAIR_BEAM& bhabhas, float part1Vx, float part1Vy, float part2Vx, float part2Vy, float e1, float e2,float ecm, float weight, MESH& mesh, int cellx, int celly,float min_z, const SWITCHES& switches, RNDM& hasard)
    {
      double bhabhan;
      float ecmratio;
      float px1, py1, pz1, en1, px2, py2, pz2, en2;
      int nbphot, numero_bhabha;
      ecmratio = ecm/switches.get_bhabha_ecmload();
      bhabhan = switches.get_bhabha_scal()*weight/(ecmratio*ecmratio);
      if (hasard.rndm()< bhabhan)
	{
	  if (	bhabhaReserve_.pick_next_bhabha(e1, e2, ecmratio, px1, py1, pz1, en1, px2, py2, pz2, en2, nbphot, numero_bhabha) )
	    {
	      boost_bhabha(part1Vx, part1Vy, part2Vx, part2Vy, e1, e2, px1, py1, pz1, en1, px2, py2, pz2, en2, nbphot, ecmratio, switches.get_do_bhabhas(), numero_bhabha);
	      //bhabhas.new_pair(mesh, cellx, celly,min_z, -1, en1, px1, py1,pz1, switches.get_bhabha_ratio(), switches.get_track_secondaries(), switches.get_store_secondaries(), hasard );  
	      //bhabhas.new_pair(mesh, cellx, celly, min_z, -1, en2, px2, py2,pz2, switches.get_bhabha_ratio(), switches.get_track_secondaries(), switches.get_store_secondaries(), hasard );
	      bhabhas.new_pair(mesh, cellx, celly, min_z, -1, en1, px1, py1,pz1, switches.get_bhabha_ratio(), switches.get_track_pairs(), switches.get_store_pairs(), hasard );  
	      bhabhas.new_pair(mesh, cellx, celly, min_z, -1, en2, px2, py2,pz2, switches.get_bhabha_ratio(), switches.get_track_pairs(), switches.get_store_pairs(), hasard );
	    }
	}
    }
  

  inline void numbers_of_loaded(int& nbhabha_ini, int& nbhabha_photon_ini) const
    {
      nbhabha_ini = nbhabha_ini_;
      nbhabha_photon_ini = nbhabha_photon_ini_;
    }
  inline void save_on_files(int do_bhabhas, string bhabha_prod, string bhphoton_prod, string bhphotons) const
    {
      bhabhaReserve_.save_on_file(bhabha_prod);
      //     bhabhaSamples_.save_on_file_pour_C(string("bhabhaIniPourC"));
      if (do_bhabhas == 2)
	{
	  bhabhaPhotonReserve_.save_on_file(bhphoton_prod);
	  boostedBhabhaPhotons_.save_on_file(bhphotons);
	}
    }
};

#endif
