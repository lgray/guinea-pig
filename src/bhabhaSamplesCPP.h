#ifndef BHABHASAMPLES_SEEN
#define BHABHASAMPLES_SEEN

#include "particlesCPP.h"
#include "pairsCPP.h"
#include "mathematicalEntities.h"

#include "abstractParticle.h"
#include "fileInputOutput.h"

#include <iostream>
#include <string>
#include <vector>

class BHABHA_PHOTON_SAMPLES : public ABSTRACT_BHABHA_PHOTON_SAMPLES
{

  std::vector<QUADRIVECTOR> bhabha_photons_;
  std::vector<int> number_bhabha_;
  int label_;
  unsigned int next_;
  public :
    
    BHABHA_PHOTON_SAMPLES() : label_(-1), next_(0) {;}
  
  
  virtual ~BHABHA_PHOTON_SAMPLES() {;}
  
  virtual inline int get_label() const  { return label_;}
  
  
  virtual inline unsigned int nb_samples() const 
    {
      if ( next_ == 0 ) return bhabha_photons_.size();
      else return next_;
    }
  
  
  virtual inline void get_parameters_for_output(unsigned int number, int& number_bhabha, double& en,double& vx,double& vy, double& vz) const
    {
      bhabha_photons_[number].trivector(vx, vy, vz);
      en = bhabha_photons_[number].energy();
      number_bhabha = number_bhabha_[number];
    }
  
  virtual inline void add_bhabha_photon(int nbhabha, double px, double py, double pz, double en)
    {
      QUADRIVECTOR bhab_phot = QUADRIVECTOR(px,py,pz,en);
      number_bhabha_.push_back(nbhabha);
      bhabha_photons_.push_back(bhab_phot);
    }
  
  
  inline void set_label(int label)  {label_ = label;}
  
/*  inline void create_bhabha_photon(int number_bhabha, double px, double py, double pz, double en)
    {
      QUADRIVECTOR bhab_phot = QUADRIVECTOR(px,py,pz,en);
      bhabha_photons_.push_back(bhab_phot);
      number_bhabha_.push_back(number_bhabha);
    }
*/
/* This function shouldn't exist - no changing the number of bhabha after loading
    inline void set_number_bhabha(int index, int num)
    {
      number_bhabha_[index] = num;
    }
*/
  inline  int load(std::string bhabhaPhotonFIleIni)
    {
      FILE_IN_OUT filin;
      filin.open_file(bhabhaPhotonFIleIni,"r");
      filin.read_bhabhaPhotonsamples(this);
      filin.close();
      return bhabha_photons_.size();
    }
  
  bool pick_next(double ecmratio, double& en,double& px,double& py,double& pz, int& found)  ; 
  
/*  inline  void set(int index, double en,double px,double py,double pz, int number_bhabha)
    {
      bhabha_photons_[index].set(px, py, pz, en);
      number_bhabha_[index] = number_bhabha;
    }
*/
  inline void save_on_file(std::string nameOfOutputFile) const 
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
    double mother1,mother2,eCM;
    unsigned int evtIndex;
    int nbphot;
  } BHABHA_INI;
  
  std::vector<BHABHA_INI> bhabha_;
  unsigned long next_;
  unsigned long prod_info_;
  
  public :
    
    BHABHASAMPLES() : next_(0), prod_info_(0) {;}
  
  virtual  ~BHABHASAMPLES() {;}
  
  
  virtual inline void get_parameters_for_output(unsigned int number, unsigned int& evtIdx, double& eCM, double& mother1_en,double&e1,double&vx1,double& vy1, double&vz1, double& mother2_en, double& e2, double& vx2, double&vy2, double&vz2, int& nbphot) const
    {
      double px, py, pz;
      
      if (number >= next_)
		{
		  std::cerr << " WARNING : BHABHASAMPLES::get_parameters_for_output: number of bhabha_prod to save is out of range number = " << number << " next_= " << next_  << std::endl;
		  return;
		}
      evtIdx = bhabha_[number].evtIndex;
      eCM = bhabha_[number].eCM;
      mother1_en = bhabha_[number].mother1;
      mother2_en = bhabha_[number].mother2;
      
      e1 = bhabha_[number].p1.energy();
      
      bhabha_[number].p1.trivector(px, py, pz);
      vx1=px/std::abs(e1);
      vy1=py/std::abs(e1);
      vz1=pz/std::abs(e1);
      
      e2=bhabha_[number].p2.energy();
      bhabha_[number].p2.trivector(px, py, pz);
      vx2=px/std::abs(e2);
      vy2=py/std::abs(e2);
      vz2=pz/std::abs(e2); 
      nbphot = bhabha_[number].nbphot;
    }
 
  virtual inline unsigned int nb_samples() const {return next_;}
  

  bool pick_next_bhabha(double e1, double e2, double ecmratio, double eCM, double& px1,double& py1,double& pz1, double& en1,double& px2,double& py2,double& pz2,double& en2, int& nbphot, unsigned int& number_bhabha);
  
  
  inline void save_on_file(std::string nameOfOutputFile) const 
    {
      FILE_IN_OUT filout;
      filout.open_file(nameOfOutputFile, "w");
      filout.save_bhabhasamples(this);
      filout.close();
    }
  
  virtual inline void add_bhabha(unsigned int evtIdx, double px1, double py1, double pz1, double e1, double px2, double py2, double pz2, double e2, int nbphot)
    {
      BHABHA_INI bhab;
      bhab.p1 = QUADRIVECTOR(px1,py1,pz1,e1);
      bhab.p2 = QUADRIVECTOR(px2,py2,pz2,e2);
      bhab.nbphot = nbphot;
      bhab.evtIndex = evtIdx;
      bhabha_.push_back(bhab);
    }
  
  inline  int load_bhabha(std::string bhabhaFIleIni)
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
  void bhabha_rotation(double theta, double phi, double& px, double& py, double& pz)
    {
      double pxin=px;
      double pyin=py;
      double pzin=pz;
      double costh = cos(theta);
      double sinth = sin(theta);
      double cosphi = cos(phi);
      double sinphi = sin(phi);
      px = costh*cosphi*pxin-sinphi*pyin+sinth*cosphi*pzin;
      py = costh*sinphi*pxin+cosphi*pyin+sinth*sinphi*pzin;
      pz = -sinth*pxin+costh*pzin;
    }
  
  int fourboost(double &en, double &px, double &py, double &pz, double beta_x, double beta_y, double beta_z)
  {
		double ein = en;
		double pxin = px;
		double pyin = py;
		double pzin = pz;
		double check_invar_mass = ein*ein - pxin*pxin - pyin*pyin - pzin*pzin;

		double betasq = beta_x*beta_x + beta_y*beta_y + beta_z*beta_z;

		if(1. - betasq < 1.e-20)
		{
			printf("Boost too far in fourboost. beta = (%f, %f, %f)\n", beta_x, beta_y, beta_z);
			px = py = pz = en = 0;
			return -1;
		}

		double gamma = 1./sqrt(1.- betasq);
		double gm1 = gamma - 1.;
		double Bxx = gm1*beta_x*beta_x/betasq;
		double Bxy = gm1*beta_x*beta_y/betasq;
		double Bxz = gm1*beta_x*beta_z/betasq;
		double Byy = gm1*beta_y*beta_y/betasq;
		double Byz = gm1*beta_y*beta_z/betasq;
		double Bzz = gm1*beta_z*beta_z/betasq;

		double enout = gamma*(ein - beta_x*pxin - beta_y*pyin - beta_z*pzin);
		double pxout = -gamma*beta_x*ein + (Bxx+1.)*pxin + Bxy*pyin + Bxz*pzin;
		double pyout = -gamma*beta_y*ein + Bxy*pxin + (Byy+1.)*pyin + Byz*pzin;
		double pzout = -gamma*beta_z*ein + Bxz*pxin + Byz*pyin + (Bzz+1.)*pzin;
		en = enout; px = pxout; py = pyout; pz = pzout;

		double end_invar_mass = enout*enout - pxout*pxout - pyout*pyout - pzout*pzout;

		if(std::abs(check_invar_mass) > 1. || std::abs(end_invar_mass) > 1. || std::abs(check_invar_mass - end_invar_mass)>1.)
		{
			printf("Start invariant mass: %f\n", check_invar_mass);
			printf("End invariant mass: %f\n", end_invar_mass);
			return -2;
		}

		return 0;
  }

/*
  void lorent_bhabha(double e1,double e2,double pz1,double pz2,double& e,double& pz)
    {
      double beta, eold, gam;
      
      beta = -(pz1 + pz2) / (e1 + e2);
      gam = 1.0 / sqrt(1.0 - beta * beta);
      eold = e;
      e = gam * (e - beta * pz);
      pz = gam * (pz - beta * eold);
    }
  
  void lorent_bhabha_back(double& e,double& pl,double beta)
    {
      double gamma,eold;
      gamma=1.0/sqrt(1.0-beta*beta);
      eold=e;
      e=gamma*(eold + beta * pl);
      pl=gamma*(pl + beta * eold);
    }
  
  void frame_change_part_of_bhabha(double partVx, double partVy, double e, double& px, double& py,double& theta, double& phi)
    {
      double cosphi,sinphi;
      px=e*partVx;
      py=e*partVy;
      
      theta=asin(sqrt(partVx*partVx+partVy*partVy));
      if(std::abs(theta)<0.00001)
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
  
  void lorent_bhabha_transformation(double e1, double e2, double pz1, double pz2,double beta_x, double beta_y, double theta, double phi, double& pxin,double& pyin,double& pzin,double& ein)
    {
      
      lorent_bhabha(e1,e2,pz1,pz2, ein,pzin);
      lorent_bhabha_back( ein,pxin,beta_x);
      lorent_bhabha_back( ein,pyin,beta_y);
      bhabha_rotation(theta,phi,pxin,pyin,pzin);
    }
*/
  void boost_bhabha(double part1Vx, double part1Vy, double part2Vx, double part2Vy,double e1, double e2, double& px1in,double& py1in,double& pz1in,double& e1in,double& px2in,double& py2in,double& pz2in,double& e2in, int nphot, double ecmratio,  int do_bhabha, int number_bhabha);
  
 public:
  
  BHABHA() : nbhabha_ini_(0),nbhabha_photon_ini_(0) {;}
  ~BHABHA() {;}
  
  inline void load_samples(int do_bhabhas, std::string bhabha_samples, std::string bhabha_photon_samples)
    {
	  std::cout << "Loading Bhabha samples\n";
      nbhabha_ini_ = bhabhaReserve_.load_bhabha(bhabha_samples);
      if (do_bhabhas > 1) 
		{
    	  std::cout << "Loading Bhabha photons\n";
		  nbhabha_photon_ini_ =  bhabhaPhotonReserve_.load(bhabha_photon_samples);
		}
    }
  
    inline void make_bhabha(PAIR_BEAM& bhabhas, double part1Vx, double part1Vy, double part2Vx, double part2Vy, double e1, double e2,double ecm, double weight, MESH& mesh, int cellx, int celly,double min_z, const SWITCHES& switches, RNDM& rndm_generator, int beamslice1, int beamslice2)
    {
      double px1, py1, pz1, en1, px2, py2, pz2, en2;
      int nbphot;
      unsigned int evtIndex;
      double ecmratio = ecm/switches.get_bhabha_ecmload();
      double bhabhan = switches.get_bhabha_scal()*weight*std::pow(double(ecmratio*ecmratio),double(-.9891));
      if (rndm_generator.rndm()< bhabhan)
        {
          if (bhabhaReserve_.pick_next_bhabha(e1, e2, ecmratio, ecm, px1, py1, pz1, en1, px2, py2, pz2, en2, nbphot, evtIndex) )
            {
              boost_bhabha(part1Vx, part1Vy, part2Vx, part2Vy, e1, e2, px1, py1, pz1, en1, px2, py2, pz2, en2, nbphot, ecmratio, switches.get_do_bhabhas(), evtIndex);
              bhabhas.new_pair(evtIndex, mesh, cellx, celly, min_z, -1, en1, px1, py1,pz1, switches.get_bhabha_ratio(), switches.get_track_pairs(), switches.get_store_pairs(), rndm_generator, beamslice1, beamslice2 );
              bhabhas.new_pair(evtIndex, mesh, cellx, celly, min_z, -1, en2, px2, py2,pz2, switches.get_bhabha_ratio(), switches.get_track_pairs(), switches.get_store_pairs(), rndm_generator, beamslice1, beamslice2  );

            }
        }
    }
  

  inline void numbers_of_loaded(int& nbhabha_ini, int& nbhabha_photon_ini) const
    {
      nbhabha_ini = nbhabha_ini_;
      nbhabha_photon_ini = nbhabha_photon_ini_;
    }
  inline void save_on_files(int do_bhabhas, std::string bhabha_prod, std::string bhphoton_prod, std::string bhphotons) const
    {
	  std::cout << "Saving Bhabha samples in " << bhabha_prod << std::endl;
      bhabhaReserve_.save_on_file(bhabha_prod);
      //     bhabhaSamples_.save_on_file_pour_C(std::string("bhabhaIniPourC"));
      if (do_bhabhas == 2)
		{
    	  std::cout << "Saving Bhabha photon samples in " << bhphoton_prod << std::endl;
		  bhabhaPhotonReserve_.save_on_file(bhphoton_prod);
    	  std::cout << "Saving boosted Bhabha photons in " << bhphotons << std::endl;
		  boostedBhabhaPhotons_.save_on_file(bhphotons);
		}
    }
};

#endif
