//#include "particlesCPP.h"
#include "physicalTools.h"
//#include "gridCPP.h"
#include "rndmCPP.h"
#include <iostream>
#include <vector>
#include <fstream>

class TRIDENT
{
 private:
  int i_equiv; //10: virtual photon spectrum up to m*m with spin correction 2: same without spin correction
  double s4,lns4;
  float one_m_x,r_phot,xmin;
  float help,e_phot,q2;
  bool virtExist,flag;

 public:
  TRIDENT();
  ~TRIDENT() {;}
  
  bool makeVirtualPhoton(float* Emother,float* e_phot,float* q2, RNDM& hasard_);
  void convertVirtualPhotons(float* Emother,vector<float> energies, vector<float>* tridents , float ups, double dz,RNDM& hasard_);

  bool pick_trident_energy(float kappa,float& energy, RNDM& hasard);

  void createTridents(float* Emother,float ups,double dz,vector<float>* electrons,vector<float>* positrons,vector<float>* virt,RNDM& hasard);

  void createTridents(float* Emother,float ups,double dz,vector<float>* electrons, vector<float>* positrons, RNDM& hasard);

  /* Field parameter as seen by the photon */
  inline float kappa(float ups,float e_phot,float Emother)
    {
      return ups*e_phot/Emother;
    }

};

