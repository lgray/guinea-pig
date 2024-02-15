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
  double one_m_x,r_phot,xmin;
  double help,e_phot,q2;
  bool virtExist,flag;

 public:
  TRIDENT();
  ~TRIDENT() {;}
  
  bool makeVirtualPhoton(double* Emother,double* e_phot,double* q2, RNDM& rndm_generator_);
  void convertVirtualPhotons(double* Emother,std::vector<double> energies, std::vector<double>* tridents , double ups, double dz,RNDM& rndm_generator_);

  bool pick_trident_energy(double kappa,double& energy, RNDM& rndm_generator);

  void createTridents(double* Emother,double ups,double dz,std::vector<double>* electrons,std::vector<double>* positrons,std::vector<double>* virt,RNDM& rndm_generator);

  void createTridents(double* Emother,double ups,double dz,std::vector<double>* electrons, std::vector<double>* positrons, RNDM& rndm_generator);

  /* Field parameter as seen by the photon */
  inline double kappa(double ups,double e_phot,double Emother)
    {
      return ups*e_phot/Emother;
    }

};

