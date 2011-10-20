#include "bhabhaSamplesCPP.h"

bool BHABHA_PHOTON_SAMPLES::pick_next(float ecmratio, float& en,float& px,float& py,float& pz, int& found) 
{
  if (next_ >= bhabha_photons_.size() ) 
    {
      cerr << " WARNING BHABHA_PHOTON_SAMPLES::pick_next() : the vector of bhabha photons samples is out  bhabha_photons_.size()= " <<  bhabha_photons_.size() << " next= " << next_ << endl;
      return false;
    }
  bhabha_photons_[next_].trivector(px, py, pz);
  en = bhabha_photons_[next_].composante4();
  px *= ecmratio;
  py *= ecmratio;
  pz *= ecmratio;
  en *= ecmratio;
  found = next_;
  next_++;
  return true;
}

bool BHABHASAMPLES::pick_next_bhabha(float e1, float e2, float ecmratio, float& px1,float& py1,float& pz1, float& en1,float& px2,float& py2,float& pz2,float& en2, int& nbphot, int& numero_bhabha) 
{
  if (next_ >= bhabha_.size() ) 
    {
      cerr << " WARNING BHABHASAMPLES::pick_next_bhabha: the vector of bhabha samples is out  bhabha_.size()= " <<  bhabha_.size() << "next= " << next_ << endl;
      return false;
    }
  if (prod_info_ == 50000)
    {
      cout << " BHABHASAMPLES::pick_next_bhabha: bhabhas produced " <<  prod_info_ << endl;
      prod_info_ = 0;
    }
  nbphot = bhabha_[next_].nbphot;
  bhabha_[next_].p1.trivector(px1, py1, pz1);
  en1 = bhabha_[next_].p1.composante4();
  bhabha_[next_].p2.trivector(px2, py2, pz2);
  en2 = bhabha_[next_].p2.composante4();
  px1 *= ecmratio;
  py1 *= ecmratio;
  pz1 *= ecmratio;
  en1 *= ecmratio;

  px2 *= ecmratio;
  py2 *= ecmratio;
  pz2 *= ecmratio;
  en2 *= ecmratio;

  bhabha_[next_].mother1 = e1;
  bhabha_[next_].mother2 = e2;
  numero_bhabha = next_;
  next_++;
  prod_info_++;
  return true;
}

void BHABHA::boost_bhabha(float part1Vx, float part1Vy, float part2Vx, float part2Vy,float e1, float e2, float& px1in,float& py1in,float& pz1in,float& e1in,float& px2in,float& py2in,float& pz2in,float& e2in, int nphot, float ecmratio,  int do_bhabha, int numero_bhabha)
{
  int k;
  float px1,py1,pz1;
  float px2,py2,pz2;
  float beta_x,beta_y;
  float theta1,theta2;
  float phi1=0.0,phi2=0.0;
  float theta,phi;
  float px,py,pz,en;
  int current_bhabha, numero_in_bhabha_file;



  frame_change_part_of_bhabha(part1Vx, part1Vy, e1, px1, py1,theta1, phi1);
  pz1 = sqrt(e1*e1-px1*px1-py1*py1);

  frame_change_part_of_bhabha(part2Vx, part2Vy, e2, px2, py2, theta2, phi2);
  pz2 = -sqrt(e2*e2-px2*px2-py2*py2);

  beta_x=(px1+px2)/(e1+e2);
  beta_y=(py1+py2)/(e1+e2);

  theta=(theta1+theta2)/2.;
  phi=0.5*(phi1+phi2);
  
  lorent_bhabha_transformation(e1, e2, pz1, pz2, beta_x, beta_y,theta,phi, px1in, py1in, pz1in, e1in);


  lorent_bhabha_transformation(e1, e2, pz1, pz2, beta_x, beta_y,theta,phi, px2in, py2in, pz2in, e2in);

  if ( pz1in < 0 && pz2in > 0 ) e1in = -e1in;
  else e2in = -e2in;


  if( do_bhabha>1)
    {
      for (k = 0; k <nphot;k ++)
	{
	  if ( bhabhaPhotonReserve_.pick_next(ecmratio, en,px,py,pz, current_bhabha) )
	    {
	      numero_in_bhabha_file = 2*numero_bhabha+1;
	      bhabhaPhotonReserve_.set_numero_bhabha(current_bhabha, numero_in_bhabha_file);
	      lorent_bhabha_transformation(e1, e2, pz1, pz2, beta_x, beta_y,theta,phi, px, py, pz, en);
	      boostedBhabhaPhotons_.create_bhabha_photon(numero_in_bhabha_file, px, py, pz, en);
	    }	
	  else
	    {
	      cerr << " GRID::boost_bhabha() : very strange, not enough photons in file, for the given number of bhabhas... " << endl;
	      exit(0);
	    }
	}
    }
}
/*end boost_bhabha*/
