#include "bhabhaSamplesCPP.h"

bool BHABHA_PHOTON_SAMPLES::pick_next(float ecmratio, float& en,float& px,float& py,float& pz, int& eventIndex)
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
  eventIndex = numero_bhabha_[next_];
  next_++;
  return true;
}

bool BHABHASAMPLES::pick_next_bhabha(float e1, float e2, float ecmratio, float eCM, float& px1,float& py1,float& pz1, float& en1,float& px2,float& py2,float& pz2,float& en2, int& nbphot, unsigned int& numero_bhabha)
{
  if (next_ >= bhabha_.size() ) 
    {
      cerr << " WARNING BHABHASAMPLES::pick_next_bhabha: the vector of bhabha samples is out  bhabha_.size()= " <<  bhabha_.size() << "next= " << next_ << endl;
      return false;
    }
  if (prod_info_ == 50000)
    {
      cout << " BHABHASAMPLES::pick_next_bhabha: bhabhas produced " <<  next_ << endl;
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

  bhabha_[next_].mother1 = fabs(e1);
  bhabha_[next_].mother2 = fabs(e2);
  bhabha_[next_].eCM = eCM;
  numero_bhabha = bhabha_[next_].evtIndex;
  next_++;
  prod_info_++;
  return true;
}

void BHABHA::boost_bhabha(float part1Vx, float part1Vy, float part2Vx, float part2Vy,float e1, float e2, float& px1in,float& py1in,float& pz1in,float& e1in,float& px2in,float& py2in,float& pz2in,float& e2in, int nphot, float ecmratio,  int do_bhabha, int bhabha_event_index)
{
// Modified by Strahinja Lukic, Nov 2011
  int k;
  float px1,py1,pz1;
  float px2,py2,pz2;
  double beta_x,beta_y,beta_z; // Velocity of the CM frame
  float prelx, prely, prelz, prelrho, prel;
  float theta,phi;
  float px,py,pz,en;
  int photon_event_index;

/*
  frame_change_part_of_bhabha(part1Vx, part1Vy, e1, px1, py1,theta1, phi1);
  pz1 = sqrt(e1*e1-px1*px1-py1*py1);

  frame_change_part_of_bhabha(part2Vx, part2Vy, e2, px2, py2, theta2, phi2);
  pz2 = -sqrt(e2*e2-px2*px2-py2*py2);
*/

  px1 = e1*part1Vx;
  py1 = e1*part1Vy;
  pz1 = e1*sqrt(double(1. - part1Vx*part1Vx - part1Vy*part1Vy));

  px2 = e2*part2Vx;
  py2 = e2*part2Vy;
  pz2 = -e2*sqrt(double(1. - part2Vx*part2Vx - part2Vy*part2Vy));
  if(pz2*pz2in<0. && pz1*pz1in<0.) cout << "Bhabha nr. " << bhabha_event_index << " is turned backwards!\n";
// Velocity of the CM frame
  beta_x=(px1+px2)/(e1+e2);
  beta_y=(py1+py2)/(e1+e2);
  beta_z=(pz1+pz2)/(e1+e2);

// Relative movement of the particles before the collision (in CM)
  if(fourboost(e1, px1, py1, pz1, beta_x, beta_y, beta_z))
  {
	  cout << "Bhabha nr. " << bhabha_event_index << ": The initial electron has anomalous invariant mass!\n";
	  cout << "Beta = (" << beta_x << ", " << beta_y << ", " << beta_z << ")\n";
  }
  if(fourboost(e2, px2, py2, pz2, beta_x, beta_y, beta_z))
  {
	  cout << "Bhabha nr. " << bhabha_event_index << ": The initial positron has anomalous invariant mass!\n";;
	  cout << "Beta = (" << beta_x << ", " << beta_y << ", " << beta_z << ")\n";
  }
  prelx = px1-px2;
  prely = py1-py2;
  prelz = pz1-pz2;
  prelrho = sqrt(prelx*prelx + prely*prely);
  prel = sqrt(prelrho*prelrho + prelz*prelz);

// Angles of the collision axis in the CM frame
  theta=asin(prelrho/prel);
  if(abs(theta)<0.00001){
    phi=0.;
  }
  else{
    float cosphi=prelx/prelrho;
    float sinphi=prely/prelrho;
    if(sinphi>0.) {phi=acos(cosphi);}
    else {phi=2*PI-acos(cosphi);}
  }

// Orient the bhabhas along the collision axis
  bhabha_rotation(theta,phi,px1in,py1in,pz1in);
  bhabha_rotation(theta,phi,px2in,py2in,pz2in);

// Go back to the lab frame
  if(fourboost(e1in, px1in, py1in, pz1in, -beta_x, -beta_y, -beta_z))
  {
	  cout << "Bhabha nr. " << bhabha_event_index << ": The final electron has anomalous invariant mass!\n";
	  cout << "Beta = (" << beta_x << ", " << beta_y << ", " << beta_z << ")\n";
  }
  if(fourboost(e2in, px2in, py2in, pz2in, -beta_x, -beta_y, -beta_z))
  {
	  cout << "Bhabha nr. " << bhabha_event_index << ": The final positron has anomalous invariant mass!\n";
	  cout << "Beta = (" << beta_x << ", " << beta_y << ", " << beta_z << ")\n";
  }
//  lorent_bhabha_transformation(e1, e2, pz1, pz2, beta_x, beta_y,theta,phi, px1in, py1in, pz1in, e1in);
//  lorent_bhabha_transformation(e1, e2, pz1, pz2, beta_x, beta_y,theta,phi, px2in, py2in, pz2in, e2in);

  if ( pz1in < 0 && pz2in > 0 ) e1in = -e1in;
  else e2in = -e2in;


  if( do_bhabha>1)
    {
      for (k = 0; k <nphot;k ++)
		{
		  if ( bhabhaPhotonReserve_.pick_next(ecmratio, en,px,py,pz, photon_event_index) )
			{
			  if(photon_event_index != bhabha_event_index)
			  {
				  cout << "Error in photon storage.\n";
				  cout << "Bhabha event index: " << bhabha_event_index << endl;
				  cout << "Photon event index: " << photon_event_index << endl;
				  if(abs(bhabha_event_index-photon_event_index)>10) exit(1);
			  }
			  bhabha_rotation(theta,phi,px,py,pz);
			  if(fourboost(en, px, py, pz, -beta_x, -beta_y, -beta_z))
			  {
				  cout << "Bhabha nr. " << bhabha_event_index << ": The photon nr. " << k << " has anomalous invariant mass!";
				  cout << "Beta = (" << beta_x << ", " << beta_y << ", " << beta_z << ")\n";
			  }
			  boostedBhabhaPhotons_.add_bhabha_photon(bhabha_event_index, px, py, pz, en);
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
