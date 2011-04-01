#include "tridentCPP.h"

#define CVELOCITY 3e+17 /* nm/s */
#define COMPTON HBAR*CVELOCITY/EMASS

/* Contact Jakob Esberg, Barbara Dalena or Daniel Schulte */

TRIDENT::TRIDENT()
{
  i_equiv=2; //10: virtual photon spectrum up to m*m with spin correction 2: same without spin correction
  s4=0;
  lns4=0;
  one_m_x=0;
  xmin=0.00001;
  r_phot=PHYSTOOLS::requiv(0,xmin,i_equiv); 
  virtExist=false;
  flag=false;
}

bool TRIDENT::makeVirtualPhoton(float* Emother,float* e_phot,float* q2,RNDM& hasard)
{
  int n_phot;
  float initialEnergy =fabs(*Emother);
  n_phot=(int)floor(r_phot);
  r_phot -= n_phot;
  if(hasard.rndm()<r_phot) 
    {
      PHYSTOOLS::mequiv(s4, lns4,xmin,initialEnergy,i_equiv,e_phot,q2,&one_m_x, hasard);    
      if(*e_phot>0.0) return true;
      else return false;
    }
  else return false;
} //jakob


/* Decides whether to create pairs from the virtual photons */
/* Modifies the energy of the particle Emother */
void TRIDENT::convertVirtualPhotons(float* Emother,vector<float> energies, vector<float>* tridents , float ups, double dz,RNDM& hasard_)
{
  double p;
  float kap;
  if(energies.size()>1) cout << "TRIDENT::createTridents. Error, tridents not able to cope with multiple photons" << endl;
  for(unsigned int i=0;i<energies.size();i++)
    {
      kap=kappa(ups,energies[i],fabs(*Emother));    
      p=PHYSTOOLS::u(kap)*ALPHA_EM*EMASS*dz/(energies[i]*COMPTON);
      if(p>1)
	{
	  cout << "WARNING, TRIDENT::convertVirtualPhotons. Probability exceeding 1" << endl;
	}
      //      else{
      if(p<0.01)
	{
	  if(p>hasard_.rndm())
	    {
	      tridents->push_back(energies[i]);
	      if (*Emother < 0.0)
		{
		  *Emother+=energies[i];
		}
	      else
		{
		  *Emother-=energies[i];
		}  
	    }
	}
      else if(1-exp(-p)>hasard_.rndm())
	{
	  tridents->push_back(energies[i]);
	  if (*Emother < 0.0)
	    {
	      *Emother+=energies[i];
	    }
	  else
	    {
	      *Emother-=energies[i];
	    }  
	  //	}
	}
    }
} //jakob

// void TRIDENT::convertVirtualPhotons(float* Emother,vector<float> energies, vector<float>* tridents , float ups, double dz,RNDM& hasard_)
// {
//   double p;
//   for(unsigned int i=0;i<energies.size();i++)
//     {
//       float kap=kappa(ups,energies[i],fabs(*Emother));    
//       p=PHYSTOOLS::u(kap)*ALPHA_EM*EMASS*dz/(energies[i]*COMPTON);
//       while(p>1)
// 	{
// 	  tridents->push_back(energies[i]);
// 	  if (*Emother < 0.0)
// 	    {
// 	      *Emother+=energies[i];
// 	    }
// 	  else
// 	    {
// 	      *Emother-=-energies[i];
// 	    }	  
// 	  p-=1;
// 	}
//       if(p>hasard_.rndm())
// 	{
// 	  tridents->push_back(energies[i]);
// 	  if (*Emother < 0.0)
// 	    {
// 	      *Emother+=energies[i];
// 	    }
// 	  else
// 	    {
// 	      *Emother-=energies[i];
// 	    }  
// 	}
//     }
// } //jakob

bool TRIDENT::pick_trident_energy(float kappa,float& energy, RNDM& hasard)
{
  double a=0.13513,eta,dxdy,x,h,hp,tmp,tmp2,y,coh;
  y=2.0*hasard.rndm()-1.0;
  eta=8.0/(3.0*kappa);
  tmp=fabs(y/(1.0+(1.0-y*y)/(2.0*sqrt(eta))));
  tmp2=1.0-tmp*tmp;
  h=-log(tmp2);
  hp=tmp/tmp2*(1.0+(1.0+y*y)/(2.0*sqrt(eta)))
    /(1.0+(1.0-y*y)/(2.0*sqrt(eta)))
    /(1.0+(1.0-y*y)/(2.0*sqrt(eta)));
  if (y>0.0){
    x=0.5*(1.0+sqrt(h/(eta+h)));
  }
  else{
    x=0.5*(1.0-sqrt(h/(eta+h)));
  }
  dxdy=hp/sqrt(h*(h+eta))*(1.0-h/(h+eta));
  coh=PHYSTOOLS::u(kappa);
  tmp=PHYSTOOLS::fcp(kappa,x)*a/coh*dxdy;
  if (hasard.rndm()<tmp){
     energy *= x;
    return true;
  } else {
    return false;
  }
}

void TRIDENT::createTridents(float* Emother,float ups,double dz,vector<float>* electrons, vector<float>* positrons,vector<float>* virt, RNDM& hasard) //dz is in nm
{  
  vector<float> energies;
  float kap;
  e_phot=0;
  q2=0;
  virtExist=makeVirtualPhoton(Emother,&e_phot,&q2,hasard);   
  if(virtExist)
    {
      energies.push_back(e_phot); 
    }
  convertVirtualPhotons(Emother,energies,electrons,ups,dz,hasard);
  if (energies.size()>1) cout << "TRIDENT::createTridents. Error, tridents not able to cope with multiple photons"<<endl;
  if (electrons->size()>1) cout << "TRIDENT::createTridents. Error, tridents not able to cope with multiple electrons"<<endl;
  for(unsigned int i=0;i<electrons->size();i++)
    {
      help=(*electrons)[i];
      flag=false;
      kap=kappa(ups,help,fabs(*Emother));
      while(!flag)
	{
	  flag=pick_trident_energy(kap,(*electrons)[i],hasard);
	}
      positrons->push_back((*electrons)[i]-help);
      virt->push_back(q2);
    }
  //  bottom: ; //to be removed
}

void TRIDENT::createTridents(float* Emother,float ups,double dz,vector<float>* electrons, vector<float>* positrons, RNDM& hasard) //dz is in nm
{  
  vector<float> energies;
  float kap;
  e_phot=0;
  q2=0;
  virtExist=makeVirtualPhoton(Emother,&e_phot,&q2,hasard);   
  if(virtExist)
    {
      energies.push_back(e_phot); 
    }
  convertVirtualPhotons(Emother,energies,electrons,ups,dz,hasard);
  if (electrons->size()>1) cout << "TRIDENT::createTridents. Error, tridents not able to cope with multiple electrons"<<endl;
  for(unsigned int i=0;i<electrons->size();i++)
    {
      help=(*electrons)[i];
      flag=false;
      kap=kappa(ups,help,fabs(*Emother));
      while(!flag)
	{
	  flag=pick_trident_energy(kap,(*electrons)[i],hasard);
	}
      positrons->push_back((*electrons)[i]-help);
      //      virt->push_back(q2);
    }
}
