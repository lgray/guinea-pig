#include "resultsCPP.h"

RESULTS::RESULTS() : switches_(NULL)
{
  int i,j;
  for (i=0;i<3;i++){
    for (j=0;j<3;j++){
      lumi[i][j]=0.0;
    }
    hadrons_ee[i]=0.0;
    hadrons_eg[i]=0.0;
    hadrons_ge[i]=0.0;
    hadrons_gg[i]=0.0;
  }
  lumi_0=0.0;
  lumi_fine=0.0;
  lumi_ee=0.0;
  lumi_ee_high=0.0;
  lumi_pp=0.0;
  lumi_eg=0.0;
  lumi_ge=0.0;
  lumi_gg=0.0;
  lumi_gg_high=0.0;
  lumi_ecm=0.0;
  lumi_ecm2=0.0;
  lumi_ecm3=0.0;
  lumi_ecm4=0.0;
  // wrong initialization commented B.Dalena 11/02/2010
  //   hadrons_ee[1]=0.0;
  //   hadrons_ee[2]=0.0;
  //   hadrons_ee[3]=0.0;
  //   hadrons_eg[1]=0.0;
  //   hadrons_eg[2]=0.0;
  //   hadrons_eg[3]=0.0;
  //   hadrons_ge[1]=0.0;
  //   hadrons_ge[2]=0.0;
  //   hadrons_ge[3]=0.0;
  //   hadrons_gg[1]=0.0;
  //   hadrons_gg[2]=0.0;
  //   hadrons_gg[3]=0.0;
  upsmax=0.0;
}



void RESULTS::bpm_signal(const BEAM& beam) 
{
  double sum_vx,sum_vy,sum_vx_2,sum_vy_2;
  int j;
  unsigned int i;
  sum_vx=0.0;
  sum_vy=0.0;
  double vx, vy;
  int n_partic = beam.particle_beam().numberOfParticles();
  if (n_partic <= 0 ) return;
  int n_slice = beam.number_of_slices();
  for (j = 0; j < n_slice; j++)
    {
      const vector<PARTICLE*>& particle = beam.getParticleVector(j);
      for (i=0;i< particle.size();i++)
	{
	  particle[i]->velocities(vx, vy);
	  sum_vx += vx;
	  sum_vy += vy;
	}
    }
  sum_vx/=(double)n_partic;
  sum_vy/=(double)n_partic;
  sum_vx_2=0.0;
  sum_vy_2=0.0;
  for (j = 0; j < n_slice; j++)
    {
      const vector<PARTICLE*>& particle = beam.getParticleVector(j);
      for (i=0;i<particle.size();i++)
	{
	  particle[i]->velocities(vx, vy);
	  sum_vx_2 += (vx-sum_vx)*(vx-sum_vx);
	  sum_vy_2 += (vy-sum_vy)*(vy-sum_vy);
	}
    }
  sum_vx_2/=(double)n_partic;
  sum_vy_2/=(double)n_partic;
  int nobeam = beam.sign_label();
  if (nobeam == 1) 
    {
      c_vx_1=sum_vx;
      c_vy_1=sum_vy;
      sig_vx_1=sqrt(sum_vx_2);
      sig_vy_1=sqrt(sum_vy_2);
    }
  else
    {
      c_vx_2=sum_vx;
      c_vy_2=sum_vy;
      sig_vx_2=sqrt(sum_vx_2);
      sig_vy_2=sqrt(sum_vy_2);
    }
}

//  routine to calculate the average and RMS angles of the beam particles after
//  the interaction
void RESULTS::bpm_signal_coherent(const BEAM& beam) 
{
  double sum_vx,sum_vy,sum_vx_2,sum_vy_2;
  int l,j,n,sum=0;
  
  unsigned int point,i;
  sum_vx=0.0;
  sum_vy=0.0;
  int n_partic = beam.numberOfParticles();
  if (n_partic <= 0 ) return;
  double vx, vy;
  int n_slice = beam.number_of_slices();
  for (j = 0; j < n_slice; j++)
    {
      const vector<PARTICLE*>& particle = beam.getParticleVector(j);
      
      for (i=0;i<particle.size();i++)
	{
	  particle[i]->velocities(vx, vy);
	  sum_vx += vx;
	  sum_vy += vy;
	}
    }
  n = beam.number_of_coherent_vectors();
  sum=n_partic;
  for (l=0; l < n; l++)
    {
      const  vector<PARTICLE*>& coherent = beam.getCoherentVector(l);

      for (point=0; point < coherent.size(); point++)
	{
	  coherent[point]->velocities(vx, vy);
	  if (coherent[point]->energy()>0.0) 
	    {
	      sum_vx += vx;
	      sum_vy += vy;
	      ++sum;
	    }
	  else 
	    {
	      sum_vx -= vx;
	      sum_vy -= vy;
	      --sum;
	    }
	}
    }
    
  sum_vx/=(double)sum;
  sum_vy/=(double)sum;
  sum_vx_2=0.0;
  sum_vy_2=0.0;
  for (j = 0; j < n_slice; j++)
    {
      const vector<PARTICLE*>& particle = beam.getParticleVector(j);
      
      for (i=0;i< particle.size();i++)
	{
	  particle[i]->velocities(vx, vy);
	  sum_vx_2 += (vx-sum_vx)*(vx-sum_vx);
	  sum_vy_2 += (vy-sum_vy)*(vy-sum_vy);
	}
    }
  sum_vx_2/=(double)n_partic;
  sum_vy_2/=(double)n_partic;
  int nobeam = beam.sign_label();
  if (nobeam == 1)
    {
      c_vx_1_coh=sum_vx;
      c_vy_1_coh=sum_vy;
      sig_vx_1=sqrt(sum_vx_2);
      sig_vy_1=sqrt(sum_vy_2);
    }
  else
    {
      c_vx_2_coh=sum_vx;
      c_vy_2_coh=sum_vy;
      sig_vx_2=sqrt(sum_vx_2);
      sig_vy_2=sqrt(sum_vy_2);
    }
}


string  RESULTS::output_flow() const 
{
  int j1,j2;
  ostringstream sortie;
  sortie << titre(string("general results"));
  sortie <<  "lumi_fine = " << lumi_fine << " m**(-2) (from charge densities) " << endl;
  sortie <<  "lumi_ee = " << lumi_ee << " m**(-2) (from beam particles collisions)" << endl;
  sortie <<  "lumi_ee_high = " << lumi_ee_high << " 1/m2 (par bunch cross. above energy ecm_min)" << endl;
  sortie <<  "lumi_pp = " << lumi_pp << " m**(-2) " << endl;
  sortie <<  "lumi_eg = " << lumi_eg << " m**(-2) (e - gamma) " << endl;
  sortie <<  "lumi_ge = " << lumi_ge << " m**(-2) ( gamma - e) " << endl;
  sortie <<  "lumi_gg = " << lumi_gg << " m**(-2) (gamma - gamma) " << endl;
  float f_rep = switches_->get_f_rep();
  float n_b = switches_->get_n_b();
  //sortie <<  "lumi_gg_high = " << lumi_gg_high*f_rep*n_b << " 1/m2 (gamma - gamma, with c.o.m energy more than gg_cut) " << endl;
  sortie <<  "lumi_gg_high = " << lumi_gg_high << " 1/m2 (gamma - gamma, with c.o.m energy more than gg_cut) " << endl;
  
  for (j1=0;j1<3;j1++)
    {
      for (j2=0;j2<3;j2++)
	{
	  sortie <<  "lumi[" << j1 << "][" << j2 << "] = " << lumi[j1][j2]*f_rep*n_b << endl;
	}
    }
  sortie <<  "upsmax= " << upsmax << " (maximal value of the beamstrahlung paramater that occured) " << endl;

  double temp1, temp2;
  const float eps=1e-6;
  if (lumi_ee > eps)
    {
      temp1=lumi_ecm/lumi_ee;
      temp2 = sqrt( max( 0.0, lumi_ecm2/max(1.0,lumi_ee)-temp1*temp1) );
    }
  else
    {
      temp1 = -1.0;
      temp2 = temp1;
    }
  sortie << "E_cm = " << temp1 << " E_cm_var = " << temp2 << endl;
  sortie << " ............................................... " << endl;

  sortie << " average and RMS angles (x or y ) of the particles of each beam after the interaction (microradians) : " << endl;

  sortie << "bpm_vx.1=" << c_vx_1*1e6  <<  ";bpm_sig_vx.1=" << sig_vx_1*1e6 << ";" << endl;
  sortie << "bpm_vy.1=" << c_vy_1*1e6  <<  ";bpm_sig_vy.1=" << sig_vy_1*1e6 << ";" << endl;
  sortie << "bpm_vx.2=" << c_vx_2*1e6  <<  ";bpm_sig_vx.2=" << sig_vx_2*1e6 << ";" << endl;
  sortie << "bpm_vy.2=" << c_vy_2*1e6  <<  ";bpm_sig_vy.2=" << sig_vy_2*1e6 << ";" << endl;
  sortie << " average and RMS angles (x or y ) including coherent particles of each beam after the interaction (microradians) : " << endl;
  sortie << "bpm_vx_coh.1=" << c_vx_1_coh*1e6 << ";" << endl; 
  sortie << "bpm_vy_coh.1=" << c_vy_1_coh*1e6 << ";" << endl;
  sortie << "bpm_vx_coh.2=" << c_vx_2_coh*1e6 << ";" << endl;
  sortie << "bpm_vy_coh.2=" << c_vy_2_coh*1e6 << ";" << endl;
  sortie << " ............................................... " << endl;
  sortie << " minimal photon-photon center of mass energies for hadronic events (GeV) : " << endl;
  sortie << "hadron_cut.1=" << 2.0 << ";"<<  "hadron_cut.2=" << 5.0 << ";"<< "hadron_cut.3=" << 10.0 << ";"<< endl;
  sortie << " nb of hadronic events per bunch crossing due to the virtual photons in ee collisions : " << endl;
  sortie << "hadrons_ee.1=" << hadrons_ee[0] << ";hadrons_ee.2=" << hadrons_ee[1] << ";hadrons_ee.3=" << hadrons_ee[2] << ";"<< endl;
  sortie << " nb of hadronic evts per bx due to e-gamma collisions " << endl;
  sortie << "hadrons_eg.1=" <<  hadrons_eg[0] << ";hadrons_eg.2=" << hadrons_eg[1] << ";hadrons_eg.3=" << hadrons_eg[2] << ";"<< endl;
  sortie << " .. due to gamma-e collisions : " << endl;
  sortie << "hadrons_ge.1=" << hadrons_ge[0] << ";hadrons_ge.2=" << hadrons_ge[1] << ";hadrons_ge.3=" << hadrons_ge[2] << ";"<< endl;
  sortie << " ... due to gamma-gamma collisions : " << endl;
  sortie << "hadrons_gg.1=" << hadrons_gg[0] << ";hadrons_gg.2=" << hadrons_gg[1] << ";hadrons_gg.3=" << hadrons_gg[2] << ";"<< endl;
  sortie << " ... due to all types of  collisions : " << endl;
  sortie << "hadrons_sum.1=" << hadrons_ee[0]+hadrons_eg[0]+hadrons_ge[0]+hadrons_gg[0] << ";hadrons_sum.2=" << hadrons_ee[1]+hadrons_eg[1]+hadrons_ge[1]+hadrons_gg[1] << ";hadrons_sum.3=" << hadrons_ee[2]+hadrons_eg[2]+hadrons_ge[2]+hadrons_gg[2] << ";" << endl;
  sortie << " ............................................... " << endl << endl;
   
  sortie <<  "lumi_ee=" << lumi_ee << ";"<< endl;   
  sortie <<  "lumi_ee_high=" << lumi_ee_high << ";"<< endl;   
  return sortie.str();
}

void PAIRS_RESULTS::set()
{
  number_=0;
  energy_=0.0;
  eproc_[0]=0.0;
  eproc_[1]=0.0;
  eproc_[2]=0.0;
  nproc_[0]=0.0;
  nproc_[1]=0.0;
  nproc_[2]=0.0;
  
  n1_=0.0;
  b1_=0.0;
  n2_=0.0;
  b2_=0.0;
  
  highptsum_ = 0.0;
  highpteng_ = 0.0;
}

// 'composante' is old 'scdn1+scdn2'
void PAIRS_RESULTS::store_full_pair(int composante, double e1,double px1,double py1,double pz1,double e2,double px2,double py2,double pz2, double wgt, bool luckyPair )
{
  update_contribution(composante,e1,px1,py1,pz1,wgt);
  update_contribution(composante,e2,px2,py2,pz2,wgt);
  eproc_[composante]+= (fabs(e1)+fabs(e2)) * wgt;
  nproc_[composante]+= 2.0*wgt;
  if (luckyPair)
    {
      number_++;
      energy_ += fabs(e1);
      number_++;
      energy_ += fabs(e2);
      ;
    }
} /* store_full_pair */

void PAIRS_RESULTS::update_contribution(int composante, double ener,double px,double py,double pz, double wgt)
{
  double pt;
  pt=sqrt(px*px+py*py);
  // particles with a transverse momentum of more than 20 MeV and 
  //an angle with respect to the beam axis of more than 150 mrad.
  if ((pt>2.0e-2)&&(atan2(pt,fabs(pz))>0.15)) 
    {
      highptsum_ += wgt;
      highpteng_ += wgt * fabs(ener);
    }
  if(pz>0.0)
    {
      n1_ += wgt;
      b1_ += wgt*fabs(ener);
    }
  else
    {
      n2_ += wgt;
      b2_ += wgt*fabs(ener);
    }
}



string PAIRS_RESULTS::output_flow() const 
{
  ostringstream sortie;
  sortie << titre(string(" results for secondaries"));
  sortie << "nb of pairs per bunch crossing : n_pairs = " << number_ << " total energy e_pairs = " << energy_ << " GeV " << endl;
  sortie << endl;
  sortie << " Breit-Wheeler process : n_BW = " << nproc_[0] <<  " e_BW = " << eproc_[0] << " GeV " << endl;
  sortie << " Bethe-Heitler process : n_BH = " << nproc_[1] <<  " e_BH = " << eproc_[1] << " GeV " << endl;
  sortie << " Landau-Lifschitz process : n_LL = " << nproc_[2] << " e_LL = " << eproc_[2] <<  " GeV " << endl;
  sortie << endl;
  sortie << " total numer (an total energy, in GeV) of the particles due to the 3 processes : " << endl;
  sortie << "n.1 = " << n1_ << " b.1 = " << b1_ << endl;
  sortie << "n.2 = " << n2_ << " b.2 = " << b2_ << endl;
  sortie << endl;
  sortie << " nb of particles (and energy in GeV) due to thr 3 previous processes with a transverse momentum of more than 20 MeV and an angle with respect to the beam axis of more than 150 mrad : " << endl;
  sortie << "n_pt=" << highptsum_<< ";e_pt=" << highpteng_ <<";" << endl;
  sortie << "n_pairs=" << number_<< ";e_pairs=" << energy_ << ";" << endl;
  return sortie.str();
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
COMPT_RESULTS::COMPT_RESULTS()
{
  number_=0;
  energy_=0.0;
  
  eproc_[0]=0.0;
  eproc_[1]=0.0;
  nproc_[0]=0.0;
  nproc_[1]=0.0;

  n1_=0.0;
  b1_=0.0;
  n2_=0.0;
  b2_=0.0;
}


// 'composante' is old 'scdn1+scdn2'
int COMPT_RESULTS::store_compt(int composante, double e,double px,double py,double pz,double wgt,RNDM& hasard)
{
    double pt;
   pt=sqrt(px*px+py*py);
    if(pz>0.0)
      {
	n1_ += wgt;
	b1_ += wgt*fabs(e);
      }
    else
      {
	n2_ += wgt;
	b2_ += wgt*fabs(e);
      }
	
    eproc_[composante] += fabs(e) * wgt;
    nproc_[composante] += wgt;
    if (wgt < hasard.rndm()) 
      {
	return 0;
      }
    number_++;
    energy_ += fabs(e);
    return 1;
} /* store_compt */


/*
  MUON_RESULTS:: MUON_RESULTS()
  {
  eproc[0]=0.0;
  eproc[1]=0.0;
  eproc[2]=0.0;
  nproc[0]=0.0;
  nproc[1]=0.0;
  nproc[2]=0.0;
  }
*/

string COHERENT_RESULTS::output_flow() const 
{
  ostringstream sortie;
  sortie << titre(string("coherent particles results"));
  
  if (do_coherent_)
    {
      sortie << "coherent.sum=" <<sum_ << ";" << endl;
      sortie << "coherent.sumeng=" << sumeng_ << ";" << endl;
      sortie << "coherent.upsmax=" << upsmax_ << ";" << endl;
      sortie <<  "coherent.probmax=" <<  probmax_ << ";" << endl;
      sortie << "coherent.count=" << count_ << ";" << endl;
      sortie << "coherent.total_energy=" << total_energy_ << ";" << endl;
    }
  else
    {

      sortie << "coherent.sum=" << -1.0 << ";" <<endl;
      sortie << "coherent.sumeng=" << -1.0 << ";" <<endl;
      sortie << "coherent.upsmax=" << -1.0 << ";" <<endl;
      sortie << "coherent.count=" << -1 << ";" <<endl;
      sortie << "coherent.total_energy=" << -1.0 << ";" <<endl;
    }
  return sortie.str();
}

string TRIDENT_RESULTS::output_flow() const 
{
  ostringstream sortie;
  sortie << titre(string("trident results"));
  
  if (do_trident_)
    {
//       sortie << "coherent.sum = " <<sum_ << endl;
//       sortie << "coherent.sumeng = " << sumeng_ << endl;
//       sortie << "coherent.upsmax = " << upsmax_ << endl;
//       sortie <<  "coherent.probmax = " <<  probmax_ << endl;
       sortie << "trident.count=" << count_ << ";" <<endl;
       sortie << "trident.total_energy=" << total_energy_ << ";" <<endl;
    }
  else
    {

//       sortie << "coherent.sum = " << -1.0 << endl;
//       sortie << "coherent.sumeng = " << -1.0 << endl;
//       sortie << "coherent.upsmax = " << -1.0 << endl;
       sortie << "trident.count=" << -1 << ";" <<endl;
       sortie << "trident.total_energy=" << -1.0 << ";" <<endl;
    }
  return sortie.str();
}
