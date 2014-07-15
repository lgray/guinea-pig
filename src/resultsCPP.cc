#include "resultsCPP.h"

#include <vector>

RESULTS::RESULTS() : switches_(NULL),minijets(0.0),eloss_1(0.0),eloss_2(0.0),
		     ephot_1(0.0),ephot_2(0.0),c_vx_1(0.0),sig_vx_1(0.0),
		     c_vy_1(0.0),sig_vy_1(0.0),c_vx_2(0.0),sig_vx_2(0.0),
		     c_vy_2(0.0),sig_vy_2(0.0),c_vx_1_coh(0.0),c_vx_2_coh(0.0),
		     c_vy_1_coh(0.0),c_vy_2_coh(0.0)
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
  upsmax=0.0;
  lumi_ee_step.reserve(2);
  lumi_ee_high_step.reserve(2);
  for(int k=0;k<2;k++) step_lumi[k]=0.0;	
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
      const std::vector<PARTICLE*>& particle = beam.getParticleVector(j);
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
      const std::vector<PARTICLE*>& particle = beam.getParticleVector(j);
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
      const std::vector<PARTICLE*>& particle = beam.getParticleVector(j);
      
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
      const  std::vector<PARTICLE*>& coherent = beam.getCoherentVector(l);

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
      const std::vector<PARTICLE*>& particle = beam.getParticleVector(j);
      
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


std::string  RESULTS::output_flow() const 
{
  int j1,j2;
  std::ostringstream out;
  out << title(std::string("general results"));
  out <<  "lumi_fine = " << lumi_fine << " m**(-2) (from charge densities) " << std::endl;
  out <<  "lumi_ee = " << lumi_ee << " m**(-2) (from beam particles collisions)" << std::endl;
  out <<  "lumi_ee_high = " << lumi_ee_high << " 1/m2 (par bunch cross. above energy ecm_min)" << std::endl;
  out <<  "lumi_pp = " << lumi_pp << " m**(-2) " << std::endl;
  out <<  "lumi_eg = " << lumi_eg << " m**(-2) (e - gamma) " << std::endl;
  out <<  "lumi_ge = " << lumi_ge << " m**(-2) ( gamma - e) " << std::endl;
  out <<  "lumi_gg = " << lumi_gg << " m**(-2) (gamma - gamma) " << std::endl;
  float f_rep = switches_->get_f_rep();
  float n_b = switches_->get_n_b();
  //out <<  "lumi_gg_high = " << lumi_gg_high*f_rep*n_b << " 1/m2 (gamma - gamma, with c.o.m energy more than gg_cut) " << std::endl;
  out <<  "lumi_gg_high = " << lumi_gg_high << " 1/m2 (gamma - gamma, with c.o.m energy more than gg_cut) " << std::endl;
  
  for (j1=0;j1<3;j1++)
    {
      for (j2=0;j2<3;j2++)
	{
	  out <<  "lumi[" << j1 << "][" << j2 << "] = " << lumi[j1][j2]*f_rep*n_b << std::endl;
	}
    }
  out <<  "upsmax= " << upsmax << " (maximal value of the beamstrahlung paramater that occured) " << std::endl;

  double temp1, temp2;
  const float eps=1e-6;
  if (lumi_ee > eps)
    {
      temp1=lumi_ecm/lumi_ee;
      temp2 = sqrt( std::max( 0.0, lumi_ecm2/std::max(1.0,lumi_ee)-temp1*temp1) );
    }
  else
    {
      temp1 = -1.0;
      temp2 = temp1;
    }
  out << "E_cm = " << temp1 << " E_cm_var = " << temp2 << std::endl;
  out << " ............................................... " << std::endl;

  out << " average and RMS angles (x or y ) of the particles of each beam after the interaction (microradians) : " << std::endl;

  out << "bpm_vx.1=" << c_vx_1*1e6  <<  ";bpm_sig_vx.1=" << sig_vx_1*1e6 << ";" << std::endl;
  out << "bpm_vy.1=" << c_vy_1*1e6  <<  ";bpm_sig_vy.1=" << sig_vy_1*1e6 << ";" << std::endl;
  out << "bpm_vx.2=" << c_vx_2*1e6  <<  ";bpm_sig_vx.2=" << sig_vx_2*1e6 << ";" << std::endl;
  out << "bpm_vy.2=" << c_vy_2*1e6  <<  ";bpm_sig_vy.2=" << sig_vy_2*1e6 << ";" << std::endl;
  out << " average and RMS angles (x or y ) including coherent particles of each beam after the interaction (microradians) : " << std::endl;
  out << "bpm_vx_coh.1=" << c_vx_1_coh*1e6 << ";" << std::endl; 
  out << "bpm_vy_coh.1=" << c_vy_1_coh*1e6 << ";" << std::endl;
  out << "bpm_vx_coh.2=" << c_vx_2_coh*1e6 << ";" << std::endl;
  out << "bpm_vy_coh.2=" << c_vy_2_coh*1e6 << ";" << std::endl;
  out << " ............................................... " << std::endl;
  out << " minimal photon-photon center of mass energies for hadronic events (GeV) : " << std::endl;
  out << "hadron_cut.1=" << 2.0 << ";"<<  "hadron_cut.2=" << 5.0 << ";"<< "hadron_cut.3=" << 10.0 << ";"<< std::endl;
  out << " nb of hadronic events per bunch crossing due to the virtual photons in ee collisions : " << std::endl;
  out << "hadrons_ee.1=" << hadrons_ee[0] << ";hadrons_ee.2=" << hadrons_ee[1] << ";hadrons_ee.3=" << hadrons_ee[2] << ";"<< std::endl;
  out << " nb of hadronic evts per bx due to e-gamma collisions " << std::endl;
  out << "hadrons_eg.1=" <<  hadrons_eg[0] << ";hadrons_eg.2=" << hadrons_eg[1] << ";hadrons_eg.3=" << hadrons_eg[2] << ";"<< std::endl;
  out << " .. due to gamma-e collisions : " << std::endl;
  out << "hadrons_ge.1=" << hadrons_ge[0] << ";hadrons_ge.2=" << hadrons_ge[1] << ";hadrons_ge.3=" << hadrons_ge[2] << ";"<< std::endl;
  out << " ... due to gamma-gamma collisions : " << std::endl;
  out << "hadrons_gg.1=" << hadrons_gg[0] << ";hadrons_gg.2=" << hadrons_gg[1] << ";hadrons_gg.3=" << hadrons_gg[2] << ";"<< std::endl;
  out << " ... due to all types of  collisions : " << std::endl;
  out << "hadrons_sum.1=" << hadrons_ee[0]+hadrons_eg[0]+hadrons_ge[0]+hadrons_gg[0] << ";hadrons_sum.2=" << hadrons_ee[1]+hadrons_eg[1]+hadrons_ge[1]+hadrons_gg[1] << ";hadrons_sum.3=" << hadrons_ee[2]+hadrons_eg[2]+hadrons_ge[2]+hadrons_gg[2] << ";" << std::endl;
  out << " ............................................... " << std::endl << std::endl;
   
  out <<  "lumi_ee=" << lumi_ee << ";"<< std::endl;   
  out <<  "lumi_ee_high=" << lumi_ee_high << ";"<< std::endl;   
  out <<  "step " <<  " " <<  " lumi total " << " " <<  " lumi peak " << std::endl;
  for (int sstep=0;sstep< (int)lumi_ee_step.size();sstep++) out << sstep+1 << "  " << lumi_ee_step[sstep]<< "  " << lumi_ee_high_step[sstep] << std::endl;
  return out.str();
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
  name_="";
}

// 'index_of_process' is old 'scdn1+scdn2'
void PAIRS_RESULTS::store_full_pair(int index_of_process, double e1,double px1,double py1,double pz1,double e2,double px2,double py2,double pz2, double wgt, bool luckyPair )
{
  update_contribution(e1,px1,py1,pz1,wgt);
  update_contribution(e2,px2,py2,pz2,wgt);
  eproc_[index_of_process]+= (fabs(e1)+fabs(e2)) * wgt;
  nproc_[index_of_process]+= 2.0*wgt;
  if (luckyPair)
    {
      number_++;
      energy_ += fabs(e1);
      number_++;
      energy_ += fabs(e2);
      ;
    }
} /* store_full_pair */

void PAIRS_RESULTS::update_contribution(double ener,double px,double py,double pz, double wgt)
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


std::string PAIRS_RESULTS::output_flow() const 
{
  std::ostringstream out;
  //PAIRS
  std::string str1 ("pairs");

  if(name_.compare("pairs")==0){
    out << title(std::string(" results for ")+(std::string)name_+std::string(" "));
    out << "nb of pairs per bunch crossing : n_pairs = " << number_ << " total energy e_pairs = " << energy_ << " GeV " << std::endl;
    out << std::endl;
    out << " Breit-Wheeler process : n_BW = " << nproc_[0] <<  " e_BW = " << eproc_[0] << " GeV " << std::endl;
    out << " Bethe-Heitler process : n_BH = " << nproc_[1] <<  " e_BH = " << eproc_[1] << " GeV " << std::endl;
    out << " Landau-Lifschitz process : n_LL = " << nproc_[2] << " e_LL = " << eproc_[2] <<  " GeV " << std::endl;
    out << std::endl;
    out << " total nb (and total energy, in GeV) of the particles due to the 3 processes : " << std::endl;
    out << "n.1 = " << n1_ << " b.1 = " << b1_ << std::endl;
    out << "n.2 = " << n2_ << " b.2 = " << b2_ << std::endl;
    out << std::endl;
    out << " nb of particles (and energy in GeV) due to thr 3 previous processes with a transverse momentum of more than 20 MeV and an angle with respect to the beam axis of more than 150 mrad : " << std::endl;
    out << "n_pt=" << highptsum_<< ";e_pt=" << highpteng_ <<";" << std::endl;
    out << "n_pairs=" << number_<< ";e_pairs=" << energy_ << ";" << std::endl;
  } else {
    // MUONS
    out << title(std::string(" results for ")+(std::string)name_+std::string(" "));
    out << "nb of pairs per bunch crossing : n_"<< name_ <<" = " << number_ << " total energy e_"<<name_<<" = " << energy_ << " GeV " << std::endl;
    out << std::endl;
    out << " Breit-Wheeler process : n_BW_"<<name_<<" = " << nproc_[0] <<  " e_BW_"<<name_<<" = " << eproc_[0] << " GeV " << std::endl;
    out << " Bethe-Heitler process : n_BH_"<<name_<<" = " << nproc_[1] <<  " e_BH_"<<name_<<" = " << eproc_[1] << " GeV " << std::endl;
    out << " Landau-Lifschitz process : n_LL_"<<name_<<" = " << nproc_[2] << " e_LL_"<<name_<<" = " << eproc_[2] <<  " GeV " << std::endl;
    out << std::endl;
    out << " total numer (an total energy, in GeV) of the particles due to the 3 processes : " << std::endl;
    out << "n.1_"<<name_<<" = " << n1_ << " b.1_"<<name_<<" = " << b1_ << std::endl;
    out << "n.2_"<<name_<<" = " << n2_ << " b.2_"<<name_<<" = " << b2_ << std::endl;
    out << std::endl;
    out << " nb of particles (and energy in GeV) due to thr 3 previous processes with a transverse momentum of more than 20 MeV and an angle with respect to the beam axis of more than 150 mrad : " << std::endl;
    out << "n_pt_"<<name_<<"=" << highptsum_<< ";e_pt_"<<name_<<"=" << highpteng_ <<";" << std::endl;
    out << "n_"<<name_<<"=" << number_<< ";e_"<<name_<<"=" << energy_ << ";" << std::endl;
  }
  return out.str();
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


// 'index_of_process' is old 'scdn1+scdn2'
int COMPT_RESULTS::store_compt(int index_of_process, double e,double /*px*/,double /*py*/,double pz,double wgt,RNDM& rndm_generator)
{
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
	
    eproc_[index_of_process] += fabs(e) * wgt;
    nproc_[index_of_process] += wgt;
    if (wgt < rndm_generator.rndm()) 
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

std::string COHERENT_RESULTS::output_flow() const 
{
  std::ostringstream out;
  out << title(std::string("coherent particles results"));
  
  if (do_coherent_)
    {
      out << "coherent.sum=" <<sum_ << ";" << std::endl;
      out << "coherent.sumeng=" << sumeng_ << ";" << std::endl;
      out << "coherent.upsmax=" << upsmax_ << ";" << std::endl;
      out <<  "coherent.probmax=" <<  probmax_ << ";" << std::endl;
      out << "coherent.count=" << count_ << ";" << std::endl;
      out << "coherent.total_energy=" << total_energy_ << ";" << std::endl;
    }
  else
    {

      out << "coherent.sum=" << -1.0 << ";" <<std::endl;
      out << "coherent.sumeng=" << -1.0 << ";" <<std::endl;
      out << "coherent.upsmax=" << -1.0 << ";" <<std::endl;
      out << "coherent.count=" << -1 << ";" <<std::endl;
      out << "coherent.total_energy=" << -1.0 << ";" <<std::endl;
    }
  return out.str();
}

std::string TRIDENT_RESULTS::output_flow() const 
{
  std::ostringstream out;
  out << title(std::string("trident results"));
  
  if (do_trident_)
    {
//       out << "coherent.sum = " <<sum_ << std::endl;
//       out << "coherent.sumeng = " << sumeng_ << std::endl;
//       out << "coherent.upsmax = " << upsmax_ << std::endl;
//       out <<  "coherent.probmax = " <<  probmax_ << std::endl;
       out << "trident.count=" << count_ << ";" <<std::endl;
       out << "trident.total_energy=" << total_energy_ << ";" <<std::endl;
    }
  else
    {

//       out << "coherent.sum = " << -1.0 << std::endl;
//       out << "coherent.sumeng = " << -1.0 << std::endl;
//       out << "coherent.upsmax = " << -1.0 << std::endl;
       out << "trident.count=" << -1 << ";" <<std::endl;
       out << "trident.total_energy=" << -1.0 << ";" <<std::endl;
    }
  return out.str();
}
