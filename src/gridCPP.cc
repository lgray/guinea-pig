#include <iostream>
#include "gridCPP.h"
#include "mathconst.h"
#include "physconst.h"
#include "mathematicalTools.h"
#include "physicalTools.h"

using namespace std;

GENERAL_GRID::GENERAL_GRID()
{
  n_cell_x_=32;
  n_cell_y_=32;
  n_cell_z_=32;
  timestep_ =10;

  min_x_ = 0.0;
  max_x_ = 0.0; 
  min_y_ = 0.0;
  max_y_ = 0.0;
  min_z_ = 0.0;
  max_z_ = 0.0;

  cut_x_ = 0.0;
  cut_y_ = 0.0;
  cut_z_ = 0.0;

  integration_method_ = 0;
  step_ = 0.0;


  rho_x_1_ = 0.0;
  rho_y_1_ = 0.0;
  rho_sum_1_ = 0.0;
  rho_x_2_ = 0.0;
  rho_y_2_ = 0.0;
  rho_sum_2_ = 0.0;

  delta_x_inverse_ = 0.0;
  delta_y_inverse_ = 0.0;

  rho_factor_ = 0.0;

}

 GENERAL_GRID::GENERAL_GRID(const GENERAL_GRID& grid)
{
  n_cell_x_ = grid.n_cell_x_;
  n_cell_y_ = grid.n_cell_y_;
  n_cell_z_ = grid.n_cell_z_;
  timestep_ = grid.timestep_;
  mesh_ = grid.mesh_;
  min_x_ = grid.min_x_;
  max_x_ = grid.max_x_; 
  min_y_ = grid.min_y_;
  max_y_ = grid.max_y_;
  min_z_ = grid.min_z_;
  max_z_ = grid.max_z_;

  cut_x_ = grid.cut_x_;
  cut_y_ = grid.cut_y_;
  cut_z_ = grid.cut_z_;
  


  integration_method_ = grid.integration_method_;
  step_ = grid.step_;
  
  
  rho_x_1_ = grid.rho_x_1_;
  rho_y_1_ = grid.rho_y_1_;
  rho_sum_1_ = grid.rho_sum_1_;
  rho_x_2_ = grid.rho_x_2_;
  rho_y_2_ = grid.rho_y_2_;
  rho_sum_2_ = grid.rho_sum_2_;
  
  delta_x_inverse_ = grid.delta_x_inverse_;
  delta_y_inverse_ = grid.delta_y_inverse_;
  
  rho_factor_ = grid.rho_factor_;
  
  slice_of_beam_[0] = SLICE_ON_GRID(grid.slice_of_beam_[0]);
  slice_of_beam_[1] = SLICE_ON_GRID(grid.slice_of_beam_[1]);

  
  champ_ = grid.champ_;
  
}

void GENERAL_GRID::interpolePotential(float xpart,float ypart, PHI_FLOAT& h_x, PHI_FLOAT& h_y, PHI_FLOAT& phi1_x, PHI_FLOAT& phi2_x, PHI_FLOAT& phi3_x, PHI_FLOAT& phi1_y,PHI_FLOAT& phi2_y, PHI_FLOAT& phi3_y, const PHI_FLOAT *phi) const
{
  int i1, i2;
  float h;
  absoluteCoordinates(xpart, ypart,h_x, h_y);

  // locate in Y
  cellLocalization(h_x, h_y, i1, i2, h);

  phiValuesY(i1, i2,h, phi, phi1_y, phi2_y, phi3_y);

  // locate in X
  cellLocalization(h_y, h_x,i2,i1,h); 


  phiValuesX(i1, i2, h, phi, phi1_x, phi2_x,phi3_x );


  h_x -= floor(h_x);
  h_y -= floor(h_y);
}

// electric field in GV/nm
TRIVECTOR GENERAL_GRID::ElectricFieldCIC(float xpart,float ypart, const PHI_FLOAT *phi) const
{
  //  int k;
  PHI_FLOAT phi1_x,phi2_x,phi3_x,phi1_y,phi2_y,phi3_y,h_x,h_y;
  TRIVECTOR Efield;
  interpolePotential(xpart, ypart, h_x, h_y, phi1_x, phi2_x, phi3_x, phi1_y, phi2_y, phi3_y, phi);

  Efield(0) = (h_x*(phi1_y-phi2_y)+(1.0-h_x)*(phi2_y-phi3_y))*delta_x_inverse_;
  Efield(1) = (h_y*(phi1_x-phi2_x)+(1.0-h_y)*(phi2_x-phi3_x))*delta_y_inverse_;

  Efield(2) = 0.0;
  // the potential seems to have been multiplied by 2 (for taking into 
  // account the magnetic field)
  Efield *= 0.5;
  return Efield;
}

void GENERAL_GRID::computeFields(int integrationMethod, float charge_sign, PHI_FLOAT *sor_parameter) 
    {
      // (see commentaries in init_grid_phys() )

      //int nn[2];
      switch (integrationMethod)
	{
	case 1:	  
	  champ_.foldfields(slice_of_beam_[0].get_rho(),slice_of_beam_[1].get_rho(),charge_sign);
	  break;
	case 2: 
	  
	  champ_.fold_fft(slice_of_beam_[0].get_rho(),slice_of_beam_[1].get_rho(),charge_sign);
	  break; 
	       
	case 3:
	  champ_.foldfronts(slice_of_beam_[0].get_rho(),slice_of_beam_[1].get_rho() , sor_parameter, charge_sign);
	}
    }


GENERAL_GRID::~GENERAL_GRID()
{
  //  if (rho1_ != NULL) delete [] rho1_;
  //  if (rho2_ != NULL) delete [] rho2_;
}



void  GENERAL_GRID::init_grid_phys (float n_particles1, float n_particles2, float cut_x,float cut_y,float cut_z, float charge_sign,   FFT_SERVER* fourier)
{
  //int i1,i2,j,j0,i;
  // PHI_FLOAT factor,phi0;
   PHI_FLOAT factor;
   //double x0,y0;
   //int nn[2];

  float offsetx, offsety;
  float deltax, deltay, deltaz;
  
  slice_of_beam_[0].update_nb_part_per_macro(n_particles1);
  slice_of_beam_[1].update_nb_part_per_macro(n_particles2);
    /////
  cut_x_ = cut_x;
  cut_y_ = cut_y;
  cut_z_ = cut_z;
  min_x_=-((float)n_cell_x_-2)/((float)n_cell_x_)*cut_x;
  max_x_=((float)n_cell_x_-2)/((float)n_cell_x_)*cut_x;
  min_y_=-((float)n_cell_y_-2)/((float)n_cell_y_)*cut_y;
  max_y_=((float)n_cell_y_-2)/((float)n_cell_y_)*cut_y;
  offsetx=0.5*(float)n_cell_x_;
  offsety=0.5*(float)n_cell_y_;
  deltax=2.0*cut_x/((float)n_cell_x_);
  deltay=2.0*cut_y/((float)n_cell_y_);
  deltaz=2.0*cut_z/((float)n_cell_z_);
  mesh_ = MESH(deltax, deltay, deltaz, offsetx, offsety);
  delta_x_inverse_ = 1.0/deltax;
  delta_y_inverse_ = 1.0/deltay;


  min_z_ = -cut_z;
  max_z_ = cut_z;
  step_=deltaz/(2.0*(float)timestep_);
  //

  factor=-charge_sign*2.0*RE/(deltaz)*EMASSeV;
  rho_factor_ = 2.0*factor;
  factor/=deltax*deltay;
  champ_.dist_init(factor, deltax, deltay,fourier);
}


GRID::GRID() : GENERAL_GRID()
{
  photon_file_ = NULL;
  hadron_file_ = NULL;
  minijets_ = NULL;
  hasard_ = NULL;
}

GRID::~GRID() 
{
  if (photon_file_ != NULL) 
    {
      photon_file_->close();
      delete photon_file_;
    }
  if (hadron_file_ != NULL) 
    {
      hadron_file_->close();
      delete hadron_file_;
    }
  if (minijets_ != NULL) delete minijets_;
}




void  GRID::lumi_init(const SWITCHES& switches)
{
  int nmax;
  if(switches.get_do_lumi()){
    if (switches.get_do_lumi()&1) 
      {
	nmax = switches.get_num_lumi();
	if (nmax < 100) nmax = 100;
	lumi_heap_ee_ = LUMI_HEAP_EE(nmax, switches.get_lumi_p(), hasard_);
      }
    if (switches.get_do_lumi()&2) 
      {
	nmax = switches.get_num_lumi_eg();
	if (nmax < 100) nmax = 100;
	lumi_heap_eg_ = LUMI_HEAP(nmax,switches.get_lumi_p_eg(), hasard_);
	lumi_heap_ge_ = LUMI_HEAP(nmax, switches.get_lumi_p_eg(), hasard_);
      }
    if (switches.get_do_lumi()&4) 
      {
	nmax = switches.get_num_lumi_gg();
	if (nmax < 100) nmax = 100;
	lumi_heap_gg_ = LUMI_HEAP( nmax,switches.get_lumi_p_gg(), hasard_);
      }
  }
}

void GRID::save_lumi_on_files(SWITCHES& switches, string lumi_ee_out, string lumi_eg_out, string lumi_ge_out, string lumi_gg_out)
{
  if(switches.get_do_lumi())
    {
      if (switches.get_do_lumi()&1) 
	{
	  lumi_heap_ee_.saveLumi(lumi_ee_out);
	}
      if (switches.get_do_lumi()&2) 
	{	
	  lumi_heap_eg_.saveLumi(lumi_eg_out);
	  lumi_heap_ge_.saveLumi(lumi_ge_out);
	}
      if (switches.get_do_lumi()&4) 
	{
	  lumi_heap_gg_.saveLumi(lumi_gg_out);
	}
    }
}

void GRID::save_tertphot_on_file(string tertphotfile)
{
  FILE_IN_OUT filout;
  filout.open_file(tertphotfile,"w");
  for(int i =0; i<tertphot_.size();i++)
    {
      filout.save_object_on_persistent_file(&tertphot_[i]);
    }
  filout.close();
}

void GRID::read(const PARAMETERS& param, int automatic)
{
  int n_n;
  float n_f;
  if (automatic != 1) 
    {
      n_n = param.readIValue("n_x");
      if (n_n > 0) n_cell_x_= n_n;
      n_n = param.readIValue("n_y");
      if (n_n > 0) n_cell_y_ = n_n;
      n_n = param.readIValue("n_z");
      if (n_n > 0) n_cell_z_ = n_n;
    }
  n_n = param.readIValue("n_t");
  if (n_n > 0) timestep_ = n_n;
  n_n = param.readIValue("n_m.1");
  //  if (n_n > 0) nb_macroparticles_[0] = n_n;
  if (n_n > 0) slice_of_beam_[0].set_macroparticles(n_n);
  n_n = param.readIValue("n_m.2");
  //  if (n_n > 0) nb_macroparticles_[1] = n_n;
  if (n_n > 0) slice_of_beam_[1].set_macroparticles(n_n);

  // l'indexation est bizarre, mais c'est comme ca dans l'original guineapig C
  n_f = param.readFValue("scale_step.1");
  //  if (n_f > 0.0) scal_step[1] = n_f;
  if (n_f > 0.0) slice_of_beam_[1].set_scal_step(n_f);  
  n_f = param.readFValue("scale_step.2");
  //  if (n_f > 0.0) scal_step[0] = n_f;
  if (n_f > 0.0) slice_of_beam_[0].set_scal_step(n_f);
}


string GRID::output_flow() const 
{
  ostringstream sortie;
  sortie << titre(string("grid parameters"));
  sortie << "n_x = " << n_cell_x_ << " n_y = " << n_cell_y_ << endl;

  sortie << "n_z = " << n_cell_z_ << " n_t = " << timestep_ << endl;

  sortie << "n_m.1 = " << slice_of_beam_[0].get_macroparticles()  << " n_m.2 = " << slice_of_beam_[1].get_macroparticles() << endl;

  sortie << "cut_x = " << cut_x_ << " nm ; cut_y = " << cut_y_ << " nm ; cut_z = " << cut_z_*1e-3 << " micrometers " << endl;

  sortie << endl;
  sortie << " ...................................................... " << endl;
  sortie << " relative amount of interacting particles that were outside the grid during one time step : " << endl;
  sortie << "beam 1 : miss = " << distribute1_.delta <<  "  beam 2 : miss = " << distribute2_.delta << endl;
  sortie << " number of interacting part. that were outside the grid during one time step : " << endl;
  sortie << "out.1=" << distribute1_.tot <<  ";out.2=" << distribute2_.tot << ";" << endl; 


 return sortie.str();
}



void GRID::init_grid_comp(int n_cell_x, int n_cell_y, int integration_method,   FFT_SERVER* fourier)
{
  n_cell_x_ = n_cell_x;
  n_cell_y_ = n_cell_y;
  integration_method_ = integration_method;
   champ_ = FIELD(integration_method, n_cell_x_, n_cell_y_, fourier);
   slice_of_beam_[0].resize(n_cell_x_, n_cell_y_);
   slice_of_beam_[1].resize(n_cell_x_, n_cell_y_);

  part_pointer1_ =  vector< list<BEAM_PARTICLE_POINTER> >(n_cell_x_*n_cell_y_);
  part_pointer2_ =  vector< list<BEAM_PARTICLE_POINTER> >(n_cell_x_*n_cell_y_);


  grid_photon1_ = vector< list<BEAM_PHOTON_POINTER> >(n_cell_x_*n_cell_y_);
  grid_photon2_ = vector< list<BEAM_PHOTON_POINTER> >(n_cell_x_*n_cell_y_);
}


/*! Routine to calculate the parameters necessary for the iterative method of
   the potential calculation sor2 */

void GRID::init_sor2 (PHI_FLOAT *parameter)
{
  PHI_FLOAT factor;
  float deltax = mesh_.get_delta_x();
  float deltay = mesh_.get_delta_y();
  //  float deltaz = mesh_.get_delta_z();

  factor=1.0/(8.0*PI*RE/mesh_.get_delta_z()*1e9/2000.0);
  parameter[0]=factor*(deltay/deltax);
  parameter[1]=parameter[0];
  parameter[2]=factor*(deltax/deltay);
  parameter[3]=parameter[2];
  parameter[4]=-2.0*(parameter[0]+parameter[2]);
  parameter[5]=(deltay*deltay*cos(PI/(float)n_cell_x_)
		+deltax*deltax*cos(PI/(float)n_cell_y_))
                        /(deltax*deltax+deltay*deltay);
}



void GRID::check_distribute(int what)
{
    float tmp;
    if (what==0){
	distribute1_.delta=0.0;
	distribute2_.delta=0.0;
	distribute1_.tot=0;
	distribute2_.tot=0;
    }
    if (what<=1){
	distribute1_.in=0; distribute1_.out=0;
	distribute2_.in=0; distribute2_.out=0;
    }
    if (what==2){
	tmp=(float)distribute1_.out
	    /((float)(distribute1_.in+distribute1_.out));
	if (tmp>distribute1_.delta) distribute1_.delta=tmp;
	tmp=(float)distribute2_.out
	    /((float)(distribute2_.in+distribute2_.out));
	if (tmp>distribute2_.delta) distribute2_.delta=tmp;
	if (distribute1_.out>distribute1_.tot){
	  distribute1_.tot=distribute1_.out;
	}
	if (distribute2_.out>distribute2_.tot){
	  distribute2_.tot=distribute2_.out;
	}
    }
    if (what==3){
      	printf("miss_1=%f;miss_2=%f;\n",distribute1_.delta,distribute2_.delta);
    }
}




/*
  Distributes the beam particles for field calculation
*/

void GRID::distribute_particles(int i_slice1,
				 int i_slice2, 
				 int electron_distribution_rho, 
				 int force_symmetric)
{
  
  slice_of_beam_[0].razRho();
  slice_of_beam_[1].razRho();
  //  razRhos(); 

  switch (electron_distribution_rho)
    {
    case 1:
      assignBeamSliceNGP(slice_of_beam_[0], i_slice1, distribute1_);
      assignBeamSliceNGP(slice_of_beam_[1], i_slice2, distribute2_);
      break;
    case 2:

      assignBeamSliceCIC(slice_of_beam_[0], i_slice1, distribute1_);
      assignBeamSliceCIC(slice_of_beam_[1], i_slice2, distribute2_);
      break;

    }

  distribute_coherent_particles(i_slice1, i_slice2, electron_distribution_rho);
  distribute_trident_particles(i_slice1, i_slice2, electron_distribution_rho);  
  if (force_symmetric) 
    {
      slice_of_beam_[0].symmetrizeCharges(n_cell_x_,n_cell_y_);
      slice_of_beam_[1].symmetrizeCharges(n_cell_x_,n_cell_y_);
    }
}




/*
  Distribute the coherent pair particles for the calculation of the fields
 */

void GRID::distribute_coherent_particles(int i_slice1,
					  int i_slice2,
					  int electron_distribution_rho)
{
  switch (electron_distribution_rho){
    case 1:


      assignCoherentBeamSliceNGP(slice_of_beam_[0], i_slice1, distribute1_);
      assignCoherentBeamSliceNGP(slice_of_beam_[1], i_slice2, distribute2_);
      break;
  case 2:

    assignCoherentBeamSliceCIC(slice_of_beam_[0], i_slice1, distribute1_);
    assignCoherentBeamSliceCIC(slice_of_beam_[1], i_slice2, distribute2_);
    break;
  }
  
}

/*
  Distribute the trident pair particles for the calculation of the fields
 */
void GRID::distribute_trident_particles(int i_slice1,int i_slice2,int electron_distribution_rho)
{
  switch (electron_distribution_rho){
  case 1:
    assignTridentBeamSliceNGP(slice_of_beam_[0], i_slice1, distribute1_);
    assignTridentBeamSliceNGP(slice_of_beam_[1], i_slice2, distribute2_);
    break;
  case 2:
    assignTridentBeamSliceCIC(slice_of_beam_[0], i_slice1, distribute1_);
    assignTridentBeamSliceCIC(slice_of_beam_[1], i_slice2, distribute2_);
    break;
  }
}

//! Distributes the particles for background calculation 

void GRID::distribute_particles_for_background(int i_slice1,
				 int i_slice2,  
				 float electron_ratio)
{
  int j,i1,i2;
  float ratio,ratio_i_1,ratio_i_2;
  const float eps=1e-5;
  
    for (i1=0;i1<n_cell_x_;i1++)
      {
	for (i2=0;i2<n_cell_y_;i2++)
	  {
	    j=i1*n_cell_y_+i2;
	    part_pointer1_[j].clear();
	    part_pointer2_[j].clear();
	  }
      }

  ratio=electron_ratio;
  if (ratio<eps) return;
  float deltax = mesh_.get_delta_x();
  float deltay = mesh_.get_delta_y();

  ratio_i_1=1e9/ratio*slice_of_beam_[0].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  ratio_i_2=1e9/ratio*slice_of_beam_[1].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  
    distributeScatter1( slice_of_beam_[0].get_beam()->getParticleVector(i_slice1), ratio, ratio_i_1, part_pointer1_);
   distributeScatter1(slice_of_beam_[1].get_beam()->getParticleVector(i_slice2), ratio, ratio_i_2, part_pointer2_);
  distribute_coherent_particles_for_background(i_slice1, i_slice2, electron_ratio);
  distribute_trident_particles_for_background(i_slice1, i_slice2, electron_ratio);
}

//
void GRID::distributeScatter1(const vector<PARTICLE*>& lesParticles, float ratio, float ratio_i, vector< list<BEAM_PARTICLE_POINTER> >& part_pointer)
{

  float xpart, ypart;
  int i1, i2;
  //  int k,j;
  int j;
  unsigned int k;
  //    const vector<PARTICLE*>& lesParticles = slice_of_beam_[i_beam-1].get_beam()->getParticleVector(i_slice);
    for (k = 0; k < lesParticles.size(); k++)
    {
      if (hasard_->rndm()<ratio)
	{
	lesParticles[k]->XYposition(xpart, ypart);
	  if (particleInGrid(xpart, ypart, i1, i2))
	    {
	      j=i1*n_cell_y_+i2;
	      part_pointer[j].push_back(BEAM_PARTICLE_POINTER(lesParticles[k], ratio_i));
	    }
	}
    }
}



// cette methode n'est pas utilisee, consomme du temps de calcul. Ne pas 
// detruire pour le moment.
void GRID::distributeScatter2(const vector<PARTICLE*>& lesParticles,float ratio, float ratio_i, vector< list<BEAM_PARTICLE_POINTER> >& part_pointer)
{

  float xpart, ypart,poids,h_x,h_y;
  int i1, i2;
  //int k,j;
  int j;
  unsigned int k;

  //  const vector<PARTICLE*>& lesParticles = slice_of_beam_[i_beam-1].get_beam()->getParticleVector(i_slice);
    for (k = 0; k < lesParticles.size(); k++)
      {
	lesParticles[k]->XYposition(xpart, ypart);
	if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2))
	  {
	    if (hasard_->rndm()<ratio)
	      {
		j=i1*n_cell_y_+i2;
		poids = (1.0-h_x)*(1.0-h_y)*ratio_i;

	      part_pointer[j].push_back(BEAM_PARTICLE_POINTER(lesParticles[k], poids));
		j=(i1+1)*n_cell_y_+i2;
		poids =  h_x*(1.0-h_y)*ratio_i;
	      part_pointer[j].push_back(BEAM_PARTICLE_POINTER(lesParticles[k], poids));
		j=i1*n_cell_y_+i2+1;
		poids =  (1.0-h_x)*h_y*ratio_i;
	      part_pointer[j].push_back(BEAM_PARTICLE_POINTER(lesParticles[k], poids));
		j=(i1+1)*n_cell_y_+i2+1;
		poids = h_x*h_y*ratio_i;
			      part_pointer[j].push_back(BEAM_PARTICLE_POINTER(lesParticles[k], poids));
	      }
	  }
      }
}


//  Distributes the particles from coherent pair creation for the background
//  calculation
 

void GRID::distribute_coherent_particles_for_background(int i_slice1,
					  int i_slice2,
					  float electron_ratio)

{

  float ratio,ratio_i_1,ratio_i_2;
  const float eps=1e-5;

  ratio=electron_ratio;
  if (ratio<eps) return;
  float deltax = mesh_.get_delta_x();
  float deltay = mesh_.get_delta_y();

  ratio_i_1=1e9/ratio*slice_of_beam_[0].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  ratio_i_2=1e9/ratio*slice_of_beam_[1].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
    distributeScatter1(slice_of_beam_[0].get_beam()->getCoherentVector(i_slice1), ratio, ratio_i_1, part_pointer1_);
        distributeScatter1(slice_of_beam_[1].get_beam()->getCoherentVector(i_slice2),ratio, ratio_i_2, part_pointer2_);
}

//  Distributes the particles from trident pair creation for the background
//  calculation
void GRID::distribute_trident_particles_for_background(int i_slice1,int i_slice2,float electron_ratio)
{
  float ratio,ratio_i_1,ratio_i_2;
  const float eps=1e-5;

  ratio=electron_ratio;
  if (ratio<eps) return;
  float deltax = mesh_.get_delta_x();
  float deltay = mesh_.get_delta_y();

  ratio_i_1=1e9/ratio*slice_of_beam_[0].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  ratio_i_2=1e9/ratio*slice_of_beam_[1].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  distributeScatter1(slice_of_beam_[0].get_beam()->getTridentVector(i_slice1), ratio, ratio_i_1, part_pointer1_);
  distributeScatter1(slice_of_beam_[1].get_beam()->getTridentVector(i_slice2), ratio, ratio_i_2, part_pointer2_);
}

//! Distributes the virtual photons 

void GRID::distribute_virtual_photons(int i_slice1,
				 int i_slice2,
				 SWITCHES* switches,
				 double s4, double lns4)
{
   const float   emass2=EMASS*EMASS;
  float xmin,r_scal;
  float ratio,ratio_i_1,ratio_i_2;
    int geom;

  r_scal=switches->get_r_scal();
  geom=switches->get_geom();
  ratio=switches->get_electron_ratio();
  float deltax = mesh_.get_delta_x();
  float deltay = mesh_.get_delta_y();
  ratio_i_1=1e9/ratio*slice_of_beam_[0].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  ratio_i_2=1e9/ratio*slice_of_beam_[1].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  clear_extra_photons();
  //   extra_photons1_.clear_extra_photons();
  //   extra_photons2_.clear_extra_photons();
  //  switch (switches->get_electron_distribution_scatter())
  //   {
  //    case 1:
  xmin=switches->get_compt_x_min()*emass2/s4;
  
  electronScatter(extra_photon_pointer1_, 1, i_slice1, xmin, s4, lns4, ratio, switches->get_ext_field(), geom, ratio_i_1, r_scal);
  
  electronScatter( extra_photon_pointer2_,  2, i_slice2, xmin, s4, lns4, ratio, switches->get_ext_field(), geom, ratio_i_2, r_scal);
  

	  
  //     break;
  // case 2:
  
  //   distributeElectronScatter2(beam1_, i_slice1, ratio, ratio_i_1, part_pointer1_, hasard);
  //   distributeElectronScatter2(beam2_, i_slice2, ratio, ratio_i_2, part_pointer2_, hasard);
  //   break;
  //   }
}

void GRID::electronScatter(vector< list<EXTRA_PHOTON_POINTER> >& extra_phot_ptr, int i_beam, int i_slice, float xmin, double s4,double lns4, float ratio, int ext_field, int geom, float ratio_i, float r_scal)
{
  unsigned int k;
  int i1, i2,j;
  int i_equiv=6;
  int n_phot, i_phot;
  float energy, r_phot;
  float e_phot,q2,one_m_x,x,y, radius, theta;
  float xVelocity, yVelocity;
  //PHYSTOOLS phys;

    const vector<PARTICLE*>& lesParticles = slice_of_beam_[i_beam-1].get_beam()->getParticleVector(i_slice);
    for (k = 0; k < lesParticles.size(); k++)
      //  for (i=beam->firstParticleOfSlice(i_slice);i<beam->firstParticleOfSlice(i_slice+1);i++)
    {
      //      energy=beam->energyOfParticle(i);
      energy=lesParticles[k]->energy();
      //r_phot=phys.requiv(lns4,xmin,i_equiv)*ratio;
      r_phot=PHYSTOOLS::requiv(lns4,xmin,i_equiv)*ratio;
      n_phot=(int)floor(r_phot);
      r_phot -= n_phot;
      if(hasard_->rndm()<r_phot) n_phot++;
      for (i_phot=0;i_phot<n_phot;i_phot++)
	{

	  //phys.mequiv(s4, lns4,xmin,energy,i_equiv,&e_phot,&q2,&one_m_x, *hasard_);
	  PHYSTOOLS::mequiv(s4, lns4,xmin,energy,i_equiv,&e_phot,&q2,&one_m_x, *hasard_);
	  //#ifdef EXT_FIELD
	  if (ext_field && (pow(q2/(EMASS*EMASS),1.5)*energy*energy < e_phot*e_phot*lesParticles[k]->getUps()) )
	    {
	      e_phot=-1.0;
	    }
	  //#endif
	  if (e_phot>0.0)
	    {
	      switch(geom)
		{
		case 0:
		  lesParticles[k]->XYposition(x, y);
		  break;
		case 1:
		  radius=HBAR*Cvelocity/sqrt(q2*one_m_x)*r_scal*1e9;
		  radius=min(radius,float(1.0e5));

		  lesParticles[k]->XYposition(x, y);
		  x += hasard_->rndm_sincos(&theta)*radius;
		  y += theta*radius;
		  break;
		case 2:
		  radius=HBAR*Cvelocity/sqrt(q2*one_m_x)*r_scal*1e9;
		  radius=min(radius,float(1e5));
		  lesParticles[k]->XYposition(x, y);
		  x += hasard_->gasdev()*radius;
		  y += hasard_->gasdev()*radius;
		  break;
		}
	      if (particleInGrid(x, y, i1, i2))
		{
		  j=i1*n_cell_y_+i2;
		  lesParticles[k]->velocities(xVelocity, yVelocity);
		  extra_phot_ptr[j].push_back(EXTRA_PHOTON_POINTER(e_phot, xVelocity, yVelocity,q2,energy,ratio_i));
		  //		  photon->store_vir_photon(e_phot,xVelocity, yVelocity,q2,energy,ratio_i,j);
		}
	      
	    }
	}
    }
}





void GRID::assignBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute)
{

  unsigned int k;
  int i1, i2;
  float xpart, ypart;
  const vector<PARTICLE*>& lesParticles = sog.get_beam()->getParticleVector(i_slice);    
  for (k = 0; k < lesParticles.size(); k++)
    {
      lesParticles[k]->XYposition(xpart, ypart);
      if (particleInGrid(xpart, ypart, i1, i2))
	{
	  sog.assignChargeToNGP( i1, i2, 1.0);
	  distribute.in++;
	}
      else   distribute.out++;		
    }
}


void GRID::assignBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute)
{
  // n_macro : nb of particles per macroparticle
  unsigned int k;
  int i1, i2;
  float h_x,h_y;
  float xpart, ypart;
const vector<PARTICLE*>& lesParticles = sog.get_beam()->getParticleVector(i_slice);
 for (k = 0; k < lesParticles.size(); k++)
	{
	lesParticles[k]->XYposition(xpart, ypart);
	  if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2))
	    {
	      sog.assignChargeToCIC( i1, i2, h_x, h_y, 1.0);
	      distribute.in++;	    
	    }
	  else
	    {
	      //	      cout << " particule no "  << k+1 << " tranche " << i_slice << " hors grille : x= " << xpart << " y= " << ypart << endl;
	      distribute.out++;
	    }
	  
	}
}


void GRID::assignCoherentBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  unsigned int particle;
  const vector<PARTICLE*>& coherent = sog.get_beam()->getCoherentVector(i_slice);
  for (particle = 0; particle < coherent.size(); particle++)
    {
      coherent[particle]->XYposition(xpart, ypart);
      if (coherent[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;
      
      if (particleInGrid(xpart, ypart,i1, i2))
	{
	  
	  sog.assignChargeToNGP( i1, i2,ch);
	  distribute.in++;
	}
      else
	{
	  distribute.out++;
	}
    }
}

void GRID::assignCoherentBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  float h_x,h_y;
  unsigned int particle;
  const vector<PARTICLE*>& coherent = sog.get_beam()->getCoherentVector(i_slice);
  for (particle = 0; particle < coherent.size(); particle++)
    {
      coherent[particle]->XYposition(xpart, ypart);
      if (coherent[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;      
      if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2))
	{
	  sog.assignChargeToCIC(i1, i2, h_x, h_y, ch);
	  distribute.in++;
	}
      else{
	distribute.out++;
      }
    }
}

void GRID::assignTridentBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  unsigned int particle;
  const vector<PARTICLE*>& trident = sog.get_beam()->getTridentVector(i_slice);
  for (particle = 0; particle < trident.size(); particle++)
    {
      trident[particle]->XYposition(xpart, ypart);
      if (trident[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;
      
      if (particleInGrid(xpart, ypart,i1, i2))
	{
	  
	  sog.assignChargeToNGP( i1, i2,ch);
	  distribute.in++;
	}
      else
	{
	  distribute.out++;
	}
    }
}

void GRID::assignTridentBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice, DISTRIBUTE& distribute)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  float h_x,h_y;
  unsigned int particle;
  const vector<PARTICLE*>& trident = sog.get_beam()->getTridentVector(i_slice);
  for (particle = 0; particle < trident.size(); particle++)
    {
      trident[particle]->XYposition(xpart, ypart);
      if (trident[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;      
      if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2))
	{
	  sog.assignChargeToCIC(i1, i2, h_x, h_y, ch);
	  distribute.in++;
	}
      else{
	distribute.out++;
      }
    }
}

/**********************************************/
/*! Routines for the generation of secondaries */
/**********************************************/

/*! This routine gives the number of equivalent photons for an electron
   with energy e above an energy fraction xmin according to spectrum number
   iflag */

// float GRID::requiv(double lns4, float xmin, int iflag)
// {
//   float help,lnxmin;

//   if (xmin>=1.0) return 0.0;
//   switch (iflag)
//     {
//     case 1: 
//       help=log(xmin);
//       return help * help * .00232460830350086;/* 1/137 / pi */
//     case 2:
//       help=log(xmin);
//       return help * help * .00232460830350086;
//     case 3:
//       return log(xmin) * -.00464921660700172 * 0.5*lns4;
//     case 4:
//       return log(xmin) * -.003951834115951462 * 0.5*lns4;
//     case 5:
//       help=lns4;
//       return 2.0*.00232461*help*help;
//     case 6:
//       help= lns4;
//       lnxmin=-log(xmin);
//       return .00232461*lnxmin*(lnxmin+help);
//     case 7:
//       help= lns4;
//       lnxmin=-log(xmin);
//       return .00232461*lnxmin*(lnxmin+help);
//     case 8:
//       help= lns4;
//       lnxmin=-log(xmin);
//       return .00232461*lnxmin*(lnxmin+help);
//     case 9:
//       help= lns4;
//       lnxmin=-log(xmin);
//       return .00232461*lnxmin*(lnxmin+help);
//     case 10:
//       help=log(xmin);
//       return (help * help + (-7.0/6.0-1/3*xmin*xmin*xmin+0.5*xmin*xmin+xmin-log(xmin)) )* .00232460830350086;
//     }
//   return 0.0;
// } /* requiv */


// /*! New version of equiv */

// void GRID::mequiv ( double s4, double lns4, float xmin,float e,int iflag,float *eph,float *q2,float *one_m_x, RNDM& hasard)
// {
//   const float emass2=EMASS*EMASS,eps=1e-5;
//   float help,q2max,q2min,lnx,x,lnxmin,z;
//   switch (iflag)
//     {
//     case 1:
//       *eph=e*powf(xmin,sqrt(1.0-hasard.rndm_equiv()));
//       *q2=0.0;
//       *one_m_x=0.0;
//       return;
//     case 2:
//       x=pow(xmin,sqrt(hasard.rndm_equiv()));
//       help=1.0-x;
//       *eph = e*x;
//       if (hasard.rndm_equiv()>(help*help+1.0)*0.5) *eph=0.0;
//       *q2=0.0;
//       *one_m_x=help;
//       return;
//     case 3:
//       x=pow(xmin,hasard.rndm_equiv());
//       help=1.0-x;
//       *eph=e*x;
//       if (hasard.rndm_equiv()>(help*help+1.0)*0.5) *eph=0.0;
//       *q2=emass2*pow(e*e/emass2,hasard.rndm_equiv());
//       *one_m_x=help;
//       return;
//     case 4:
//       *eph=pow(xmin,hasard.rndm_equiv());
//       help=1.0-*eph;
//       *eph*=e;
//       if (hasard.rndm_equiv()>(help*help+1.0)*0.5) *eph=0.0;
//       *q2=emass2*pow(e*e/emass2,hasard.rndm_equiv());
//       *one_m_x=help;
//       return;
//     case 5:
// 	if(hasard.rndm_equiv()<0.5){
// 	    lnx=-sqrt(hasard.rndm_equiv())*lns4;
// 	    x=exp(lnx);
// 	    q2min=x*x*emass2;
// 	    q2max=emass2;
// 	}
// 	else{
// 	    lnx=-hasard.rndm_equiv()*lns4;
// 	    x=exp(lnx);
// 	    q2min=emass2;
// 	    q2max=s4;
// 	}
// 	if((1.0+(1.0-x)*(1.0-x))*0.5<hasard.rndm_equiv()){
// 	    *eph=0.0;
// 	    *q2=0.0;
// 	}
// 	else{
// 	    *eph=e*x;
// 	    *q2=q2min*pow(q2max/q2min,hasard.rndm_equiv());
// 	}
// /*	if (*q2*(1.0-x)<x*x*emass2) *eph=0.0;*/
// /*	if (*q2<x*x*emass2/(1.0-x)*exp(1.0/(1.0+0.5*x*x/(1.0-x)))) *eph=0.0;*/
// /*	if (hasard.rndm_equiv()>(log(*q2*(1.0-x)/(x*x*emass2))
// 			  -2.0*(1.0-x)/(1.0+(1.0-x)*(1.0-x)))
// 	    /log(*q2*(1.0-x)/(x*x*emass2)))
// 	    *eph=0.0;*/
// 	*one_m_x=1.0-x;
// 	return;
//     case 6:
//         lnxmin=-log(xmin);
// 	if(hasard.rndm_equiv()<lnxmin/(lnxmin+lns4)){
// 	    lnx=-sqrt(hasard.rndm_equiv())*lnxmin;
// 	    x=exp(lnx);
// 	    q2min=x*x*emass2;
// 	    q2max=emass2;
// 	}
// 	else{
// 	    lnx=-hasard.rndm_equiv()*lnxmin;
// 	    x=exp(lnx);
// 	    q2min=emass2;
// 	    q2max=s4;
// 	}
// 	if((1.0+(1.0-x)*(1.0-x))*0.5<hasard.rndm_equiv()){
// 	    *eph=0.0;
// 	    *q2=0.0;
// 	}
// 	else{
// 	    *eph=e*x;
// 	    *q2=q2min*pow(q2max/q2min,hasard.rndm_equiv());
// 	}
// 	if (*q2*(1.0-x)<x*x*emass2) *eph=0.0;
// 	*one_m_x=1.0-x;
// 	return;
//     case 7:
//         lnxmin=-log(xmin);
// 	if(hasard.rndm_equiv()<lnxmin/(lnxmin+lns4)){
// 	    lnx=-sqrt(hasard.rndm_equiv())*lnxmin;
// 	    x=exp(lnx);
// 	    q2min=x*x*emass2;
// 	    q2max=emass2;
// 	}
// 	else{
// 	    lnx=-hasard.rndm_equiv()*lnxmin;
// 	    x=exp(lnx);
// 	    q2min=emass2;
// 	    q2max=s4;
// 	}
// 	if((1.0+(1.0-x)*(1.0-x))*0.5<hasard.rndm_equiv()){
// 	    *eph=0.0;
// 	    *q2=0.0;
// 	}
// 	else{
// 	    *eph=e*x;
// 	    *q2=q2min*pow(q2max/q2min,hasard.rndm_equiv());
// 	}
// 	*one_m_x=1.0-x;
// 	return;
//     case 8:
//         lnxmin=-log(xmin);
// 	if(hasard.rndm_equiv()<lnxmin/(lnxmin+lns4)){
// 	    lnx=-sqrt(hasard.rndm_equiv())*lnxmin;
// 	    x=exp(lnx);
// 	    q2min=x*x*emass2;
// 	    q2max=emass2;
// 	}
// 	else{
// 	    lnx=-hasard.rndm_equiv()*lnxmin;
// 	    x=exp(lnx);
// 	    q2min=emass2;
// 	    q2max=s4;
// 	}
//         *eph=e*x;
//         *q2=q2min*pow(q2max/q2min,hasard.rndm_equiv());
// 	*one_m_x=1.0-x;
// 	return;
//     case 9:
//         lnxmin=-log(xmin);
// 	if(hasard.rndm_equiv()*(lnxmin+lns4)<lnxmin){
// 	    lnx=-sqrt(hasard.rndm_equiv())*lnxmin;
// 	    x=exp(lnx);
// 	    q2min=x*x*emass2;
// 	    q2max=emass2;
// 	}
// 	else{
// 	    lnx=-hasard.rndm_equiv()*lnxmin;
// 	    x=exp(lnx);
// 	    q2min=emass2;
// 	    q2max=s4;
// 	}
//         *eph=e*x;
//         z=q2min*pow(q2max/q2min,hasard.rndm_equiv());
// 	*q2=z-x*x*emass2;
// 	if (2.0*(1.0-x)* *q2+x*x*z<hasard.rndm_equiv()*2.0*z){
// 	  *q2=0.0;
// 	  *eph=0.0;
// 	}
// 	if (*q2*e*e<x*x*emass2*emass2) {
// 	  *q2=0.0;
// 	  *eph=0.0;
// 	}
// 	else{
// 	  if (1.0-x>eps){
// 	    *q2/=(1.0-x);
// 	  }
// 	  else{
// 	    *eph=0.0;
// 	    *q2=0.0;
// 	  }
// 	}
// 	*one_m_x=1.0-x;
// 	return;
//     case 10:
//       double spin= .00232460830350086*(-7.0/6.0-1/3*xmin*xmin*xmin+0.5*xmin*xmin+xmin-log(xmin)); //total spin part of spectrum
//       double eps = 0.00001; // how close it should be to the random number
//       if (spin/(spin+requiv(lns4,xmin,iflag))>hasard.rndm_equiv() )
// 	{
// 	  double yk=hasard.rndm_equiv()*spin/.00232460830350086;
// 	  double y1=yk;
// 	  double x0;
// 	  double x1=exp(-y1);;
// 	  help=0;
// 	  while(fabs( (-7.0/6.0-1.0/3.0*x1*x1*x1+0.5*x1*x1+x1-log(x1) -yk))>eps )
// 	    {	  
// 	      x0=exp(-y1);
// 	      y1= yk- (-7.0/6.0-1.0/3.0*x0*x0*x0+0.5*x0*x0+x0-log(x0))+help;
// 	      x1=exp(-y1);
// 	      help=y1;
// 	      //	      printf("yk= %e  y1= %e x0= %e  x1= %e \n" ,yk,y1,x0,x1);
// 	    }
// 	  *eph = e*x1;
// 	  *q2=0.0;
// 	  *one_m_x=0;
// 	}
//       else
// 	{
// 	  x=pow(xmin,sqrt(hasard.rndm_equiv()));
// 	  help=1.0-x;
// 	  *eph = e*x;
// 	  if (hasard.rndm_equiv()>(help*help+1.0)*0.5) *eph=0.0;
// 	  *q2=0.0;
// 	  *one_m_x=0;
// 	}
//       return;
//     }
// } /* mequiv */


/*! This routine calculates the luminosity from the collision of the two
slices i1 and i2 of grid. The dimensions of the grids have to be the same. */

void GRID::step_lumi(float min_z, PAIR_BEAM& secondaries,int time_counter, SWITCHES& switches)
{
  float sum=0.0;
  int i1,i2,j;
  list<BEAM_PARTICLE_POINTER>::const_iterator pointer1,pointer2;

  const PHI_FLOAT* rho1 = slice_of_beam_[0].get_rho();
  const PHI_FLOAT* rho2 = slice_of_beam_[1].get_rho();
  for (i1=0;i1 < n_cell_x_;i1++)
    {
      for (i2=0;i2 < n_cell_y_;i2++)
	{
	  j=i1*n_cell_y_+i2;
	  sum += rho1[j]*rho2[j];
	  for (pointer1=part_pointer1_[j].begin(); pointer1 != part_pointer1_[j].end(); pointer1++)
	    {
	      for( pointer2 = part_pointer2_[j].begin(); pointer2 != part_pointer2_[j].end(); pointer2++)
		{
		  collide_ee(i1, i2, min_z, *pointer1, *pointer2,switches, secondaries,time_counter);
		}
	    }
	}
    }

  results_.add_lumi_fine( sum*1e18/(mesh_.get_delta_x()*mesh_.get_delta_y()*timestep_) );	    
}

void GRID::collide_ee(int cellx, int celly,float min_z, const BEAM_PARTICLE_POINTER& pointer1, const BEAM_PARTICLE_POINTER& pointer2, SWITCHES& switches,PAIR_BEAM& secondaries, int time_counter)
{
  //  JET_FLOAT bhabhan;
  float help,ecm,e1,e2;
  //float ecmratio;
  //int nbphot, numero_bhabha;
  float weight = pointer1.weight() * pointer2.weight();
  int j1=0,j2=1;
  int bmt_precession = switches.get_bmt_precession();

  e1 =  fabs(pointer1.energy());
  e2 =  fabs(pointer2.energy());
  if (switches.get_do_espread())
    {
      e1 = spread_energy(e1, switches.get_which_espread1(), switches.get_espread1(), *hasard_);
      e2 = spread_energy(e2, switches.get_which_espread2(), switches.get_espread2(), *hasard_);
    }

  if (switches.get_do_isr()) isr2(e1,e2,&e1,&e2, *hasard_);
  
  ecm=sqrt(4.0*e1*e2);

  float energie1 = pointer1.energy();
  float energie2 = pointer2.energy();
  if (switches.get_do_lumi()&1)
    {
      float p1Vx, p1Vy, p2Vx, p2Vy;
      pointer1.velocities(p1Vx, p1Vy);
      pointer2.velocities(p2Vx, p2Vy);
      if ( bmt_precession ) 
      //      if (SPIN)
	{
	  lumi_heap_ee_.lumi_store_ee(mesh_, cellx, celly,min_z, energie1, p1Vx, p1Vy, energie2, p2Vx, p2Vy, weight, pointer1.getSpin(),pointer2.getSpin() , time_counter);
	}
      else lumi_heap_ee_.lumi_store_ee(mesh_, cellx, celly,min_z, energie1, p1Vx, p1Vy, energie2, p2Vx, p2Vy, weight,time_counter);
    }

  if (energie1 < 0.0) j1=1;
  if (energie2 < 0.0) j2=0;
  results_.add_lumi(j1, j2, weight);

  if (energie1*energie2 > 0.0)
    {
      if (switches.get_do_cross())
	{
	  cross_->cross_add(e1,e2,weight);     
	}
      results_.add_lumi_ee(weight);

      if ( bmt_precession ) 
	//     if (SPIN)
	{
	  float spin1= pointer1.getSpin()(2);
	  float spin2= pointer2.getSpin()(2);
	  help=0.5*(1.0+spin1*spin2);
	  results_.add_lumi_pp(weight*help);
	}

      if (ecm > switches.get_ecm_min()) results_.add_lumi_ee_high(weight);
      help= weight*ecm;
      results_.add_lumi_ecm(help);
      results_.add_lumi_ecm2(help*ecm);
    
    }
  
  // ici les bhabha de Cecile

  if (switches.get_do_bhabhas())
    {
      float part1Vx, part1Vy, part2Vx, part2Vy;
      pointer1.velocities(part1Vx, part1Vy);
      pointer2.velocities(part2Vx, part2Vy);
      bhabhas_.make_bhabha(secondaries, part1Vx, part1Vy, part2Vx, part2Vy, e1, e2, ecm, weight, mesh_, cellx, celly,min_z, switches, *hasard_);
    }

  if (switches.get_do_jets())
    {
      minijets_->mkjll_(secondaries.get_pair_parameters(),pointer1.energy(),pointer2.energy(), weight, switches, *hasard_);

    }
  if (switches.get_do_prod()==1)
    {

      cerr << " GRID::collide_ee : do_prod non programme " << endl;
      exit(0);
    }  
}


void GRID::isr2(float e1,float e2,float *e1p,float *e2p, RNDM& hasard)
{
    double c=2.0*ALPHA_EM/PI,emass2=EMASS*EMASS;
    double x,s,beta,corr,tmp;

    s=4.0*e1*e2;
    beta=c*(log(s/emass2)-1.0);
    do{
	x=pow(1.0-hasard.rndm(),2.0/beta);
	tmp=0.5*beta*pow(x,0.5*beta-1.0)*(1.0+0.375*beta);
	corr=(tmp-0.25*beta*(2.0-x))/tmp;
    }
    while(hasard.rndm()>corr);
    *e1p=e1*(1.0-x);
    do{
	x=pow(1.0-hasard.rndm(),2.0/beta);
	tmp=0.5*beta*pow(x,0.5*beta-1.0)*(1.0+0.375*beta);
	corr=(tmp-0.25*beta*(2.0-x))/tmp;
    }
    while(hasard.rndm()>corr);
    *e2p=e2*(1.0-x);
}


/*! This routine moves the particles one timestep. */

void GRID::move_particles(const vector<GENERAL_GRID*>& grids,  int i_beam, int i_slice,int interpolation,int do_beamstrahlung,int do_trident, int sokolov, float emin,int do_prod, int extra_grids, float charge_sign, int bmt_rotate)
{
  float vx,vy,xpos,ypos;

  // ultrarelivistic particles 
  //      Bfield *= 1.0/Cvelocity;
  // champ electrique : GV/nm
  // pour le champ magnetique, on envoie le produit c.B
  // aussi en GV/nm
  //       Efield = EBfield;
  //       // ultrarelivistic particles  
  //       //      Bfield *= 1.0/Cvelocity;
  //       // champ electrique : GV/nm
  //       // pour le champ magnetique, on envoie le produit c.B
  //       // aussi en GV/nm

  unsigned int k,i;
  const PHI_FLOAT *phi;
  //  BEAM* beam;
  if (i_beam==1) 
    {
      //     beam = slice_of_beam_[0].get_beam();
      //      phi=phi2_;
      phi= champ_.get_phi(2);
    }
  else 
    {
      //     beam = slice_of_beam_[1].get_beam();
      //     phi=phi1_;
      phi= champ_.get_phi(1);
    }
  const vector<PARTICLE*>& lesParticles = slice_of_beam_[i_beam-1].get_beam()->getParticleVector(i_slice);
  switch (interpolation)
    {
    case 1:
      cerr << " GRID::move_particles interpolation = 1, a traiter " << endl;
      exit(0);
      break;
    case 2:
      {
	TRIVECTOR EBfield;
	TRIDVECTOR Efield, Bfield;
	float dzOnRadius,oldener;
	////////////////////////////////////////////////
	// en cas de do_beamstrahlung, voir l'implementation du cas  (do_prod>1)&&(i_beam==1) 
	///////////////////////////////////////////////
	PARTICLE* particuleCourante;
	for (k = 0; k < lesParticles.size(); k++)
	  {
	    particuleCourante = lesParticles[k];
	    oldener=particuleCourante->getEnergy();
	    EBfield = EBfieldOnParticle( particuleCourante, grids, phi, i_beam, extra_grids);


	    dzOnRadius = particuleCourante->advanceDueToEBfield(EBfield, step_,slice_of_beam_[i_beam-1].get_scal_step());


	    Efield = EBfield;
	    Bfield = TRIDVECTOR(EBfield(1), -EBfield(0), 0.0);      
	    if (bmt_rotate)
	      {
		particuleCourante->rotateBMT(Efield, Bfield, charge_sign, step_);
		if (do_beamstrahlung)
		  { 
		    vector<float> photonEnergies;
		    results_.updateUpsmax(particuleCourante->getUps());
		    if (sokolov) 
		      {
			particuleCourante->beamstrahlungSokolov(Efield, Bfield, dzOnRadius, emin, charge_sign,photonEnergies, *hasard_);
		      }
		    else 
		      {
			particuleCourante->beamstrahlung(Efield, Bfield, dzOnRadius,emin,photonEnergies, *hasard_);
		      }
		    registerPhotons(photonEnergies, *(particuleCourante), i_beam, i_slice);    
		  }
	      }
	    else
	      {
		if(do_beamstrahlung)
		  {
		    vector<float> photonEnergies;
		    results_.updateUpsmax(particuleCourante->getUps());
		    // 		    particleBeamstrahlung(particuleCourante,Efield, Bfield, dzOnRadius, emin, photonEnergies);
		    particuleCourante->beamstrahlung(Efield, Bfield, dzOnRadius, emin, photonEnergies, *hasard_);
		    registerPhotons(photonEnergies, *(particuleCourante), i_beam, i_slice);    
		  }
	      }
	    if(do_trident)
	      {
		vector<float> electrons,positrons,virt;
		if(do_beamstrahlung)
		  {
		    particuleCourante->setUps(particuleCourante->getUps()*particuleCourante->getEnergy()/oldener);
		  }
		float oldnyr=particuleCourante->getEnergy();
		particuleCourante->createTridents(step_*slice_of_beam_[i_beam-1].get_scal_step(),&electrons,&positrons,&virt,*hasard_);
		for(i=0;i<electrons.size();i++)
		  {
		    particuleCourante->XYposition(xpos,ypos);
		    particuleCourante->velocities(vx,vy);
		    store_trident_particle(*slice_of_beam_[i_beam-1].get_beam(),electrons[i],vx,vy,xpos,ypos,i_slice);
		    store_trident_particle(*slice_of_beam_[i_beam-1].get_beam(),positrons[i],vx,vy,xpos,ypos,i_slice);
		  }
	      }
	  }
	break;
      }
    default:
      cerr << " GRID::move_particles : unknown value of interpolation type :  " << interpolation << endl;
      exit(0);
    }
}


void GRID::move_pairs(const vector<GENERAL_GRID*>& grids, PAIR_BEAM& pairs, int i_slice, double d_eps_1, double d_eps_2, int extra_grids, float charge_sign_0, RNDM& hasard)
{
  float stepLocal,d;
  int n_pair_steps;
  double mass=pairs.get_pair_parameters().get_mass();

  vector<PAIR_PARTICLE>& les_paires = pairs.get_pairs(i_slice);
  //  list<PAIR_PARTICLE>::iterator itr;
  unsigned int k;
  //  for (itr = les_paires.begin(); itr !=  les_paires.end(); itr++)
  for (k = 0; k <  les_paires.size(); k++)
    {      
      //      d=sqrt((rho_sum_1_*d_eps_1 + rho_sum_2_*d_eps_2)/fabs(itr->energy()));
      d=sqrt( (rho_sum_1_*d_eps_1 + rho_sum_2_*d_eps_2)/les_paires[k].unsigned_energy() );
      n_pair_steps=(int)d+1;
      pairs.add_pair_steps(n_pair_steps);
      stepLocal=step_/(float)n_pair_steps;
      step_pair_1(grids,les_paires[k],mass,stepLocal,n_pair_steps, extra_grids, charge_sign_0, hasard);
    }
}

void GRID::move_pairs_tertphot(const vector<GENERAL_GRID*>& grids, PAIR_BEAM& pairs, int i_slice, double d_eps_1, double d_eps_2, int extra_grids, float charge_sign_0, RNDM& hasard)
{
  float stepLocal,d;
  int n_pair_steps;
  double mass=pairs.get_pair_parameters().get_mass();

  vector<PAIR_PARTICLE>& les_paires = pairs.get_pairs(i_slice);
  //  list<PAIR_PARTICLE>::iterator itr;
  unsigned int k;
  //  for (itr = les_paires.begin(); itr !=  les_paires.end(); itr++)
  for (k = 0; k <  les_paires.size(); k++)
    {      
      //      d=sqrt((rho_sum_1_*d_eps_1 + rho_sum_2_*d_eps_2)/fabs(itr->energy()));
      d=sqrt( (rho_sum_1_*d_eps_1 + rho_sum_2_*d_eps_2)/les_paires[k].unsigned_energy() );
      n_pair_steps=(int)d+1;
      pairs.add_pair_steps(n_pair_steps);
      stepLocal=step_/(float)n_pair_steps;
      step_pair_1_tertphot(grids,les_paires[k],mass,stepLocal,n_pair_steps, extra_grids, charge_sign_0, hasard);
    }
}

void GRID::deltaVelocityFromFieldCIC(float xpart,float ypart, float energy, PHI_FLOAT *phi, float pasDeTemps, float& deltavx, float& deltavy)
{
  PHI_FLOAT phi1_x,phi2_x,phi3_x,phi1_y,phi2_y,phi3_y,h_x,h_y;

  interpolePotential(xpart, ypart, h_x, h_y, phi1_x, phi2_x, phi3_x, phi1_y, phi2_y, phi3_y, phi);

  deltavx = (h_x*(phi1_y-phi2_y)+(1.0-h_x)*(phi2_y-phi3_y))*delta_x_inverse_/energy;
  deltavy = (h_y*(phi1_x-phi2_x)+(1.0-h_y)*(phi2_x-phi3_x))*delta_y_inverse_/energy;
  //  cout << " :deltaVelocityFromFieldCIC Ex/g= " << deltavx << endl;
      deltavx = -deltavx*pasDeTemps;
      deltavy = -deltavy*pasDeTemps;
}


TRIVECTOR GRID::electric_field_out_of_main_grid(const vector<GENERAL_GRID*>& grids,int beam,PHI_FLOAT x,PHI_FLOAT y, int extra_grids) const
{
  int k;
  PHI_FLOAT tmp,dx,dy;
  const PHI_FLOAT *phi;
  // float aux;
  float ax, ay;
  TRIVECTOR EBfield;
  ax = 0.0;
  ay = 0.0;
  for(k = 0;k <= extra_grids;k++)
    {
      if (grids[k]->coordinatesInGridRange(x,y))
	{
	  if (beam == 1 ) phi = grids[k]->get_phi2();
	  else phi = grids[k]->get_phi1();

	  // electric field in V/m
	  EBfield = grids[k]->ElectricFieldCIC(x, y, phi);
	  return EBfield;
	}
    }
  if (extra_grids<2) 
    {
      EBfield.setComponents(0.0, 0.0, 0.0);
      return EBfield;
    }
  // particle is not in largest grid 
  if (beam==1)
    {
      dx=(x-grids[0]->get_rho_x_2());
      dy=(y-grids[0]->get_rho_y_2());
      tmp=grids[0]->get_rho_factor()*grids[0]->get_rho_sum_2()/(dx*dx+dy*dy);
      ax=dx*tmp;
      ay=dy*tmp;
    }
  else
    {
      dx=(x-grids[0]->get_rho_x_1());
      dy=(y-grids[0]->get_rho_y_1());
      tmp=grids[0]->get_rho_factor()*grids[0]->get_rho_sum_1()/(dx*dx+dy*dy);
      ax=dx*tmp;
      ay=dy*tmp;
    }

      EBfield.setComponents(0.5*ax, 0.5*ay, 0.0);
      return EBfield;
}


void GRID::registerPhotons(const vector<float>& photonEnergies, PARTICLE& particle, int i_beam, int i_slice)
{
  unsigned int k;
  BEAM* beam = slice_of_beam_[i_beam-1].get_beam();
  //  int  i_beam = beam.label();
  //  cout << " GRID::createPhotons : nb photons recuperes " << photonEnergies.size() << endl;
  for (k=0; k < photonEnergies.size() ;k++)
    {
      
      const PHOTON& foton = beam->new_photon(photonEnergies[k],particle, i_slice);

      //     energy -= photonEnergies[k];
      //      const PHOTON& foton = beam->new_photon(photonEnergies[k],particle, i_slice);
//       if (photon_file_ != NULL) 
// 	{
// 	  if (i_beam != 1) 
// 	    {
// 	      PHOTON aux(foton);
// 	      float energyAux = aux.energy();
// 	      aux.setEnergy(-energyAux);
// 	      photon_file_->save_object_on_persistent_file(&aux);
// 	    }
// 	  else photon_file_->save_object_on_persistent_file(&foton);	  
// 	}
    }
}


void GRID::beamstrahlungSingleCoherentParticle(PARTICLE* particle, TRIVECTOR EBfield, float dzOnRadius, const vector<GENERAL_GRID*>& grids, int i_beam,  int i_slice, const PHI_FLOAT *phi, float pasDeTemps, float emin,int do_prod, int extra_grids, float charge_sign)
{	
  float upsilon = particle->getUps();

  float energy = 0.0;
  float initialEnergy = fabs(particle->energy());
  results_.updateUpsmax(upsilon);

  vector<float> photonEnergies;
  TRIDVECTOR ev1, ev2, ev3; 
  energy = PHYSTOOLS::synrad_no_spin_flip(upsilon, initialEnergy, dzOnRadius,photonEnergies, *hasard_);
  registerPhotons(photonEnergies, *(particle), i_beam, i_slice);
  if (energy < emin)
    {
      cout << "PARTICLE:: e_low2: " << energy << endl;
      energy = emin;
    }
  if ((do_prod>1)&&( i_beam ==1))
    {
      cerr << " GRID:: do_prod a implanter "  << endl;
      exit(0);
    }
  if (particle->energy() < 0.0)
    {
      particle->setEnergy( -energy);
    }
  else
    {
      particle->setEnergy(energy);
    }	  
    
}

// // transverse Lorents force in eV/m, with respect to a particle trajectory
// // data E and B, in eV/m
// TRIDVECTOR GRID::transverse_Lorentz_force(PARTICLE& part, TRIDVECTOR E, TRIDVECTOR B)
// {
//   TRIDVECTOR  Ftrans;
//   double xcomp, ycomp, zcomp;
//   float vx,vy;
//   float dvx, dvy;
//   // the longitudinal betaz is assumed to be equal to 1
//   part.velocities(vx,vy);
//   dvx = (double)vx;
//   dvy = (double)vy;
//   // projection of E on the momentum
//   double ElongProjection = dvx*E(0) + dvy*E(1) + E(2);
//  xcomp = E(0) - dvx*ElongProjection + dvy*B(2) - B(1);
//  ycomp = E(1) - dvy*ElongProjection + B(0) - dvx*B(2);
//  zcomp = E(3) - ElongProjection     + dvx*B(1) - dvy*B(0);
//  Ftrans.setComponents(xcomp, ycomp, zcomp);
//  return Ftrans;
// }
void GRID::move_coherent_particles(const vector<GENERAL_GRID*>& grids,  int i_beam, int i_slice, int interpolation, int do_beamstrahlung, float emin, int do_prod, int extra_grids,float charge_sign)
{
  unsigned int k; 
  //  float initialEnergyForLoss;
  const PHI_FLOAT *phi;

  //  BEAM* beam;

  if (i_beam==1) 
    {
      //     beam = slice_of_beam_[0].get_beam();
      //      phi=phi2_;  
      phi = champ_.get_phi(2);
    }
  else
    { 
      //     beam = slice_of_beam_[1].get_beam();
      //     phi=phi1_;
      phi = champ_.get_phi(1);
    }
  const vector<PARTICLE*>& coherent = slice_of_beam_[i_beam-1].get_beam()->getCoherentVector(i_slice);

  switch (interpolation)
    {
    case 1:
      cerr << " GRID::move_coherent_particles interpolation = 1, a traiter " << endl;
      exit(0);
      break;
    case 2:
      {
	for (k = 0; k < coherent.size(); k++)
	  {

	    TRIVECTOR EBfield = EBfieldOnParticle(coherent[k], grids, phi, i_beam, extra_grids);
	    float dzOnRadius = coherent[k]->advanceDueToEBfield(EBfield, step_,  slice_of_beam_[i_beam-1].get_scal_step());
	    // (no spin rotation for coherent particles, nor sokolov ternov spin flip)
	    if(do_beamstrahlung)
	      {
		beamstrahlungSingleCoherentParticle(coherent[k], EBfield, dzOnRadius,grids, i_beam, i_slice, phi, step_, emin, do_prod, extra_grids, charge_sign);	
	      }  
	  }      
	break;
      }
    default:
      cerr << " GRID::move_coherent_particles : unknown value of interpolation type :  " << interpolation << endl;
      exit(0);
    }
}

void GRID::move_trident_particles(const vector<GENERAL_GRID*>& grids,  int i_beam, int i_slice, int interpolation, int do_beamstrahlung, float emin, int do_prod, int extra_grids,float charge_sign)
{
  unsigned int k; 
  const PHI_FLOAT *phi;
  if (i_beam==1) 
    {
      phi = champ_.get_phi(2);
    }
  else
    { 
      phi = champ_.get_phi(1);
    }
  const vector<PARTICLE*>& trident = slice_of_beam_[i_beam-1].get_beam()->getTridentVector(i_slice);

  switch (interpolation)
    {
    case 1:
      cerr << " GRID::move_trident_particles interpolation = 1, a traiter " << endl;
      exit(0);
      break;
    case 2:
      {
	for (k = 0; k < trident.size(); k++)
	  {
	    TRIVECTOR EBfield = EBfieldOnParticle(trident[k], grids, phi, i_beam, extra_grids);
	    float dzOnRadius = trident[k]->advanceDueToEBfield(EBfield, step_, slice_of_beam_[i_beam-1].get_scal_step());
	    /* This routine includer transverse momentum rotation due to non-transversality if the fields */
	    //	    float dzOnRadius = trident[k]->advanceTridentDueToEBfield(EBfield, step_, slice_of_beam_[i_beam-1].get_scal_step());
	    // (no spin rotation for coherent particles, nor sokolov ternov spin flip)
	    if(do_beamstrahlung)
	      {
		beamstrahlungSingleCoherentParticle(trident[k], EBfield, dzOnRadius,grids, i_beam, i_slice, phi, step_, emin, do_prod, extra_grids, charge_sign);	
	      }  
	  }      
	break;
      }
    default:
      cerr << " GRID::move_trident_particles : unknown value of interpolation type :  " << interpolation << endl;
      exit(0);
    }
}


/* Distributes the photons for background calculation */

void GRID::distribute_photons(int slice_1,int slice_2, float photon_ratio, RNDM& hasard)
{
  //int k;
  float ratio,ratio_i_1,ratio_i_2;
  const float eps=1e-5;

  //  cout << " GRID::distribute_photons " << endl;

  ratio = photon_ratio;
  clear_photon_pointer();
 
  if (ratio<eps) return;

  float deltax = mesh_.get_delta_x();
  float deltay = mesh_.get_delta_y();

  ratio_i_1=1e9/ratio*slice_of_beam_[0].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);
  ratio_i_2=1e9/ratio*slice_of_beam_[1].get_nb_part_per_macro()/sqrt(deltax*deltay*timestep_);

  distributePhotonInBeam(grid_photon1_, 1, slice_1, ratio, ratio_i_1, hasard);
  distributePhotonInBeam(grid_photon2_, 2, slice_2, ratio, ratio_i_2, hasard);

}

//void GRID::distributePhotonInBeam(   vector< list<BEAM_PARTICLE_POINTER> >& grid_photon, BEAM* beam, int slice,float ratio, float ratio_i, RNDM& hasard)
void GRID::distributePhotonInBeam(   vector< list<BEAM_PHOTON_POINTER> >& grid_photon,int i_beam, int slice,float ratio, float ratio_i, RNDM& hasard)
{
  float xphot, yphot;
  int i1, i2, j;
   

  float aux;

  const vector<PHOTON>& lesPhotons =  slice_of_beam_[i_beam-1].get_beam()->getPhotonVector(slice);
  unsigned int k;
  for(k=0; k < lesPhotons.size(); k++)
    {
      aux = hasard.rndm();
      if (aux<ratio)
	{
	  lesPhotons[k].XYposition(xphot, yphot);

	  if (particleInGrid(xphot, yphot, i1, i2))
	    {
	      j = i1*n_cell_y_ + i2;
	      grid_photon[j].push_back(BEAM_PHOTON_POINTER(&lesPhotons[k],ratio_i ));

		
	    }
	}
    }
}

/*! This routine calculates the electron-photon, positron-photon and
   photon-photon luminosities */

void GRID::photon_lumi(float min_z, SWITCHES& switches, PAIR_BEAM& secondaries, PAIR_BEAM& muons, RNDM& hasard)
{

  int i1,i2,j,n_x,n_y;
  //  list<BEAM_PARTICLE_POINTER>::iterator photon_pointer, photon_pointer1,photon_pointer2;
  list<BEAM_PHOTON_POINTER>::const_iterator photon_pointer; 
  list<BEAM_PHOTON_POINTER>::const_iterator photon_pointer1,photon_pointer2;
  list<BEAM_PARTICLE_POINTER>::const_iterator particle_pointer;
  n_x = n_cell_x_;
  n_y = n_cell_y_;

  const PAIR_PARAMETER& pair_parameter = secondaries.get_pair_parameters();

  for (i1=0;i1<n_x;i1++)
    {
      for (i2=0;i2<n_y;i2++)
	{
	  j=i1*n_y+i2;
	  for (photon_pointer=grid_photon1_[j].begin();photon_pointer != grid_photon1_[j].end(); photon_pointer++)
	    {
	      for (particle_pointer = part_pointer2_[j].begin(); particle_pointer != part_pointer2_[j].end();  particle_pointer++)
		{
		  collide_ge(i1, i2, min_z,*photon_pointer,*particle_pointer, secondaries, switches, pair_parameter, hasard);
		}
	    }
	  
	  for (photon_pointer=grid_photon2_[j].begin(); photon_pointer!= grid_photon2_[j].end(); photon_pointer++)
	    {
	      for (particle_pointer= part_pointer1_[j].begin(); particle_pointer != part_pointer1_[j].end(); particle_pointer++)
		{
		  collide_eg( i1, i2, min_z,*particle_pointer, *photon_pointer, secondaries, switches,pair_parameter,hasard);		  
		}
	    }
	  
	  
	  double returnCrossSectionToAdd[3];
	  double tempor[3];
	  int k;
	  for (k=0; k<3;k++) {
	    tempor[k]=0.0;
	    returnCrossSectionToAdd[k] = 0.0;
	  }
	  
	  for ( photon_pointer1 = grid_photon1_[j].begin(); photon_pointer1 != grid_photon1_[j].end(); photon_pointer1++)
	    {
	      for( photon_pointer2 = grid_photon2_[j].begin(); photon_pointer2 != grid_photon2_[j].end(); photon_pointer2++)
		{
		  collide_gg(i1, i2, min_z,*photon_pointer1, *photon_pointer2, secondaries, muons, switches, hasard, returnCrossSectionToAdd);
		  for (k=0; k<3;k++) tempor[k] += returnCrossSectionToAdd[k];
		}
	      
	    }	  
	  results_. cumulate_hadrons_gg(tempor);
	}
    }
}


void GRID::collide_ge(int cellx, int celly,float min_z, const BEAM_PHOTON_POINTER& photon_pointer, const BEAM_PARTICLE_POINTER& particle_pointer, PAIR_BEAM& secondaries, SWITCHES& switches, const PAIR_PARAMETER& pair_parameter, RNDM& hasard)
{
   float photonEnergy = photon_pointer.energy();
  if (photonEnergy <= 0.0) return; 

  float particleEnergy = particle_pointer.energy(); 
  int j2=1;
  if (particleEnergy <0.0 ) j2=0;

   float weight = photon_pointer.weight()*particle_pointer.weight(); 

  results_.add_lumi(2, j2, weight);

  if (switches.get_do_lumi()&2)
    {
      lumi_heap_ge_.lumi_store(mesh_, cellx, celly,min_z, photonEnergy,particleEnergy, weight);
    }

  results_.add_lumi_ge(weight);

  if (switches.get_do_jets())
    minijets_->mkjbh1_(pair_parameter, photonEnergy,particleEnergy,weight, switches, hasard);
  // ancien appel de collid_ge_2
  collide_compton(0, cellx, celly,min_z, photon_pointer, particle_pointer, secondaries, switches,hasard, 1);		    
}

void GRID::collide_eg(int cellx, int celly,float min_z, const BEAM_PARTICLE_POINTER& particle_pointer, const BEAM_PHOTON_POINTER& photon_pointer, PAIR_BEAM& secondaries, SWITCHES& switches, const PAIR_PARAMETER& pair_parameter, RNDM& hasard)
{
  float photonEnergy = photon_pointer.energy();
  if (photonEnergy <= 0.0) return;

  float particleEnergy = particle_pointer.energy();  
  int j1=0;
  if (particleEnergy < 0.0) j1=1;
 
  float weight = photon_pointer.weight()*particle_pointer.weight();

  results_.add_lumi(j1, 2, weight);

  if (switches.get_do_lumi()&2)
    {
      lumi_heap_eg_.lumi_store(mesh_, cellx, celly,min_z,particleEnergy,photonEnergy, weight);
    }

  results_.add_lumi_eg(weight);

  if (switches.get_do_jets())
    minijets_->mkjbh2_(pair_parameter, particleEnergy,photonEnergy,weight, switches, hasard);
  // ancien appel de collide_eg_2 
  collide_compton(0,cellx, celly, min_z, photon_pointer, particle_pointer, secondaries, switches, hasard, 2);		    
}

void GRID::collide_gg(int cellx, int celly, float min_z, const BEAM_PHOTON_POINTER& photon_pointer1, const BEAM_PHOTON_POINTER& photon_pointer2, PAIR_BEAM& secondaries ,PAIR_BEAM& muons, SWITCHES& switches, RNDM& hasard, double *returnCrossSectionToAdd)
{
  JET_FLOAT ecm;
  //double beta_x,beta_y;
  // float particleVx1,particleVy1,particleVx2,particleVy2 ;
  float photonEnergy1 = photon_pointer1.energy();
  float photonEnergy2 = photon_pointer2.energy();
  float weight  = photon_pointer1.weight()*photon_pointer2.weight();

  ecm=sqrt(4.0*photonEnergy1*photonEnergy2);

  results_.add_lumi(2, 2, weight);
  if (ecm<=0.0) {
    return;
  }

  if (switches.get_do_lumi()&4){

    lumi_heap_gg_.lumi_store(mesh_, cellx, celly,min_z,photonEnergy1, photonEnergy2, weight);
  }
  results_.add_lumi_gg(weight);
  if (ecm*ecm > switches.get_gg_smin())
    {
      results_.add_lumi_gg_high(weight);
    }
  if (switches.get_do_jets())
    minijets_->mkjbw_(photonEnergy1,photonEnergy2,weight, switches, hasard );

  collide_gg_XX(0,cellx, celly, min_z,photon_pointer1, photon_pointer2, switches, secondaries, muons, hasard, returnCrossSectionToAdd);
}

//void GRID::collide_gg_XX(int index_of_process, int cellx, int celly,float min_z, const PARTICLE_POINTER& photon1, const PARTICLE_POINTER& photon2, SWITCHES& switches, PAIR_BEAM& secondaries, RNDM& hasard, double returnCrossSectionToAdd[3])
void GRID::collide_gg_XX(int index_of_process, int cellx, int celly,float min_z, const PHOTON_POINTER& photon1, const PHOTON_POINTER& photon2, SWITCHES& switches, PAIR_BEAM& secondaries, PAIR_BEAM& muons, RNDM& hasard, double *returnCrossSectionToAdd)
{
  double beta_x,beta_y;
  float particle1Vx,particle1Vy,particle2Vx,particle2Vy ;
  float phot1_q2, phot1_eorg;
  float phot2_q2, phot2_eorg;
  float weight;
  float energy1 = photon1.energy();
  float energy2 = photon2.energy();
  if (energy1*energy2 <= 0.0) return;

  photon1.velocities(particle1Vx,particle1Vy);
  photon2.velocities(particle2Vx,particle2Vy);
  beta_x=(particle1Vx*energy1 + particle2Vx*energy2)/(energy1+energy2);
  beta_y=(particle1Vy*energy1 + particle2Vy*energy2)/(energy1+energy2);

  phot1_q2 =  photon1.q2();
  phot2_q2 =  photon2.q2();
  phot1_eorg = photon1.eorg();
  phot2_eorg = photon2.eorg();

  weight  = photon1.weight()*photon2.weight();

  if (switches.get_do_pairs())
    {
      secondaries. make_pair_bw(mesh_, cellx, celly,min_z,index_of_process, energy1,phot1_q2,phot1_eorg,energy2,phot2_q2,phot2_eorg, weight,beta_x,beta_y, switches,hasard);
    }
  if (switches.get_do_muons())
    {
      muons.make_muon(mesh_, cellx, celly,min_z, index_of_process, energy1,phot1_q2,energy2,phot2_q2,weight,beta_x,beta_y, switches, hasard);
    }
  if (switches.get_do_hadrons()) make_hadrons_gg2(min_z, energy1,phot1_q2,energy2,phot2_q2,weight, switches, returnCrossSectionToAdd, hasard);
}

/**********************************************************************/
/* Routines for the creation and storage of hadrons                   */
/**********************************************************************/

void GRID::make_hadrons_gg2(float min_z, float energy1,float q2_1,float energy2,float q2_2,float lumi, SWITCHES& switches, double *returnCrossSectionToAdd, RNDM& hasard)
{
  double s,h;
  double cross_section=0.0;
  int num,k;
  double eps=0.0808,mu=0.4525;

  for (k=0; k<3; k++) returnCrossSectionToAdd[k] = 0.0;

  s=energy1*energy2;
  h=max(1.0,pow(s/100.0,0.43));
  if ((q2_1>h)||(q2_2>h)) {
    return;
  }
  s*=4.0;
  if (s<4.0) return;
  switch (switches.get_do_hadrons())
    {
    case 1:
      cross_section=(211.0*pow(s,eps)+297.0*pow(s,-mu))*1e-37*lumi;
      break;
    case 2:
      cross_section=lumi*200e-37*(1.0+6.3e-3*pow(log(s),2.1)+1.96*pow(s,-0.37));
      break;
    case 3:
      cross_section=(211.0*pow(s,eps)+215.0*pow(s,-mu))*1e-37*lumi;
      break;
    case 4:
      cross_section=(-0.99244+0.0989203*log(s+22865.9))*1e-34*lumi;
      break;
    }
  if (hadron_file_ != NULL)
    {
      h=cross_section*switches.get_hadron_ratio();
      num=(int)floor(h);
      h-=num;
      if(h>hasard.rndm()) num++;
      saveOnHadronFile(num, energy1, energy2,min_z, hasard );
  }
  returnCrossSectionToAdd[0] = cross_section;
  if (s<25.0) return;
  returnCrossSectionToAdd[1] = cross_section;
  if (s<100.0) return;
  returnCrossSectionToAdd[2] = cross_section;
}

/*! This routine calculates the background using the virtual photons with the
   apropriate impact parameter. In the moment it just does pair creation. */

void GRID::photon_lumi_2(float min_z,SWITCHES& switches, PAIR_BEAM& secondaries, PAIR_BEAM& muons, RNDM& hasard)
{
  int i1,i2,j,n_x,n_y;
  list<BEAM_PHOTON_POINTER>::iterator photon_pointer;
  //double beta_x,beta_y,e1,e2;
  n_x = n_cell_x_;
  n_y = n_cell_y_;
  for (i1=0;i1<n_x;i1++)
    {
      for (i2=0;i2<n_y;i2++)
	{
	  j=i1*n_y+i2;

	  double returnCrossSectionToAdd[3];
	  double tempor[3];
	  int k;
	  for (k=0; k<3;k++) {
	    tempor[k]=0.0;
	    returnCrossSectionToAdd[k]=0.0;
	  }


	  list<EXTRA_PHOTON_POINTER>::const_iterator extr_phot;
	  //float extr_phot_e, extr_phot_weight, extr_phot_q2, extr_phot_eorg;
	  for (photon_pointer=grid_photon1_[j].begin(); photon_pointer != grid_photon1_[j].end(); photon_pointer++)
	    {
	      const list<EXTRA_PHOTON_POINTER>& extr_phot_list = extra_photon_pointer2_[j];
	      
	      for (extr_phot = extr_phot_list.begin() ; extr_phot!=extr_phot_list.end(); extr_phot++)
		{

		  collide_gg_XX(1, i1,i2,min_z, *photon_pointer, *extr_phot,switches, secondaries, muons, hasard, returnCrossSectionToAdd);
		  for (k=0; k<3;k++) tempor[k] += returnCrossSectionToAdd[k];
		}
	    }

	  results_. cumulate_hadrons_ge(tempor);
	  for (k=0; k<3;k++) {
	    tempor[k]=0.0;
	    returnCrossSectionToAdd[k]=0.0;
	  }
	  for( photon_pointer=grid_photon2_[j].begin(); photon_pointer!= grid_photon2_[j].end(); photon_pointer++)
	    {
	      const list<EXTRA_PHOTON_POINTER>& extr_phot_list = extra_photon_pointer1_[j];
	      for ( extr_phot = extr_phot_list.begin(); extr_phot != extr_phot_list.end(); extr_phot++)
		{
		  collide_gg_XX(1,i1, i2,min_z,  *extr_phot, *photon_pointer,switches,secondaries, muons, hasard, returnCrossSectionToAdd);
		  for (k=0; k<3;k++) tempor[k] += returnCrossSectionToAdd[k];
		}
	    }
	  results_. cumulate_hadrons_eg(tempor);
	  for (k=0; k<3;k++) {
	    tempor[k]=0.0;
	    returnCrossSectionToAdd[k]=0.0;
	  }

	  const list<EXTRA_PHOTON_POINTER>& extr_phot1_list = extra_photon_pointer1_[j];
	  list<EXTRA_PHOTON_POINTER>::const_iterator extr_phot1;
	  list<EXTRA_PHOTON_POINTER>::const_iterator extr_phot2;
	  for( extr_phot1 = extr_phot1_list.begin(); extr_phot1!=extr_phot1_list.end(); extr_phot1++)
	    {

	      const list<EXTRA_PHOTON_POINTER>& extr_phot2_list = extra_photon_pointer2_[j];

	      for ( extr_phot2 = extr_phot2_list.begin(); extr_phot2!=extr_phot2_list.end(); extr_phot2++)
		{

		  collide_gg_XX(2, i1, i2,min_z, *extr_phot1, *extr_phot2, switches,secondaries, muons, hasard, returnCrossSectionToAdd);
		  for (k=0; k<3;k++) tempor[k] += returnCrossSectionToAdd[k];
		}
	    }
	  results_. cumulate_hadrons_ee(tempor);
	}
    }
}


void GRID::photon_lumi_3(float min_z,SWITCHES& switches, PAIR_BEAM& secondaries,  RNDM& hasard)
{
  int i1,i2,j,n_x,n_y;

  list<BEAM_PARTICLE_POINTER>::iterator particle_pointer;

    list<EXTRA_PHOTON_POINTER>::const_iterator extr_phot1, extr_phot2;

  n_x = n_cell_x_;
  n_y = n_cell_y_;
  for (i1=0;i1<n_x;i1++)
    {
      for (i2=0;i2<n_y;i2++)
	{
	  j=i1*n_y+i2;
	  	  const list<EXTRA_PHOTON_POINTER>& extr_phot1_list = extra_photon_pointer1_[j];

	  for( extr_phot1 = extr_phot1_list.begin(); extr_phot1!=extr_phot1_list.end(); extr_phot1++)
	    {
	      for ( particle_pointer = part_pointer2_[j].begin(); particle_pointer != part_pointer2_[j].end(); particle_pointer++)
		{
		  collide_compton(1,i1, i2,min_z, *extr_phot1, *particle_pointer, secondaries, switches, hasard, 1);
		}
	    }

	  	  const list<EXTRA_PHOTON_POINTER>& extr_phot2_list = extra_photon_pointer2_[j];

	  for (extr_phot2 = extr_phot2_list.begin(); extr_phot2!=extr_phot2_list.end(); extr_phot2++)
	    {
	      for (particle_pointer = part_pointer1_[j].begin(); particle_pointer != part_pointer1_[j].end(); particle_pointer++)
		{

		  collide_compton(1,i1, i2,min_z, *extr_phot2, *particle_pointer, secondaries, switches, hasard,2);
		}
	    }
	}
    }
}


void GRID::move_photons2(BEAM& beam,int ibeam,int i_slice, RNDM& hasard)
{
  float ax,ay,wgt;
  PHI_FLOAT phi1_x,phi2_x,phi3_x,phi1_y,phi2_y,phi3_y,h_x,h_y;
  float radius_i;
  const PHI_FLOAT *phi;
    PHI_FLOAT upsilon,upsilon0,tmp;
  PHI_FLOAT x,y;
  int i, nd;


  float scal_step_local;
  //  scal_step_local = scal_step[ibeam-1];
  scal_step_local = slice_of_beam_[ibeam-1].get_scal_step();
  if (ibeam==1){
    //    phi=phi2_;
    phi = champ_.get_phi(2);
    wgt= slice_of_beam_[0].get_nb_part_per_macro();
  }
  else{
    //   phi=phi1_;
    phi = champ_.get_phi(1);
    wgt= slice_of_beam_[1].get_nb_part_per_macro();
  }
  vector<PHOTON>& lesPhotons =  beam.getPhotonVector(i_slice);
  unsigned int k;
  for(k=0; k < lesPhotons.size(); k++)
    {
      lesPhotons[k].XYposition(x, y);
      if ( coordinatesInGridRange(x, y) )
	{ 
	  interpolePotential(x, y, h_x, h_y, phi1_x, phi2_x, phi3_x, phi1_y, phi2_y, phi3_y, phi);
	  ax = (h_x*(phi1_y-phi2_y)+(1.0-h_x)*(phi2_y-phi3_y))*delta_x_inverse_;
	  ay = (h_y*(phi1_x-phi2_x)+(1.0-h_y)*(phi2_x-phi3_x))*delta_y_inverse_;
	}
      else
	{
	  ax = 0.0;
	  ay = 0.0;
	}
      float facteur = step_*scal_step_local;
      //      photItr->move(facteur);
      lesPhotons[k].advancePosition(facteur);

      radius_i = sqrt(ax*ax+ay*ay)*1e9/scal_step_local;
      
      upsilon0=LAMBDA_BAR/EMASS*radius_i;
      //     upsilon=upsilon0*photItr->get_energy()/EMASS;
      upsilon=upsilon0*lesPhotons[k].energy()/EMASS;
      
      if (upsilon>1e-10)
	{
	  tmp=ALPHA_EM*ALPHA_EM/RE*step_*1e-9*scal_step_local*upsilon0*0.23*exp(-8.0/(3.0*upsilon));
	  tmp*=pow(1.0+0.22*upsilon,-1.0/3.0);
	  coherent_results_.updateProbmax(tmp);
	  tmp*=1.36;
	  if (hasard.rndm()<tmp) {
	    if (tmp<0.1) {
	      if (coherent_generate(beam,i_slice,upsilon,lesPhotons[k], hasard)) {
		coherent_results_.updateSumeng(wgt*lesPhotons[k].energy());
		lesPhotons[k].razEnergy();//photon->energy=0.0;
		tmp=1.0;
		//nyes++;
	      }
	      else {
		tmp=0.0;
		//nno++;
	      }
	    }
	    else {
	      nd=(int)(tmp/0.1+1.0);
	      int ds=tmp/nd;
	      tmp=0.0;
	      for (i=0;i<nd;++i) {
		if (hasard.rndm()<ds) {
		  if (coherent_generate(beam,i_slice,upsilon,lesPhotons[k], hasard)) {
		    coherent_results_.updateSumeng(wgt*lesPhotons[k].energy());
		    lesPhotons[k].razEnergy();//photon->energy=0.0;
		    tmp=1.0;
		    //nyes++;
		    i=nd;
		  }
		}
	      }
	    }
	  }
	  else {
	    tmp=0.0;
	  }
	  /*
	    if (tmp>1.0)
	    {
	    //	      coherent_results.updateSumeng(wgt*photItr->get_energy());
	    coherent_results_.updateSumeng(wgt*lesPhotons[k].energy());
	    coherent_generate(beam,i_slice,upsilon,lesPhotons[k], hasard);
	    //	      photItr->razEnergy();
	    lesPhotons[k].razEnergy();
	    tmp=1.0;
	    }
	    else
	    {
	    if (hasard.rndm()<tmp*1.36) 
	    {
	    coherent_generate(beam,i_slice,upsilon,lesPhotons[k], hasard);
	    //		  photItr->razEnergy();
	    lesPhotons[k].razEnergy();
	    }
	    }
	    // end test
	  */ 
	  tmp*=wgt;
	  //	  coherent_results.updateSums(tmp, photItr->get_energy(), upsilon0);
	  coherent_results_.updateSums(tmp, lesPhotons[k].energy(), upsilon0);
	}
      //      photItr++;
    }
}




//void GRID::coherent_generate(BEAM& beam,int i_slice, float upsilon,PHOTON& phot, RNDM& hasard)
int GRID::coherent_generate(BEAM& beam,int i_slice, float upsilon,PHOTON& phot, RNDM& hasard)
{
  float x,y;
  float vx,vy;
  float energy;
  energy =  phot.energy();
  // number adjustment to be done for total cross section
  if (pick_coherent_energy(upsilon,energy, hasard))
    {
      phot.XYposition(x, y);
      phot.velocities(vx,vy);
      store_coherent_particle(beam,energy,vx,vy,x,y,i_slice);
      store_coherent_particle(beam,-(phot.energy()-energy),vx,vy,x,y,
		     i_slice);
      return 1;
    }
  return 0;
}

bool GRID::pick_coherent_energy(float ups,float& energy, RNDM& hasard)
{
  float a=0.13513,eta,dxdy,x,h,hp,tmp,tmp2,y,coh;
  y=2.0*hasard.rndm()-1.0;
  eta=8.0/(3.0*ups);
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
  coh=PHYSTOOLS::u(ups);
  tmp=PHYSTOOLS::fcp(ups,x)*a/coh*dxdy;
  if (hasard.rndm()<tmp){
     energy *= x;
    return true;
  } else {
    return false;
  }
}


// void GRID::step_pair_1X(const vector<GENERAL_GRID*>& grids, PAIR_PARTICLE& paire,float step,int nbSteps, int extra_grids, float charge_sign_0, RNDM& hasard)
// {
//   const float eps=1e-35,emass2=EMASS*EMASS;
//   float x,y,z,vx,vy,vz,e_inv2,e_inv,step_2,step_q,vold2,scal,thetamax;
//   float energie;
//   float ex,ey,bx,by,b_norm,b_norm_i,theta,a1,a2,a3,vb,eng;
//   float ph[1000],vx0,vy0,vz0,eng0;
//   int nph;
//   int i,icharge,j;
  
//   paire.get_parameters(x,y,z,vx,vy,vz,energie);
//   if (energie >0.0)
//     {
//       eng = energie;
//       icharge=0;
//     }
//   else
//     {
//       eng=-energie;
//       icharge=1;
//     }
//   e_inv=1.0/eng;
//   e_inv2=e_inv*e_inv;
//   step_2=0.5*step;
//   step_q=step*step;
  
//   // initial half step 
  
//   field_pair(grids,x,y,&ex,&ey,&bx,&by, extra_grids, charge_sign_0);
//   if (icharge)
//     {
//       ex=-ex;
//       ey=-ey;
//       bx=-bx;
//       by=-by;
//     }
//   b_norm=sqrt(bx*bx+by*by);
//   b_norm_i=1.0/max(b_norm,eps);
//   bx*=b_norm_i;
//   by*=b_norm_i;
//   vb=vx*by-vy*bx;
// #ifdef PAIR_SYN
//   vx0=vx;
//   vy0=vy;
//   vz0=vz;
// #endif
//   theta=0.25*(vz*vz*(bx*bx+by*by)+vb*vb)*b_norm*b_norm*e_inv2*step_q;
//   a3=0.5*theta;
//   a1=1.0-a3;
//   theta=sqrt(theta);
//   thetamax=2.0*theta;
//   a2=theta*vz;
//   a3*=vx*bx+vy*by;
//   vz=vz*a1+theta*vb;
//   vx=vx*a1-a2*by+a3*bx;
//   vy=vy*a1+a2*bx+a3*by;
//   vold2=vx*vx+vy*vy+vz*vz;
//   vx+=step_2*ex*e_inv;
//   vy+=step_2*ey*e_inv;
// #ifdef SCALE_ENERGY
//   scal=sqrt((vold2*eng*eng+emass2)/((vx*vx+vy*vy+vz*vz)*eng*eng+emass2));
//   vx*=scal;
//   vy*=scal;
//   vz*=scal;
// #ifdef PAIR_SYN
//   vx0*=scal;
//   vy0*=scal;
//   vz0*=scal;
// #endif
//   eng/=scal;
//   e_inv=1.0/eng;
//   e_inv2=e_inv*e_inv;
// #endif
// #ifdef PAIR_SYN
//   synrad(eng,
// 	 sqrt((vx-vx0)*(vx-vx0)+(vy-vy0)*(vy-vy0)+(vz-vz0)*(vz-vz0))
// 	 /(0.5*step*1e-9),
// 	 step*1e-9,ph,&nph, hasard);
//   if (nph>0) {
//     eng0=eng;
//     for (j=0;j<nph;j++){
//       eng-=ph[j];
//     }
//     scal=sqrt(((eng*eng-emass2)*eng0*eng0)
// 	      /((eng0*eng0-emass2)*eng*eng));
//     vx*=scal;
//     vy*=scal;
//     vz*=scal;
//   }
// #endif
//   /* loop over steps */
//   for(i=1;i<nbSteps;i++)
//     {
//       x+=vx*step;
//       y+=vy*step;
//       z+=vz*step;
//       field_pair(grids,x,y,&ex,&ey,&bx,&by, extra_grids, charge_sign_0);
//       if (icharge){
// 	ex=-ex;
// 	ey=-ey;
// 	bx=-bx;
// 	by=-by;
//       }
      
// #ifdef PAIR_SYN
//       vx0=vx;
//       vy0=vy;
//       vz0=vz;
// #endif
//       /* scd new */
//       vold2=vx*vx+vy*vy+vz*vz;
//       vx+=step_2*ex*e_inv;
//       vy+=step_2*ey*e_inv;
// #ifdef SCALE_ENERGY
//       scal=sqrt((vold2*eng*eng+emass2)/((vx*vx+vy*vy+vz*vz)*eng*eng+emass2));
//       vx*=scal;
//       vy*=scal;
//       vz*=scal;
// #ifdef PAIR_SYN
//       vx0*=scal;
//       vy0*=scal;
//       vz0*=scal;
// #endif
//       eng/=scal;
//       e_inv=1.0/eng;
//       e_inv2=e_inv*e_inv;
// #endif
//       b_norm=sqrt(bx*bx+by*by);
//       b_norm_i=1.0/max(b_norm,eps);
//       bx*=b_norm_i;
//       by*=b_norm_i;
//       vb=vx*by-vy*bx;
//       /*vb=0.0;*/
//       theta=(vz*vz*(bx*bx+by*by)+vb*vb)*b_norm*b_norm*e_inv2*step_q;
//       a3=0.5*theta;
//       a1=1.0-a3;
//       theta=sqrt(theta);
//       thetamax=max(thetamax,theta);
//       a2=theta*vz;
//       a3*=vx*bx+vy*by;
//       /*a3=0.0;
// 	a1=1.0;*/
//       vz=vz*a1+theta*vb;
//       vx=vx*a1-a2*by+a3*bx;
//       vy=vy*a1+a2*bx+a3*by;
      
//       /* scd new */
//       vold2=vx*vx+vy*vy+vz*vz;
//       vx+=step_2*ex*e_inv;
//       vy+=step_2*ey*e_inv;
// #ifdef SCALE_ENERGY
//       scal=sqrt((vold2*eng*eng+emass2)/((vx*vx+vy*vy+vz*vz)*eng*eng+emass2));
//       vx*=scal;
//       vy*=scal;
//       vz*=scal;
// #ifdef PAIR_SYN
//       vx0*=scal;
//       vy0*=scal;
//       vz0*=scal;
// #endif
//       eng/=scal;
//       e_inv=1.0/eng;
//       e_inv2=e_inv*e_inv;
// #endif
// #ifdef PAIR_SYN
//       /* new for test */
//       synrad(eng,
// 	     sqrt((vx-vx0)*(vx-vx0)+(vy-vy0)*(vy-vy0)+(vz-vz0)*(vz-vz0))/step*1e9, step*1e-9,ph,&nph, hasard);
//       if (nph>0) {
// 	eng0=eng;
// 	for (j=0;j<nph;j++){
// 	  eng-=ph[j];
// 	}
// 	scal=sqrt(((eng*eng-emass2)*eng0*eng0)
// 		  /((eng0*eng0-emass2)*eng*eng));
// 	vx*=scal;
// 	vy*=scal;
// 	vz*=scal;
//       }
// #endif
//     }
//   /* last half step */
//   x+=vx*step;
//   y+=vy*step;
//   z+=vz*step;
//   field_pair(grids,x,y,&ex,&ey,&bx,&by, extra_grids, charge_sign_0);
//   if (icharge){
//     ex=-ex;
//     ey=-ey;
//     bx=-bx;
//     by=-by;
//   }
  
//   /* scd new */
//   vold2=vx*vx+vy*vy+vz*vz;
//   vx+=step_2*ex*e_inv;
//   vy+=step_2*ey*e_inv;
// #ifdef SCALE_ENERGY
//   scal=sqrt((vold2*eng*eng+emass2)/((vx*vx+vy*vy+vz*vz)*eng*eng+emass2));
//   vx*=scal;
//   vy*=scal;
//   vz*=scal;
//   eng/=scal;
//   e_inv=1.0/eng;
//   e_inv2=e_inv*e_inv;
// #endif
//   b_norm=sqrt(bx*bx+by*by);
//   b_norm_i=1.0/max(b_norm,eps);
//   bx*=b_norm_i;
//   by*=b_norm_i;
//   vb=vx*by-vy*bx;
//   /*vb=0.0;*/
//   theta=0.25*(vz*vz*(bx*bx+by*by)+vb*vb)*b_norm*b_norm*e_inv2*step_q;
//   a3=0.5*theta;
//   a1=1.0-a3;
//   theta=sqrt(theta);
//   thetamax=max(thetamax,float(2.0)*theta);
//   a2=theta*vz;
//   a3*=vx*bx+vy*by;
//   /*a3=0.0;
//     a1=1.0;*/
//   vz=vz*a1+theta*vb;
//   vx=vx*a1-a2*by+a3*bx;
//   vy=vy*a1+a2*bx+a3*by;
// #ifdef SCALE_ENERGY
//   scal=sqrt((eng*eng-emass2)/((eng*eng)*(vx*vx+vy*vy+vz*vz)));
//   if (fabs(scal-1.0)>0.01)
//     printf("> %g %g %g %g %g\n",eng,eng-fabs(energie),thetamax,
// 	   sqrt(vx*vx+vy*vy+vz*vz),scal);
// #else
//   scal=1.0;
// #endif
  
//   if (icharge) energie = -eng;
//   else energie =eng; 
//   paire.set(x,y,z,vx*scal,vy*scal,vz*scal,energie);
// }
void GRID::step_pair_1(const vector<GENERAL_GRID*>& grids, PAIR_PARTICLE& paire, double mass,float step,int nbSteps, int extra_grids, float charge_sign_0, RNDM& hasard)
{
  float thetamax;
  float ex,ey,bx,by;
  int i;
  // initial half step 
 
  // on recupere les champs E et B
  field_pair(paire, grids,ex,ey,bx,by, extra_grids, charge_sign_0);
  thetamax = 2.0*paire.apply_initial_half_step_fields(step, mass, ex,ey, bx, by, hasard);
  /* loop over steps */
  for(i=1;i<nbSteps;i++)
    {
      paire.advancePosition(step);

      // on recupere les champs E et B 
      field_pair(paire, grids,ex,ey,bx,by, extra_grids, charge_sign_0);    
      thetamax = max (thetamax, paire.apply_full_step_fields(step, mass, ex,ey, bx, by, hasard)); 
    }
  /* last half step */
  paire.advancePosition(step);
  field_pair(paire, grids,ex,ey,bx,by, extra_grids, charge_sign_0);
  thetamax = max( thetamax, (float)2.0*paire.apply_final_half_step_fields(step, mass, ex,ey, bx, by, thetamax,hasard));
  if (	!paire.last_rescaling_ok() )
    {
      cerr << " GRID::step_pair_1() : " << " thetamax " << thetamax << endl;
    }	  
}

void GRID::step_pair_1_tertphot(const vector<GENERAL_GRID*>& grids, PAIR_PARTICLE& paire, double mass,float step,int nbSteps, int extra_grids, float charge_sign_0, RNDM& hasard)
{
  float thetamax;
  float ex,ey,bx,by;
  int i;
  vector<float> photon_e;
  // initial half step 
 
  // on recupere les champs E et B
  field_pair(paire, grids,ex,ey,bx,by, extra_grids, charge_sign_0);
  thetamax = 2.0*paire.apply_initial_half_step_fields(step, mass, ex,ey, bx, by, &photon_e, hasard);
  /* loop over steps */
  for(i=1;i<nbSteps;i++)
    {
      paire.advancePosition(step);

      // on recupere les champs E et B 
      field_pair(paire, grids,ex,ey,bx,by, extra_grids, charge_sign_0);    
      thetamax = max (thetamax, paire.apply_full_step_fields(step, mass, ex,ey, bx, by, &photon_e, hasard)); 
    }
  /* last half step */
  paire.advancePosition(step);
  field_pair(paire, grids,ex,ey,bx,by, extra_grids, charge_sign_0);
  thetamax = max( thetamax, (float)2.0*paire.apply_final_half_step_fields(step, mass, ex,ey, bx, by, thetamax,hasard));
  for(i=0;i<photon_e.size();i++)
    {
      tertphot_.push_back(TERTPHOTON(photon_e[i],paire));
    }
  if (	!paire.last_rescaling_ok() )
    {
      cerr << " GRID::step_pair_1_tertphot() : " << " thetamax " << thetamax << endl;
    }	  
}

void GRID::apply_magnetic_field_on_pair(float fac_theta, float step_q, float e_inv2, float bx, float by, float& vx,float& vy, float& vz, float& theta) const
{
  const float eps=1e-35;
  float b_norm, b_norm_i,a1,a2,a3,vb;

  // on normalise B
  b_norm=sqrt(bx*bx+by*by);
  b_norm_i=1.0/max(b_norm,eps);
  bx *= b_norm_i;
  by *= b_norm_i;

  // v x B
  vb = vx*by-vy*bx;


  theta=fac_theta*(vz*vz*(bx*bx+by*by)+vb*vb)*b_norm*b_norm*e_inv2*step_q;
  a3=0.5*theta;

  // 1 - theta**2/2
  a1=1.0-a3;

  // theta
  theta=sqrt(theta);
  a2=theta*vz;

  // a3 = (theta**/2).(v.B)  [B norme]
  a3*=vx*bx+vy*by;

  // application du champ magnetique (??)
  vx=vx*a1-a2*by+a3*bx;
  vy=vy*a1+a2*bx+a3*by;
  vz=vz*a1+theta*vb;
}


int GRID::field_pair(const PAIR_PARTICLE& paire, const vector<GENERAL_GRID*>& grids,float& ex,float& ey, float& bx,float& by, int extra_grids, float charge_sign_0)
{
  PHI_FLOAT ax_1,ay_1,ax_2,ay_2;
  PHI_FLOAT tmp,dx,dy;
  float Ex1, Ey1,Ex2, Ey2;
  float x,y;
  int i;
  TRIVECTOR Efield1, Efield2;

  paire.XYposition(x, y);
  for(i=0;i<=extra_grids;i++)
    {
      if (grids[i]->coordinatesInGridRange(x,y))
	{
	  Efield1 = grids[i]->ElectricFieldCIC(x, y, grids[i]->get_phi1());
	  Efield2 = grids[i]->ElectricFieldCIC(x, y, grids[i]->get_phi2());

	  Efield1 *= -charge_sign_0;

	  Ex1 = Efield1(0);
	  Ey1 = Efield1(1);
	  Ex2 = Efield2(0);
	  Ey2 = Efield2(1);
	  ex = -(Ex2-Ex1);
	  ey = -(Ey2-Ey1);
	  by =  (Ex2 + Ex1);
	  bx = -(Ey2 + Ey1);

	  return i;
	}
    }
/* particle is not in largest grid */
  dx=(x-grids[0]->get_rho_x_1());
  dy=(y-grids[0]->get_rho_y_1());
  tmp=grids[0]->get_rho_factor()*grids[0]->get_rho_sum_1()/(dx*dx+dy*dy);
  tmp*=-charge_sign_0;

/*  tmp=grids[0].rho_factor*grids[0].rho_sum_1/(dx*dx+dy*dy);*/
  ax_1=dx*tmp;
  ay_1=dy*tmp;
  dx=(x-grids[0]->get_rho_x_2());
  dy=(y-grids[0]->get_rho_y_2());
  tmp=grids[0]->get_rho_factor()*grids[0]->get_rho_sum_2()/(dx*dx+dy*dy);
  ax_2=dx*tmp;
  ay_2=dy*tmp;

  ex = -0.5*(ax_2-ax_1);
  ey = -0.5*(ay_2-ay_1);
  by = 0.5*(ax_2+ax_1);
  bx = -0.5*(ay_2+ay_1);
  return -1;
}



EXTRA_GRID::~EXTRA_GRID()
{
}

//same as distribute_particles0 but without counting particles as being on or off the grid 

void EXTRA_GRID::distribute_particles(int i_slice1,
				   int i_slice2,
				   int electron_distribution_rho, 
				   int force_symmetric)
{

  slice_of_beam_[0].razRho();
  slice_of_beam_[1].razRho();

  //  razRhos(); 
  switch (electron_distribution_rho)
    {
    case 1:

//       assignBeamSliceNGP(beam1_, i_slice1, rho1_, nb_part_per_macro_[0]);
//       assignBeamSliceNGP(beam2_, i_slice2, rho2_, nb_part_per_macro_[1]);
      assignBeamSliceNGP(slice_of_beam_[0], i_slice1);
      assignBeamSliceNGP(slice_of_beam_[1], i_slice2);

      break;
    case 2:

//        assignBeamSliceCIC(beam1_, i_slice1, rho1_, nb_part_per_macro_[0]);
//       assignBeamSliceCIC(beam2_, i_slice2, rho2_, nb_part_per_macro_[1]);
      assignBeamSliceCIC(slice_of_beam_[0], i_slice1);
      assignBeamSliceCIC(slice_of_beam_[1], i_slice2);
      break;
    }
  distribute_coherent_particles(i_slice1, i_slice2, electron_distribution_rho);
  distribute_trident_particles(i_slice1, i_slice2, electron_distribution_rho);
  if (force_symmetric)
    {
      slice_of_beam_[0].symmetrizeCharges(n_cell_x_,n_cell_y_);
      slice_of_beam_[1].symmetrizeCharges(n_cell_x_,n_cell_y_);
    }

    // symmetrizeCharges();
}


//  Distribute the coherent pair particles for the calculation of the fields
//  but does not count if they are on the grid or not (this routine is used for
//  the larger grids)
 

void EXTRA_GRID::distribute_coherent_particles(int i_slice1,
					    int i_slice2,
					    int electron_distribution_rho)
{
  switch (electron_distribution_rho)
    {
    case 1:


//       assignCoherentBeamSliceNGP(beam1_, i_slice1, rho1_, nb_part_per_macro_[0]);
//       assignCoherentBeamSliceNGP(beam2_, i_slice2, rho2_, nb_part_per_macro_[1]);
      assignCoherentBeamSliceNGP(slice_of_beam_[0], i_slice1);
      assignCoherentBeamSliceNGP(slice_of_beam_[1], i_slice2);
      break;
    case 2:


      assignCoherentBeamSliceCIC(slice_of_beam_[0], i_slice1);
      assignCoherentBeamSliceCIC(slice_of_beam_[1], i_slice2);
      break;
    }
}

void EXTRA_GRID::distribute_trident_particles(int i_slice1,
					    int i_slice2,
					    int electron_distribution_rho)
{
  switch (electron_distribution_rho)
    {
    case 1:
      assignTridentBeamSliceNGP(slice_of_beam_[0], i_slice1);
      assignTridentBeamSliceNGP(slice_of_beam_[1], i_slice2);
      break;
    case 2:
      assignTridentBeamSliceCIC(slice_of_beam_[0], i_slice1);
      assignTridentBeamSliceCIC(slice_of_beam_[1], i_slice2);
      break;
    }
}

// void EXTRA_GRID::assignCoherentBeamSliceNGP(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro)
// {
//   int i1, i2;
//   float xpart, ypart;  
//   float ch;
//   unsigned int particle;
//   const vector<PARTICLE*>& coherent = beam->getCoherentVector(i_slice);
//   for (particle = 0; particle < coherent.size(); particle++)
//     {
//       coherent[particle]->XYposition(xpart, ypart);
//       if (coherent[particle]->energy() < 0.0) ch=-1.0;
//       else ch=1.0;
//       if (particleInGrid(xpart, ypart,i1, i2)) assignChargeToNGP(rho, i1, i2, n_macro*ch);	    
//     }
// }
void EXTRA_GRID::assignCoherentBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  unsigned int particle;
  const vector<PARTICLE*>& coherent = sog.get_beam()->getCoherentVector(i_slice);
  for (particle = 0; particle < coherent.size(); particle++)
    {
      coherent[particle]->XYposition(xpart, ypart);
      if (coherent[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;
      if (particleInGrid(xpart, ypart,i1, i2)) sog.assignChargeToNGP(i1, i2,ch);	    
    }
}
void EXTRA_GRID::assignTridentBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  unsigned int particle;
  const vector<PARTICLE*>& trident = sog.get_beam()->getTridentVector(i_slice);
  for (particle = 0; particle < trident.size(); particle++)
    {
      trident[particle]->XYposition(xpart, ypart);
      if (trident[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;
      if (particleInGrid(xpart, ypart,i1, i2)) sog.assignChargeToNGP(i1, i2,ch);	    
    }
}
// void EXTRA_GRID::assignBeamSliceNGP(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro)
// {
//   int k, i1, i2;
//   float xpart, ypart;
//     vector<PARTICLE*>& lesParticles = beam->getParticleVector(i_slice);
//     for (k = 0; k < lesParticles.size(); k++)
//    {
// 	lesParticles[k]->XYposition(xpart, ypart);
//       if (particleInGrid(xpart, ypart, i1, i2)) assignChargeToNGP(rho, i1, i2, n_macro);
//     }
// }
void EXTRA_GRID::assignBeamSliceNGP(SLICE_ON_GRID& sog, int i_slice)
{
  unsigned int k;
  int i1, i2;
  float xpart, ypart;
    const vector<PARTICLE*>& lesParticles = sog.get_beam()->getParticleVector(i_slice);
    for (k = 0; k < lesParticles.size(); k++)
   {
	lesParticles[k]->XYposition(xpart, ypart);
      if (particleInGrid(xpart, ypart, i1, i2)) sog.assignChargeToNGP(i1, i2,1.0);
    }
}
// void EXTRA_GRID::assignBeamSliceCIC(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro)
// {
//   int k, i1, i2;
//   float h_x,h_y;
//   float xpart, ypart;
//     vector<PARTICLE*>& lesParticles = beam->getParticleVector(i_slice);
//     for (k = 0; k < lesParticles.size(); k++)
// 	{
// 	lesParticles[k]->XYposition(xpart, ypart);
// 	  if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2)) assignChargeToCIC(rho, i1, i2, h_x, h_y, n_macro);	    	  
// 	}
// }
void EXTRA_GRID::assignBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice)
{
  unsigned int k;
  int i1, i2;
  float h_x,h_y;
  float xpart, ypart;
  const vector<PARTICLE*>& lesParticles = sog.get_beam()->getParticleVector(i_slice);
  for (k = 0; k < lesParticles.size(); k++)
    {
      lesParticles[k]->XYposition(xpart, ypart);
      if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2)) sog.assignChargeToCIC(i1, i2, h_x, h_y, 1.0);	    	  
    }
}
// void EXTRA_GRID::assignCoherentBeamSliceCIC(BEAM *beam, int i_slice, PHI_FLOAT *rho, float n_macro)
// {
//   int i1, i2;
//   float xpart, ypart;  
//   float ch;
//   float h_x,h_y;
//   unsigned int particle;
//   const vector<PARTICLE*>& coherent = beam->getCoherentVector(i_slice);
//   for (particle = 0; particle < coherent.size(); particle++)
//     {
//       coherent[particle]->XYposition(xpart, ypart);
//       if (coherent[particle]->energy() < 0.0) ch=-1.0;
//       else ch=1.0;
//       if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2)) assignChargeToCIC(rho, i1, i2, h_x, h_y, n_macro*ch);	
//     }
// }
void EXTRA_GRID::assignCoherentBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  float h_x,h_y;
  unsigned int particle;
  const vector<PARTICLE*>& coherent = sog.get_beam()->getCoherentVector(i_slice);
  for (particle = 0; particle < coherent.size(); particle++)
    {
      coherent[particle]->XYposition(xpart, ypart);
      if (coherent[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;
      if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2)) sog.assignChargeToCIC(i1, i2, h_x, h_y, ch);	
    }
}

void EXTRA_GRID::assignTridentBeamSliceCIC(SLICE_ON_GRID& sog, int i_slice)
{
  int i1, i2;
  float xpart, ypart;  
  float ch;
  float h_x,h_y;
  unsigned int particle;
  const vector<PARTICLE*>& trident = sog.get_beam()->getTridentVector(i_slice);
  for (particle = 0; particle < trident.size(); particle++)
    {
      trident[particle]->XYposition(xpart, ypart);
      if (trident[particle]->energy() < 0.0) ch=-1.0;
      else ch=1.0;
      if (particleInGrid(xpart, ypart, h_x, h_y, i1, i2)) sog.assignChargeToCIC(i1, i2, h_x, h_y, ch);	
    }
}
