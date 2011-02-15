#include <iostream>
#include <cstdlib>
#include <sstream>
#include "switchesCPP.h"
#include "physconst.h"

using namespace std;
SWITCHES::SWITCHES() 
{
  electron_distribution_rho=2;
  //    electron_distribution_scatter=1;
  electron_ratio=1.0;
  do_lumi=0;
  num_lumi=10000;
  do_cross=0;
  lumi_p=1e-23;
  do_photons_[0]=0;
  do_photons_[1]=0;
  write_photons=0;
  photon_distribution=1;
  photon_ratio=1.0;
  do_hadrons=0;
  store_hadrons=0;
  hadron_ratio=1e4;
  do_jets=0;
  jet_store=0;
  jet_pstar=3.2;
  jet_ratio=1e5;
  do_pairs=0;
  load_event=0;
  //track_secondaries=0;
  track_pairs=0;
  pair_ratio=1.0;
  bhabha_ratio=1.0;
  pair_ecut=EMASS;
  integration_method=2;
  extra_grids=0;
  time_order=2;
  interpolation=2;
  adjust=0;
  geom=1;
  r_scal=1.0;
  jet_pythia=0;
  jet_select=0;
  pair_q2=1;
  load_photon=0;
  load_beam=0;
  cuts_from_loaded_beam = 0;
  bmt_precession_ = 0;
  ST_spin_flip_ = 0;
  automatic_grid_sizing = 0;
  emin=1.0;
  charge_sign=-1.0;
  do_beamstrahlung_=1;
  store_beam=0;
  do_cross_gg=0;
  force_symmetric=0;
  do_isr=1;
  do_espread=0;
  do_coherent=0;
  do_trident=0;
  gg_smin=4.0*150.0*150.0;
  twobeam=0;
}


// void SWITCHES::lectureFirstBeamParameters(const PARAMETERS& param)
// {
//    n_b= param.readDValue("n_b");
//    f_rep=param.readDValue("f_rep");
//    charge_sign=param.readFValue("charge_sign");

//   if(charge_sign>0.0){
//     charge_sign_0=1.0;
//   }
//   if(charge_sign<0.0){
//     charge_sign_0=-1.0;
//   }

// }

void SWITCHES::read(const PARAMETERS& param)
{

   n_b= param.readDValue("n_b");
   f_rep=param.readDValue("f_rep");
   charge_sign=param.readFValue("charge_sign");

  if(charge_sign>0.0){
    charge_sign_0=1.0;
  }
  if(charge_sign<0.0){
    charge_sign_0=-1.0;
  }


  integration_method = param.readIValue("integration_method");
  silent= param.readIValue("silent");


  extra_grids = param.readIValue("grids")-1;


  if (extra_grids<0) extra_grids=0;

  load_photon = param.readIValue("load_photons");
  write_photons = param.readIValue("store_photons");
  do_photons_[0] = param.readIValue("do_photons.1");

  do_photons_[1] = param.readIValue("do_photons.2");
  
  do_beamstrahlung_ = param.readIValue("do_eloss");

  ecm_min = param.readFValue("ecm_min"); 

  float aux = param.readFValue("ecm_min_gg");
  gg_smin=4*aux*aux;

  do_hadrons = param.readIValue("do_hadrons");

  store_hadrons = param.readIValue("store_hadrons");

  hadron_ratio = param.readFValue("hadron_ratio");

  do_jets = param.readIValue("do_jets");

  do_pairs = param.readIValue("do_pairs");

  beam_pair = param.readIValue("beam_pair");

  //store_secondaries = param.readIValue("store_secondaries");
  store_pairs = param.readIValue("store_pairs");

  do_muons = param.readIValue("do_muons");

  do_coherent = param.readIValue("do_coherent");

  do_trident = param.readIValue("do_trident");

  emin = param.readFValue("emin");

  //track_secondaries = param.readIValue("track_secondaries");
  track_pairs = param.readIValue("track_pairs");

  pair_ecut = param.readFValue("pair_ecut");

  pair_step = param.readFValue("pair_step");

  electron_ratio = param.readFValue("electron_ratio");

  do_lumi = param.readIValue("do_lumi");

  bhabha_scal = param.readFValue("bhabha_scal");

  bhabha_ecmload = param.readFValue("bhabha_ecmload");

  do_bhabhas = param.readIValue("do_bhabhas");

  rndm_save = param.readIValue("rndm_save");

  rndm_load = param.readIValue("rndm_load");

  rndm_seed = param.readIValue("rndm_seed");

  do_lumi_ee_2 = param.readIValue("do_lumi_ee_2");

  do_size_log = param.readIValue("do_size_log");

  lumi_ee_2_n = param.readIValue("lumi_ee_2_n");

  lumi_ee_2_xmin=param.readFValue("lumi_ee_2_min"); 

  lumi_ee_2_xmax=param.readFValue("lumi_ee_2_max"); 

  do_lumi_eg_2=param.readIValue("do_lumi_eg_2");

  lumi_eg_2_n=param.readIValue("lumi_eg_2_n");

  lumi_eg_2_xmin=param.readFValue("lumi_eg_2_min"); 

  lumi_eg_2_xmax=param.readFValue("lumi_eg_2_max"); 

  do_lumi_ge_2=param.readIValue("do_lumi_ge_2");

  lumi_ge_2_n=param.readIValue("lumi_ge_2_n");

  lumi_ge_2_xmin=param.readFValue("lumi_ge_2_min"); 

  lumi_ge_2_xmax=param.readFValue("lumi_ge_2_max"); 

  do_lumi_gg_2=param.readIValue("do_lumi_gg_2");

  lumi_gg_2_n=param.readIValue("lumi_gg_2_n");

  lumi_gg_2_xmin=param.readFValue("lumi_gg_2_min"); 

  lumi_gg_2_xmax=param.readFValue("lumi_gg_2_max"); 

  do_cross=param.readIValue("do_cross");

  do_prod=param.readIValue("do_prod");

  load_event=param.readIValue("load_events");

  prod_e=param.readFValue("prod_e");

  prod_scal=param.readFValue("prod_scal");

  do_compt=param.readIValue("do_compt");

  do_compt_phot=param.readIValue("do_compt_phot");

  compt_x_min=param.readFValue("compt_x_min");

  compt_scale=param.readFValue("compt_scale");

  compt_emax=param.readFValue("compt_emax");

  do_isr=param.readIValue("do_isr");

  do_espread=param.readIValue("do_espread");

  espread1=param.readFValue("espread.1");

  which_espread1=param.readIValue("which_espread.1");

  espread2=param.readFValue("espread.2");

  which_espread2=param.readIValue("which_espread.2");

  num_lumi=param.readIValue("num_lumi");

  num_lumi_eg=param.readIValue("num_lumi_eg");

  num_lumi_gg=param.readIValue("num_lumi_gg");

  lumi_p=param.readFValue("lumi_p");

  lumi_p_eg=param.readFValue("lumi_p_eg");

  lumi_p_gg=param.readFValue("lumi_p_gg");

  photon_ratio=param.readFValue("photon_ratio");

  pair_ratio=param.readFValue("pair_ratio");

  bhabha_ratio=param.readFValue("bhabha_ratio");

  jet_ratio=param.readFValue("jet_ratio");

  jet_pstar=param.readFValue("jet_ptmin");

  jet_store=param.readIValue("store_jets");

  jet_select=param.readIValue("jet_log");

  geom=param.readIValue("beam_size");

  r_scal=param.readFValue("beam_size_scale"); 

  ext_field=param.readIValue("ext_field");

//   if (!EXT_FIELD) 
//     {
//       if (ext_field){
// 	//	fprintf(stderr,"EXT_FIELD = false (LesDifines.h) \n");
// 	fprintf(stderr,"Cannot use flag ext_field\n");
// 	//	fprintf(stderr,"Please recompile with EXT_FIELD = true \n");
// 	exit(1);
//       }
      //    }
  pair_q2=param.readIValue("pair_q2");

  store_beam=param.readIValue("store_beam");

  load_beam=param.readIValue("load_beam");

  cuts_from_loaded_beam = param.readIValue("cuts_from_loaded_beam");

  bmt_precession_   = param.readIValue("bmt_precession");
  ST_spin_flip_   = param.readIValue("ST_spin_flip");

  automatic_grid_sizing = param.readIValue("automatic_grid_sizing");

  load_photon=param.readIValue("load_photons");

  jet_pythia=param.readIValue("jet_pythia");

  force_symmetric=param.readIValue("force_symmetric");

  charge_symmetric=param.readIValue("charge_symmetric");

  beam_vx_min=param.readFValue("beam_vx_min");

  beam_vx_max=param.readFValue("beam_vx_max");

  beam_vx_interval=param.readIValue("beam_vx_interval");

  beam_vy_min=param.readFValue("beam_vy_min");

  beam_vy_max=param.readFValue("beam_vy_max");

  beam_vy_interval=param.readIValue("beam_vy_interval");

  do_dump=param.readIValue("do_dump");

  dump_step=param.readIValue("dump_step");

  dump_particle=param.readIValue("dump_particle");
  if ( dump_step <= 0 ) dump_step = 1;
  if ( dump_particle <= 0 ) dump_particle = 1;
  check_consistency();
}

void SWITCHES::check_consistency() const
{
  if ( do_bhabhas && do_pairs)
    {
      cerr << " do_pairs= " << do_pairs << " do_bhabhas = " << do_bhabhas << endl;
      cerr << " ERROR : it is not allowed to have do_bhabhas= 1 together with do_pairs = 1 " << endl;
      exit(0);
    } 
  if ( do_bhabhas && do_compt)
    {
      cerr << " do_pairs= " << do_pairs << " do_compt = " << do_compt << endl;
      cerr << " ERROR : it is not allowed to have do_bhabhas= 1 together with do_compt = 1 " << endl;
      exit(0);
    } 
  //  if (store_secondaries && !track_secondaries)
  if (store_pairs && !track_pairs)
    {
      //cerr << " store_secondaries = " << store_secondaries << " track_secondaries = " << track_secondaries << endl;
      //cerr << " WARNING : it is not very consistent to store secondaries without tracking them! " << endl;
      cerr << " store_pairs = " << store_pairs << " track_pairs = " << track_pairs << endl;
      cerr << " WARNING : it is not very consistent to store pairs without tracking them! " << endl;
    }
  if (cuts_from_loaded_beam > 0  && load_beam == 0) 
    {
      cerr << " WARNING : the switch cuts_from_loaded_beam is without effect with load_beam = 0 " << endl;
    }
}

void SWITCHES::lectureTWOBEAM(const PARAMETERS& param)
{
  //  VALUE value;
  //  recuperer_variable("twobeam",&value);
  //  twobeam=CONTENTS(value);
  twobeam=param.readIValue("twobeam");
}

void SWITCHES::lectureCharge_sign_2(const PARAMETERS& param)
{
  //  VALUE value;
  //  recuperer_variable("charge_sign_2",&value);
  //  charge_sign_2=CONTENTS(value);
  charge_sign_2=param.readFValue("charge_sign_2");
}

string SWITCHES::output_flow() const 
{
  ostringstream sortie;
  sortie << titre(string("SWITCHES : "));
  sortie << "charge_sign = " << charge_sign << endl;
  sortie << "bmt_precession = " << bmt_precession_ << endl;
  sortie << "ST_spin_flip = " << ST_spin_flip_ << endl;
  sortie << "automatic_grid_sizing = " << automatic_grid_sizing << endl;
  sortie << "integration_method = " << integration_method << " force_symmetric = " << force_symmetric << endl;
  sortie <<  "rndm_load = " << rndm_load << "rndm_save = " << rndm_save << " rndm_seed = " << rndm_seed << endl;
  sortie << "do_photons.1 = " << do_photons_[0] << " do_photons.2 = " << do_photons_[1] << endl;
  sortie << "write_photons = " << write_photons << endl;
  sortie << "do_comp = " << do_compt << " do_prod = " << do_prod << endl;
  sortie << "electron_ratio = " << electron_ratio << endl;
  sortie << "compt_x_min = " << compt_x_min << " compt_emax = " << compt_emax << " GeV ; compt_scale = " << compt_scale << endl;
  sortie << "do_lumi = " << do_lumi << "num_lumi = " << num_lumi << " lumi_p = " << lumi_p << endl;
  sortie << "do_cross = " << do_cross << " do_isr = "<< do_isr << " do_espread = " << do_espread << endl;
  sortie << "photon_ratio = " << photon_ratio << endl;
  sortie << "do_hadrons = " << do_hadrons << " store_hadrons = " << store_hadrons << " hadron_ratio = " << hadron_ratio << endl;
  sortie << "do_jets = " << do_jets << " store_jets = " << jet_store << endl;
  sortie << "do_pairs = " << do_pairs << " do_muons = " << do_muons << " load_events = " << load_event << endl;
  //sortie << "track_secondaries = " << track_secondaries << " pair_step = " << pair_step << endl;
  //sortie << "store_secondaries = " << store_secondaries << endl;
  sortie << "track_pairs = " << track_pairs << " pair_step = " << pair_step << endl;
  sortie << "store_pairs = " << store_pairs << endl;
  sortie << "bhabha_scal = " << bhabha_scal << " bhabha_ecmload = " << bhabha_ecmload << " GeV " <<  endl;
  sortie << "do_coherent = " << do_coherent << endl;
  sortie << "do_trident = " << do_trident << endl;
  sortie << "emin = " << emin << endl;
  sortie << "grids = " << extra_grids+1 << endl;
  sortie << "pair_ecut = " << pair_ecut << " GeV " << endl;
  sortie << "pair_ratio = " << pair_ratio << endl;
  sortie << "bhabha_ratio = " << bhabha_ratio << endl;
  sortie << "pair_q2 = " << pair_q2 << endl;
  sortie << "beam_pair = " << beam_pair << endl;
  sortie << "jet_ratio = " << jet_ratio << endl;
  sortie << "jet_ptmin = " << jet_pstar << endl;
  sortie << "jet_pythia = " << jet_pythia << endl;
  sortie << "jet_log = " << jet_select << endl;
  sortie << "beam_size = " << geom << " beam_size_scale = " << r_scal << " ext_field = " << ext_field << endl;
  sortie << "espread.1 = " << espread1 << " which_espread.1 = " << which_espread1 << " espread.2 = " << espread2 << " which_espread.2 = " << which_espread2 << endl;
  sortie << "f_rep = " << f_rep << " n_b = " << n_b << endl;

  return sortie.str();
}


