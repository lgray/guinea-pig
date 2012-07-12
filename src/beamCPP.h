#ifndef BEAM_SEEN
#define BEAM_SEEN


#include "beamParametersCPP.h"
#include "particleBeamCPP.h"
#include "mathematicalEntities.h"
#include "fileInputOutput.h"

class  BEAM : public ABSTRACT_BEAM, public ABSTRACT_IO_CLASS
{
  PARTICLE_BEAM particle_beam_;
  PHOTON_BEAM photon_beam_;
  BEAM_PARAMETERS* beam_parameters_;
  int beam_label_;
  int sign_label_;
  FILE_IN_OUT* size_log_file_;
  
 public:

  BEAM() : beam_parameters_(NULL), sign_label_(0), size_log_file_(NULL) {;};

  ~BEAM()
    {
      if (size_log_file_ != NULL)
	{
	  size_log_file_->close();
	  delete size_log_file_;
	}
    }

  inline virtual string output_flow() const
    {
      ostringstream out;
      double esum;
      double sumphot;
      int numphot;
      ostringstream aux;
      aux <<  beam_label_ ;
      string start("beam");
      start.append(aux.str());
      out << title(start);
      out << particle_beam_.output_flow();
      esum = particle_beam_.meanLostEnergy(beam_parameters_->ebeam());
      out << "average energy loss of the beam particles (GeV) : de = " <<  esum << endl;

      out << " ---photon statistics : " << endl;
      out << " initial : " << endl;
      double  nb_part = (double)particle_beam_.numberOfParticles();
      PHOTON_COUNTER pc = photon_beam_.get_photon_counter();
      string photstat;
      photstat+=" initial average photon energy (GeV) phot-e";
      photstat+=aux.str();
      photstat+=" : ";
      out << photstat <<  pc.getSum()/max(1,pc.getNumber()) << endl;
      photstat.erase();
      photstat+=" initial number of phot. per tracked macropart.";
      photstat+=aux.str();
      photstat+=" : ";
      out << photstat << pc.getNumber()/nb_part << endl;
      out << " final : " << endl;
      photon_beam_.photon_info(sumphot,numphot);
      photstat.erase();
      photstat+=" final average photon energy (GeV) phot-e";
      photstat+=aux.str();
      photstat+=" : ";
      out << photstat << sumphot/max(1,numphot) <<endl;
      photstat.erase();
      photstat+=" final number of phot. per tracked macropart.";
      photstat+=aux.str();
      photstat+=" : ";
      out <<  photstat << numphot/nb_part << endl << endl;
      // final quantities for post processing
      if ( beam_label_ == 1) {
	out << "phot-e1=" << sumphot/max(1,numphot) << ";n_phot1=" <<  numphot/nb_part << ";" << endl;
	out << "de1=" << esum <<  ";" << endl;
      } else if ( beam_label_ == 2 ) {
	out << "phot-e2=" << sumphot/max(1,numphot) << ";n_phot2=" <<  numphot/nb_part << ";" << endl;
	out << "de2=" << esum <<  ";" << endl;
      }
      return out.str();
    } 

  inline const PARTICLE_BEAM& particle_beam() const { return particle_beam_;}
  inline const PHOTON_COUNTER& get_photon_counter() const { return photon_beam_.get_photon_counter();}
  
  inline void connect_parameters(BEAM_PARAMETERS* beamParameters)
    {
      beam_parameters_ = beamParameters;
      beam_label_ = beam_parameters_->label();
      if (beam_label_ == 1) sign_label_ = 1;
      else if (beam_label_ == 2) sign_label_ = -1;
      else
	{
	  cerr << " BEAM:: connect_parameters WARNING : bad beam label = " << beam_label_ << endl;
	}
    }

  inline void transverse_sigmas(float& sigX, float& sigY) const 
    {
      particle_beam_.transverse_sigmas(sigX, sigY);
    }
  
  inline float sigma_xyz_from_data(char dir) const
    {
      if (dir == 'x') return beam_parameters_->sigma_x();
      else if (dir == 'y') return beam_parameters_->sigma_y();
      else if (dir == 'z') return beam_parameters_->sigma_z();
      else 
	{
	  cerr << "BEAM::sigma_xyz : ERROR unknown direction for sigma : " << dir << endl;
	  exit(0);
	}
    }
  
  inline int label() const { return beam_label_;}
  inline int sign_label() const { return sign_label_;}
  inline void write_size_init(string filename)
    {
      ostringstream nom;
      int label = beam_parameters_->label();
      nom << filename << label << ".dat" << ends;
      string nomfic = nom.str();
      if (size_log_file_ != NULL)
	{
	  cerr << " BEAM::write_size_init : WARNING : only one output beamsize.dat file is allowed " << endl;
	  exit(0);
	}
      size_log_file_ = new FILE_IN_OUT();
      size_log_file_->open_file(nomfic, "w");
      nom.seekp(0);
      nom << string("#slice x_rms y_rms xmin xmax xmean ymin ymax ymean") << endl;
      string out = nom.str();
      size_log_file_->write_line(out);
    }

  
  void write_size(int first_slice, int last_slice);
  
  
  inline int numberOfParticles() const {return particle_beam_.numberOfParticles();}
  
  
  inline void make_beam(int n_slice, int bmt_rotate, RNDM* rndm_generator) 
    {
      photon_beam_ = PHOTON_BEAM(n_slice);
      particle_beam_ = PARTICLE_BEAM(n_slice,bmt_rotate,beam_parameters_->get_polar(), rndm_generator);
    }
  
  
  inline void newCoherent(int slice,float x,float y, float vx,float vy, float energy)
    {
      particle_beam_.newCoherent(slice, x, y, vx, vy, energy);
    }
  
  inline void newTrident(int slice,float x,float y, float vx,float vy, float energy)
    {
      particle_beam_.newTrident(slice, x, y, vx, vy, energy);
    }

  inline void rotate_particles()
    {
      particle_beam_.rotate_particles(beam_parameters_->phi_angle());
    }
  
  inline void set_angle_particles(float delta_z)
    {
      particle_beam_.set_angle_particles(beam_parameters_->x_angle(), beam_parameters_->y_angle(), delta_z);
    }
  
  inline void set_particles_offset()
    {
      particle_beam_.set_particles_offset( beam_parameters_->offset_x(), 
					   beam_parameters_->offset_y(), 
					   beam_parameters_->waist_x(),
					   beam_parameters_->waist_y());
    }
  
  
  inline void adjustToMeanPosition()
    {
      particle_beam_.adjustToMeanPosition();
    }
  
  inline void backstep (int beam, float max_z, float step,  int timestep, float scal_step[2])
    {
      particle_beam_.backstep(beam, beam_parameters_->trav_focus(), max_z, step,  timestep, scal_step);
    }
  
  inline void backstep2 (int nbeam, float max_z, float step,  int timestep, float scal_step[2])
    {
      particle_beam_.backstep2(nbeam, max_z, step, timestep, scal_step);
    }
  
  inline unsigned int load_particles(BEAM_FROM_FILE*& bff, float emin, float zmin, float deltaz, float sigx, float sigy, float sigz) 
    {
      return particle_beam_.load_particles(bff, emin, zmin, deltaz, sigx, sigy, sigz);
    }
  
  
  inline void load_photons(string filename, int type_of_beam, float delta_z, float max_z, int n_cell_z)
    {
      photon_beam_.load_photons(filename, type_of_beam, delta_z, max_z, n_cell_z);
    }
  
  
  inline void init_particles(int nbpart, float delta_z, int charge_symmetric)
    {
      particle_beam_.init_particles (nbpart,  beam_parameters_->sigma_x(),
				     beam_parameters_->sigma_y(), 
				     beam_parameters_->sigma_z(),
				     beam_parameters_->dist_x(),
				     beam_parameters_->dist_z(),
				     beam_parameters_->em_x(), 
				     beam_parameters_->em_y(),
				     delta_z,
				     beam_parameters_->ebeam(),
				     charge_symmetric);
    }
  
  inline float get_ebeam() const { return beam_parameters_->ebeam();};
  
  inline void dump_beam(string name, int istep, int every_particle, int timestep, float step, float max_z)
    {
      particle_beam_.dump_beam(name, istep, every_particle, timestep, step, max_z,  sign_label_);
    }
  
  inline void dump_photons(string name,int istep, int every_particle,int timestep, float step, float max_z)
  {
    photon_beam_.dump_photons(name, istep,every_particle, timestep, step, max_z,sign_label_);
  }
 
  
  inline void meanPositionOfSlice(int slice, float& x,float& y) const
    {
      particle_beam_.meanPositionOfSlice(slice,x, y);
    }
  inline void beamXyRms(float& x0, float& y0,  float& sigmax, float& sigmay) const
    {
      particle_beam_.beamXyRms(x0, y0, sigmax, sigmay);
    }
  
  inline void beamZRms(float& z0,  float& sigmaz) const
    {
      particle_beam_.beamZRms(z0, sigmaz);
    }
  
  inline void emittances(float& emittx, float& emitty) const
    {
      particle_beam_.emittances(emittx, emitty);
    }
  
  inline float gamma() const { return particle_beam_.gamma();}
  
  virtual inline int numberOfParticlesOfSlice(int slice) const {return particle_beam_.numberOfParticlesOfSlice(slice);};

  inline void get_named_slices_vector(named_int_vector& vec, string name) const 
    {
      int k;
      vec.clear();
      vec.put_name(name);
      for (k=0; k< number_of_slices(); k++)
	{
	  vec.add_component(particle_beam_.numberOfParticlesOfSlice(k));
	}
      
    }
  
  inline  vector<PARTICLE*>& getCoherentVector(int slice) {return  particle_beam_.getCoherentVector(slice);}
  inline  vector<PARTICLE*>& getTridentVector(int slice) {return  particle_beam_.getTridentVector(slice);}
  inline  const vector<PARTICLE*>& getCoherentVector(int slice) const {return  particle_beam_.getCoherentVector(slice);}
  inline  const vector<PARTICLE*>& getTridentVector(int slice) const {return  particle_beam_.getTridentVector(slice);}
 
  inline const vector<PARTICLE*>& getParticleVector(int slice) const {return particle_beam_.getParticleVector(slice);}
  inline  vector<PARTICLE*>& getParticleVector(int slice)  {return particle_beam_.getParticleVector(slice);}
  
  inline int number_of_coherent_vectors() const {return particle_beam_.number_of_coherent_vectors();}
  virtual   inline int number_of_slices() const { return particle_beam_.number_of_slices();}
  
  virtual inline int sizeOfPhotonSlice(int slice) const
    {
      return photon_beam_.sizeOfSlice(slice);
    }
  
  inline vector<PHOTON>& getPhotonVector(int slice) { return photon_beam_.getPhotonVector(slice);}
  inline const vector<PHOTON>& getPhotonVector(int slice) const { return photon_beam_.getPhotonVector(slice);}
  
  
  inline int store_beam(string name) const 
    {
      return particle_beam_.store_beam(name);
    }
  inline void store_coherent_beam(string name) const 
    {
      particle_beam_.store_coherent_beam(name);
    }
  inline void store_trident_beam(string name) const 
    {
      particle_beam_.store_trident_beam(name);
    }


  inline void ang_dis(unsigned int n_bin, vector< vector<float> >& bin ) const {particle_beam_.ang_dis(n_bin, bin);} 
  
  inline const PARTICLE_BEAM& particle_beam() { return particle_beam_;}
  
  inline const PHOTON& new_photon(float energy, PARTICLE& particle, int slice) 
    {
      
      return photon_beam_.new_photon(energy, particle, slice);
      
    }

  inline void move_photon_beam(int i_slice, float delta)
    {
      photon_beam_.move_photons(i_slice, delta);
    }
  
  inline void photon_info(double& sum, int& number) const
    {
      photon_beam_.photon_info(sum,number);
    }
  
};




#endif
