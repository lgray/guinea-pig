#ifndef FILEOUTPUTASCII_SEEN
#define FILEOUTPUTASCII_SEEN

#include "IfileInputOutput.h"
#include <iostream>
#include <fstream>
#include <iomanip>

class FILE_IN_OUT_ASCII : public IFILE_IN_OUT
{
      ifstream infile_;
      ofstream outfile_;

      public :

	FILE_IN_OUT_ASCII()  {;}
      virtual ~FILE_IN_OUT_ASCII() {;}
      
      virtual void open_file(string name, const char* mode)
	{  
	  switch (mode[0])
	    {
	      // lecture 
	    case 'r':
	      infile_.open(name.c_str(), ios::in);
	      if (!infile_) 
		{
		  cerr << " error opening input stream " << name << endl;
		  exit(0);
		}
	      break;
	      // ecriture
	    case 'w':
	      
	      outfile_.open(name.c_str(), ios::out);
	      if (!outfile_) 
		{
		  cerr << " error opening output stream " << name << endl;
		  exit(0);
		}
	      break;
	    default:
	      cerr << " FILE_IN_OUT_ASCII::open_file: unknown open mode : " << mode[0] << endl;
	      exit(0);
	    }
	}

      virtual void set_header(string head)
	{  
	  outfile_ << head << endl;
	}
      
      virtual void close()
	{
	  infile_.close();
	  outfile_.close();
	}
 
      virtual bool read_line(string& ligne) 
	{
	  return getline(infile_, ligne);
	}
      
      virtual void write_line(string& ligne) 
	{
	  outfile_ << ligne << endl;
	}
      virtual void set_jet_header(float ptmin, float sqrt_s)
	{
	  outfile_ << ptmin << " " << sqrt_s << endl;
	}
      
      virtual void  save_jet(float energy1, float energy2, int process) 
	{   
	  outfile_ << energy1 << " "  << energy2 << " "  << process << endl;
	}
      virtual void  save_jet(float energy1, float energy2, float pz1, float pz2, float pt, int process) 
	{   
	  outfile_ << energy1 << " "  << energy2 << " "  << pz1 << " " << pz2 << " " << pt << " " << process << endl;
	}
      
      virtual void  save_hadron(float energy1, float energy2, float z) 
	{   
	  outfile_ << energy1 << " "  << energy2 << " "  << z << endl;
	}

      /*  virtual void  save_photon(const ABSTRACT_PARTICLE& part, int no_beam)  */
      /*   {    */
      /*     float  energy, xpos, ypos, zpos, vx,vy, dummy; */
      /*     part.get_parameters_for_output(energy, xpos, ypos, zpos, vx,vy, dummy, dummy, dummy); */
      /*     if (no_beam != 1) energy = -energy; */
      /*     outfile_ << energy << " "  << vx << " "  << vy << endl; */
      /*   } */

      virtual void save_compton_photon(float y, float px, float py)
	{
	  outfile_ << y << " " << px << " " << py << endl;
	}
      
      /* virtual void  save_particle(const ABSTRACT_PARTICLE& part)  */
      /*   {    */
      /*     float  energy, xpos, ypos, zpos, vx,vy, spinx, spiny, spinz ; */
      
      /*     part.get_parameters_for_output(energy, xpos, ypos, zpos, vx,vy,spinx, spiny, spinz ); */
      /*     outfile_ << energy << " "  << xpos << " "  << ypos <<" "  <<  zpos << " "  << vx << " "  << vy << " "  << spinx << " "  << spiny << " "  << spinz << endl; */
      /*   } */
      
      virtual bool  read_particle(PARTICLE_INTERFACE& part) 
	{   
	  int k;
	  // float  energy, xpos, ypos, zpos, vx,vy ;
	  float readValue[9];
	  for ( k =0; k < 9; k++) readValue[k] = 0.0;
	  bool test = false;
	  string slu;
	  int number = 0;
	  if ( getline(infile_, slu) )
	    {
	      bool goodline = false;
	      istringstream ss(slu);
	      float aux;
	      while ( ss >> aux && number < 9) 
		{
		  goodline = true;
		  readValue[number] = aux;
		  number++;
		}
	      part.init_from_input_file(readValue[0],readValue[1],readValue[2],readValue[3], readValue[4], readValue[5], readValue[6], readValue[7], readValue[8] );
	      if ( !goodline) cerr << " read_particle : reading failure " << endl;
	      test = true;
	    }
	  return test;
	}
      
      virtual bool  read_photon(PARTICLE_INTERFACE& part) 
	{   
	  float  energy, xpos, ypos, zpos, vx,vy, hel ;
	  float dumy = 0.0;
	  float dumz = 0.0;
	  bool test;
	  if ( infile_ >> energy >> xpos >> ypos >> zpos >> vx >> vy >> hel) 
	    {
	      part.init_from_input_file(energy, xpos, ypos, zpos, vx, vy, hel, dumy,dumz);
	      test= true;
	    }
	  else test =  false;
	  return test;
	}
      
      /* virtual void save_pair_particle(const ABSTRACT_PAIR_PARTICLE& pair_part)  */
      /*   { */
      /*     float x,y,z,vx,vy,vz, energy; */
      /*     pair_part.get_parametersX(x,y,z,vx,vy,vz,energy); */
      /*     //   outfile_ << energy << " "  << vx << " "  << vy << " "  << vz << " "  << x << " "  << y << " "  << z << " "  << pair_part.get_process() <<  " "  << pair_part.get_label() << endl; */
      /*        outfile_ << energy << " "  << vx << " "  << vy << " "  << vz <<  " "  << pair_part.get_process() <<  " "  << pair_part.get_label() << endl; */
      /*   } */
 
      virtual void save_bhabhasamples(const ABSTRACT_BHABHASAMPLES* const bhabhas)  
	{
	  
	  unsigned int k, evtIndex;
	  float eCM, mother1_en, e1, vx1, vy1, vz1, mother2_en, e2, vx2, vy2, vz2;
	  int nbphot;
	  for (k=0; k < bhabhas->nb_samples(); k++)
	    {
	      bhabhas->get_parameters_for_output(k, evtIndex, eCM, mother1_en, e1, vx1, vy1, vz1, mother2_en, e2, vx2, vy2, vz2, nbphot);
	      outfile_ << setw(8) << evtIndex << " " << setw(10) << fixed << setprecision(4) << eCM << " ";
	      outfile_ << setw(10) << setprecision(4) << mother1_en << " " << setw(10) << e1 << " "  << setprecision(6) << setw(10) << vx1 << " "  << setw(10) << vy1 << " "  << setw(10) << vz1 << " ";
	      outfile_ << setw(10) << setprecision(4) << mother2_en << " " << setw(10) << e2 << " "  << setprecision(6) << setw(10) << vx2 << " "  << setw(10) << vy2 << " "  << setw(10) << vz2 << " ";
	      outfile_ << setw(4) << nbphot << endl;
	    }
	}
      
      
  virtual bool read_bhabhasamples(ABSTRACT_BHABHASAMPLES* const bhabhas)
	{
	  float px1, py1, pz1, e1, px2, py2, pz2, e2;
	  unsigned int evtIdx, nbphot;
	  int i=0;
	  while (  infile_ >> evtIdx >> px1 >> py1 >> pz1 >> e1 >> px2 >> py2 >> pz2 >> e2 >> nbphot)
	    {
/*		  cout << "Read a new bhabha\n";
		  cout << "P1 = (" << px1 << ", " << py1 << ", " << pz1 << ", " << e1 << ")\n";
		  cout << "P2 = (" << px2 << ", " << py2 << ", " << pz2 << ", " << e2 << ")\n";
		  cout << "m1**2 = " << e1*e1 - px1*px1 - py1*py1 - pz1*pz1 << endl;
		  cout << "m2**2 = " << e2*e2 - px2*px2 - py2*py2 - pz2*pz2 << endl;
*/	      bhabhas->add_bhabha(evtIdx, px1, py1, pz1, e1, px2, py2, pz2, e2, nbphot);
	      i++;
	    }
	  return true;
	}
      
      virtual void save_bhabhaPhotonSamples(const ABSTRACT_BHABHA_PHOTON_SAMPLES*  bhabhaPhot)  
	{
	  
	  unsigned int k;
	  float en, vx, vy, vz;
	  int evtIndex;
	  for (k=0; k < bhabhaPhot->nb_samples(); k++)
	    {
	      bhabhaPhot->get_parameters_for_output(k,evtIndex, en, vx, vy, vz);
	      outfile_ << setw(8) << k+1 << " " << setw(8) << evtIndex << " " << setw(10) << fixed << setprecision(4) << en << " "  << setw(10) << vx << " "  << setw(10) << vy << " "  << setw(10) << vz << endl;
	    }
	}
      
    virtual bool read_bhabhaPhotonsamples(ABSTRACT_BHABHA_PHOTON_SAMPLES* const bhabhasPhoton)
	{
      int evtIdx;
	  float px, py, pz, en;
	  //int nbphot;
	  while (infile_ >> evtIdx >> px >> py >> pz >> en)
	    {
	      bhabhasPhoton->add_bhabha_photon(evtIdx, px, py, pz, en);
	    }
	  return true;
	}
      
      
      virtual int read_pythia_file(int& logx, int& logy, vector<double>& x, vector<double>& y)
	{
	  int k;
	  int nentries;
	  x.clear();
	  y.clear();
	  if (infile_ >> nentries >> logx >> logy )
	    {
	      x.resize(nentries);
	      y.resize(nentries);     
	      for (k = 0; k < nentries;k++)
		{
		  if ( !(infile_ >> x[k] >> y[k]) )
		    {
		      cerr << " FILE_IN_OUT_ASCII::read_pythia_file : the list of pythia values is interrupted " << endl;
		      exit(0);
		    }
		}
	    }
	  else
	    {
	      cerr << " FILE_IN_OUT_ASCII::read_pythia_file error reading pythia file " << endl;
	      exit(0);
	    }
	  return nentries; 
	}
      
      virtual bool read_cross(ABSTRACT_CROSS_DATA* const crossIni)
	{
	  int k;
	  int n, nval, dummy1, dummy2;
	  float energy;
	  float* values = NULL;
	  bool test = false;
	  if (infile_ >> n >> nval >> dummy1 >> dummy2)
	    {
	      crossIni->resize(n, nval);
	      values = new float[nval];
	      while (  infile_ >> energy )
		{
		  for (k=0; k< nval; k++)
		    {
		      if (!(infile_ >> values[k]))
			{
			  cerr << " FILE_IN_OUT_ASCII::read_cross : the list of cross values is interrupted " << endl;
			  exit(0);
			}
		    }
		  crossIni->add_data(energy, values);
		}
	      if (values != NULL) delete [] values;
	      test = true;
	    }
	  else
	    {
	      cerr << " FILE_IN_OUT_ASCII::read_cross error reading cross values from file " << endl;
	      exit(0);
	    }
	  return test; 
	}


      virtual void save_lumi_heap(const ABSTRACT_LUMI_HEAP* const lumi_heap) 
	{
	  outfile_ << lumi_heap->persistent_flow() << endl;
	}
      
      virtual void save_object_on_output_listing(const ABSTRACT_IO_CLASS* const obj)
	{
	  outfile_ << obj->output_flow();
	}
      
      virtual void save_object_on_persistent_file(const ABSTRACT_IO_CLASS* const obj)
  	{
  	  outfile_ << obj->persistent_flow() << endl;
  	}

      virtual void save_object_on_persistent_file(const int eventIndex, const ABSTRACT_IO_CLASS* const obj)
  	{
  	  outfile_ << eventIndex << " " << obj->persistent_flow() << endl;
  	}

};



#endif
