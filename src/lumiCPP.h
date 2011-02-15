#ifndef LUMI_SEEN
#define LUMI_SEEN
#include <vector>
//#include "pairsCPP.h"
#include "rndmCPP.h"
#include "meshCPP.h"
#include "abstractParticle.h"
#include "fileInputOutput.h"
#include "mathematicalEntities.h"
#include "abstractIOclass.h"
using namespace std;


class LUMI_PAIR : public ABSTRACT_IO_CLASS
{
 protected:
  
  float e1_,e2_,x_,y_,z_;
  
  
  public :
    
    LUMI_PAIR() {;}
  
  LUMI_PAIR(float energy1, float energy2)
    {
      e1_ = energy1; 
      e2_ = energy2; 
    }
  
  ~LUMI_PAIR() {;}
  
  
  inline void random_position(const MESH& mesh, int cellx, int celly,float min_z, RNDM& hasard)
    {
      mesh.guess_position_in_cell(cellx, celly,min_z, x_, y_, z_, hasard);
    }
   
  
  virtual inline string name_of_class() const 
    {
      return string("LUMI_PAIR");
    }
  
  virtual string output_flow() const 
    {
      ostringstream sortie;
      sortie << " LUMI_PAIR:: no data for output file " << endl;
      return sortie.str();
    }
  
  virtual string persistent_flow() const
    {
      ostringstream sortie;
      sortie << e1_ << " " << e2_ << " " << x_ << " " << y_ << " " << z_*1e-3 << " " ;
      return sortie.str();
    }
  
  inline void get_parameters_for_output(float& e1,float& e2,float& x,float& y,float& z) const
    {
      e1 = e1_;
      e2 = e2_;
      x = x_;
      y = y_;
      z = z_*1e-3;
    }
  
  
  inline void get_impulsion_parameters(float& vx1,float& vy1,float& vx2,float& vy2, int& t) const {;}
  
};

class LUMI_PAIR_EE  : public LUMI_PAIR

{

protected :

  float vx1_,vy1_,vx2_,vy2_;
  TRIDVECTOR spin1_, spin2_;
  int t_;
  
  
 public:
  
  
  LUMI_PAIR_EE() {;}
  
  LUMI_PAIR_EE(float energy1, float energy2, float p1Vx, float p1Vy, float p2Vx, float p2Vy,int time_counter) : LUMI_PAIR(energy1, energy2)
    {
     vx1_ = p1Vx;
     vx2_ = p2Vx;
     vy1_ = p1Vy;
     vy2_ = p2Vy;
     t_ = time_counter;
    }

  ~LUMI_PAIR_EE() {;}

  inline  void set_spins(const TRIDVECTOR& s1, const TRIDVECTOR& s2 )
    {
      spin1_ = s1;
      spin2_ = s2;
    }

  void get_spins(float& sx1, float& sy1, float& sz1,float& sx2, float& sy2, float& sz2 ) const
    {
      sx1 = spin1_(0);
      sy1 =  spin1_(1);
      sz1 =  spin1_(2);
      
      sx2 =  spin2_(0);
      sy2 =  spin2_(1);
      sz2 =   spin2_(2);
    }
  

  inline void get_impulsion_parameters(float& vx1,float& vy1,float& vx2,float& vy2, int& t) const
    {
      vx1 = vx1_;
      vy1 = vy1_;
      vx2 =vx2_;
      vy2 = vy2_;
      t = t_;
    }
  
  virtual inline string name_of_class() const 
    {
      return string("LUMI_PAIR_EE");
    }
  
  
  virtual string persistent_flow() const
    {
      ostringstream sortie;
      
      sortie << LUMI_PAIR::persistent_flow();
      sortie << " " << t_ << " " << vx1_ << " " << vy1_ << " " << vx2_ << " " << vy2_ << " " <<  spin1_(0) << " " << spin1_(1) << " " << spin1_(2) << " " << spin2_(0) << " " << spin2_(1) << " " <<  spin2_(2) << " ";
      
      return sortie.str();
    }
};

class GENERAL_LUMI_HEAP : public ABSTRACT_LUMI_HEAP
{
  protected : 
    
    RNDM* hasard_;
  
  int nmax_;
  int nb_pairs_;
  float p_;
  
  public : 
    GENERAL_LUMI_HEAP()  { hasard_ = NULL;}
  
  
  GENERAL_LUMI_HEAP(int nmax, float p, RNDM* hasard)
    {
      hasard_ = hasard;
      nmax_ = nmax;
      p_ = p;
      nb_pairs_ = 0;
    }
  
  inline bool random_ok() const
    {
      return hasard_ != NULL;
    }
  
  
  virtual inline int nb_pairs() const {return nb_pairs_;}
  
  
  inline int numberToStore(float store) const
    {
      int nstore;
      nstore=(int)store;
      store -= nstore;
      if (hasard_->rndm() <= store) nstore++;
      return nstore;
    }
  
  int stack_vector(int nstore, vector<int>& selectionned_indices);
  
  
};

class LUMI_HEAP : public GENERAL_LUMI_HEAP
{
  protected :
    
    vector<LUMI_PAIR> data_;
  
 public:
  
  LUMI_HEAP(): GENERAL_LUMI_HEAP() {;};
  
  LUMI_HEAP( int nmax, float p, RNDM* hasard) : GENERAL_LUMI_HEAP( nmax,p, hasard)
    {
      data_.reserve(nmax_);
    }
  
  virtual ~LUMI_HEAP() {;}

  void lumi_store(const MESH& mesh, int cellx, int celly,float min_z, float energy1, float energy2, float weight);
  
  
  
  virtual inline void get_parameters_for_output(int numero, float& e1,float& e2,float& x,float& y,float& z) const
   {
     data_[numero].get_parameters_for_output(e1,e2,x,y,z);
   }
  
  virtual inline void get_parameters_for_output(int numero, float& e1,float& e2,float& x,float& y,float& z, float& vx1,float& vy1,float& vx2,float& vy2, float& sx1, float& sy1, float& sz1, float& sx2, float& sy2, float& sz2, int& t) const {;}
  

  virtual inline string name_of_class() const 
    {
      return string("LUMI_HEAP");
    }

  
  virtual string persistent_flow() const
    {
      //unsigned long int k;
      int k;
      ostringstream sortie;
      vector<unsigned long int> ordre;
      //unsigned long int npairs  =  data_.size();
      int npairs = (int) data_.size();
      // vefif a supprimer 
      if ( npairs != nb_pairs_) 
	{
	  cerr << " LUMI_HEAP_EE problem avec npombre de paires npairs = " << npairs << " nb_pairs_ = " << nb_pairs_ << endl;
	  exit(0);
	}
      
      hasard_->getShuffledIntegerSequence(npairs, ordre);
      for (k = 0; k < npairs; k++)
	{
	  sortie << data_[ordre[k] - 1 ].persistent_flow() << " " << ordre[k] << endl;
	}
      return sortie.str();
    }
  
  
  
  inline void saveLumi(string nameOfOutputFile) const
    {
      FILE_IN_OUT filout;
      filout.open_file(nameOfOutputFile, "w");
      filout.save_lumi_heap(this);
      filout.close();
    }
  
  int  prepare_store(float weight);
  
  
};

class LUMI_HEAP_EE : public GENERAL_LUMI_HEAP
{
 protected: 
  
  vector<LUMI_PAIR_EE> data_;
  
  public :
    
    LUMI_HEAP_EE() : GENERAL_LUMI_HEAP() {;}
  
  LUMI_HEAP_EE( int nmax, float p, RNDM* hasard) : GENERAL_LUMI_HEAP( nmax,p, hasard)
    {
      data_.reserve(nmax_);
    }
  
  
  virtual inline void get_parameters_for_output(int numero, float& e1,float& e2,float& x,float& y,float& z) const {;}
  
  virtual inline void get_parameters_for_output(int numero, float& e1,float& e2,float& x,float& y,float& z, float& vx1,float& vy1,float& vx2,float& vy2, float& sx1, float& sy1, float& sz1, float& sx2, float& sy2, float& sz2,int& t) const
    {
      data_[numero].get_parameters_for_output(e1,e2,x,y,z);
      data_[numero].get_impulsion_parameters(vx1,vy1,vx2,vy2,t);
      data_[numero].get_spins(sx1, sy1, sz1, sx2, sy2, sz2);
      /*      cout << " get_parameters_for_output les spins :  " << endl; */
      /*      cout << sx1 << " " << sy1 << " " << sz1 << endl; */
      /*      cout << sx2 << " " << sy2 << " " << sz2 << endl; */
    }
  
  virtual inline string name_of_class() const 
    {
  return string("LUMI_HEAP_EE");
    }
  

  virtual string persistent_flow() const
    {
      //unsigned long int k;
      int k;
      ostringstream sortie;
      vector<unsigned long int> ordre;
      //unsigned long int npairs  =  data_.size();
      int npairs = (int) data_.size();
      
      // vefif a supprimer 
      if ( npairs != nb_pairs_) 
	{
	  cerr << " LUMI_HEAP_EE problem avec npombre de paires npairs = " << npairs << " nb_pairs_ = " << nb_pairs_ << endl;
	  exit(0);
	}
      hasard_->getShuffledIntegerSequence(npairs, ordre);
      for (k = 0; k < npairs; k++)
	{
	  sortie << data_[ ordre[k] - 1 ].persistent_flow() << " " << ordre[k] << endl;
	}
      return sortie.str();
    }
  
  
  inline void saveLumi(string nameOfOutputFile) const
    {
      FILE_IN_OUT filout;
      filout.open_file(nameOfOutputFile, "w");
      //    filout.save_lumi_heap_full(this);
      filout.save_lumi_heap(this);
      filout.close();
    }
  
  
  void lumi_store_ee(const MESH& mesh, int cellx, int celly,float min_z, float energy1, float p1Vx, float p1Vy, float energy2, float p2Vx, float p2Vy,float weight,int time_counter);
  
  void lumi_store_ee(const MESH& mesh, int cellx, int celly,float min_z, float energy1, float p1Vx, float p1Vy, float energy2, float p2Vx, float p2Vy,float weight, const TRIDVECTOR& s1, const TRIDVECTOR& s2, int time_counter); 

  int  prepare_store(float weight);
  
  
};



#endif
