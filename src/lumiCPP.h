#ifndef LUMI_SEEN
#define LUMI_SEEN

#include <sstream>
#include <string>
#include <vector>
//#include "pairsCPP.h"
#include "rndmCPP.h"
#include "meshCPP.h"
#include "abstractParticle.h"
#include "fileInputOutput.h"
#include "mathematicalEntities.h"
#include "abstractIOclass.h"

class LUMI_PAIR : public ABSTRACT_IO_CLASS
{
 protected:
  
  double e1_,e2_,x_,y_,z_;
  
 public :
    
 LUMI_PAIR() : e1_(0.0),e2_(0.0),x_(0.0),y_(0.0),z_(0.0) {;}
  
 LUMI_PAIR(double energy1, double energy2) : x_(0.0),y_(0.0),z_(0.0)
    {
      e1_ = energy1; 
      e2_ = energy2; 
    }
  
  ~LUMI_PAIR() {;}
  
  inline void random_position(const MESH& mesh, int cellx, int celly,double min_z, RNDM& rndm_generator)
    {
      mesh.guess_position_in_cell(cellx, celly,min_z, x_, y_, z_, rndm_generator);
    }
   
  virtual std::string output_flow() const 
    {
      std::ostringstream out;
      out << " LUMI_PAIR:: no data for output file " << std::endl;
      return out.str();
    }
  
  virtual std::string persistent_flow() const
  {
    std::ostringstream out;
    out << e1_ << " " << e2_ << " " << x_ << " " << y_ << " " << z_*1e-3 << " " ;
    return out.str();
  }
  
  inline void get_parameters_for_output(double& e1,double& e2,double& x,double& y,double& z) const
    {
      e1 = e1_;
      e2 = e2_;
      x = x_;
      y = y_;
      z = z_*1e-3;
    }
  
  
  inline void get_impulsion_parameters(double& /*vx1*/,double& /*vy1*/,double& /*vx2*/,double& /*vy2*/, int& /*t*/) const {;}
  
};

class LUMI_PAIR_EE  : public LUMI_PAIR

{

protected :

  double vx1_,vy1_,vx2_,vy2_;
  TRIDVECTOR spin1_, spin2_;
  int t_;
  
  
 public:
  
  
 LUMI_PAIR_EE() : vx1_(0.0),vy1_(0.0),vx2_(0.0),vy2_(0.0),t_(0) {;}
  
  LUMI_PAIR_EE(double energy1, double energy2, double p1Vx, double p1Vy, double p2Vx, double p2Vy,int time_counter) : LUMI_PAIR(energy1, energy2)
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

  void get_spins(double& sx1, double& sy1, double& sz1,double& sx2, double& sy2, double& sz2 ) const
    {
      sx1 = spin1_(0);
      sy1 =  spin1_(1);
      sz1 =  spin1_(2);
      
      sx2 =  spin2_(0);
      sy2 =  spin2_(1);
      sz2 =   spin2_(2);
    }
  

  inline void get_impulsion_parameters(double& vx1,double& vy1,double& vx2,double& vy2, int& t) const
    {
      vx1 = vx1_;
      vy1 = vy1_;
      vx2 =vx2_;
      vy2 = vy2_;
      t = t_;
    }
  
  virtual std::string persistent_flow() const
    {
      std::ostringstream out;
      
      out << LUMI_PAIR::persistent_flow();
      out << " " << t_ << " " << vx1_ << " " << vy1_ << " " << vx2_ << " " << vy2_ << " " <<  spin1_(0) << " " << spin1_(1) << " " << spin1_(2) << " " << spin2_(0) << " " << spin2_(1) << " " <<  spin2_(2) << " ";
      
      return out.str();
    }
};

class GENERAL_LUMI_HEAP : public ABSTRACT_LUMI_HEAP
{
  protected : 
    
    RNDM* rndm_generator_;
  
  int nmax_;
  int nb_pairs_;
  double p_;
  
  public : 
 GENERAL_LUMI_HEAP() : nmax_(0),nb_pairs_(0),p_(0.0) { rndm_generator_ = NULL;}
  
  
  GENERAL_LUMI_HEAP(int nmax, double p, RNDM* rndm_generator)
    {
      rndm_generator_ = rndm_generator;
      nmax_ = nmax;
      p_ = p;
      nb_pairs_ = 0;
    }
  
  inline bool random_ok() const
    {
      return rndm_generator_ != NULL;
    }
  
  
  virtual inline int nb_pairs() const {return nb_pairs_;}
  
  
  inline int numberToStore(double store) const
    {
      int nstore;
      nstore=(int)store;
      store -= nstore;
      if (rndm_generator_->rndm() <= store) nstore++;
      return nstore;
    }
  
  int stack_vector(int nstore, std::vector<int>& selected_indices);
  
  
};

class LUMI_HEAP : public GENERAL_LUMI_HEAP
{
  protected :
    
    std::vector<LUMI_PAIR> data_;
  
 public:
  
  LUMI_HEAP(): GENERAL_LUMI_HEAP() {;};
  
  LUMI_HEAP( int nmax, double p, RNDM* rndm_generator) : GENERAL_LUMI_HEAP( nmax,p, rndm_generator)
    {
      data_.reserve(nmax_);
    }
  
  virtual ~LUMI_HEAP() {;}

  void lumi_store(const MESH& mesh, int cellx, int celly,double min_z, double energy1, double energy2, double weight);
  
  
  
  virtual inline void get_parameters_for_output(unsigned int number, double& e1,double& e2,double& x,double& y,double& z) const
   {
     data_[number].get_parameters_for_output(e1,e2,x,y,z);
   }
  
  virtual inline void get_parameters_for_output(unsigned int /*number*/, double& /*e1*/,double& /*e2*/,double& /*x*/,double& /*y*/,double& /*z*/, double& /*vx1*/,double& /*vy1*/,double& /*vx2*/,double& /*vy2*/, double& /*sx1*/, double& /*sy1*/, double& /*sz1*/, double& /*sx2*/, double& /*sy2*/, double& /*sz2*/, int& /*t*/) const {;}
  
  virtual std::string persistent_flow() const
    {
      //unsigned long int k;
      int k;
      std::ostringstream out;
      std::vector<unsigned long int> order;
      //unsigned long int npairs  =  data_.size();
      int npairs = (int) data_.size();
      // check 
      if ( npairs != nb_pairs_) 
	{
	  std::cerr << " LUMI_HEAP_EE problem with the number of pairs: npairs = " << npairs << " nb_pairs_ = " << nb_pairs_ << std::endl;
	  exit(0);
	}
      
      rndm_generator_->getShuffledIntegerSequence(npairs, order);
      for (k = 0; k < npairs; k++)
	{
	  out << data_[order[k] - 1 ].persistent_flow() << " " << order[k] << std::endl;
	}
      return out.str();
    }
  
  
  
  inline void saveLumi(std::string nameOfOutputFile) const
    {
      FILE_IN_OUT filout;
      filout.open_file(nameOfOutputFile, "w");
      filout.save_lumi_heap(this);
      filout.close();
    }
  
  int  prepare_store(double weight);
  
  
};

class LUMI_HEAP_EE : public GENERAL_LUMI_HEAP
{
 protected: 
  
  std::vector<LUMI_PAIR_EE> data_;
  
  public :
    
    LUMI_HEAP_EE() : GENERAL_LUMI_HEAP() {;}
  
  LUMI_HEAP_EE( int nmax, double p, RNDM* rndm_generator) : GENERAL_LUMI_HEAP( nmax,p, rndm_generator)
    {
      data_.reserve(nmax_);
    }
  
  
  virtual inline void get_parameters_for_output(unsigned int /*number*/, double& /*e1*/,double& /*e2*/,double& /*x*/,double& /*y*/,double& /*z*/) const {;}
  
  virtual inline void get_parameters_for_output(unsigned int number, double& e1,double& e2,double& x,double& y,double& z, double& vx1,double& vy1,double& vx2,double& vy2, double& sx1, double& sy1, double& sz1, double& sx2, double& sy2, double& sz2,int& t) const
    {
      data_[number].get_parameters_for_output(e1,e2,x,y,z);
      data_[number].get_impulsion_parameters(vx1,vy1,vx2,vy2,t);
      data_[number].get_spins(sx1, sy1, sz1, sx2, sy2, sz2);
      /*      std::cout << " get_parameters_for_output, the spins :  " << std::endl; */
      /*      std::cout << sx1 << " " << sy1 << " " << sz1 << std::endl; */
      /*      std::cout << sx2 << " " << sy2 << " " << sz2 << std::endl; */
    }
  
  virtual std::string persistent_flow() const
    {
      //unsigned long int k;
      int k;
      std::ostringstream out;
      std::vector<unsigned long int> order;
      //unsigned long int npairs  =  data_.size();
      int npairs = (int) data_.size();
      
      // check 
      if ( npairs != nb_pairs_) 
	{
	  std::cerr << " LUMI_HEAP_EE problem with the number of pairs: npairs = " << npairs << " nb_pairs_ = " << nb_pairs_ << std::endl;
	  exit(0);
	}
      rndm_generator_->getShuffledIntegerSequence(npairs, order);
      for (k = 0; k < npairs; k++)
	{
	  out << data_[ order[k] - 1 ].persistent_flow() << " " << order[k] << std::endl;
	}
      return out.str();
    }
  
  
  inline void saveLumi(std::string nameOfOutputFile) const
    {
      FILE_IN_OUT filout;
      filout.open_file(nameOfOutputFile, "w");
      //    filout.save_lumi_heap_full(this);
      filout.save_lumi_heap(this);
      filout.close();
    }
  
  void lumi_store_ee(const MESH& mesh, int cellx, int celly,double min_z, double energy1, double p1Vx, double p1Vy, double energy2, double p2Vx, double p2Vy,double weight,int time_counter);
  
  void lumi_store_ee(const MESH& mesh, int cellx, int celly,double min_z, double energy1, double p1Vx, double p1Vy, double energy2, double p2Vx, double p2Vy,double weight, const TRIDVECTOR& s1, const TRIDVECTOR& s2, int time_counter); 

  int  prepare_store(double weight);
  
  
};

#endif
