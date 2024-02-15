#include "lumiCPP.h"

int GENERAL_LUMI_HEAP::stack_vector(int nstore, std::vector<int>& selected_indices)
  {
     int k;
     double store,scal;
     scal=1.0/(double)((int)((nstore+nb_pairs_)/(nmax_))+1);
     p_ *= scal;      
     for (k = 0; k< nb_pairs_; k++)
	{
	  if(rndm_generator_->rndm()<scal)
	    {
              selected_indices.push_back(k);
	    }
	}
      store = nstore*scal;
      nstore = numberToStore(store);
      return nstore;
}

 
void LUMI_HEAP::lumi_store(const MESH& mesh, int cellx, int celly,double min_z, double energy1, double energy2, double weight) 
{
  int k;
  int nstore = prepare_store(weight);
  for (k=0; k < nstore; k++)
    {
      data_.push_back(LUMI_PAIR(energy1, energy2));
      data_.back().random_position(mesh,cellx, celly, min_z,*rndm_generator_);
      nb_pairs_++;
    }
}

int  LUMI_HEAP::prepare_store(double weight)
{
  unsigned int j;
  int nstore;
  double store;
  //double scal;
  store=weight*p_;
  nstore = numberToStore(store);
  if(nstore + nb_pairs_ > nmax_ )
    {
      std::vector<int> selected_indices;
      nstore = stack_vector(nstore, selected_indices);
      nb_pairs_ = 0;
      for (j=0; j<selected_indices.size(); j++)
          {
	      data_[nb_pairs_++] = data_[selected_indices[j]];
          }
      data_.resize(nb_pairs_);
      // at the end of this method: certrin elements of the rndm_generator are
      // stored in the elements of data_ and they are transferred to the first
      // elements of vector. nstore is the number of stored elements.
    }
  return nstore;
}

int  LUMI_HEAP_EE::prepare_store(double weight)
{
  unsigned int j;
  int nstore;
  double store;
  //double scal;
  store=weight*p_;
  nstore = numberToStore(store);
  if(nstore + nb_pairs_ > nmax_ )
    {
      std::vector<int> selected_indices;
      // recuperate the pair indices to store
      nstore = stack_vector(nstore, selected_indices);
      nb_pairs_ = 0;
      for (j=0; j<selected_indices.size(); j++)
	{
	  data_[nb_pairs_++] = data_[selected_indices[j]];
	}
      data_.resize(nb_pairs_);
      // at the end of this method: certrin elements of the rndm_generator are
      // stored in the elements of data_ and they are transferred to the first
      // elements of vector. nstore is the number of stored elements.
    }
  return nstore;
}


void LUMI_HEAP_EE::lumi_store_ee(const MESH& mesh, int cellx, int celly,double min_z, double energy1, double p1Vx, double p1Vy, double energy2, double p2Vx, double p2Vy,double weight,int time_counter) 
{
  int k;
  int nstore = prepare_store(weight);
  for (k=0; k < nstore; k++)
    {
      data_.push_back(LUMI_PAIR_EE(energy1, energy2, p1Vx,p1Vy,p2Vx, p2Vy,time_counter));
      data_.back().random_position(mesh,cellx, celly, min_z,*rndm_generator_);
      nb_pairs_++;
    }
}

void LUMI_HEAP_EE::lumi_store_ee(const MESH& mesh, int cellx, int celly,double min_z, double energy1, double p1Vx, double p1Vy, double energy2, double p2Vx, double p2Vy,double weight, const TRIDVECTOR& s1, const TRIDVECTOR& s2, int time_counter) 
{
  int k;
  int nstore = prepare_store(weight);
  for (k=0; k < nstore; k++)
    {
      data_.push_back(LUMI_PAIR_EE(energy1, energy2, p1Vx,p1Vy,p2Vx, p2Vy,time_counter));
      data_.back().random_position(mesh,cellx, celly, min_z,*rndm_generator_);
      data_.back().set_spins( s1,s2);
      nb_pairs_++;
    }
}
