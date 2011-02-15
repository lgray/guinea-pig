#include "lumiCPP.h"

int GENERAL_LUMI_HEAP::stack_vector(int nstore, vector<int>& selectionned_indices)
  {
     int k;
     float store,scal;
     scal=1.0/(float)((int)((nstore+nb_pairs_)/(nmax_))+1);
     p_ *= scal;      
     for (k = 0; k< nb_pairs_; k++)
	{
	  if(hasard_->rndm()<scal)
	    {
              selectionned_indices.push_back(k);
	    }
	}
      store = nstore*scal;
      nstore = numberToStore(store);
      return nstore;
}

 
void LUMI_HEAP::lumi_store(const MESH& mesh, int cellx, int celly,float min_z, float energy1, float energy2, float weight) 
{
  int k;
  int nstore = prepare_store(weight);
  for (k=0; k < nstore; k++)
    {
      data_.push_back(LUMI_PAIR(energy1, energy2));
      data_.back().random_position(mesh,cellx, celly, min_z,*hasard_);
      nb_pairs_++;
    }
}

int  LUMI_HEAP::prepare_store(float weight)
{
  unsigned int j;
  int nstore;
  float store;
  //float scal;
  store=weight*p_;
  nstore = numberToStore(store);
  if(nstore + nb_pairs_ > nmax_ )
    {
      vector<int> selectionned_indices;
      nstore = stack_vector(nstore, selectionned_indices);
      nb_pairs_ = 0;
      for (j=0; j<selectionned_indices.size(); j++)
          {
	      data_[nb_pairs_++] = data_[selectionned_indices[j]];
          }
      data_.resize(nb_pairs_);
      // a la sortie de cette boucle : on a conserve, au hasard certains 
      // elements de data_ et on les a translates dans les premiers 
      // elements du vecteur. n_ est le nombre des elements conserves
      // index_ aussi (?)
    }
  return nstore;
}

int  LUMI_HEAP_EE::prepare_store(float weight)
{
  unsigned int j;
  int nstore;
  float store;
  //float scal;
  store=weight*p_;
  nstore = numberToStore(store);
  if(nstore + nb_pairs_ > nmax_ )
    {
      vector<int> selectionned_indices;
      // recuperer les indices de paires a conserver
      nstore = stack_vector(nstore, selectionned_indices);
      nb_pairs_ = 0;
      for (j=0; j<selectionned_indices.size(); j++)
	{
	  data_[nb_pairs_++] = data_[selectionned_indices[j]];
	}
      data_.resize(nb_pairs_);
      // a la sortie de cette boucle : on a conserve, au hasard certains 
      // elements de data_ et on les a translates dans les premiers 
      // elements du vecteur. n_ est le nombre des elements conserves
      // index_ aussi (?)
    }
  return nstore;
}


void LUMI_HEAP_EE::lumi_store_ee(const MESH& mesh, int cellx, int celly,float min_z, float energy1, float p1Vx, float p1Vy, float energy2, float p2Vx, float p2Vy,float weight,int time_counter) 
{
  int k;
  int nstore = prepare_store(weight);
  for (k=0; k < nstore; k++)
    {
      data_.push_back(LUMI_PAIR_EE(energy1, energy2, p1Vx,p1Vy,p2Vx, p2Vy,time_counter));
      data_.back().random_position(mesh,cellx, celly, min_z,*hasard_);
      nb_pairs_++;
    }
}

void LUMI_HEAP_EE::lumi_store_ee(const MESH& mesh, int cellx, int celly,float min_z, float energy1, float p1Vx, float p1Vy, float energy2, float p2Vx, float p2Vy,float weight, const TRIDVECTOR& s1, const TRIDVECTOR& s2, int time_counter) 
{
  int k;
  int nstore = prepare_store(weight);
  for (k=0; k < nstore; k++)
    {
      data_.push_back(LUMI_PAIR_EE(energy1, energy2, p1Vx,p1Vy,p2Vx, p2Vy,time_counter));
      data_.back().random_position(mesh,cellx, celly, min_z,*hasard_);
      data_.back().set_spins( s1,s2);
      nb_pairs_++;
    }
}
