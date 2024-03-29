#ifndef RNDM_SEEN
#define RNDM_SEEN

#include <cstdio>   //   exit function 
#include <cstdlib>  //   exit function

#include "define.h"

#define RNDM_EPS 6e-8

#include <vector>

/** 
    Random number generator class
    Implements several random number generators.
    Currently number 7 is used in all cases
 */

class RNDM
{

  int m1_;
  int ia1_;
  int ic1_;
  double rm1_;
  int m2_;
  int ia2_;
  int ic2_;
  double rm2_;
  int m3_;
  int ia3_;
  int ic3_;
  int iff_;
  int ix1_;
  int ix2_;
  int ix3_;
  double rm3_;

  int iset_;
  double v1_,v2_;

  /** counter to count how often rndm7 is called */ 
  int counterRndm7_;
  unsigned long rndm7_a_coeff_, rndm7_c_coeff_;
  double rndm7_modulo_dividing_factor_;

  struct
  {
    int i;
  } rndm0_store;
  
  struct
  {
    long i,p,is[32];
  } rndm1_store;
  
  struct
  {
    int i1,i2,p,is[32];
  } rndm2_store;
  
  struct
  {
    double u[97],c,cd,cm;
    int i,j;
  } rndm5_store;
  
  struct
  {
    int i;
  } rndm6_store;
  
  struct
  {
    unsigned long i;
  } rndm7_store;
  
  struct
  {
    double r[97];
    int ix1,ix2,ix3;
  } rndm8_store;
  
  
  /** random number generators
      only rndm7 is used at the moment
   */  
  void rndmst0(int i);
  double rndm0();
  void rndmst1(int i);
  double rndm1();
  void rndmst2(int i);
  double rndm2();
  void rndmst5(int na1,int na2,int na3, int nb1);
  double rndm5();
  void rndmst6(int i);
  double rndm6();

  void rndmst7(unsigned long i);
  bool rndm_test7(unsigned long a_coeff, unsigned long c_coeff, const unsigned long* test_values);
  double rndm7();

  void rndmst8(int idummy);
  double rndm8();
  double expdev();
  double exp_dev();
  


  /* choice of random number generators 5 and 7 are good */
 public:
  /** constructor with seed for random number generator 7 */
  RNDM(unsigned long seed = 1);
  
  /** gaussian distributed random number*/
  double gasdev();
  /** save random number status */
  int rndm_save();
  /** load random number status from file rndm.save */
  int rndm_load();
  
  /** method draws random number theta between [0,2*Pi] and return value is sin(theta), and *c is cos(theta)  */
  double rndm_sincos(double *c);

  /** get a vector containing the sequence of integers from 1 to maxInt 
      randomly shuffled */
  void getShuffledIntegerSequence(int maxInt, std::vector<unsigned long int>& vect);
  
  /** public random number generators */
  inline double rndm() {return rndm7();}
  inline double rndm_hadron(){return rndm7();}
  inline double rndm_synrad() {return rndm7();}
  inline double rndm_equiv() {return rndm7();}
  inline double rndm_jet() {return rndm7();}
  inline double rndm_pairs() {return rndm7();}
  
};
#endif
