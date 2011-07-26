#ifndef RNDM_SEEN
#define RNDM_SEEN

#include "LesDefines.h"
#include "stdio.h"   //   exit function 
#include "stdlib.h"  //   exit function

#define RNDM_EPS 6e-8

#include <vector>
using namespace std;



class RNDM
{

  int m1_;
  int ia1_;
  int ic1_;
  float rm1_;
  int m2_;
  int ia2_;
  int ic2_;
  float rm2_;
  int m3_;
  int ia3_;
  int ic3_;
  int iff_;
  int ix1_;
  int ix2_;
  int ix3_;
  float rm3_;

  int iset_;
  float v1_,v2_;

  int compteurRndm7_;
  unsigned long rndm7_a_coeff_, rndm7_c_coeff_;
  float rndm7_modulo_dividing_factor_;

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
  int in1,in2,is[55];
} rndm3_store;

 struct
{
  float u[97],c,cd,cm;
  int i,j;
} rndm5_store;

 struct
{
  int i;
} rndm6_store;

 struct
{
  //  unsigned int i;
  unsigned long i;
  float scal;
  unsigned int n;
} rndm7_store;

 struct
{
  int i;
  float r[97];
  int ix1,ix2,ix3;
} rndm8_store;




void rndmst0(int i);
float rndm0();
void rndmst1(int i);
float rndm1();
void rndmst2(int i);
float rndm2();
void rndmst3(int i);
float rndm3();
void rndmst5(int na1,int na2,int na3, int nb1);
float rndm5();
void rndmst6(int i);
float rndm6();
//void rndmst7(int i);
void rndmst7(unsigned long i);
bool rndm_test7(unsigned long a_coeff, unsigned long c_coeff, const unsigned long* test_values);

float rndm7();
void rndmst8(int idummy);
float rndm8();
float expdev();
float exp_dev();



/* choice of random number generators 5 and 7 are good */
public:

//void rndmst();
RNDM(unsigned long seed = 1);


float gasdev();
int rndm_save();
int rndm_load();

float rndm_sincos(float *c);

void getShuffledIntegerSequence(int maxInt, vector<unsigned long int>& vect);

inline float rndm() {return rndm7();}
inline float rndm_hadron(){return rndm7();}
inline float rndm_synrad() {return rndm7();}
inline float rndm_equiv() {return rndm7();}
inline float rndm_jet() {return rndm7();}
inline float rndm_pairs() {return rndm7();}

};
#endif
