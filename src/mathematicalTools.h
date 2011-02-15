#ifndef MATHEMATICALTOOLS_SEEN
#define  MATHEMATICALTOOLS_SEEN
#include <iostream>
#include "besselCPP.h"
#include "typeDefs.h"
#include "rndmCPP.h"
#include "physconst.h"

class TOOLS 
{

  public : 

    /* This routine takes an electron with energy e and longitudinal momentum
       pz in the center of mass frame of the two photons e1 and e2 and boosts the
       momentum to the laborytory frame. */
    
    /* inline static void lorent_pair(float e1,float e2,double& e,double& pz) */
    /* { */
    /*     double beta, eold, gam; */
    
    /*     beta = -((double)e1 - (double)e2) / ((double)e1 + (double)e2); */
    /*     gam = (double)1.0 / sqrt((double)1.0 - beta * beta); */
    /*     eold = e; */
    /*     e = gam * (e - beta * pz); */
    /*     pz = gam * (pz - beta * eold); */
    /* } */
    
    /* inline static void lorent_pair(float e1,float e2, float& e, float& pz) */
    /*   { */
    /*     double ed = (double)e; */
    /*     double pzd = (double)pz; */
    /*     lorent_pair(e1, e2, ed, pzd); */
    /*     e = (float)ed; */
    /*     pz = (float)pzd; */
    /*   } */
    
    /* inline static  void lorent(double& e,double& pl,double beta) */
    /* { */
    /* /\* +beta because want to transform back *\/ */
    /*     double gamma,eold; */
    /*     gamma=1.0/sqrt(1.0-beta*beta); */
    /*     eold = e; */
    /*     e = gamma*(eold + beta * pl); */
    /*     pl = gamma*(pl + beta * eold); */
    /* } */
    
    /* inline static  void lorent(float& e, float& pl, float beta) */
    /* { */
    /*     double ed = (double)e; */
    /*     double pld = (double)pl; */
    /*     lorent(ed, pld, (double)beta); */
    /*     e = (float)ed; */
    /*     pl = (float)pld; */
    /* } */
    
    /* /\* Differential probability for producing a pair * sqrt(3)*PI *\/ */
    
    /* inline static float fcp(float ups,float x) */
    /* { */
    /*   float eta; */
    /*   eta=2.0/(3.0*ups*x*(1.0-x)); */
    /*   return bessel.ki13(eta)+(x/(1.0-x)+(1.0-x)/x)*bessel.k23(eta); */
    /* } */
    
    /* /\* total probability of producing a pair *\/ */
    
    /* inline static float u(float ups) */
    /* { */
    /*   return 0.23*ups*exp(-8.0/(3.0*ups))*pow(1.0+0.22*ups,-1.0/3.0); */
    /* } */
    
    inline static int verifPuissanceDe2(int aTester )
    {
      int 	ntest=0;
      int nh=  aTester;
      while (nh%2 == 0) 
	{
	  nh = nh >> 1;
	  ntest++;
	}
      if (nh>1)
	{
	  ntest = 0;
	}
      return ntest;
    }
  
  inline static int nearest_power_of_2(float x)
    {
      int value = (int)x;
      int index = 2;
      while (index < value) index *= 2;
      int index2 = index/2; 
      return (value - index2) < index - value ? index2 : index;
    }
  
  inline static void equal_newton(double (*fint)(double),double (*fdiff)(double),
				  double xmin,double xmax,double y,double& x)
    {
      //    double eps=1e-6,tiny=1e-20;
      double eps=1e-6;
      //    double ytry,xtry,dfdx;
      double ytry,xtry;
      int i=0;
      xtry = x;
      ytry=(*fint)(xtry);
      while (fabs(ytry-y)>(fabs(ytry)+fabs(y))*eps
	     && (xmax-xmin)>eps) {
	i++;
	xtry-=(ytry-y)/(*fdiff)(xtry);
	if ((xtry>=xmax)||(xtry<=xmin)) {
	  xtry=0.5*(xmax+xmin);
	}
	ytry=(*fint)(xtry);
	if(ytry<y) {
	  xmin=xtry;
	}
	else {
	  xmax=xtry;
	}
      }
      x = xtry;
    }
  
  
  inline static void equal_newton(double (*fint)(double, double),double (*fdiff)(double, double), double xmin,double xmax,double y, double scdgam2i, double& x)
    {
      //    double eps=1e-6,tiny=1e-20;
      double eps=1e-6;
      //    double ytry,xtry,dfdx;
      double ytry,xtry;
      int i=0;
      xtry = x;
      ytry=(*fint)(xtry, scdgam2i);
      while (fabs(ytry-y)>(fabs(ytry)+fabs(y))*eps
	     && (xmax-xmin)>eps) {
	i++;
	xtry-=(ytry-y)/(*fdiff)(xtry, scdgam2i);
	if ((xtry>=xmax)||(xtry<=xmin)) {
	  xtry=0.5*(xmax+xmin);
	}
	ytry=(*fint)(xtry, scdgam2i);
	if(ytry<y) {
	  xmin=xtry;
	}
	else {
	  xmax=xtry;
	}
      }
      x = xtry;
    }
  
  
  /*  inline static  JET_FLOAT pair_ang_d(JET_FLOAT x, double scdgam2i) */
  /* { */
  /*     JET_FLOAT t,u; */
  /*     t=-(1.0-x); */
  /*     u=-(1.0+x); */
  /*     return t/u+u/t-2.0*scdgam2i*(1.0/t+1.0/u)-scdgam2i*scdgam2i*(1.0/t+1.0/u)*(1.0/t+1.0/u); */
  /* } */
  
  /* inline static  JET_FLOAT pair_ang_i(JET_FLOAT x, double scdgam2i) */
  /* { */
  /*   return 2.0*((1.0+scdgam2i*(1.0-0.5*scdgam2i))*log((1.0+x)/(1.0-x))- (1.0+scdgam2i*scdgam2i/(1.0-x*x))*x); */
  /* } */
  
  /* inline static JET_FLOAT hcd_q1q2_q1q2(JET_FLOAT x) */
  /* { */
  /*     JET_FLOAT t,u; */
  /*     t=0.5*(1.0-x); */
  /*     u=0.5*(1.0+x); */
  /*     return 4.0/9.0*(1.0+u*u)/(t*t); */
  /* } */
  
  /* inline static JET_FLOAT hci_q1q2_q1q2(JET_FLOAT x) */
  /* { */
  /*     return 4.0/9.0*(x+8.0/(1.0-x)+4.0*log(1.0-x)); */
  /* } */
  
  /* inline static JET_FLOAT hcd_qqb(JET_FLOAT x) */
  /* { */
  /*     JET_FLOAT t,u; */
  /*     t=1.0-x; */
  /*     u=1.0+x; */
  /*     return (t*t+u*u)/(t*u); */
  /* } */
  
  /* inline static JET_FLOAT hci_qqb(JET_FLOAT x) */
  /* { */
  /*     return 2.0*(log((1.0+x)/(1.0-x))-x); */
  /* } */
  
  /* /\* sigma=t/s+s/t *\/ */
  
  /* inline static JET_FLOAT hcd_qph(JET_FLOAT x) */
  /* { */
  /*     JET_FLOAT t,s; */
  /*     t=1.0-x; */
  /*     s=2.0; */
  /*     return (s*s+t*t)/(s*t); */
  /* } */
  
  /* inline static JET_FLOAT hci_qph(JET_FLOAT x) */
  /* { */
  /*     return -2.0*log(1.0-x)+0.5*x*(1.0-0.5*x); */
  /* } */
  
  
  /* inline static JET_FLOAT hcd_q1q1_q1q1(JET_FLOAT x)  */
  /* { */
  /*     JET_FLOAT t,u; */
  /*     t=0.5*(1.0-x); */
  /*     u=0.5*(1.0+x); */
  /*     return 4.0/9.0*((1.0+u*u)/(t*t)+(1.0+t*t)/(u*u))-8.0/(27.0*u*t); */
  /* } */
  
  /* inline static JET_FLOAT hci_q1q1_q1q1(JET_FLOAT x) */
  /* { */
  /*     return 8.0/9.0*(x*(1.0+8.0/(1.0-x*x))+8.0/3.0*log((1.0-x)/(1.0+x))); */
  /* } */
  
  /* inline static JET_FLOAT hcd_q1q1b_q2q2b(JET_FLOAT x)  */
  /* { */
  /*     JET_FLOAT t,u; */
  /*     t=0.5*(1.0-x); */
  /*     u=0.5*(1.0+x); */
  /*     return 4.0/9.0*(t*t+u*u); */
  /* } */

  /* inline static  JET_FLOAT hci_q1q1b_q2q2b(JET_FLOAT x) */
  /* { */
  /*     return 1.0/27.0*x*(3.0+x*x); */
  /* } */
  
  /* inline static JET_FLOAT hcd_q1q1b_q1q1b(JET_FLOAT x)  */
  /* { */
  /*     JET_FLOAT t,u; */
  /*     t=-0.5*(1.0-x); */
  /*     u=-0.5*(1.0+x); */
  /*     return 4.0/9.0*((1.0+u*u)/(t*t)+(t*t+u*u))-8.0*u*u/(27.0*t); */
  /* } */
  
  /* inline static JET_FLOAT hci_q1q1b_q1q1b(JET_FLOAT x) */
  /* { */
  /*     return 2.0/9.0*(x*(1.0-x*(1.0-x))+16.0/(1.0-x)+16.0/3.0*log(1.0-x)); */
  /* } */
  
  /* inline static JET_FLOAT hcd_q1q1b_gg(JET_FLOAT x)  */
  /* { */
  /*     JET_FLOAT t,u; */
  /*     t=0.5*(1.0-x); */
  /*     u=0.5*(1.0+x); */
  /*     return 32.0/27.0*(u/t+t/u)-8.0/3.0*(t*t+u*u); */
  /* } */
  
  /* inline static JET_FLOAT  hci_q1q1b_gg(JET_FLOAT x) */
  /* { */
  /*     return 4.0/9.0*(x*(-25.0/3.0-x*x)+16.0/3.0*log((1.0+x)/(1.0-x))); */
  /* } */
  
  
  /* inline static JET_FLOAT hcd_gg_q1q1b(JET_FLOAT x)  */
  /* { */
  /*     JET_FLOAT t,u; */
  /*     t=0.5*(1.0-x); */
  /*     u=0.5*(1.0+x); */
  /*     return 1.0/6.0*(u/t+t/u)-3.0/8.0*(t*t+u*u); */
  /* } */
  
  /* inline static JET_FLOAT hci_gg_q1q1b(JET_FLOAT x) */
  /* { */
  /*     return -25.0/48.0*x+log((1.0+x)/(1.0-x))/3.0-x*x*x/16.0; */
  /* } */
  
  /* inline static  JET_FLOAT hcd_qg_gq(JET_FLOAT x)  */
  /* { */
  /*     JET_FLOAT t,u; */
  /*     t=-0.5*(1.0-x); */
  /*     u=-0.5*(1.0+x); */
  /*     return -4.0/9.0*(1.0+u*u)/u+(u*u+1.0)/(t*t); */
  /* } */
  
  /* inline static  JET_FLOAT hci_qg_gq(JET_FLOAT x) */
  /* { */
  /* /\*    return 8.0/9.0*log(0.5*(1.0+x))+11.0/9.0*x+x*x/9.0+8.0/(1.0-x) */
  /* 	+4.0*log(1.0-x);*\/ */
  /*     return 8.0/9.0*log(1.0+x)+11.0/9.0*x+x*x/9.0+8.0/(1.0-x) */
  /* 	+4.0*log(1.0-x); */
  /* } */
  
  /* inline static  JET_FLOAT hcd_gg_gg(JET_FLOAT x) */
  /* { */
  /*     JET_FLOAT t,u; */
  /*     t=-0.5*(1.0-x); */
  /*     u=-0.5*(1.0+x); */
  /*     return 9.0/2.0*(3.0-u*t-u/(t*t)-t/(u*u)); */
  /* } */
  
  /* inline static  JET_FLOAT hci_gg_gg(JET_FLOAT x)  */
  /* { */
  /*     return 99.0/8.0*x+3.0/8.0*x*x*x-9.0*log((1.0+x)/(1.0-x))+36.0*x/(1.0-x*x); */
  /* } */
  
  
  /* static void mkit(double gam2i,double& c, RNDM& hasard); */
  
  /*  inline static float synrad_p0(float eng,float radius_i,float dz)  */
  /*  { */
  /*    // 4**(1/3)*(1/137)/(Gamma(4/3)*emass)  */
  /*    //  const double pconst=25.4e0; */
  /*    // cf p21 daniel's thesis and ref 21 */
  /*    return PCONST*eng*dz*radius_i; */
  /*  } */
  
  /*  // integral K(5/3,u) du  over x < u < infinity */
  /*  inline static void  synradKi53andK23(double x, double& Ki53, double& K23); */
  
  /* inline static double synradK13(double x)  */
  /* { */
  /*   if (x <= 165.) return K13(x); */
  /*   else return 0.0; */
  /* } */
  
  // MODIFIED BESSEL FUNCTION K(1/3,X) (from Yokoya, corresponds to L=1 in CAIN)
  static double K13(double x);
  static double Ki13(double x);
  
  static void rotmat(double angle, double axis[3], double r[3][3]);
  //????????????????????????? not in this class?????????????????????
  //static double BformFunction(double ups);
  
};

#endif
