#include "mathematicalTools.h"

//double PAIR_TOOLS::scdgam2i_ = 0.0;
// void TOOLS::mkit(double gam2i,double& c, RNDM& hasard)
// {
//    JET_FLOAT x,sigma0,beta,y;
//     beta=sqrt((double)1.0-gam2i);
//     //   scdgam2i_ = gam2i;
//     sigma0= pair_ang_i(beta, gam2i);
//     y=(2.0*hasard.rndm_pairs()-1.0)*sigma0;
//     if (y<=0.0)
//       {
// 	x=-0.5*beta;
// 	TOOLS::equal_newton(&TOOLS::pair_ang_i,&TOOLS::pair_ang_d,-beta,0.0,y,gam2i, x);

//       }
//     else
//       {
// 	x=0.5*beta;
// 	TOOLS::equal_newton(&TOOLS::pair_ang_i,&TOOLS::pair_ang_d,0.0,beta,y,gam2i, x);
//       }
//     c= x/beta;
// }

// void TOOLS::synrad (float eng, double* e2, double* e3, double* polar, int sokolov,float radius_i,float dz,vector<float>&  photon,int* number, RNDM& hasard) 
// {
//   int n,i,j=0;
//   float tmp;
//   float photener = 0.0;
//   // corresponds to A const p.21 in DS's thesis
//   tmp=synrad_p0(eng,radius_i,dz);
//   n=(int)(tmp*10.0)+1;

//   dz/=(double)n;
//   for (i=0;i<n;i++)
//     {
//       //      int aux = synrad_0(eng,e2, e3, polar, sokolov, radius_i,dz,photon+j, hasard);
//       int aux = synrad_0(eng,e2, e3, polar, sokolov, radius_i,dz,&photener, hasard);
//       if (aux) 
// 	{
// 	  photon.push_back(photener);
// 	  radius_i *= eng/(eng-photon[j]);
// 	  if (photon[j]<=0.0) 
// 	    {
// 	      cerr << "warning TOOLS::synrad " << photon[j] << " " << eng << " " << j << " " << n << endl;
// 	    }
	  
// 	  eng -= photon[j];
// 	  j++;
// 	  // the limit of 1000 corresponds to the dimension of the input 
// 	  // array photon : method to be modified
// 	  // 
// 	  if (j>=1000)
// 	    {
// 	      cerr << " TOOLS::synrad too many photons (>= 1000) produced by one particle, j= " << j << endl;
// 	      exit(-1);
// 	    }
// 	}
      
//     }
//   *number=j;
// }
// // adopt the variables change of CAIN (Yokoya, user's guide 235.1 p. 111))
// // eng : energie en GeV
// // dz : metres 
// // avec spin flip
// int TOOLS::synrad_0 (float eng, double* e2, double* e3, double* polar, int sokolov,float radius_i,float dz,float* photonEnergy, RNDM& hasard)
// {
//   int j,k;
//   double x, s2, s3;
//   double fu0, fusp;
//   double p0,p1,v1,v3,g;
//   double fK13, fKi13, fKi53, fK23;
//   j=0;
//   if (eng<=0.0)
//     {
//       cerr << "Initial particle energy below zero : " << eng << endl;
//       return 1;
//     }

//   // upsilon_bar = h_bar. omegac / E
//   double upsilon_bar = CCRIT*eng*eng*radius_i;
//   double upsilon = 0.6666666667 * upsilon_bar;
//   double gamma = eng/EMASS;
//   double factor =  pow ( (1.0 + 0.5 * upsilon_bar ), 0.33333333333);
//   p0 = CONST1 * dz * gamma * radius_i / factor ; 
//   //  cout << " p0 = " << p0 << " a l'ancienne " << synrad_p0(eng,radius_i,dz) << endl;
//   if (hasard.rndm_synrad()>p0) return 0;
//   p1=hasard.rndm_synrad();
//   while((v1=hasard.rndm_synrad())==0.0) ; /* v1!= 0.0 */
//   v3 = v1*v1*v1;
//   double xden = 1.0 - v3 + 0.5 * upsilon_bar * ( 1.0 + v3 * v3 );
//   x = upsilon_bar * v3 / xden;
//   double x1 = 1.0 - x;
//   double z = x/(upsilon_bar * x1);
//   synradKi53andK23(z, fKi53, fK23);
//   double F00 = fKi53 + (x*x/x1) * fK23;
//   double F00star;
//   if ( sokolov == 1 )
//     {
//       s2 = polar[0]*e2[0] + polar[1]*e2[1] + polar[2]*e2[2];
//       fK13 = synradK13(z);
//       F00star = F00 - s2 * x * fK13; 
//     }
//   else F00star = F00;
//   double dxdy = 3.0 * v1 * v1 * (  upsilon_bar + x * (1.0 - upsilon_bar * v3 ) ) /xden;
//   g  = F00 * dxdy * factor /(CONST0  * upsilon );
//   if ( p1 < g) 
//     {
//       if (sokolov) 
// 	{
// 	  s3 = polar[0]*e3[0] + polar[1]*e3[1] + polar[2]*e3[2];
// 	  fKi13 = Ki13(z);
// 	  for (k=0; k < 3; k++)
// 	    {
// 	      polar[k] = ( F00*polar[k] - x/(1.0-x) * fK13 * e2[k] - x*x/(1.0-x) * ( s3*e3[k]* fKi13 + (polar[k] - s3 * e3[k]  ) * fK23 ) ) / F00star;
// 	    }
// 	}
//       *photonEnergy = eng * x;
//       return 1;
//     }
//   else 
//     {
//       if (sokolov) 
// 	{
// 	  s2 = polar[0]*e2[0] + polar[1]*e2[1] + polar[2]*e2[2];
// 	  fu0 = 1.0 - CONST3 * gamma*radius_i * dz * fradu0(upsilon);
// 	  fusp = CONST3 * gamma*radius_i * dz * fradsp(upsilon);
// 	  double c0 = fu0 + fusp * s2;
// 	  for (k=0; k<3; k++)
// 	    {
// 	      polar[k] = ( polar[k] * fu0 + fusp * e2[k] ) / c0;
// 	    }
// 	  double sum = polar[0]*polar[0] + polar[1]*polar[1] + polar[2]*polar[2];
// 	  if (sum > 1.0)
// 	    {
// 	      sum = 1.0 / sqrt(sum);
// 	      for ( k=0; k < 3; k++) polar[k] *= sum;
// 	    }
// 	}
//       *photonEnergy = 0.0;
//       return 0;
//     }
// }

void TOOLS::rotmat(double angle, double axis[3], double r[3][3])
  // from CAIN
{
  double a0, a1, a2, a3, si;
  a0 = cos(angle/2.);
  si = sin(angle/2.);
  a1 = si*axis[0];
  a2 = si*axis[1];
  a3 = si*axis[2];
  r[0][0] = 1. - 2.0*(a2*a2 + a3*a3);
  r[1][0] = 2.0*(a1*a2 - a3*a0);
  r[2][0] = 2.0*(a1*a3 + a2*a0);
  r[0][1] = 2.0*(a1*a2 + a3*a0);
  r[1][1] = 1. - 2.0*(a3*a3 + a1*a1);
  r[2][1] = 2.0*(a2*a3 - a1*a0);
  r[0][2] = 2.0*(a1*a3 - a2*a0);
  r[1][2] = 2.0*(a2*a3 + a1*a0);
  r[2][2] = 1. - 2.0*(a1*a1 + a2*a2);
}

// double TOOLS::BformFunction(double ups)
// {
//   // function giving the coefficient for interpolating the anomalous 
//   // magnetic moment, as a function of upsilonfollowing Baier

// //   a = alpha/(2*pi)*F(U) 
// //             2             x*dx                   x      t**3 
// //      F(U)= --- integral -------- * integral sin[---(t + ----)]*dt
// //             U           (1+x)**3                 U       3     
// //      F(0)=1

//   double mu;
//   double ups0 = 0.6125;
//   double acoef[8] = { 12.190054896 , 24.159295870, -37.341016656 , -190.56332408 , -267.06921477, -80.540475512, -95.539356489, 246.29207314};

//   double bcoef[8] = { 0.51911469393, 0.75556983816, -0.98346938317, 0.31707943752, -1.5451047715, 0.67601308567, -0.061924565451, -0.23548134968 };

//   double u2, u3,u4,logups;

//   if (ups <= 0.0) mu = 1.0;
//   else
//     {
//       if ( ups <= ups0) 
// 	{
// 	  logups = log(ups);
// 	  mu = 1 + ups*ups*(acoef[0]*logups+ acoef[1] + 
// 				ups*(acoef[2]*logups + acoef[3] + 
// 				     ups*(acoef[4]*logups + acoef[5] +
// 					  ups*(acoef[6]*logups + 
// 					       acoef[7]))));
// 	}
//       else
// 	{
// 	  u2 = 1./(ups*ups);
// 	  u3 = pow(u2, 0.3333333333333333);
// 	  u4 = u3*u3;
// 	  logups = log(ups);
// 	  mu = bcoef[0]*u3 + bcoef[1]*u4 + (bcoef[2]*logups + bcoef[3])*u2 + 
// 	    (bcoef[4]*u3 + bcoef[5]*u4 + (bcoef[6]*logups + bcoef[7])*u2)*u2;
// 	}
//     }
//   mu *= ANOMEL;
//   return mu;
// }

//  // integral K(5/3,u).du  over x < u < infinity
// void  TOOLS::synradKi53andK23(double x, double& Ki53, double& K23)
//  {
//    if (x <= 1.54 )
//      {
//        double x1 = pow(x, 0.6666666667);
//        double x2 = x1*x1;
//        double x1inv = 1.0/x1;
//        Ki53 = ( ( 0.03363782042 * x1 - 0.1134842702 ) * x2 +
//                 0.3944669710 ) * x2  -
// 	 1.812672515 + 2.1495282415344901 * x1inv;

//        K23 = ( ( ( 0.04122878139 * x1 - 0.1494040962 ) * x2 +
// 		 0.7862616059 ) * x1 - 1.258219373 ) * x1 + 1.074647298 * x1inv;

//      }
//    else
//      if ( x <= 4.48 )
//        {
// 	 Ki53 = ( ( 0.09120815010 * x  -1.229105693) * x + 4.442223505 ) /
// 	   ( ( ( x - 0.6903991322) * x + 5.651947051 ) * x - 0.9691386396 );
// 	 K23 = ( ( 0.08194471311 * x - 1.112728296 ) * x + 4.052334415 ) /
//            ( ( ( x - 0.6524469236 ) * x + 6.1800441958 ) * x - 0.4915880600 );
//        }
//      else
//        if ( x <= 1.65 )
// 	 {
// 	   double c = exp(-x) / sqrt(x);
// 	   Ki53 = c * ( 2.187014852 + 1.253535946 * x ) / ( 0.9949036186 + x );
// 	   K23  = c * ( 0.6120387636 + 1.253322122 * x ) / ( 0.3915531539 + x );
// 	 }
//        else
// 	 {
// 	   Ki53 = 0.0;
// 	   K23  = 0.0;
// 	 }
//  }

// MODIFIED BESSEL FUNCTION K(1/3,X) (from Yokoya, corresponds to L=1 in CAIN)
double TOOLS::K13(double x) 
{
  const double A[4] = { 1.687570834 , -1.611154585 , 0.611515182 , -0.256982212 };
  const double B[4] = { 1.253273433 ,  0.671444784 , 0.604094190 ,  0.010909437 };
  const double x0 = 0.546;
  double resu = 0.0;
  double x13, x2,y;
  if ( x <= 0.0 ) 
    {
      cerr << " mathematicalTools::K13 : invalid arguments : x = " << x << endl;
      return 0.0;
    }
  if ( x <= x0 ) 
    {
      x13 = pow( x, 0.3333333333);
      x2 = x*x;
      resu = (A[0]+A[2]*x2)/x13+(A[1]+A[3]*x2)*x13;
    }
  else
    {
      y = 1.0/x;
      resu = sqrt(y) * ( B[0] + B[1] * y ) / ( 1.0 + y * ( B[2] + y * B[3] ) );
      if ( x >= 130. ) resu = 0.0;
      else resu *= exp(-x);
    }
  return resu;
}

// Integral of  MODIFIED BESSEL FUNCTION K(1/3,X) (from Yokoya, corresponds 
// to L=1 in CAIN)
// integral K(1/3,u).du  over x < u < infinity
double TOOLS::Ki13(double x) 
{
  double x2, x23, y;
  const double A[5] = {1.8136145, -2.5278090, 1.1979670, -0.20490413, 0.058671692 };
  const double B[4] = {1.2531864, 1.6215470, 1.8548547, 0.28146211 };
  const double x0 = 1.2777;
  double resu = 0.0;
  if ( x < 0.0 ) 
    {
      cerr << " mathematicalTools::Ki13 : invalid arguments : x = " << x << endl;
      return 0.0;
    }
  if ( x < x0 ) 
    {
      if (x == 0.0 ) resu = A[0];
      else
	{
	  x23 = pow( x, 0.666666667);
	  x2 = x*x;
	  resu = A[0] + x23 * (A[1] + A[3] * x2 + x23 * (A[2] + A[4] * x2));
	}
    }
  else
    {
      y = 1.0/x;
      resu = sqrt(y) * ( B[0] + B[1] * y ) / ( 1.0 + y * ( B[2] + y * B[3] ) );
      if ( x >= 130. ) resu = 0.0;
      else resu *= exp(-x);
    }
  return resu;
}

// //  FRADU0 = (total radiation rate in quantum theory)/(that in classical)
// //  as a function of Upsilon.
// //   Accuracy:  Max.relative error < 5.05D-6  (for any Upsilon>=0)
// // (from yokoya)
// double TOOLS::fradu0(double ups)
// {
//   double x;
//   const double A[6] = {8.0503959320, 10.9756518968, 1.4194297054, 8.9730711310, 15.8341489137, 4.1313700056 };
//   const double B[6] = {1.0118692717, 2.9309973207, 1.6930111582, 2.8972660432, 2.2315495296, 1.7387035105 };
//   const double ups0 = 1.3261;
//   double resu;
//   if (ups <= ups0) resu = ( ( ( A[2] * ups +A[1] ) * ups +A[0] ) * ups +1.0 ) / 
// 		     ( ( ( A[5] * ups +A[4] ) * ups +A[3] ) * ups +1.0);
//   else
//     {
//       x = pow(1.0/ups, 0.3333333333);
//       resu = ( ( B[2] * x +B[1] ) * x +B[0] ) * x / ( ( ( B[5] * x + B[4] ) * x + B[3] ) * x + 1.0);
//     }
//   return resu;
// }

//   Integral of the spin-dependent term of radiation,
//   normalized by the classical total rate of radiation.
//                  2
//   FRADSP(Y)= ------ * Integral x*dx*K(1/3,(2/3/Y)*x/(1-x)) 
//               5*Pi*Y                      (over 0<x<1)
//   Relative error:  < 0.936E-5 for any Y.
// (from yokoya)
// double TOOLS::fradsp(double ups)
// {
//   double x;
//   const double A[7] = { .299997194, 2.62120367, .386895664, 15.6610660, 57.6471880, 26.4479482, .543635258 };
//   const double B[6] = {9.91434629e-02, 4.92000917e-04, 1.45454865, 1.78219471e-03, 3.90751942e-01, 1.88294532e-01 };
//   const double ups0 = 1.06237;
//   double resu = 0.0;
//   if (ups <= ups0 )
//     {
//       resu = ups * ( A[0] + ups * ( A[1] + ups * A[2] ) ) / 
// 	( 1.0 + ups * ( A[3] + ups * ( A[4] + ups  * ( A[5] + ups * A[6] ) ) ) );
//     }
//   else
//     {
//       x = pow(1.0/ups, 0.3333333333);
//       resu = x*x * B[0] / ( 1.0 + x * ( B[1] + x * ( B[2] + x * ( B[3] + x* ( B[4] + x * B[5] ) ) ) ) );
//     }
//   return resu;
// }
