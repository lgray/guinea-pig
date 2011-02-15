#ifndef LESDEFINES_SEEN
#define LESDEFINES_SEEN



//#ifndef LONGPHOT
//#define LONGPHOT
//#endif

#ifndef SHORTCUT_FFT
#define SHORTCUT_FFT
#endif



/* this is necessary only for some workstations */
#ifndef CLK_TCK 
#define CLK_TCK 1000000
#endif


#ifndef SCALE_ENERGY
#define SCALE_ENERGY
#endif

#ifndef PAIR_SYN
#define PAIR_SYN
#endif

const int ZPOS=1;
//#ifndef ZPOS
//#define ZPOS
//#endif

#ifndef XYPOS
#define XYPOS
#endif

//const int SPIN = 0;
//const int EXT_FIELD = 1;

//#ifndef EXT_FIELD
//#define EXT_FIELD
//#endif

const int MCROSS = 1;
const int  AVERCROSS = 0;

#define NCROSSMAX 201
#define NVALMAX 21


#ifndef MOD 
#define MOD %
#endif


#endif
