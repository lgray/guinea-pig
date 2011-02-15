#include "guineapigCPP.h"
#include "particleBeamCPP.h"

VAR_HEAP var_heap;
INPUT_MAISON input;

int mainGuineapig(char *argv[])
{
   GUINEA guinee(argv[1]);
      guinee.run(argv[2],argv[3]);
    return 0;
}

           
  
int main (int argc,char *argv[])
{

  if(argc!=4) 
    {
      if (argc>4) 
	{
	  cout << "Too many arguments for guinea" << endl;
	}
      else 
	{
	  cout << "Not enough arguments for guinea" << endl;
	  cout << " principal , argc= " << argc << endl;
	}
      cout << "Usage : guinea accelerator parameter_set output_file " << endl;
      return 1;
    }
  mainGuineapig(argv);
  return 0;
}
