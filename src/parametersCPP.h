#ifndef PARAMETERS_SEEN
#define PARAMETERS_SEEN

#include <iostream>

#include <string>
extern "C"
{
#include "readData.h"
}

#include "option_args.h"

class PARAMETERS
{

  MEMORY_ACCOUNT m_account_;
  DATEI datei_;


  inline int readIntegerValue(std::string nom) const
  {
    VALUE value;
    get_named_variable( const_cast<char*>(nom.c_str()), &value);
    return (int)CONTENTS(value);
  }
  inline float readFloatValue(std::string nom) const
  {
    VALUE value;
    get_named_variable( const_cast<char*>(nom.c_str()), &value);
    return (float)CONTENTS(value);
  }

  inline double readDoubleValue(std::string nom) const
  {
    VALUE value;
    get_named_variable( const_cast<char*>(nom.c_str()), &value);
    return (double)CONTENTS(value);
  }


  inline int fileOpen(DATEI *datei,std::string name, std::string type) const
  {
    return  file_open(datei,const_cast<char*>(name.c_str()),const_cast<char*>(type.c_str()));
  }


  inline int fileFindWord(DATEI *datei, std::string word) const
  {
    return  file_find_word(datei,   const_cast<char*>(word.c_str()));
  }

  inline int fileNextWord(DATEI *datei, std::string word) const
  {
    return  file_next_word(datei, const_cast<char*>(word.c_str()));
  }

  inline int fileReadBraces(DATEI *datei,std::string begin, std::string end, char* buff,int n_max) const
  {
    return file_read_braces(datei,const_cast<char*>(begin.c_str()),const_cast<char*>(end.c_str()),buff,n_max);
    //  return file_read_braces(datei,begin, end,buff,n_max);
  }

 public: 

  inline void read_accelerator(char *name) 
  {
    int n=max_buffer2;
    char buffer[max_buffer2];
  
    def_acc(&m_account_);

    do
      {
	do
	  {
	    if (fileFindWord(&datei_,"$ACCELERATOR::")==0) 
	      {
		std::cerr << "Error: Accelerator not found with name "<< name << std::endl;
		return;
	      }
	  }
	while(fileNextWord(&datei_,name)==0);
      }
    while(!fileReadBraces(&datei_,"{","}",buffer,n)) ;
    input_values(buffer, &m_account_);
  }

  inline void read_parameters(char *par ) 
  {
    int n=max_buffer2;
    char buffer[max_buffer2];
  
    do
      {
	do
	  {
	    if (fileFindWord(&datei_,"$PARAMETERS::")==0) 
	      {
		std::cerr << "Error: Parameters  not found with name " << par << std::endl;
		return;
	      }
	  }
	while(file_next_word(&datei_,par)==0);
      }
    while(!fileReadBraces(&datei_,"{","}",buffer,n)) ;
    def_param(&m_account_);
    input_values(buffer, &m_account_);
  }



  PARAMETERS() 
    {
      init_memory_account(&m_account_,0);
      int testFile =0;
      testFile = fileOpen(&datei_,get_acc_filename(),"r");
  
      if (!testFile) 
	{
	  std::cerr << " PARAMETERS:: error in opening file acc.dat " << std::endl;
	  exit(0);
	}
      init_named_variable(200, &m_account_);
    }

  ~PARAMETERS() {  file_close(&datei_);}


  void setDoubleValue(std::string nom, double d) const
  {
    VALUE value;
    double_to_value(d,&value);
    set_named_variable( const_cast<char*>(nom.c_str()), value);  
  }


  inline double readDValue(std::string nom) const
  {
    return readDoubleValue(nom);
  }

  inline float readFValue(std::string nom) const
  {
    return readFloatValue(nom);
  }

  inline int readIValue(std::string nom) const
  {
    return readIntegerValue(nom);
  }


};

#endif
