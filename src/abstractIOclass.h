#ifndef ABSTRACTIOCLASS_SEEN
#define ABSTRACTIOCLASS_SEEN

class  ABSTRACT_IO_CLASS
{

 protected :
 inline string titre(string str) const 
   
{
  ostringstream sortie;
  
  sortie << "                                         " << endl;
  sortie << " *********************************************************************** " << endl;
  sortie << " ------------- " << str << " ------------ " << endl;
  sortie << "             " << endl;
  return sortie.str();
}

 public :
   
   ABSTRACT_IO_CLASS() {;}
 virtual  ~ABSTRACT_IO_CLASS() {;}
 
 virtual string name_of_class() const = 0;
 virtual string output_flow() const =0;
 virtual string persistent_flow() const {return "";}

};


#endif
