#ifndef ABSTRACTIOCLASS_SEEN
#define ABSTRACTIOCLASS_SEEN

class  ABSTRACT_IO_CLASS
{

 protected :
 inline string title(string str) const 
   
{
  ostringstream out;
  
  out << "                                         " << endl;
  out << " *********************************************************************** " << endl;
  out << " ------------- " << str << " ------------ " << endl;
  out << "             " << endl;
  return out.str();
}

 public :
   
 ABSTRACT_IO_CLASS() {;}
 virtual  ~ABSTRACT_IO_CLASS() {;}
 
 virtual string output_flow() const =0;
 virtual string persistent_flow() const {return "";}

};


#endif
