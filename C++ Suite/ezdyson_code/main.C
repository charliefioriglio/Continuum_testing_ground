#include "dyson_main.h"

#include <cstdlib>
#include <string>
#include <iostream>
#include <time.h>

// Returns a string with the current time:
std::string GetTime(){ time_t rawtime;  time(&rawtime);  return  ctime(&rawtime); }

int main(int argc, char *argv[])
{
  if ((argc-1)!=1) {
    std::cout << argv[0] << ": only one argument is required, <input.xml>.\n" << std::flush;
    exit(1);
  }

  std::cout << "Job \"" << argv[0] << ' ' << argv[1] << "\" has been started: " << GetTime() << '\n' << std::flush;

 #if 0
  //FIXIT: 
  simpleXMLparser xmlF;
  xmlF.assignFile(argv[1]); 

  std::string job;
  job=xmlF.reset().node("root").value("job");
#endif
  
  std::string job="dyson";
  bool done = false;
  if (job == "dyson")
    done = dyson_main(argv[1]);
  
  if ( !done )
    std::cout << "Method \"" << job <<"\" is unknown, or it has been failed. \n";

  //FIXIT: will do it later
  // copy <input.xml> in the output:
  //xmlF.printInputFile(); 
  
  std::cout << '\n' << "Job \"" << argv[0] << ' ' << argv[1] << "\" has been finished: " <<  GetTime() << '\n' << std::flush ;

  return EXIT_SUCCESS;
}




