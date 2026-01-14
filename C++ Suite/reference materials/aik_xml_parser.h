#ifndef _aik_xml_parser_h
#define _aik_xml_parser_h

/*! \file aik_xml_parser.h
\brief A simple class for XML parsing (AIK-2020)
*/
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>

class xml_node_info
{
 public:
  std::string name;
  std::streampos begin_sec;
  std::streampos end_sec;

  xml_node_info() {}
 xml_node_info(const std::string& name_,  std::streampos begin_, std::streampos end_) :
  name(name_),  begin_sec(begin_), end_sec(end_) {}
 xml_node_info(const xml_node_info& other) : name(other.name), begin_sec(other.begin_sec), end_sec(other.end_sec) {}
  xml_node_info& operator=(const xml_node_info& other) {
    
    if(this!=&other) {
      name=other.name;
      begin_sec=other.begin_sec;
      end_sec=other.end_sec;
    }
    return *this;
  }

  /*
  bool operator==(const xml_node_info& other ) const {
    return (name==other.name &&  begin_sec==other.begin_sec && end_sec==other.end_sec);
    }*/

  void print(std::ostream& out=std::cout) const;
  
};


//Keeps content of one xml node 
class xml_node
{
 private:

  std::string xml_node_name;
  std::istream& xml_file; //reference to the stream
  
  std::streampos begin_sec;
  std::streampos end_sec;

  //TOC to keep track of the structure
  std::vector<xml_node_info> xml_toc;
  void make_toc();
  void print_toc(std::ostream& out=std::cout) const;
  
 public:
  xml_node(std::string name, std::ifstream& xml_file);
  xml_node(const xml_node& other);
  //xml_node(xml_node&& other); //Not needed
  //For recursive use: create a node by finding a subnode #subnode_number in the node
  //nodes counted from 0 ... nnodes-1
  xml_node(const xml_node& other,std::string subnode_name, std::size_t subnode_number);
  xml_node& operator=(const xml_node& other);
  
  //For recursive use: find how many subnodes are in the node
  std::size_t find_subnode(std::string name); 
  
  //Access functions to read content, e.g., read values such as number_of_atoms="4", units="au", etc
  //FIXIT: how to handle correctly  the case when field is not found
  std::string read_string_value(std::string field) const;
  std::size_t read_int_value(std::string field) const;
  bool read_bool_value(std::string field) const;
  double read_double_value(std::string field) const;
  bool read_flag_value(std::string field);

  //Read the entrure contetnt of the node as a string
  std::string read_node_string_value() const;
  //Read the entrure contetnt of the sub-node as a string, i.e., content of <ip> 1.0 </ip> returns '<ip> 1.0'
  std::string read_node_string_value(std::string field) const;
  //Read the contetnt of the sub-node as a double, , i.e., content of <ip> 1.0 </ip> returns 1.0
  double read_node_double_value(std::string field) const;
  //Read the contetnt of the node as a double, , i.e., content of <ip> 1.0 </ip> returns 1.0
  double read_node_double_value() const;

  bool is_valid_state() const;
  void print(std::ostream& out);
  void print_status(std::ostream& out);
};


/*! Helper class to extract what we need from the string: 
slightly modified code from the old XML parser */
class My_istringstream
{
  std::istringstream iStr;
  
  //Need to take this out of the class  
  bool ifLetterOrNumber(char Ch) {
    bool return_bool;
    
    if ( ((int(Ch)<=int('Z')) and (int(Ch)>=int('A')))  
	 or ((int(Ch)<=int('z')) and (int(Ch)>=int('a'))) 
	 or ((int(Ch)<=int('9')) and (int(Ch)>=int('0')))  )
      return_bool=true;
    else 
      return_bool=false;
    
    return return_bool;
  };
  
 public:
  //My_istringstream(const char* str){ iStr.str(str); iStr.clear(); };
  //AIK: I think this is safer then above. But it does not make a difference
  My_istringstream(const std::string& str){ iStr.str(str); iStr.clear();
    
    //std::cout << "My_istringstream(const std::string& str)   str=" << iStr.str() << std::endl;
  };
  //My_istringstream(const std::string str) { iStr.str(str); iStr.clear(); };
  My_istringstream(const My_istringstream& other){ iStr.str(other.iStr.str()); iStr.clear();
    //std::cout << "My_istringstream(const std::string& other)   str=" << iStr.str() << std::endl;
  };
  
  std::string str() { return iStr.str(); };
  bool fail() { return iStr.fail(); } //|| iStr.eof(); }; //AIK added eof
  
  int getNextInt() {
    int next;
    iStr>>next;
    return next;
  };

  double getNextDouble() {
    double next;
    iStr>>next;
    return next;
  };

  bool getNextWord(std::string& next) {
    char tmp_char;
    
    //skip spaces (non-letters/numbers):
    tmp_char=' ';
    while ( not(ifLetterOrNumber(tmp_char)) and not(iStr.fail()) and not(iStr.eof())) {
      iStr.get(tmp_char);
      //std::cout << "tmpcha=" << tmp_char << std::endl;
    }
    
#if 0
    if ( !ifLetterOrNumber(tmp_char) )
	std::cout << "tmpcha=" << tmp_char << std::endl;
      if (iStr.fail()) {
	std::cout << "FAIL" << std::endl;
	std::cout << "tmpcha=" << tmp_char << std::endl;
      }
#endif      
      
      //read the word:
      next="";
      if (not(iStr.fail())) // and not(iStr.eof()))
	do {
	  //std::cout << "Next= "<< next << std::endl;
	  next+=tmp_char;
	  //std::cout << "New next= "<< next << std::endl;
	  iStr.get(tmp_char);
	  //std::cout << "tmpcha=" << tmp_char << std::endl;
	} while(ifLetterOrNumber(tmp_char) && not(iStr.eof()));
      
      //std::cout << "NEXT=" << next << std::endl;
      return !(iStr.fail()); //returns true for normal execution
    };
};

#endif
