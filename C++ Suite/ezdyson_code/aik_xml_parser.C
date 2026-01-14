#include "aik_xml_parser.h"

/*! Finds word in istream from specified postion, returns posion relative to the begining of the stream
 right at the begining of the target 
*/
std::streampos find_in_stream(std::streampos start_at,  std::streampos end_at,
			      std::istream& xml_file, std::string target) {

  std::streampos begin_sec=-1, previous_step=start_at;
  xml_file.clear();
  xml_file.seekg(start_at);

  //std::cout <<"Looking for " << target << " between position=" << start_at << " and position=" << end_at << std::endl;
  
  std::string word;
  bool found=false;
  bool stop=false;
  
  while (xml_file >> word && !found && !stop) {
    if(word==target) {
      //std::cout << "Word found:" << target << std::endl;
      found=true;
      begin_sec = previous_step;
    }
    stop = (-1==end_at) ? false : (xml_file.tellg() >= end_at);
    previous_step=xml_file.tellg();
    //std::cout << "Prev step " << previous_step << " Stop=" << stop << std::endl;
  } 
  return begin_sec;
}

/*! Finds sequence (which may be part of the word or a sepate word) in istream from specified postion, 
returns posion of the word containing the sequence relative to the begining of the stream.
If the target equal to the word, then result should be identical to find_seq_in_stream
*/
std::streampos find_seq_in_stream(std::streampos start_at,  std::streampos end_at,
			      std::istream& xml_file, std::string target) {

  std::streampos begin_sec=-1, previous_step=start_at;
  xml_file.clear();
  xml_file.seekg(start_at);

  //std::cout <<"Looking for " << target << " between position=" << start_at << " and position=" << end_at << std::endl;
  
  std::string word;
  bool found=false;
  bool stop=false;
  
  while (xml_file >> word && !found && !stop) {

    size_t pos_match=word.rfind(target);

    if(pos_match < word.length()) {
      //std::cout << "Word found:" << target << std::endl;
      found=true;
      begin_sec = previous_step;
    }
    stop = (-1==end_at) ? false : (xml_file.tellg() >= end_at);
    previous_step=xml_file.tellg(); //This is now position of the word in which the sequense is found
    //std::cout << "Prev step " << previous_step << " Stop=" << stop << std::endl;
  } 
  return begin_sec;
}

//! Read string between two specified positions
bool read_string(std::streampos begin_sec, std::streampos end_sec, std::istream& xml_file,
		 std::string& result) {

  if(begin_sec == -1 || end_sec == -1)
  //if(begin_sec == -1)
    return false;
  else {
    std::size_t length=end_sec-begin_sec;
    xml_file.clear();
    xml_file.seekg(begin_sec,std::ios::beg);
    char *buffer=new char [length];
    xml_file.read(buffer,length);
    result=std::string(buffer);
    delete [] buffer;
    return true;
  }
}

/*! Finds the first instance of '<target' or '<target>', returns stream position where tag begins.
  Assuming only one type is present, will get confused otherwise. */
std::streampos find_head_in_stream(std::streampos start_at,  std::streampos end_at, std::istream& xml_file, std::string name) {

  std::string b_name1="<"+name, b_name2="<"+name+">";
  //std::cout << "Looking for " << b_name1 << " or " << b_name2 <<  " and " << e_name << std::endl;
  //std::streampos current=xml_file.tellg();
  //std::cout << "Current position= " << current << std::endl;
  std::streampos begin_sec = find_in_stream(start_at,end_at,xml_file,b_name1);
  if ( begin_sec == -1) //not found, try bname_2
    begin_sec = find_in_stream(start_at,end_at,xml_file,b_name2);
  return begin_sec;
}

//! Finds the first instance of '</target>', should be separated from whatever before it by space or newline (will fix later)
std::streampos find_tail_in_stream(std::streampos start_at,  std::streampos end_at, std::istream& xml_file, std::string name) {

  std::string e_name="</"+name+">";
  std::streampos tail_pos=find_in_stream(start_at,end_at,xml_file,e_name);
  //std::streampos tail_pos2=find_in_stream(tail_pos,end_at,xml_file,">"); //Does not work!!!

  xml_file.clear();
  xml_file.seekg(tail_pos);

  std::string word;
  xml_file >> word;
  tail_pos=xml_file.tellg();

  return tail_pos;
}

/* Finds first xml subsection in the stream, returns true if found, puts information in node_info */ 
bool find_section_in_stream(std::streampos start_at,  std::streampos end_at, std::istream& xml_file, xml_node_info& node_info) {

  bool found = false;
  if ( !(start_at < end_at))
    return false;

  std::streampos current_pos=start_at;
  xml_file.clear();
  xml_file.seekg(current_pos);

  std::streampos head_found=-1;
  std::string word;

  //Find first instance 
  while( !found && current_pos < end_at) {

    head_found=find_seq_in_stream(current_pos,end_at,xml_file,"<"); 
    current_pos=xml_file.tellg();
    
    if(head_found < end_at && (head_found >-1)) {

      xml_file.seekg(head_found);
      xml_file >> word;
      current_pos=xml_file.tellg();
    
      //Identify comments and closing tags and skip them
      if(!(word[1]=='!' || word[1]=='/')) 
	found=true;
    }
  }

  //std::cout << "Head found at " << head_found << std::endl;

  if(found) {
    
    //std::cout << "Found first < at " << head_found << " in word" << word << std::endl;
    size_t adjust= (word[word.length()-1]=='>') ? 1 : 0;
    std::string sec_name (word,1,word.length()-1-adjust);
    //std::cout << "Sec name = " << sec_name << std::endl;
      
    //Find end: 
    std::streampos tail_found=find_tail_in_stream(head_found,end_at,xml_file,sec_name); 
    if (!(tail_found <= end_at) || (tail_found == -1)) {
      std::cout << "Closing tag for section " <<sec_name << " not fond\n";
      exit(1);
    }
      
    node_info=xml_node_info(sec_name,head_found,tail_found);
    //std::cout << "\n\n XML info: ";
    //node_info.print();
  }

  return found;
}


void xml_node_info::print(std::ostream& out) const {
  
  
  out << std::setw(20) <<std::left  << name   <<std::right << std::setw(10) << begin_sec << std::setw(10) << end_sec ; 
}

/* I am for now assuming that all tags have the same structure: <nodename attributes> element content</nodename>. I will not deal with empty-element tags, such as <nodename ..... />. These to be replaced 
by the construct avbove. That means that the scripts should be modified accordingly. */
xml_node::xml_node(std::string name, std::ifstream& file) : xml_node_name(name), xml_file(file),
							    begin_sec(-1), end_sec(-1) {

  std::streampos b_line;
  std::streampos e_line;

  std::streampos current=xml_file.tellg();
  //std::cout << "Current position= " << current << std::endl;

  begin_sec = find_head_in_stream(current,-1,xml_file,xml_node_name);
  if (begin_sec != -1) 
    end_sec = find_tail_in_stream(current,-1,xml_file,xml_node_name);
  
  //std::cout << "Node " << name << " initialized in state=" << is_valid_state() << std::endl;
  //std::cout << "Current pos in xml_file=" << xml_file.tellg() << std::endl;

  make_toc();
  //print_status(std::cout);
  //print_toc(std::cout);
}

xml_node::xml_node(const xml_node& other) :
  xml_node_name(other.xml_node_name), xml_file(other.xml_file),begin_sec(other.begin_sec), end_sec(other.end_sec) {

  std::cout << "Using copy-constructor" << std::endl; //not used
  xml_toc=other.xml_toc;
}

xml_node::xml_node(const xml_node& parent,std::string subnode_name, std::size_t subnode_number) :
  xml_node_name(subnode_name), xml_file(parent.xml_file),begin_sec(-1), end_sec(-1) {

  if ( parent.is_valid_state() ) {

    //std::cout << "Find subnode "<< subnode_name << " in node " << parent.xml_node_name << std::endl;
    //std::cout << "Parent TOC:" << std::endl;
    //parent.print_toc();
    
    std::size_t nodes_found=0, found_at=0;
    
    for(size_t i=0; i < parent.xml_toc.size() && nodes_found < subnode_number+1;  i++, found_at++)
      nodes_found+= (parent.xml_toc[i].name==subnode_name) ? 1 : 0;
    found_at--;
    
    //std::cout << "Nodes_found = " << nodes_found << " Found at =" << found_at << "Node name=" << parent.xml_toc[found_at].name << std::endl;

    if(subnode_number<nodes_found) {

      //std::cout << "Found subnode.... " << std::endl;
      begin_sec=parent.xml_toc[found_at].begin_sec;
      end_sec=parent.xml_toc[found_at].end_sec;
      make_toc();
    }
    
    //print_status(std::cout);
    //print_toc(std::cout);
  }
    
}


xml_node& xml_node::operator=(const xml_node& other) {

  std::cout << "Use operator=" << std::endl;
  
  if(this!=&other) {
    xml_node_name=other.xml_node_name;
    begin_sec=other.begin_sec;
    end_sec=other.end_sec;
    xml_toc=other.xml_toc;
  }
  return *this;
}

void xml_node::make_toc() {

  if(xml_toc.size()) {
    std::cout << "xml_node::make_toc: table already exists" << std::endl;
    exit(1);
  }

  //std::cout << "Start making TOC ....  TOC size="  << xml_toc.size() << "\n";
  //Fill TOC here
  
  if ( is_valid_state() ) {

    xml_node_info master_sec, tmp_node;
    find_section_in_stream(begin_sec,end_sec,xml_file,master_sec);
    //master_sec.print(); //This is current section 
    //std::cout << std::endl;
      
    xml_file.clear();
    xml_file.seekg(master_sec.begin_sec);
    std::string tmp_word;
    xml_file >> tmp_word; //advance into the current xml section
    std::streampos current_pos=xml_file.tellg();
    //Now figure out how to chop off the closing master section piece 
    std::streampos new_end=find_in_stream(current_pos, master_sec.end_sec,xml_file,"</"+master_sec.name+">");
    //std::cout << "Current pos =" << current_pos << "End pos = " << new_end << std::endl;

    while ( (current_pos>-1) && (current_pos < new_end) && find_section_in_stream(current_pos, new_end, xml_file, tmp_node) ) {
      xml_toc.push_back(tmp_node);
      xml_file.seekg(tmp_node.end_sec);
      current_pos=xml_file.tellg();
    }
  }
    
  //  print_toc();
}


void xml_node::print_toc(std::ostream& out) const {
  
  out << "TOC (size=" << xml_toc.size() << ")" << " for node " << xml_node_name<< ":\n";
  if ( xml_toc.size() < 25) {
   
    for(size_t i=0; i<xml_toc.size(); i++) {
      xml_toc[i].print();
      out << std::endl;
    }
    out << std::endl;
  }
}

std::size_t xml_node::find_subnode(std::string subname) {

  //std::cout << "Find subnode "<< subname << " in node " << xml_node_name << std::endl;
  std::size_t nodes_found=0;
  for(size_t i=0; i < xml_toc.size();  i++)
    nodes_found+= (xml_toc[i].name==subname) ? 1 : 0; 
  
  
#if 0
  //Old code, does not  follow hiearchy
  std::size_t nodes_found=0;
  
  if ( is_valid_state() ) {
    
    std::streampos current_pos=begin_sec;
    xml_file.clear();
    xml_file.seekg(current_pos);
      
    while( !xml_file.eof() && xml_file.tellg() < end_sec) {

      current_pos=find_head_in_stream(xml_file.tellg(),end_sec,xml_file,subname);
      if(current_pos > 0 ) 
	nodes_found++;
      //std::cout << "current_pos= "<< current_pos << " tellg=" << xml_file.tellg() << std::endl;
    }
  }
#endif  
  return nodes_found;
}


//! Access functions to read content: reads string appearing in quotes after field
std::string xml_node::read_string_value(std::string field) const  { 

  std::string current_string;
  std::string return_string;
  read_string(begin_sec,end_sec,xml_file,current_string);

  std::size_t length=0;
  std::size_t found=current_string.find(field);

  if (!(found == std::string::npos)) {

    std::size_t  found_first_quotes=current_string.find("\"",found+1);
    std::size_t  found_second_quotes=current_string.find("\"",found_first_quotes+1);
    length=found_second_quotes-found_first_quotes-1;

    //std::cout << "Quotes found at: "<< found_first_quotes << " and " << found_second_quotes << std::endl;
    //std::cout << "Length= " << length;
  
    return_string.assign(current_string,found_first_quotes+1,length);
    //std::cout << "  Return string= " << return_string << std::endl;
  }
  
  return return_string;
}

std::size_t xml_node::read_int_value(std::string field) const  { 

  return std::stoi(read_string_value(field));
}

bool xml_node::read_bool_value(std::string field) const  { 

  return (read_string_value(field)=="true");
}


bool xml_node::read_flag_value(std::string field) { 

  bool flag_value=false;
  if(find_subnode(field)) 
    flag_value=xml_node(*this,field,1).read_bool_value("flag");
  
  return flag_value;
}

double xml_node::read_double_value(std::string field) const  { 

  return std::stod(read_string_value(field));
}

//Read the entrure contetnt of the node as a string, sans closing tag
std::string xml_node::read_node_string_value() const {

  std::string current_string;
  read_string(begin_sec,end_sec,xml_file,current_string);
  return current_string;
}

//Read the entire contetnt of the sub-node as a string
std::string xml_node::read_node_string_value(std::string field) const {

  return xml_node(*this,field,0).read_node_string_value();
}

//Read the contetnt of the sub-node as a double
double xml_node::read_node_double_value(std::string field) const {

  return xml_node(*this,field,0).read_node_double_value();

  /*  My_istringstream iStr(read_node_string_value(field));
  std::string tmp;
  iStr.getNextWord(tmp);
  
  return iStr.getNextDouble();*/ 
}

//Read the contetnt of the node as a double
double xml_node::read_node_double_value() const {

  std::istringstream node_string(read_node_string_value());
  
  std::string tmp;
  size_t pos;
  bool found=false;

  do {
    node_string >> tmp;
    pos=tmp.rfind(">");
    found=pos<tmp.length();
    //std::cout << "Word= " << tmp << " Pos=" << pos << " found=" << found << std::endl;
  } while( !found && !node_string.fail());

  if(found==false) {

    std::cout << "xml_node::read_node_double_value(): closing > not found" << std::endl;
    std::cout << " xml_node::read_node_double_value(): String value "  << node_string.str() << std::endl;
    exit(1);
  }
    
  //std::cout << " xml_node::read_node_double_value() "  << node_string.str() << tmp << std::endl;

  double next;
  node_string >> next;
  
  return next; 
}

void xml_node::print(std::ostream& out) {

  out << "Node "<<  xml_node_name << std::endl;
  print_toc(out);
  //std::string current_string;
  //read_string(begin_sec,end_sec,xml_file,current_string);
  //out <<  current_string << std::endl;
  out << read_node_string_value() << std::endl;
}

void xml_node::print_status(std::ostream& out) {

  out << "Node "<<  xml_node_name << std::endl;
  out <<  "Start at=" << begin_sec << " and end at=" << end_sec << std::endl;
  out << "State=" << is_valid_state() << std::endl;
}

bool xml_node::is_valid_state() const {

  if ( (begin_sec==-1) || (end_sec==-1) || (begin_sec > end_sec) )
    //if ( (begin_sec==-1) || (begin_sec > end_sec) )
    return false;
  else
    return true;
}
  
