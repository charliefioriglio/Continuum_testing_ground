#include "readwrite.h"

void ReadCharges(int natoms, double *charges, const char* xmlFileName) {

  std::ifstream xml_file(xmlFileName);
  xml_file.clear();
  xml_file.seekg(0, std::ios::beg);
  xml_node node_free_electron("free_electron",xml_file);

  if ( node_free_electron.find_subnode("charges") ) {

    xml_node node_charges(node_free_electron, "charges",0);
    std::stringstream tmp_iStr(node_charges.read_string_value("text"));
    for(int i=0;  i<natoms; i++)
      tmp_iStr >> charges[i];
  }
  else {
    double radfn=node_free_electron.read_double_value("charge_of_ionized_core");
    for(int i=0;  i<natoms; i++)
      charges[i]=radfn;
  }
  
  std::cout << "Charge of ionized core for multi-center calculation:\n";
  for(int i=0;  i<natoms; i++)
    std::cout << charges[i] << "  ";
  std::cout << std::endl;  

  xml_file.close();
}

//! Which MOs in AO basis to read from QChem molden_format output?
void MONumbers(int *monumbers, const char* xmlFileName, int nomos)
{
  std::ifstream xml_file(xmlFileName);
  xml_file.clear();
  xml_file.seekg(0, std::ios::beg);
  xml_node node_job_parameters("job_parameters",xml_file);

  std::istringstream tmp_iStr;
  
  if (nomos>0) {
    
    tmp_iStr.str(node_job_parameters.read_string_value("MOs_to_plot") );
    for (int i=0; i<nomos; i++) {
      tmp_iStr >> monumbers[i];
      check(tmp_iStr.fail(),"MONumbers: invalid value="+tmp_iStr.str());
    }

    std::cout << nomos << " orbitals to plot with the following numbers: ";
    for (int i=0; i<nomos; i++)
      std::cout << monumbers[i] << ' ';
    std::cout << '\n';
  } 
  else
    std::cout << "No orbitals were requested to plot.\n";

  xml_file.close();
}


//!
void ReadOrbitalFromFile(Orbital& orb, int nb, const char* xmlFileName)
{
  std::ifstream xml_file(xmlFileName); 
  xml_file.clear();
  xml_file.seekg(0, std::ios::beg);
  xml_node node_jobparams("job_parameters",xml_file);

  //FIXIT: unrestricted won't work, but it was not working before either ...
  //  bool ifrestr = !( xmlF.getBoolValue("unrestricted") );

/*! SG: Here I will override the above and set ifrestr = true always
 * This means regardless of what is in the input we will always read only two
 * dyson orbitals (left-right) and not four (left/alpha-left/beta-right/alpha-right/beta).
 * This is because in ground to IP cases we only ionize one electron (alpha)
 * and so beta Dyson orbitals are always 0.
 */
  bool ifrestr = true;

  int transno = node_jobparams.read_int_value("Dyson_MO_transitions");

  double ln, rn;

  //! Read dyson Norm
  xml_node dyson_mos("dyson_molecular_orbitals",xml_file);
    
  //transno = 1,2... but nodes (DMOs) counted 0,1, 2 ... 
  //not sure if it will work for unrestricted.... 
  std::size_t mo_position_l=(transno*2-1)-1;
  std::size_t mo_position_r=(ifrestr) ? transno*2-1 : transno*2; 
  
  //  std::cout << "Transno =" << transno << " DMO number for left=" <<  mo_position_l << " DMO number for right=" <<  mo_position_r << std::endl;
  
  xml_node dyson_left(dyson_mos,"DMO", mo_position_l);
  //dyson_left.print(std::cout);
  ln=dyson_left.read_double_value("norm");
  //std::cout << "LN= " << ln << std::endl;

  xml_node dyson_right(dyson_mos,"DMO", mo_position_r);
  rn=dyson_right.read_double_value("norm");
  //std::cout << "RN= " << rn << std::endl;
  
  double *tmpR=new double[nb];
  double *tmpL=new double[nb];

  {
    std::stringstream tmp_iStr(dyson_left.read_string_value("text"));
    for (int j=0; j<nb; j++)
      tmp_iStr >> tmpL[j];
  }
  {
    std::stringstream tmp_iStr(dyson_right.read_string_value("text"));
    for (int j=0; j<nb; j++)
      tmp_iStr >> tmpR[j];
  }
  
  orb=Orbital(nb,ifrestr,transno,ln,rn,tmpL,tmpR);
  delete [] tmpR;
  delete [] tmpL;

  orb.Print();
  xml_file.close();
}

//! Read MO to plot from qchem.
void ReadMOFromFile(Orbital &mo, int nb, int monumber, const char* xmlFileName)
{
  std::ifstream xml_file(xmlFileName);
  xml_file.clear();
  xml_file.seekg(0, std::ios::beg);
  xml_node node_job_parameters("job_parameters",xml_file);

  bool ifrestr = !( node_job_parameters.read_bool_value("unrestricted"));
  
  double *tmpA= new double[nb];
  double *tmpB= new double[nb];
 
  //! Read MOs in AO basis coeffs from qchem.
  xml_node node_molecular_orbitals("molecular_orbitals",xml_file);
  //Numbering of sections 0,1, ....
  xml_node node_alpha_mo(node_molecular_orbitals,"alpha_MOs",0);
  xml_node node_alpha_mo_to_plot(node_alpha_mo,"MO",monumber-1); 
  
  {
    std::stringstream tmp_iStr(node_alpha_mo_to_plot.read_string_value("text"));
    for (int k=0; k<nb; k++)
      tmp_iStr >> tmpA[k];
  }
  
  if (ifrestr==false) {

    xml_node node_beta_mo(node_molecular_orbitals,"beta_MOs",0);
    xml_node node_beta_mo_to_plot(node_beta_mo,"MO",monumber-1); 
    std::stringstream tmp_iStr(node_beta_mo_to_plot.read_string_value("text"));
    for (int k=0; k<nb; k++)
      tmp_iStr >> tmpB[k];
  }

  mo=Orbital(nb,ifrestr,monumber,1.0,1.0,tmpA,tmpB);
  delete [] tmpA;
  delete [] tmpB;

  mo.Print();
  xml_file.close();
}


//!Choose whether to compute overlap integral.
//SG:This is still in testing mode and not fully implemented.
Overlap ReadOverlap(const char* xmlFileName)
{
  std::string tmpStr;
  Overlap overlap=NO;

  std::ifstream xml_file(xmlFileName);
  xml_node node_root("root",xml_file);

  if ( node_root.find_subnode("overlap") ) {  

    xml_node node_overlap(node_root,"overlap",0);
    tmpStr=node_overlap.read_string_value("calc_overlap");

    if (tmpStr=="both")
      overlap=BOTH;
    else if (tmpStr=="only")
      overlap=ONLY;
    else if (tmpStr=="no")
      overlap=NO;
    else
      error("ReadOverlap: Invalid value "+tmpStr);
  }
  else
    overlap=NO;
  
  fprintf(stdout,"\nCalculation of overlaps: %d\n",overlap); fflush(stdout);
  xml_file.close();
  
  return overlap;
}


//!Read info about averaging over molec orientations in lab frame.
//SG: Now we only give a choice to the user to use avg or not
MolAvg ReadMolecOrientAvg(const char* xmlFileName)
{
  std::string tmpStr;
  MolAvg molavg;

  std::ifstream xml_file(xmlFileName);
  xml_node node_root("root",xml_file);
  if ( node_root.find_subnode("averaging") ) {

    xml_node node_averaging(node_root,"averaging",0);
    
    tmpStr=node_averaging.read_string_value("method");

    if (tmpStr=="noavg")
      molavg=NOAVG;
    else if (tmpStr=="avg")
      molavg=AVG;
    else if (tmpStr=="num")
      molavg=NUM;
    else
      error("ReadMolecOrientAvg: unrecognized value "+tmpStr);
  }
  else
    molavg=AVG;

  fprintf(stdout,"\nAveraging scheme: %d\n",molavg); fflush(stdout);
  xml_file.close();
  
  return molavg;
}


//!x,y,z lab frame components of ionization laser polarization.
void ReadRIonz(double &ie, double *rioniz, const char* xmlFileName)
{
  std::ifstream xml_file(xmlFileName);
  xml_node node_laser("laser",xml_file);

  // read and convert to a.u.
  ie = node_laser.read_double_value("ionization_energy") / HAR_TO_EV;

  // SG: read laser polarization (now only applicable for fixed molecular orientation)
  xml_node node_laser_pol(node_laser,"laser_polarization",0);
  //node_laser_pol.print(std::cout);
  
  rioniz[0]=node_laser_pol.read_double_value("x=");
  rioniz[1]=node_laser_pol.read_double_value("y=");
  rioniz[2]=node_laser_pol.read_double_value("z=");

  std::cout << "IE=" << ie*HAR_TO_EV <<" eV;  Polarization=(" << rioniz[0] << ',' << rioniz[1] << ',' << rioniz[2] << ")\n";
  xml_file.close();
}

//! Print
void PrintOrbitalsScanZ(Orbital *mostoplot, int nomos, XYZGrid& labgrid)
{
  FILE *file_orbs=fopen("mosplot.dat","w");
  double *gridptr_z=labgrid.GetGridPtr(Z);
  double xyz[XYZ], Ldys_value, Rdys_value;

  xyz[X]=xyz[Y]=0.0;

  fprintf(file_orbs, "Z,A RDys_val RD*Z \n");

  for(int i=0; i<labgrid.NPoints(Z); i++)
    {
      xyz[Z]=gridptr_z[i];
      fprintf(file_orbs,"%lE",xyz[Z]*BOHR_TO_ANGS);
      for (int j=0; j<nomos; j++)
        {
         mostoplot[j].CalcDyson_xyz(Ldys_value,Rdys_value,0.0,0.0,xyz[Z]);
         fprintf(file_orbs," %lE %lE   ", Rdys_value, Rdys_value*xyz[Z]);
        }
      fprintf(file_orbs,"\n");
    }

  fclose(file_orbs);
}



