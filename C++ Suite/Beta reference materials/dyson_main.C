#include "dyson_main.h"

//! Driver for computing overlap between the dyson orbital and free wave
void CalcDysonFreeElOverlap(double charge_radfn, XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid);
//! Driver for computing xsections without partial waves, with numeric averaging
void CalcXSecWithNumericAveraging(double IE, double degeneracyfactor, XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid);
//! Driver for computing xsections with partial waves, without numeric averaging (analytic or no average)
void CalcXSecWithoutNumericAveraging(MolAvg  molavg, double IE, double degeneracyfactor, double charge_radfn, double *RIonz, XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid);
//! Driver for multi-center xsections with partial waves, without numeric averaging (analytic or no average)
void CalcMultiCenter(MolAvg  molavg, double IE, double degeneracyfactor, double *charge_radfn, double *RIonz, XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid);
//! Driver for multi-center xsections with partial waves, without numeric averaging (analytic or no average)
void CalcMultiCenterWithCoherences(MolAvg  molavg, double IE, double degeneracyfactor, double *charge_radfn, double *RIonz, XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid);

//! Final printing
void PrintXSecAndPad(KLMPoints& klmgrid, double *xsectot, PAD& totalpad, double IE);

void CalcOverapMatrixTest(const XYZGrid& grid);
void CalcLowdinTransf(const XYZGrid& grid);
  
bool dyson_main(const char* xmlFileName)
{
  std::cout << "Starting ezDyson version 5.0\n\n";
  //std::cout << "Reading the input file\n\n";

  time_t start,end;
  time (&start);
  
  //! Read atoms and AO basis information from QChem output.
  printf("Starting AO BASIS: Reading the atom and basis set information \n\n"); fflush(stdout);
  TheAOBasis().IniAOBasis(xmlFileName); 
  printf("AO BASIS done\n\n"); fflush(stdout);

  //! Initialize/read info about Dyson MOs.
  printf("Reading Dyson Orbitals\n\n"); fflush(stdout);
  Orbital dysorb;
  ReadOrbitalFromFile(dysorb, TheAOBasis().NBasis(), xmlFileName);
  printf("Dyson Orbital done\n\n"); fflush(stdout);

  std::ifstream xml_file(xmlFileName);
  xml_file.clear();
  xml_file.seekg(0, std::ios::beg);
  xml_node node_free_electron("free_electron",xml_file);
  xml_file.clear();
  xml_file.seekg(0, std::ios::beg);
  xml_node node_job_parameters("job_parameters",xml_file);

  //!Read lab frame x,y,z grid info.
  XYZGrid labgrid(xmlFileName);
  labgrid.PrintGridInfo(stdout);
  //labgrid.Print("labgrid.chk");

  //! Read MOs to plot and plot them
  int nomos = node_job_parameters.read_int_value("number_of_MOs_to_plot");
  if (nomos) {
    printf("Plotting of MOs requested.\n\n"); fflush(stdout);
    int *monumbers = new int[nomos];
    MONumbers(monumbers,xmlFileName,nomos);
    Orbital *mostoplot = new Orbital[nomos];
    for (int i=0; i<nomos; i++)
      ReadMOFromFile(mostoplot[i],TheAOBasis().NBasis(),monumbers[i],xmlFileName);
    delete [] monumbers;
    //! Plot MOs and Z*MOs after labframe center shift.  
    
    PrintOrbitalsScanZ(mostoplot,nomos,labgrid);
    delete [] mostoplot;
    printf("Plotting MOs done.\n\n"); fflush(stdout);
  }
  //end plotting MOs

  //! Calc center of (un)averaged DysonMO(left one).
  double dyscenter[XYZ], dnorm;
  dnorm=dysorb.GetOrbitalNormAndCenter(dyscenter,labgrid);
  std::cout << "Dyson orbital centroid: " << dyscenter[X] << " " << dyscenter[Y] << " " << dyscenter[Z] << std::endl;
  
  bool shif_to_dyson_cenroid=true;
  if ( node_free_electron.find_subnode("wave_origin")) {

      xml_node node_wave_origin(node_free_electron,"wave_origin",0);
      if(!node_wave_origin.read_bool_value("dyson_centroid")) {
	std::stringstream tmp_iStr(node_wave_origin.read_string_value("text"));
  // SG: modified below so input can conveniently be entered in Angstroms.
        for (int j=0; j<XYZ; j++){
          tmp_iStr >> dyscenter[j];
          dyscenter[j]*=ANGS_TO_BOHR;
        }
      }
  }
  std::cout <<  std::setw(8) << std::setprecision(6) << std::fixed;
  std::cout << "Origin of the expansion: " << dyscenter[X] << " " << dyscenter[Y] << " " << dyscenter[Z] << std::endl;
  TheAOBasis().ShiftGeom(dyscenter);

  
  //! Read info for generating the k,l,m points of the spherical waves.
  KLMPoints klmgrid(xmlFileName);
  klmgrid.PrintGridInfo();

  Overlap overlap_type = ReadOverlap(xmlFileName);

  if (overlap_type == BOTH || overlap_type == ONLY) {

    //Radial Function from input
    double radfn;
    radfn=node_free_electron.read_double_value("charge_of_ionized_core");
    std::cout << "Charge of ionized core for overlap calculation= " << radfn <<'\n';
    CalcDysonFreeElOverlap(radfn,labgrid,dysorb,klmgrid);
  }

  //! Read averaging method from input.
  MolAvg molavg = ReadMolecOrientAvg(xmlFileName);
      
  double spinfactor, orbfactor;
  spinfactor=node_job_parameters.read_double_value("spin_degeneracy");
  orbfactor=node_job_parameters.read_double_value("orbital_degeneracy");
  std::cout << "\nOther input parameters:\n";
  std::cout << "Spin degeneracy=" << spinfactor << ". Orbital degeneracy="<<orbfactor <<"\n";
  
  //! Read direction of ionization laser and ionization energy (read in eV, convert to au)
  double RIonz[XYZ], IE;
  ReadRIonz(IE, RIonz, xmlFileName);
      
  //Now decide about multicenter verus single-center
  bool do_multicenter=node_free_electron.read_bool_value("multicenter");
 
  if(!do_multicenter)  {  

    std::cout << "Single-center expansion of the free-electron state is requested" <<'\n';
    double radfn;
    radfn=node_free_electron.read_double_value("charge_of_ionized_core");
    std::cout << "Charge of ionized core for single-center calculation= " << radfn <<'\n';

    if (overlap_type == NO || overlap_type == BOTH) {
      
      //Numeric averaging, does not use partial wave. Only for plane wave.
      if (molavg == NUM)
	CalcXSecWithNumericAveraging(IE,spinfactor*orbfactor,labgrid,dysorb,klmgrid);
      //Analytic averaging or no averaging. Use partial wave expansion. Available for plane and Coulomb wave.
      else if (molavg==AVG || molavg==NOAVG)
	CalcXSecWithoutNumericAveraging(molavg,IE,spinfactor*orbfactor,radfn,RIonz,labgrid,dysorb,klmgrid);
    }
  }
  else { //multicenter treatment 
    std::cout << "Multi-center expansion of the free-electron state is requested" <<'\n';
    //Can only do analytic averaging or no averaging
    check(molavg==NUM, "WARNING: Numeric averaging not yet implemented for multi-center treatment");
    int natoms=TheAOBasis().NAtoms();
    double *charges=new double[natoms];
    ReadCharges(natoms,charges,xmlFileName);
    bool do_coherences=node_free_electron.read_bool_value("coherences");
    if(!do_coherences)
      CalcMultiCenter(molavg,IE,spinfactor*orbfactor,charges,RIonz,labgrid,dysorb,klmgrid);
    else
      CalcMultiCenterWithCoherences(molavg,IE,spinfactor*orbfactor,charges,RIonz,labgrid,dysorb,klmgrid);
    delete [] charges;
  }

  time (&end);
  double dif = difftime (end,start);
  printf("\nJob complete. Thank you for using ezDyson!\n");
  printf("Job time: %lf seconds\n",dif);
  xml_file.close();
  
  return true;
}

//Comment out to compile on fluffy
#if 0
//This is stupid function, for debugging. Very slow!
void CalcOverapMatrixTest(const XYZGrid& grid) {
  
  //! Initialize fake Dyson orbital for calculations of the overlap:
  int nbas=TheAOBasis().NBasis();

  arma::Col<double> lorb_1(nbas), rorb_1(nbas), lorb_2(nbas), rorb_2(nbas);
  //"Left" orbitals are fake, just set them to zero and do not touch 
  lorb_1.zeros(); lorb_2.zeros();

  arma::Mat<double> overlap(nbas,nbas);
  
  for(int i=0; i<nbas; i++)
    {
      rorb_1.zeros(); 
      rorb_1[i]=1.0; 
      for(int j=i; j<nbas; j++)
	{
	  rorb_2.zeros();
	  rorb_2[j]=1.0;
	  Orbital dysorb1(nbas,true,0,1.0,1.0,rorb_1.memptr(),lorb_1.memptr());
	  Orbital dysorb2(nbas,true,0,1.0,1.0,rorb_2.memptr(),lorb_2.memptr());
	  double srr, sll;
	  dysorb1.GetOverlap(dysorb2,grid,srr,sll);
	  overlap(i,j)=overlap(j,i)=srr;
	}
    }

  //Now can print overlap and work with it
  printf("Overlap matrix:\n");
  overlap.print();
  //Can diagonalize using diag function....
}
#endif


/*! Driver for computing overlap between the dyson orbital and free wave
  Executed when (overlap_type == BOTH || overlap_type == ONLY) 
 */
void CalcDysonFreeElOverlap(double charge_radfn, XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid) {

  NumEikr allEikr;
  CklmCoeff allCklm; 
  std::cout << "Computing overlap between normalized Dyson orbitals and free electron" << std::endl;
            
  allCklm.IniCklmCoeff_Overlap(labgrid,dysorb,charge_radfn,klmgrid);
  //Now computing total overlap
  int nkv = klmgrid.NKV();
  double *overlap = new double[nkv];
  memset(overlap,0,sizeof(double)*nkv);
  for (int k=0; k<nkv; k++) //LOOP BY k
    overlap[k]=allCklm.GetSumCklm(k);
         
  //Assuming our PW is normalized as 4pi 
  //and that the rest of coefficients computed correctly
  double norm=1./(4.0*M_PI);
  
  printf("\n Overlap |<phi^d|Psi_el>|^2 as function of electron energy:\n");
  printf("\nE_k,eV        overlap     overlap/(4pi)\n");
  for (int k=0; k<nkv; k++) {
    double kwave=klmgrid.GetKV(k);
    double ene=kwave*kwave/2.;  //energy above ionzaton/detach threshold
    
    printf("%lf", ene*HAR_TO_EV);
    printf(" %13.6lf   %13.6lf\n",overlap[k], overlap[k]*norm);
  }
  delete [] overlap;
}


void CalcXSecWithNumericAveraging(double IE, double degeneracyfactor, XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid) {

  NumEikr allEikr;
  AngleGrid anggrid;
  //AIK: I assume this is PW-only
  allEikr.IniNumEikr(labgrid,dysorb,anggrid,klmgrid);
  
  double norm=(1.)/(2*M_PI*C_AU);
  double dyson_norm=dysorb.GetLRNorm();
  //Calculate the conversion factors for the xsections
  double au_to_cm_sq = pow(BOHR_TO_ANGS*pow(10.,-8),2);
  printf("\nConversion factors for x-sections:\n");
  printf("Xsec(cm^-2) = Xsec(au^-2) * %10lE\n", au_to_cm_sq);
  printf("Xsec(au^-2) = Xsec(cm^-2) * %10lE\n", 1./au_to_cm_sq);
  
  printf("\n\nPADs:\n");
  printf("NOTE: PADs are only correct in molecular frame for Z-polarized light\n");
  printf("____________________________________________________________\n");
  printf("E=IE+E_k,eV   Sigma_par     Sigma_perp    Sigma_tot     beta\n");

  int nkv = klmgrid.NKV();
  for (int k=0; k<nkv; k++) {
    double kwave=klmgrid.GetKV(k);
    double ene=kwave*kwave/2. + IE;  //energy of laser, in AU
    
    //Bethe salpter.
    //double scale=norm*ene*kwave*dyson_norm*spinfactor*orbfactor;
    double scale=norm*ene*kwave*dyson_norm*degeneracyfactor;
    
    printf("%lf", ene*HAR_TO_EV);
    double sigma_par=allEikr.GetCPar(k)*scale;
    double sigma_perp=allEikr.GetCPerp(k)*scale;
    double sigma_tot=(sigma_par+2.*sigma_perp)*4.0*M_PI/3.;
    double beta=2.*(sigma_par-sigma_perp)/(sigma_par+2.*sigma_perp);
    
    printf(" %13.6lf %13.6lf %13.6lf %13.6lf\n",sigma_par,sigma_perp,sigma_tot,beta);
  }
}

void CalcXSecWithoutNumericAveraging(MolAvg  molavg, double IE, double degeneracyfactor, double charge_radfn, double *RIonz, XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid) {      

  //! Calc avg |Cklm|^2.
  CklmCoeff allCklm; //this is essentially photoelectron matrix element, but without the sum over l and m
  
  //Algorithm: go over each center, compute |Cklm|^2 and add then to the master allCklm
  if (NOAVG==molavg)
    allCklm.IniCklmCoeff(labgrid,dysorb,RIonz,charge_radfn,klmgrid);
  else // if (AVG==molavg)
    allCklm.IniCklmCoeff(labgrid,dysorb,charge_radfn,klmgrid);

  int nkv = klmgrid.NKV();
  int lmax = klmgrid.LMax();
  int ntheta = 101;
  PAD totalpad(ntheta,nkv);
  totalpad.CalcPad(allCklm);

  //FIXIT: modify to use PrintXSecAndPad(KLMPoints& klmgrid, double *xsectot, PAD& totalpad, double IE)

  double *xsectot = new double[nkv];
  double *xsec1 = new double[nkv];
  memset(xsec1,0,sizeof(double)*nkv);
  memset(xsectot,0,sizeof(double)*nkv);
 
  for (int k=0; k<nkv; k++) //LOOP BY k
    xsec1[k]=allCklm.GetSumCklm(k);
  
  //! Bethe-Salpeter factor of 4*pi^2/c devided by factor from Sakurai of (2*pi). See updated derivation
  double norm=(8*M_PI*M_PI*4)/(3*C_AU);
  
  printf("\n\nTotal cross-sections (in bohr^2) vs E(eV)\n");
  printf("See the ezDyson manual for equations and details of the normalization\n");

  double dyson_norm=dysorb.GetLRNorm();
  
//  printf("[E*2*Pi/(c*k) * |Cklm|^2]\n");
  printf("E=IE+E_k,eV   xsec,a.u.\n");
  for (int k=0; k<nkv; k++) {
    double kwave=klmgrid.GetKV(k);
    double ene=kwave*kwave/2. + IE;  //energy of laser, in AU
    
    //double scale=norm*kwave*spinfactor*orbfactor;
    double scale=norm*kwave*degeneracyfactor;
    xsec1[k] *= scale*dyson_norm;
    xsectot[k] = xsec1[k];
    
    printf("%lf", ene*HAR_TO_EV);
    printf(" %13.6lf\n",xsectot[k]*ene);
  }
  printf("\n\nTotal cross-section / E (in a.u.) vs. Ek (in eV):\n");
  printf("To be used with xsecFCF script for incorporating FCFs\n");
  printf("E_k,eV        xsec/E,a.u. \n");
  
  for (int k=0; k<nkv; k++) {
    double kwave=klmgrid.GetKV(k);
    double ene=kwave*kwave/2.;  //energy of the electron
    
    printf("%lf", ene*HAR_TO_EV);
    printf("%14.6lf\n",xsectot[k]);
  }
  //Calculate the conversion factors for the xsections
  double au_to_cm_sq = pow(BOHR_TO_ANGS*pow(10.,-8),2);
  printf("\nConversion factors for x-sections:\n");
  printf("Xsec(cm^-2) = Xsec(au^-2) * %10lE\n", au_to_cm_sq);
  printf("Xsec(au^-2) = Xsec(cm^-2) * %10lE\n", 1./au_to_cm_sq);
  
  printf("\n\nPADs:\n");
  printf("NOTE: PADs are only correct in molecular frame for Z-polarized light\n");
  printf("____________________________________________________________\n");
  printf("E=IE+E_k,eV   Sigma_par     Sigma_perp    Sigma_tot     beta\n");
  for (int k=0; k<nkv; k++) {
    double kwave=klmgrid.GetKV(k);
    double ene=kwave*kwave/2. + IE;  //energy of laser, in AU
    
    //Bethe salpter.
    //double scale=norm*ene*kwave*dyson_norm*spinfactor*orbfactor;
    double scale=norm*ene*kwave*dyson_norm*degeneracyfactor;
    
    printf("%lf", ene*HAR_TO_EV);
    double sigma_par=totalpad.XSec_par(k)*scale;
    double sigma_perp=totalpad.XSec_perp(k)*scale;
    double sigma_tot=(sigma_par+2.*sigma_perp)*4.0*M_PI/3.;
    double beta=2.*(sigma_par-sigma_perp)/(sigma_par+2.*sigma_perp);
    
    printf(" %13.6lf %13.6lf %13.6lf %13.6lf\n",sigma_par,sigma_perp,sigma_tot,beta);
  }
  
  //Write PADs on disk
  totalpad.Print(klmgrid);
      
  delete [] xsec1;
  delete [] xsectot;
}
  

//Here we can do incoherent multi-center expansion
void CalcMultiCenter(MolAvg  molavg, double IE, double degeneracyfactor, double *charge_radfn, double *RIonz, XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid) {      
  
  std::cout << "Calculating MC expansion without coherences" <<  std::endl;
  //! Bethe-Salpeter factor of 4*pi^2/c devided by factor from Sakurai of (2*pi). See updated derivation
  double norm=(8*M_PI*M_PI*4)/(3*C_AU);
  
  int nkv = klmgrid.NKV();
  int lmax = klmgrid.LMax();
  int ntheta = 101;
  PAD totalpad(ntheta,nkv);
  
  double *xsectot = new double[nkv];
  memset(xsectot,0,sizeof(double)*nkv);
  
  //We loop over centers and aggregate contributions from multiple centers to the xsec and to the pad
  AOBasis *pbasis=&(TheAOBasis());
  int natoms=pbasis->NAtoms();

  //Initialize Lowdin transformation here, before using for each center transformation
  pbasis->CalcLowdinTransf(labgrid);

  double total_normsq=0.0;
  double orig_normsq=dysorb.GetLRNorm();
    
  for(int nc=0; nc<natoms; nc++) {

    std::cout << "Computing XSec for atom=" << nc << std::endl;
    Orbital local_dyson(dysorb,nc);
    std::cout << "Local Dyson orbital:\n";  
    local_dyson.Print();
    double dyson_norm=local_dyson.GetLRNorm();
    total_normsq+=dyson_norm;
    
    double thresh = 1e-6;
    if (dyson_norm > thresh ) {  //only bother if norm is large enough
      
      //Shift origin to atom nc
      pbasis->ShiftGeomToAtom(nc);
      
      //! Calc avg |Cklm|^2 for the given center
      CklmCoeff allCklm; 
      
      if (NOAVG==molavg)
	allCklm.IniCklmCoeff(labgrid,local_dyson,RIonz,charge_radfn[nc],klmgrid);
      else // if (AVG==molavg)
	allCklm.IniCklmCoeff(labgrid,local_dyson,charge_radfn[nc],klmgrid);
      
      
      //Scale factor for cross-sections, incorporating all non-energy dependent coefficients
      double scale=norm*dyson_norm*degeneracyfactor; 
      totalpad.CalcPad(allCklm, scale);
      
      double *xsec1 = new double[nkv];
      memset(xsec1,0,sizeof(double)*nkv);
      for (int k=0; k<nkv; k++) //LOOP BY k
	xsec1[k]=allCklm.GetSumCklm(k)*scale;
      
      //Now add these to the total xsec
      for (int k=0; k<nkv; k++) {
	double kwave=klmgrid.GetKV(k);
	double ene=kwave*kwave/2. + IE;  //energy of laser, in AU
	xsectot[k] += xsec1[k]*kwave; //this is energy factor
      }
      delete [] xsec1;
    }
  }

  std::cout << "\nMulticeneter Dyson LRnorm= " << orig_normsq << ";  sum of local LRnorms=" << total_normsq << " Local/Orig=" <<  total_normsq/orig_normsq << std::endl;  
  PrintXSecAndPad(klmgrid, xsectot,totalpad,IE);
  
  delete [] xsectot;
}

//Here we can do incoherent multi-center expansion
void CalcMultiCenterWithCoherences(MolAvg  molavg, double IE, double degeneracyfactor, double *charge_radfn, double *RIonz, XYZGrid& labgrid, Orbital& dysorb, KLMPoints& klmgrid) {      

  std::cout << "Calculating MC expansion with coherences" <<  std::endl;
  //! Bethe-Salpeter factor of 4*pi^2/c devided by factor from Sakurai of (2*pi). See updated derivation
  double norm=(8*M_PI*M_PI*4)/(3*C_AU);
  
  int nkv = klmgrid.NKV();
  int lmax = klmgrid.LMax();
  int ntheta = 101;
  PAD totalpad(ntheta,nkv);
  
  double *xsectot = new double[nkv];
  memset(xsectot,0,sizeof(double)*nkv);
  
  //We loop over centers and aggregate contributions from multiple centers to the xsec and to the pad
  AOBasis *pbasis=&(TheAOBasis());
  int natoms=pbasis->NAtoms();

  //Initialize Lowdin transformation here, before using for each center transformation
  pbasis->CalcLowdinTransf(labgrid);
  
  //! Calc avg |Cklm|^2 for the given center
  CklmCoeff *allCklm=new CklmCoeff[natoms]; 
  double *do_norms=new double[natoms];
  double thresh = 1e-6;

  double total_normsq=0.0;
  double orig_normsq=dysorb.GetLRNorm();
    
  for(int nc=0; nc<natoms; nc++) {
    
    std::cout << "Computing Cklm for atom=" << nc << std::endl;
    Orbital local_dyson(dysorb,nc);
    std::cout << "Local Dyson orbital:\n";  
    local_dyson.Print();
    do_norms[nc]=local_dyson.GetLRNorm();
    total_normsq+=do_norms[nc];
    
    if (do_norms[nc] > thresh ) {  //only bother if norm is large enough
      
      //Shift origin to atom nc
      pbasis->ShiftGeomToAtom(nc);
      
      if (NOAVG==molavg)
	allCklm[nc].IniCklmCoeff(labgrid,local_dyson,RIonz,charge_radfn[nc],klmgrid);
      else // if (AVG==molavg)
	allCklm[nc].IniCklmCoeff(labgrid,local_dyson,charge_radfn[nc],klmgrid);
    }
  }

  for(int nc1=0; nc1<natoms; nc1++) {

    double coord1[XYZ];
    pbasis->GetGeom(nc1,coord1);
    
    for(int nc2=0; nc2<natoms; nc2++) {

      double coord2[XYZ];
      pbasis->GetGeom(nc2,coord2);

      //Now need to assemble these into PAD and cross sections
      double scale1=norm*do_norms[nc1]*degeneracyfactor;
      double scale2=norm*do_norms[nc2]*degeneracyfactor;

      if(scale1*scale2>thresh) { // &&  nc1==nc2) {

	double *xsec1 = new double[nkv];
	memset(xsec1,0,sizeof(double)*nkv);
	totalpad.CalcMCPadWithCoherences(klmgrid,allCklm[nc1], allCklm[nc2], coord1, coord2, scale1, scale2, xsec1);
	//std::cout << "TEMP PAD from nc1="<< nc1 << " and nc2="<<nc2 <<std::endl;
	//PrintXSecAndPad(klmgrid, xsectot,totalpad,IE);  //Nans for last energy point for H
	for (int k=0; k<nkv; k++) {
	  double kwave=klmgrid.GetKV(k);
	  xsectot[k] += xsec1[k]*kwave; //this is energy factor
	}
	delete [] xsec1;
      }
    }
  }

  std::cout << "\nMulticeneter Dyson LRnorm= " << orig_normsq << ";  sum of local LRnorms=" << total_normsq << " Local/Orig=" <<  total_normsq/orig_normsq << std::endl;
  
  PrintXSecAndPad(klmgrid, xsectot,totalpad,IE);

  delete [] xsectot;
  delete [] do_norms;
  delete [] allCklm;
}


//! Final print of the result
void PrintXSecAndPad(KLMPoints& klmgrid, double *xsectot, PAD& totalpad, double IE) {

  //Now print the final result
  printf("\n\nTotal cross-sections (in bohr^2) vs E(eV)\n");
  printf("See the ezDyson manual for equations and details of the normalization\n");
  printf("E=IE+E_k,eV   xsec,a.u.\n");

  int nkv = klmgrid.NKV();

  for (int k=0; k<nkv; k++) {
    double kwave=klmgrid.GetKV(k);
    double ene=kwave*kwave/2. + IE;  //energy of laser, in AU
    printf("%lf", ene*HAR_TO_EV);
    printf(" %13.6lf\n",xsectot[k]*ene);
  }
  printf("\n\nTotal cross-section / E (in atomic units) vs. Ek (in eV):\n");
  printf("To be used with xsecFCF script for incorporating FCFs\n");
  printf("E_k,eV        xsec/E,a.u. \n");
  for (int k=0; k<nkv; k++) {
    double kwave=klmgrid.GetKV(k);
    double ene=kwave*kwave/2.;  //energy of the electron
    
    printf("%lf", ene*HAR_TO_EV);
    printf("%14.6lf\n",xsectot[k]);
  }

  //Calculate the conversion factors for the xsections
  double au_to_cm_sq = pow(BOHR_TO_ANGS*pow(10.,-8),2);
  double au_to_mb = pow(BOHR_TO_ANGS*pow(10.,1),2);
  printf("\nCross-sections above are in units of bohr^-2\n");
  printf("Useful conversion factors for x-sections:\n");
  printf("Xsec(cm^-2) = Xsec(bohr^-2) * %10lE\n", au_to_cm_sq);
  printf("Xsec(bohr^-2) = Xsec(cm^-2) * %10lE\n", 1./au_to_cm_sq);
  printf("Xsec(Mb) = Xsec(bohr^-2) * %10lE\n", au_to_mb);
  printf("Xsec(bohr^-2) = Xsec(Mb) * %10lE\n", 1./au_to_mb);
 
  printf("\n\nAnisotropy information:\n");
  printf("NOTE: PADs and beta are only correct in molecular frame for Z-polarized light\n");
  printf("____________________________________________________________\n");
  printf("E=IE+E_k,eV   Sigma_par     Sigma_perp    Sigma_tot     beta\n");
  for (int k=0; k<nkv; k++) {
    double kwave=klmgrid.GetKV(k);
    double ene=kwave*kwave/2. + IE;  //energy of laser, in AU
    
    //Bethe salpter.
    //double scale=norm*ene*kwave*dyson_norm*degeneracyfactor; //Put these into PAD when computing it
    double scale=ene*kwave; 
      
    printf("%lf", ene*HAR_TO_EV);
    double sigma_par=totalpad.XSec_par(k)*scale;
    double sigma_perp=totalpad.XSec_perp(k)*scale;
    double sigma_tot=(sigma_par+2.*sigma_perp)*4.0*M_PI/3.;
    double beta=2.*(sigma_par-sigma_perp)/(sigma_par+2.*sigma_perp);
    
    printf(" %13.6lf %13.6lf %13.6lf %13.6lf\n",sigma_par,sigma_perp,sigma_tot,beta);
  }
  
  //Write PADs on disk
  totalpad.Print(klmgrid);
 
}

