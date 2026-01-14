#include "aobasis.h"
#include "aik_xml_parser.h"

//For printing
#include <iostream>
#include <iomanip>
//#include <ios>


void AOBasis::Alloc(int nb, int na) 
{
  check(!(nb>0),"AOBasis::Alloc() : Invalid nbasis\n");
  check(!(na>0),"AOBasis::Alloc() : Invalid natoms\n");

  Free();

  nbasis=nb;
  natoms=na;
  basis=new Gauss[nbasis];
  xgeom=new double[natoms];
  ygeom=new double[natoms];
  zgeom=new double[natoms];
  gtos_atom=new int[natoms];
  is_lowdin_ready=false;
}

void AOBasis::Free()
{
  if(IfAlloc())
   {
     delete [] basis;
     delete [] xgeom;;
     delete [] ygeom;
     delete [] zgeom;
     delete [] gtos_atom;
     nbasis=0;
     natoms=0;
     is_lowdin_ready=false;
   }
}

AOBasis::~AOBasis() 
{ 
  Free();
}

void AOBasis::ShiftGeom(double *newcenter)
{
  printf("\nShifting molecular geometry to the new center\n");
  fflush(stdout);
  for(int i=0; i<natoms; i++)
    {
      xgeom[i]-=newcenter[X];
      ygeom[i]-=newcenter[Y];
      zgeom[i]-=newcenter[Z];
    }
  printf("\nNew molecular geometry is:\n");
  printf(" atom         X             Y             Z\n");
  for (int i=0; i<natoms; i++)
    printf("%4d %13lf %13lf %13lf\n",i+1,xgeom[i]/ANGS_TO_BOHR,ygeom[i]/ANGS_TO_BOHR,zgeom[i]/ANGS_TO_BOHR);
}

void AOBasis::ShiftGeomToAtom(int atom) {

  check(!(atom<natoms),"AOBasis::ShiftGeomToAtom(): Invalid atom="+atom);
  double new_origin[XYZ];

  new_origin[X]=xgeom[atom];
  new_origin[Y]=ygeom[atom];
  new_origin[Z]=zgeom[atom];

  ShiftGeom(new_origin);
}

void AOBasis::GetGeom(int at, double *coord) const{

  coord[X]=xgeom[at];
  coord[Y]=ygeom[at];
  coord[Z]=zgeom[at];
}

void AOBasis::PrintAOBasis(FILE *fout)
{
  // Print purecart.//
  fprintf(fout,"Reading purecart: %d%d%d%d\n",purecart[0],purecart[1],purecart[2],purecart[3]);

  //Print standard geom.
  fprintf(fout,"\nStandard geometry\n");
  fprintf(fout," atom         X             Y             Z\n");
  for (int i=0; i<natoms; i++)
    fprintf(fout,"%4d %13lf %13lf %13lf\n",i+1,xgeom[i]/ANGS_TO_BOHR,ygeom[i]/ANGS_TO_BOHR,zgeom[i]/ANGS_TO_BOHR);
  
  // Print AO basis information.
  fprintf(fout,"\nAO basis\n");
  fprintf(fout," atom  AO    NctG ctrG  L    ia   ja   ka       exp           coef        coeff_norm\n");
  for (int n=0,k=0,i=0; i<nbasis; i++)
    {
      for (int j=0; j<basis[i].NContr(); j++)
	fprintf(fout,"%5d%5d%5d%5d%5d%5d%5d%5d%15lE%15lE%15lE\n",n+1,i+1,j+1,basis[i].NContr(),basis[i].AngMom(),
		basis[i].Ia(),basis[i].Ja(),basis[i].Ka(),basis[i].Alpha(j),basis[i].Coeff(j),basis[i].CoeffNorm(j));
      basis[i].CheckAngMom();
      k++;
      if (k == gtos_atom[n])
	{
	  k=0;
         n++;
	}
     }   

  // Print no contracted GTO's per atom
  fprintf(fout,"\n atom   no_basis_fns\n");
  for (int n=0; n<natoms; n++)
    fprintf(fout,"%4d %10d\n",n+1,gtos_atom[n]);
  fprintf(fout,"-----------------------------------------------------------------------------------\n");
}

void AOBasis::CalcOverlapMatrix(const XYZGrid& grid, arma::mat& smatrix) const
{

  smatrix.set_size(nbasis,nbasis);
  smatrix.zeros();
  arma::vec ao_value(nbasis);

  double dV=grid.DXYZ(X)*grid.DXYZ(Y)*grid.DXYZ(Z);
  double xp, yp, zp;
    
  for(int p=0; p<grid.NXYZ(); p++)
    {
      grid.GetPoint(p,xp,yp,zp,X,Y,Z);
      
      //At each grid point get the value of each basis function and compute overlap
  
      //In which order basis functions should be ordered? Judging from Orbital class,
      //(CalcDyson_xyz), they are ordered first for atom=0, then for atom=1, etc
      //n counts atoms, k bf on the atom and i total basis functions
      for(int n=0,k=0,i=0; i<nbasis; i++,k++) {
	if(k==GetGTOsAtom(n)) {
	  n++;
	  k=0;
	}
	ao_value[i] = CalcGaussValue(i,n,xp,yp,zp);
      }

      for(int i=0; i<nbasis; i++)
	for(int j=i; j<nbasis; j++)
	  smatrix(i,j)+=ao_value[i]*ao_value[j];
    }
  
  for(int i=0; i<nbasis; i++)
    for(int j=i; j<nbasis; j++)
      smatrix(i,j)*=dV;

  for(int i=0; i<nbasis; i++)
    for(int j=i+1; j<nbasis; j++)
      smatrix(j,i)=smatrix(i,j);

//Now can print overlap and work with it
//std::cout << "Nbasis=" << nbasis << std::endl;
//std::cout << std::setprecision(3) << std::setw(4);
//smatrix.print("Overlap matrix from AObasis:");
}

//WRONG!
//X_mat: symm_AO(row)=X_min orig_AO(col), M_min_mat: orig_AO(row)=X symm_AO(col)
void AOBasis::CalcLowdinTransf(const XYZGrid& grid) {
  //! Initialize fake Dyson orbital for calculations of the overlap:

  int nbas=NBasis();
  
  std::cout << "\nCalcLowdinTransf()"  <<  std::endl;
  std::cout << "Nbasis=" << nbas << std::endl; 

  arma::mat overlap(nbas,nbas);
  CalcOverlapMatrix(grid,overlap); 

  std::cout <<  std::setw(8) << std::setprecision(5) << std::fixed;
  //std::cout.setf(std::ios::fixed);
  //std::cout.precision(5);
  //overlap.raw_print("S-matrix:");
  overlap.diag().print("Diagonal of the S-matrix:");
  S_mat=overlap; //do not really need it, only for debug purposes

  //Now can compute Lowdin's matrices X and Xminus

  //Now do symmetric orthogonalization
  arma::vec s_eval;
  arma::mat s_evec;
  arma::eig_sym(s_eval,s_evec,overlap);
  s_eval.print("S-matrix eigenvalues:");

  arma::mat s_eval_mat = arma::diagmat(s_eval);
  arma::mat s_inv_sqrt = arma::sqrt(arma::inv(s_eval_mat));
  arma::mat s_sqrt = arma::sqrt(s_eval_mat);
  
  X_mat = s_evec*s_inv_sqrt*s_evec.t(); //S^-1/2
  X_min_mat = s_evec*s_sqrt*s_evec.t(); //S^1/2
   
  //Now these are ready for use!
  is_lowdin_ready=true;
  std::cout << "CalcLowdinTransf(): Done\n"  <<  std::endl;
}


AOBasis& TheAOBasis()
{
  static AOBasis theAOBasis;
  //if(!theAOBasis)
  //theAOBasis=new AOBasis;
  
  //Note there is no check if theAOBasis was initialized
  return theAOBasis;
}

// ======================================================================
//  AOBasis::IniAOBasis()
// ======================================================================

void AOBasis::IniAOBasis(const char* xmlFileName)
{
  static bool if_ini=0;
  check(if_ini,"AOBasis::IniAOBasis() : Already initialized\n");
  if_ini=1;

  std::ifstream xml_file(xmlFileName); 
  xml_node node_geometry("geometry",xml_file);
  int na = node_geometry.read_int_value("n_of_atoms");

  //rewind just in case geometry is placed after basis
  xml_file.clear();
  xml_file.seekg(0, std::ios::beg);
  
  xml_node node_basis("basis",xml_file);
  int nb = node_basis.read_int_value("n_of_basis_functions");
  
  Alloc(nb,na);

  std::string Str;

  //-------------------- GetAtomicXYZ --------------------
  My_istringstream geom_iStr(node_geometry.read_string_value("text"));

  for (int i=0; i<natoms; i++) {
    
      //discard an atomic name:
      geom_iStr.getNextWord(Str);
      //get coordinates & convert to a.u.
      xgeom[i] = geom_iStr.getNextDouble()*ANGS_TO_BOHR;
      ygeom[i] = geom_iStr.getNextDouble()*ANGS_TO_BOHR;
      zgeom[i] = geom_iStr.getNextDouble()*ANGS_TO_BOHR;
      check(geom_iStr.fail(), "IniAOBasis: Error. Wrong format in geometry: ["+geom_iStr.str()+"]\n");
    }

  //-------------------- ReadPureCart --------------------
  int tempStr = node_basis.read_int_value("purecart");

  std::ostringstream convert;
  convert << tempStr;
  Str = convert.str();

  check( (Str.size()>4), "IniAOBasis: Wrong format of \"purecart\": "+Str+". Purecart should be of no more than four digits and each one should be either 1 or 2.\n");

  int n=0; //purecrt index 0..3;
  for (int i=0; i<Str.size(); i++, n++)
    if (Str[i]=='1')
      purecart[i]=1;
    else if (Str[i]=='2')
      purecart[i]=2;
    else {
      std::ostringstream msg;
      //std::string msg;
      msg << "IniAOBasis: Wrong format: \"purecart["<<i+1<<"]\"="<<Str[i]<<" in \"" << Str <<"\"\n"
	      <<"Purecart should be of no more than four digits and each one should be either 1 or 2.\n";
      error(msg.str());
    }
 
  //fill the rest of purecart (default length=4) with the defalt values {1,1,1,1};
  for (int i=n; i<4; i++) {
    
      if (i==3) purecart[i]=1;
      else purecart[i]=1;
  }
  
  std::string AOordering;
  int AOorder = 0;
  AOordering=node_basis.read_string_value("AO_ordering");

  if (AOordering=="Q-Chem")
    AOorder = 1; // Q-Chem AO order
  else if (AOordering=="Molden")
    AOorder = 2; // Molden AO order
  else 
    error("IniAOBasis: AO_ordering must be 1 (for Q-Chem) or 2 (for Molden).\n");

  //-------------------- ReadAOBasis --------------------
  //! number of gaussians in a contracted b.fn.
  int nc; 
  //! amplitude
  double sc; 
  //! arrays of exponents and coefficients of the gaussians in each contraction
  double *exps=NULL, *coefs=NULL, *pcoefs=NULL; 
  
  //! total number of AO b.fns. loaded for the last contracted b.fn.
  int gtos_loaded; 
  //! total number of AO b.fns. on the current atom
  int gtos_on_current_atom; 
  //! total number of AO b.fns. (for the whole molecule)
  int gtos_count=0; 

  //! if a basis set on the current atom is over
  bool ifEndOfBasisOnAtom;

  std::cout << "Number of atoms: " << natoms << std::endl;
  //! for each atom:
  for (int atom_number=0; atom_number<natoms; atom_number++)
    {
      //! read the basis on the n-th atom
      xml_node node_atom(node_basis, "atom", atom_number);
      My_istringstream basis_iStr(node_atom.read_string_value("text") );

      //! read and discard the atom name:
      basis_iStr.getNextWord(Str);
 
      //! read and discard the coef.
      basis_iStr.getNextDouble();

      ifEndOfBasisOnAtom=false;
      gtos_on_current_atom=0;

      //! read first orbital type to Str:
      basis_iStr.getNextWord(Str);
 
     //! load b.fns. on the current atom (**** is the end)
      while ( not(ifEndOfBasisOnAtom) and not(basis_iStr.fail()) )
	{
	  //! read the number of contractions:
	  nc = basis_iStr.getNextInt();
	  sc = basis_iStr.getNextDouble();
	  check(basis_iStr.fail(),"IniAOBasis: Error in ["+basis_iStr.str()+"]\n");

	  exps=new double[nc]; // exponents
	  coefs=new double[nc]; // coefficients
	  pcoefs=new double[nc]; // coefficients for SP

	  //! load AO b.fns.
	  if (Str=="S") {
	    
	      for (int j=0; j<nc; j++)
		{
		  exps[j] = basis_iStr.getNextDouble();
		  coefs[j]= basis_iStr.getNextDouble();
		}
	      basis[gtos_count].IniGauss(0,nc,0,0,0,0,0,0,0,0,0,exps,coefs,sc);
	      gtos_loaded=1;
	    }

	  else if (Str=="P")
	    {
	      for (int j=0; j<nc; j++)
		{
		  exps[j] = basis_iStr.getNextDouble();
		  coefs[j]= basis_iStr.getNextDouble();
		}

	      basis[gtos_count].IniGauss(1,nc,1,0,0,0,0,0,0,0,0,exps,coefs,sc);
	      basis[gtos_count+1].IniGauss(1,nc,0,1,0,0,0,0,0,0,0,exps,coefs,sc);
	      basis[gtos_count+2].IniGauss(1,nc,0,0,1,0,0,0,0,0,0,exps,coefs,sc);
	      gtos_loaded=3;
	    }

	  else if (Str=="D")
	    {
	      for (int j=0; j<nc; j++)
		{
		  exps[j] = basis_iStr.getNextDouble();
		  coefs[j]= basis_iStr.getNextDouble();
		}

	      if (purecart[3] == 1 && AOorder == 1)
		{
		  basis[gtos_count].IniGauss(2,nc,1,1,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+1].IniGauss(2,nc,0,1,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+2].IniGauss(2,nc,0,0,0,0,2,0,0,0,0,exps,coefs,sc);                 
		  basis[gtos_count+3].IniGauss(2,nc,1,0,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+4].IniGauss(2,nc,0,0,0,2,0,0,0,0,0,exps,coefs,sc);
		  gtos_loaded=5;
		}
	      if (purecart[3] == 2 && AOorder == 1)
		{
		  basis[gtos_count].IniGauss(2,nc,2,0,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+1].IniGauss(2,nc,1,1,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+2].IniGauss(2,nc,0,2,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+3].IniGauss(2,nc,1,0,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+4].IniGauss(2,nc,0,1,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+5].IniGauss(2,nc,0,0,2,0,0,0,0,0,0,exps,coefs,sc);
		  gtos_loaded=6;
		}
              if (purecart[3] == 1 && AOorder == 2)
                {
                  basis[gtos_count].IniGauss(2,nc,0,0,0,0,2,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(2,nc,1,0,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(2,nc,0,1,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(2,nc,0,0,0,2,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(2,nc,1,1,0,0,0,0,0,0,0,exps,coefs,sc);
                  gtos_loaded=5;
                }
              if (purecart[3] == 2 && AOorder == 2)
                {
                  basis[gtos_count].IniGauss(2,nc,2,0,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(2,nc,0,2,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(2,nc,0,0,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(2,nc,1,1,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(2,nc,1,0,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(2,nc,0,1,1,0,0,0,0,0,0,exps,coefs,sc);
                  gtos_loaded=6;
                }

	    }
	
	  else if (Str=="F")
	    {
	      for (int j=0; j<nc; j++)
		{
		  exps[j] = basis_iStr.getNextDouble();
		  coefs[j]= basis_iStr.getNextDouble();
		}

	      if (purecart[2] == 1 && AOorder == 1)
		{
		  basis[gtos_count].IniGauss(3,nc,0,0,0,0,0,3,1,0,0,exps,coefs,sc);
		  basis[gtos_count+1].IniGauss(3,nc,1,1,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+2].IniGauss(3,nc,0,0,0,0,0,3,3,0,0,exps,coefs,sc);
		  basis[gtos_count+3].IniGauss(3,nc,0,0,0,0,0,3,4,0,0,exps,coefs,sc);
		  basis[gtos_count+4].IniGauss(3,nc,0,0,0,0,0,3,5,0,0,exps,coefs,sc);
		  basis[gtos_count+5].IniGauss(3,nc,0,0,0,0,0,3,6,0,0,exps,coefs,sc);
		  basis[gtos_count+6].IniGauss(3,nc,0,0,0,0,0,3,7,0,0,exps,coefs,sc);
		  gtos_loaded=7;
		}
	      if (purecart[2] == 2 && AOorder == 1)
		{
		  basis[gtos_count].IniGauss(3,nc,3,0,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+1].IniGauss(3,nc,2,1,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+2].IniGauss(3,nc,1,2,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+3].IniGauss(3,nc,0,3,0,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+4].IniGauss(3,nc,2,0,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+5].IniGauss(3,nc,1,1,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+6].IniGauss(3,nc,0,2,1,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+7].IniGauss(3,nc,1,0,2,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+8].IniGauss(3,nc,0,1,2,0,0,0,0,0,0,exps,coefs,sc);
		  basis[gtos_count+9].IniGauss(3,nc,0,0,3,0,0,0,0,0,0,exps,coefs,sc);
		  gtos_loaded=10;
		}
              if (purecart[2] == 1 && AOorder == 2)
                {
                  basis[gtos_count].IniGauss(3,nc,0,0,0,0,0,3,4,0,0,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(3,nc,0,0,0,0,0,3,5,0,0,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(3,nc,0,0,0,0,0,3,3,0,0,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(3,nc,0,0,0,0,0,3,6,0,0,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(3,nc,1,1,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(3,nc,0,0,0,0,0,3,7,0,0,exps,coefs,sc);
                  basis[gtos_count+6].IniGauss(3,nc,0,0,0,0,0,3,1,0,0,exps,coefs,sc);
                  gtos_loaded=7;
                }
              if (purecart[2] == 2 && AOorder == 2)
                {
                  basis[gtos_count].IniGauss(3,nc,3,0,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(3,nc,0,3,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(3,nc,0,0,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(3,nc,1,2,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(3,nc,2,1,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(3,nc,2,0,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+6].IniGauss(3,nc,1,0,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+7].IniGauss(3,nc,0,1,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+8].IniGauss(3,nc,0,2,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+9].IniGauss(3,nc,1,1,1,0,0,0,0,0,0,exps,coefs,sc);
                  gtos_loaded=10;
                }
	    }

          else if (Str=="G")
            {
              for (int j=0; j<nc; j++)
                {
                  exps[j] = basis_iStr.getNextDouble();
                  coefs[j]= basis_iStr.getNextDouble();
                }

              if (purecart[1] == 1 && AOorder == 1)
                {
                  basis[gtos_count].IniGauss(4,nc,0,0,0,0,0,0,0,4,1,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(4,nc,0,0,0,0,0,0,0,4,2,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(4,nc,0,0,0,0,0,0,0,4,3,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(4,nc,0,0,0,0,0,0,0,4,4,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(4,nc,0,0,0,0,0,0,0,4,5,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(4,nc,0,0,0,0,0,0,0,4,6,exps,coefs,sc);
                  basis[gtos_count+6].IniGauss(4,nc,0,0,0,0,0,0,0,4,7,exps,coefs,sc);
                  basis[gtos_count+7].IniGauss(4,nc,0,0,0,0,0,0,0,4,8,exps,coefs,sc);
                  basis[gtos_count+8].IniGauss(4,nc,0,0,0,0,0,0,0,4,9,exps,coefs,sc);
                  gtos_loaded=9;
                }
              if (purecart[1] == 2 && AOorder == 1)
                {
                  basis[gtos_count].IniGauss(4,nc,4,0,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(4,nc,3,1,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(4,nc,2,2,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(4,nc,1,0,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(4,nc,0,4,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(4,nc,3,0,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+6].IniGauss(4,nc,2,1,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+7].IniGauss(4,nc,1,2,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+8].IniGauss(4,nc,0,3,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+9].IniGauss(4,nc,2,0,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+10].IniGauss(4,nc,1,1,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+11].IniGauss(4,nc,0,2,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+12].IniGauss(4,nc,1,0,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+13].IniGauss(4,nc,0,1,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+14].IniGauss(4,nc,0,0,4,0,0,0,0,0,0,exps,coefs,sc);
                  gtos_loaded=15;
                }
              if (purecart[1] == 1 && AOorder == 2)
                {
                  basis[gtos_count].IniGauss(4,nc,0,0,0,0,0,0,0,4,5,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(4,nc,0,0,0,0,0,0,0,4,6,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(4,nc,0,0,0,0,0,0,0,4,4,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(4,nc,0,0,0,0,0,0,0,4,7,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(4,nc,0,0,0,0,0,0,0,4,3,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(4,nc,0,0,0,0,0,0,0,4,8,exps,coefs,sc);
                  basis[gtos_count+6].IniGauss(4,nc,0,0,0,0,0,0,0,4,2,exps,coefs,sc);
                  basis[gtos_count+7].IniGauss(4,nc,0,0,0,0,0,0,0,4,9,exps,coefs,sc);
                  basis[gtos_count+8].IniGauss(4,nc,0,0,0,0,0,0,0,4,1,exps,coefs,sc);
                  gtos_loaded=9;
                }
              if (purecart[1] == 2 && AOorder == 2)
                {
                  basis[gtos_count].IniGauss(4,nc,4,0,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+1].IniGauss(4,nc,0,4,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+2].IniGauss(4,nc,0,0,4,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+3].IniGauss(4,nc,3,1,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+4].IniGauss(4,nc,3,0,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+5].IniGauss(4,nc,1,0,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+6].IniGauss(4,nc,0,3,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+7].IniGauss(4,nc,1,0,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+8].IniGauss(4,nc,0,1,3,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+9].IniGauss(4,nc,2,2,0,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+10].IniGauss(4,nc,2,0,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+11].IniGauss(4,nc,0,2,2,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+12].IniGauss(4,nc,2,1,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+13].IniGauss(4,nc,1,2,1,0,0,0,0,0,0,exps,coefs,sc);
                  basis[gtos_count+14].IniGauss(4,nc,1,1,2,0,0,0,0,0,0,exps,coefs,sc);
                  gtos_loaded=15;
                }

            }
	  else if (Str=="SP")
	    {
	      for (int j=0; j<nc; j++)
		{
		  exps[j]  = basis_iStr.getNextDouble();
		  coefs[j] = basis_iStr.getNextDouble();
		  pcoefs[j]= basis_iStr.getNextDouble();
		}

	      basis[gtos_count].IniGauss(0,nc,0,0,0,0,0,0,0,0,0,exps,coefs,sc);
	      basis[gtos_count+1].IniGauss(1,nc,1,0,0,0,0,0,0,0,0,exps,pcoefs,sc);
	      basis[gtos_count+2].IniGauss(1,nc,0,1,0,0,0,0,0,0,0,exps,pcoefs,sc);
	      basis[gtos_count+3].IniGauss(1,nc,0,0,1,0,0,0,0,0,0,exps,pcoefs,sc);
	      gtos_loaded=4;
	    }
	  else
	    error("IniAOBasis: Unknown l of atomic orbital: \""+ Str+"\"\n");

	  gtos_on_current_atom+=gtos_loaded;
	  gtos_count+=gtos_loaded;

	  delete [] exps; exps=NULL;
	  delete [] coefs; coefs=NULL;
	  delete [] pcoefs; pcoefs=NULL;

	  //! check if the end of the basis for the atom (marked by "****")
	  basis_iStr.getNextWord(Str);
	  if ( Str=="****" )
	    ifEndOfBasisOnAtom=true;
	  //else: Str is the next orbital type
	}
      gtos_atom[atom_number] = gtos_on_current_atom; 
    }

  //-----------------------------------------------------

  xml_file.close();
  PrintAOBasis();
}




// ======================================================================



