#! /usr/bin/env python

"""
/*! \file make_xml.py  
\this scripts creates an input XML file for ezDyson\n

Requires: Q-Chem 5.1+ output file. If you need a script compatible with an older version of Q-Chem, please contact 
sgozem@gsu.edu or krylov@usc.edu.

The script copies the following information from a Q-Chem output to an XML file:
(1) Molecular geometry from the last Q-Chem job in the output;
(2) Atomic basis set on every atom of the molecule in Q-Chem format;
(3) Two Dyson MOs (right and left) for each transition,
(4) If the purecart keyword (pure or cartesian atomic orbitals) is found in the Q-Chem job then purecart is set. Otherwise a manual input is required.

Script assigns default values to: job_parameters, laser polarization and energy, lab_xyz_grid, k_grid, radial_function, degeneracies, and averaging parameters.
Please check and change them accordingly before running ezDyson.

*/
"""

import sys

print ('')

if len(sys.argv) != 3:
    print ('To create a Dyson XML input from a Q-Chem output please type:\n\
    \"make_xml.py <filename.xml> <Q-Chem_dyson_job.out>\" \n')
    sys.exit()

print ('This script will create an XML input for the ezDyson program')

# check that there is one and only one dyson job in the Q-Chem output and Q-Chem is ver 3.2 or above;
n_dyson_jobs=0
if_correct_version="false"
if_trans_prop="false"
qchemF=open(sys.argv[2],'r')
line=qchemF.readline()
while line:
    if (line.lower().find("cc_do_dyson")!=-1):
        n_dyson_jobs+=1
    if (line.lower().find("cc_trans_prop")!=-1):
        if_trans_prop="true"
    if (line.find("Mol. Phys. 113")!=-1):
        if_correct_version="true"
    line=qchemF.readline()
qchemF.close()
print ('Dyson jobs in Q-Chem output: '+str(n_dyson_jobs))
if (n_dyson_jobs!=1):
    print ('Error: Q-Chem output should contain one and only one Dyson job\n')
    sys.exit()
if (if_trans_prop=="false"):
    print ('Error: Q-Chem input did not contain "cc_trans_prop"\n')
    sys.exit()
if (if_correct_version=="false"):
    print ('Error: Q-Chem version should be 4.3 or above\n')
    sys.exit()


# create an xml file
xmlF=open(sys.argv[1],'w')
xmlF.write("""<?xml version="1.0" encoding="ISO-8859-1"?>\n\n""")
xmlF.write("""<root\n  job = "dyson"\n  >\n\n""")


# check if last job in Q-Chem output is restricted/unrestricted (as in the last HF job)
# load molecular geometry (from the last job), count dyson transitions, get transition types and norms,
# read number of basis set functions, try to find purecart
qchemF=open(sys.argv[2],'r')
line=qchemF.readline()
dyson_transitions=0
d_transition=[]
d_norm=[]
n_of_basis_functions=0
tmp_str="RESTRICTED"
if_purecart_found="false"
if_unrestricted="false"
while line:
#    if (line.find('A unrestricted SCF calculation will be')!=-1):
#        if_unrestricted="true"
#        tmp_str="UNRESTRICTED"
    # if geometry found -- load the geometry (overwrite for evry q-chem job, until the last one)
    if (line.find('Standard Nuclear Orientation')>=0):
        geometry=""
        n_of_atoms=0
        qchemF.readline()
        qchemF.readline()
        line=qchemF.readline()
        while line.find('----') == -1:
            n_of_atoms+=1
            geometry+=line[5:]
            line=qchemF.readline()
    if (line.find('Reference -- EOM-IP-CCSD state')!=-1 or line.find('Reference  transition')!=-1) or (line.find('EOM-EE-CCSD state')!=-1) or (line.find('EOM-SF-CCSD state')!=-1):
        d_transition.append("["+line[:-11]+"]")
        d_transition.append("["+line[:-11]+"]")
#        if (if_unrestricted=="true"):
#            d_transition.append("["+line+"]")
#            d_transition.append("["+line+"]")
        qchemF.readline()
        qchemF.readline()
        line=qchemF.readline()
        if (line.find('Dyson orbital')!=-1):
            dyson_transitions+=1
        # find two DMO's norms (or four for unrestricted) (but no DMOs)
        while line.find('End transition') == -1:
            if (line.find('Dyson orbital norm is')!=-1):
                d_norm.append(line[-9:-1])
            line=qchemF.readline()
    # read number of basis functions
    if (line.find('basis functions\n')!=-1) and (n_of_basis_functions==0):
        n_of_basis_functions=line[line.rfind('shells and')+11:line.find('basis functions')-1]
    # try to find purecart
    if (line.lower().find('purecart')!=-1):
        if_purecart_found="true"
        pure_cart=line[9:-1].split()[0]
    line=qchemF.readline()
qchemF.close()

print ('Last Q-Chem job in \"'+sys.argv[2]+'\" was identified as '+tmp_str+'')
print ('Dyson MO transitions found: '+str(dyson_transitions))
if (if_purecart_found=="true"):
    print ('purecart was set to: \"'+pure_cart+'\" (please check if it is correct)')
else:
    print ('No purecart value was found in Q-Chem output. Please set manually (see xml file).\n')
    pure_cart="Please set purecart manually. Q-Chem default value for the Pople and Dunning basis sets are purecart =2222 and =1111 respectively. The order is hgfd, and 1 is pure (spherical) 2 is cartesian representation."

# write job parameters to the XML file
xmlF.write("""<geometry
 n_of_atoms=\""""+str(n_of_atoms)+"""\"
 text = \"
"""+geometry+"""
  \" >
</geometry>

<free_electron
   l_max = \"5\"
   charge_of_ionized_core = \"0\" multicenter=\"false\" coherences=\"false\" >
  <k_grid n_points=\"11\" min=\"0.1\" max=\"10.1\" >
  </k_grid>
  <wave_origin dyson_centroid=\"true\"
    text = " 0.0000    0.0000     0.0000" >
  </wave_origin>
  <charges> text = \" 0.0  0.0  0.0 \"
  </charges>
</free_electron>

<averaging
   method= \"avg\"
   method_possible_values=\"noavg, avg, num\" >
</averaging>

<laser
   ionization_energy = \"1.0\" >
   <laser_polarization x=\"0.0\" y=\"0.0\" z=\"1.0\" >
   </laser_polarization>
</laser>

<lab_xyz_grid>
 <axis n_points=\"201\" min=\"-10.0\" max=\"10.0\" >
 </axis>
</lab_xyz_grid>

<job_parameters
   unrestricted = \"false\"
   Dyson_MO_transitions = \"1\"
   spin_degeneracy = \"2\"
   orbital_degeneracy = \"1\"
   number_of_MOs_to_plot=\"0\"
   MOs_to_plot = \"\" >
</job_parameters> \n\n""")

# read the basis set
qchemF=open(sys.argv[2],'r')
line=qchemF.readline()
while line.find('Basis set in general basis input format') == -1:
    line=qchemF.readline()
qchemF.readline()
qchemF.readline()
basis_set=[]
basis_on_current_atom=""
for atom in range(n_of_atoms):
    line=qchemF.readline()
    while line.find('****') == -1:
        basis_on_current_atom+=line
        line=qchemF.readline()
    # add "****" also
    basis_on_current_atom+=line
    basis_set.append(basis_on_current_atom)
    basis_on_current_atom=""
qchemF.close()
    
# copy general basis
xmlF.write("""<basis
  n_of_basis_functions=\""""+n_of_basis_functions+"""\"
  AO_ordering = "Q-Chem"
  purecart=\""""+pure_cart+"""\">\n""")
for i in range(n_of_atoms):
    xmlF.write("""  <atom\n  text = \"\n""")
    xmlF.write(basis_set[i])
    xmlF.write("""       \" >\n</atom>\n""")
xmlF.write("""</basis>\n\n""")

# write DMOs the XML file
xmlF.write("""<!-- DMOs and MOs BELOW ARE FROM THE \""""+sys.argv[2]+"""\" Q-CHEM OUTPUT -->\n""")
xmlF.write("""<dyson_molecular_orbitals>\n""")

if if_unrestricted=="true":
    dmo_comments=["dyson left alpha", "dyson right alpha"]
else:
    dmo_comments=["dyson right", "dyson left"]

qchemF=open(sys.argv[2],'r')
norm_trans_iter=0
for dyson_transition in range(dyson_transitions):
    for dmo_number in range(len(dmo_comments)):
        #for restricted case only! should be *4 for unrestricted
        norm_trans_iter=dyson_transition*2+dmo_number
#        print 'dyson_tranition equals '+str(dyson_transition)
#        print 'dmo_number, in its turn, equals '+str(dmo_number)
#        print 'norm_trans_iter is '+str(norm_trans_iter)
#        Next 2 lines added to avoid ezDyson crashing for molecules with Cs symmetry with A' and A" irreps
        d_transition[norm_trans_iter]=d_transition[norm_trans_iter].replace("A\'","A1");
        d_transition[norm_trans_iter]=d_transition[norm_trans_iter].replace("A\"","A2");
        xmlF.write("""  <DMO   norm=\""""+d_norm[norm_trans_iter]+"""\"    transition=\""""+d_transition[norm_trans_iter]+"""\"    comment=\""""
                   +dmo_comments[dmo_number]+"""\"\n    text=\"\n""")
        line=qchemF.readline()
        while (line.find('Decomposition over AOs')==-1) and line:
            line=qchemF.readline()
        if not(line):
            print ('\n ERROR: less than expected (less than '+ str( len(dmo_comments)*dyson_transitions ) +') DMOs were found in \"'+sys.argv[2]+'\"\n')
            sys.exit(1)
        line=qchemF.readline()
        while (line.find('*****')==-1):
            xmlF.write(line)
            line=qchemF.readline()
        xmlF.write("""    \" >\n </DMO> \n\n""")
qchemF.close()

xmlF.write("""</dyson_molecular_orbitals>\n\n""")

# Count number of MOs in the q-chem output to the XML file (if any):
number_MOs=0
qchemF=open(sys.argv[2],'r')
line=qchemF.readline()
while line:
    if (line.find('ALPHAMO[')!=-1):
        number_MOs+=1
    line=qchemF.readline()
qchemF.close()


# Copy MOs from the Q-Chem output to the XML file (if any):
xmlF.write("""<molecular_orbitals total_number=\""""+ str(number_MOs)+"""\">\n""")

if (number_MOs!=0):
    qchemF=open(sys.argv[2],'r')

    #copy alpha MOs
    xmlF.write("""  <alpha_MOs>\n""")
    #find start of MOs section:
    line=qchemF.readline()
    while (line.find('ALPHAMO[')==-1):
        line=qchemF.readline()
    number_MOs=0
    while (line.find('$end')==-1):
        number_MOs+=1
        xmlF.write("""    <MO                          type="alpha" number=\""""+str(number_MOs)+"""\"   text=\"\n""")
        line=qchemF.readline()
        while (line.find('*****')==-1):
            xmlF.write(line)
            line=qchemF.readline()
        line=qchemF.readline()
        xmlF.write("""    \" >\n </MO> \n\n""")
    xmlF.write("""  </alpha_MOs>\n\n""")

    #copy beta MOs
    xmlF.write("""  <beta_MOs>\n""")
    #find start of MOs section:
    line=qchemF.readline()
    while (line.find('BETAMO[')==-1):
        line=qchemF.readline()
    number_MOs=0
    while (line.find('$end')==-1):
        number_MOs+=1
        xmlF.write("""    <MO                          type="beta" number=\""""+str(number_MOs)+"""\"   text=\"\n""")
        line=qchemF.readline()
        while (line.find('*****')==-1):
            xmlF.write(line)
            line=qchemF.readline()
        line=qchemF.readline()
        xmlF.write("""    \" >\n </MO> \n\n""")
    xmlF.write("""  </beta_MOs>\n\n""")

    qchemF.close()
xmlF.write("""</molecular_orbitals>\n\n""")


# write a copy of the restored input from the Q-Chem output:
xmlF.write("""<!-- Q-CHEM INPUT IS RESTORED FROM THE \""""+sys.argv[2]+"""\" Q-CHEM OUTPUT -->\n""")
xmlF.write("""<qchem_input>\n""")
separator_string="--------------------------------------------------------------"
xmlF.write(separator_string+"\n")

qchemF=open(sys.argv[2],'r')
number_of_inputs=0

line=qchemF.readline()
while line:
    while (line.find('User input:')==-1) and line:
        line=qchemF.readline()
    #if there is still an input left (not eof):
    if line:
        number_of_inputs+=1
        if (number_of_inputs>1):
            xmlF.write("""@@@\n\n""")
            #find "----------..."line (real mark of the input start)
            while (line.find(separator_string)==-1):
                line=qchemF.readline()
            line=qchemF.readline() #and skip "User input:" line
        line=qchemF.readline() #skip "----------..."line after User input:
        line=qchemF.readline() 
        while (line.find(separator_string)==-1):
            xmlF.write(line)
            line=qchemF.readline()
        xmlF.write("\n")
        line=qchemF.readline()
qchemF.close()
xmlF.write(separator_string+"\n")
xmlF.write("""</qchem_input>\n\n""")


xmlF.write("""</root>\n""")
print ('\nThe following keywords have been assigned default values, please update accordingly (see ezDyson manual for more information):')
print ('Free electron parameters, averaging method, ionization energy, laser polarization.')
print ('Please check the spin and orbital degeneracies and Dyson_MO_transition requested. ALL Dyson orbitals are included, starting from high spin reference to IP/EA.')
print ('Order: Ref-IP1, Ref-IP2, ..., SF1-IP1, SF1-IP2, ..., SF2-IP1, SF2-IP2, ...; EA1-Ref, EA2-Ref, ..., EA1-SF1, EA2-SF1, ..., EA1-SF2, EA2-SF2, ...')
print ('xyz grid parameters should be modified if the total Norm as integrated on the grid is far from 1.00.')
print ('\n'+sys.argv[1]+' has been successfully created\n')


