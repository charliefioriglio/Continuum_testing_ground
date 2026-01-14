#ifndef _dyson_main_
#define _dyson_main_

/*! \file dyson_main.h
\brief  the "main" for dyson job
*/

#include "cklm.h"
#include "eikr.h"
#include "orbital.h"
#include "aobasis.h"
#include "xyzgrid.h"
#include "klmgrid.h"
//#include "radialfns.h"
#include "readwrite.h"
#include "ylm.h"
#include "sph.h"
#include "pad.h"
#include "rotnmatr.h"
#include "anglegrid.h"

#include "aik_xml_parser.h"

#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
//For printing
#include <iostream>
#include <iomanip>

//Armadillo include, need for matrices
#include <armadillo>

//! main dyson code
bool dyson_main(const char* xmlFileName);

#endif
