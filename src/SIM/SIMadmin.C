// $Id$
//==============================================================================
//!
//! \file SIMadmin.C
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Administration base class for FEM simulators.
//!
//==============================================================================

#include "SIMadmin.h"
#include "ProcessAdm.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"
#include <fstream>

#if defined(__MINGW32__) || defined(__MINGW64__)
#define strcasestr strstr
#endif

int SIMadmin::msgLevel = 2;


SIMadmin::SIMadmin (const char* heading) : opt(myOpts), adm(nullptr)
{
#if defined(HAS_PETSC) || defined(HAVE_MPI)
  adm   = new ProcessAdm(true);
  myPid = adm->getProcId();
  nProc = adm->getNoProcs();
#else
  myPid = 0;
  nProc = 1;
#endif

  myOpts = IFEM::getOptions(); // Initialize options from command-line arguments

  if (heading) myHeading = heading;
}


SIMadmin::SIMadmin (SIMadmin& anotherSIM) : opt(anotherSIM.myOpts), adm(nullptr)
{
  if (anotherSIM.adm)
    adm = new ProcessAdm(*anotherSIM.adm); //TODO: Consider sharing the adm too
  myPid = anotherSIM.myPid;
  nProc = anotherSIM.nProc;
}


SIMadmin::~SIMadmin ()
{
  if (adm) delete adm;
}


void SIMadmin::printHeading (int& subStep) const
{
  if (myHeading.empty())
    return;

  size_t n = myHeading.find_last_of('\n');
  if (n+1 < myHeading.size()) n = myHeading.size()-n;
  IFEM::cout <<"\n"<< ++subStep <<". "<< myHeading <<"\n";
  for (size_t i = 0; i < 3+n && i < 3+myHeading.size(); i++) IFEM::cout <<'=';
  if (subStep > 9) IFEM::cout <<'=';
  IFEM::cout << std::endl;
}


bool SIMadmin::read (const char* fileName)
{
  static int substep = 0;
  this->printHeading(substep);

  bool result;
  if (strcasestr(fileName,".xinp"))
    result = this->readXML(fileName);
  else
    result = this->readFlat(fileName);

  // Let command-line options override settings on the input file
  IFEM::applyCommandLineOptions(opt);

  return result;
}


bool SIMadmin::readFlat (const char* fileName)
{
  std::ifstream is(fileName);
  if (!is)
  {
    std::cerr <<"\n *** SIMadmin::read: Failure opening input file \""
              << fileName <<"\"."<< std::endl;
    return false;
  }

  IFEM::cout <<"\nReading input file "<< fileName << std::endl;

  char* keyWord = nullptr;
  while (is.good() && (keyWord = utl::readLine(is)))
    if (!this->parse(keyWord,is))
    {
      std::cerr <<" *** SIMadmin::read: Failure occured while parsing \""
                << keyWord <<"\"."<< std::endl;
      return false;
    }

  IFEM::cout <<"\nReading input file succeeded."<< std::endl;

  return true;
}


bool SIMadmin::parse (char*, std::istream&)
{
  std::cerr <<" *** SIMadmin::parse(char*,std::istream&):"
            <<" The flat file format is depreciated.\n"
            <<"     Use the XML format instead."<< std::endl;
  return false;
}


const ProcessAdm& SIMadmin::getProcessAdm () const
{
  if (adm) return *adm;

  // This is a serial build - no parallel administrator needed, return a dummy
#ifdef SP_DEBUG
  std::cerr <<"  ** SIMadmin::getProcessAdm: Process administrator requested"
            <<" in a serial build. Check logic.."<< std::endl;
#endif
  static ProcessAdm dummy;
  return dummy;
}


utl::LogStream& SIMadmin::getLogStream () const
{
  if (adm) return adm->cout;

  // When no parallel administrator, the log stream defaults to std::cout
  static utl::LogStream sout(std::cout);
  return sout;
}


int SIMadmin::getProcId () const
{
  return adm ? adm->getProcId() : 0;
}
