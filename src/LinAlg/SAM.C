// $Id$
//==============================================================================
//!
//! \file SAM.C
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of FE matrices into system matrices.
//!
//==============================================================================

#include "SAM.h"
#include "SystemMatrix.h"
#include <numeric>
#include <iomanip>

#ifdef USE_F77SAM
#if defined(_WIN32)
#define elmeq_  ELMEQ
#define syseq_  SYSEQ
#define addev2_ ADDEV2
#define expand_ EXPAND
#elif defined(_AIX)
#define elmeq_  elmeq
#define syseq_  syseq
#define addev2_ addev2
#define expand_ expand
#endif

extern "C" {
//! \brief Finds the matrix of element equation numbers \a MEEN for an element.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void elmeq_(const int* madof, const int* mnpc, const int* mpmceq,
	    const int* meqn, const int& nenod, int* meen,
	    int& nedof, int& neslv, int &neprd);

//! \brief Determines the control matrix \a MEQN.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void syseq_(const int* msc, const int* mpmceq, const int* mmceq,
	    const int& lpu, int* mpar, int* meqn, int& ierr);

//! \brief Adds an element vector \a eS into the system vector \a sysRHS.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void addev2_(const Real* eS, const Real* ttcc, const int* mpar,
	     const int* madof, const int* meqn, const int* mpmnpc,
	     const int* mmnpc, const int* mpmceq, const int* mmceq,
	     const int& iel, const int& nedof, const int& lpu,
	     Real* sysRHS, int* work, int& ierr);

//! \brief Expands and rearranges a system vector from equation to DOF-order.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void expand_(const Real* solVec, const Real* ttcc,
	     const int* mpmceq, const int* mmceq, const int* meqn,
	     const Real& ff, const Real& fs, const int& ndof, const int& neq,
	     Real* sv);
}
#endif


SAM::SAM () : nnod(mpar[0]), nel(mpar[1]), ndof(mpar[2]),
	      nspdof(mpar[5]), nceq(mpar[6]), neq(mpar[10]),
	      nmmnpc(mpar[14]), nmmceq(mpar[15])
{
  // Initialize the parameters array to zero
  memset(mpar,0,sizeof(mpar));

  // Initialize the array pointers
  mpmnpc = nullptr;
  mmnpc  = nullptr;
  madof  = nullptr;
  msc    = nullptr;
  mpmceq = nullptr;
  mmceq  = nullptr;
  ttcc   = nullptr;
  minex  = nullptr;
  meqn   = nullptr;
}


SAM::~SAM ()
{
  // Deallocate all dynamically allocated arrays
  delete[] mpmnpc;
  delete[] mmnpc;
  delete[] madof;
  delete[] msc;
  delete[] mpmceq;
  delete[] mmceq;
  delete[] ttcc;
  delete[] minex;
  delete[] meqn;
}


std::pair<int,int> SAM::getNodeAndLocalDof (int idof) const
{
  for (int n = 1; n <= nnod; n++)
    if (madof[n] > idof)
      return std::make_pair(minex ? minex[n-1] : n, idof-madof[n-1]+1);

  return std::make_pair(0,0);
}


//! \brief Global stream operator printing a \a (node,localDOF) pair.

std::ostream& operator<< (std::ostream& os, const std::pair<int,int>& p)
{
  return os <<'('<< p.first <<','<< p.second <<')';
}


void SAM::printCEQ (std::ostream& os, int iceq) const
{
  int ip = mpmceq[iceq++]-1;
  std::pair<int,int> nodeDof = this->getNodeAndLocalDof(mmceq[ip]);
  std::string indent(11+log10(nodeDof.first),' ');
  os << nodeDof <<" =";
  if (ttcc[ip] || ip+2 >= mpmceq[iceq])
    os <<" "<< ttcc[ip];
  for (int jp = ip+1; jp < mpmceq[iceq]-1; jp++)
  {
    if ((jp-ip)%10 == 1 && jp > ip+1) os <<"\n"<< indent;
    if (ttcc[ip] || jp > ip+1) os <<" +";
    os <<" "<< ttcc[jp] <<"*"<< this->getNodeAndLocalDof(mmceq[jp]);
  }
}


void SAM::print (std::ostream& os) const
{
  os <<"\n\nSAM::mpar: "<< mpar[0];
  for (int i = 1; i < 30; i++)
    os << ((i%10) ? " " : "\n           ") << mpar[i];
  os << std::endl;

  if (mmnpc && mpmnpc && nel > 0)
  {
    os <<"\n\nElement --> Nodes";
    for (int e = 0; e < nel; e++)
      if (mpmnpc[e+1] > mpmnpc[e])
      {
        os <<'\n'<< std::setw(4) << e+1 <<":";
        for (int i = mpmnpc[e]; i < mpmnpc[e+1]; i++)
          if (mmnpc[i-1] > 0)
            os <<' '<< (minex ? minex[mmnpc[i-1]-1] : mmnpc[i-1]);
          else if (mmnpc[i-1] < 0)
            os <<' '<< (minex ? -minex[-mmnpc[i-1]-1] : mmnpc[i-1]);
      }
    os << std::endl;
  }

  if (ttcc && mmceq && mpmceq && nceq > 0)
  {
    os <<"\n\nConstraint Equations";
    for (int i = 0; i < nceq; i++)
    {
      os <<'\n'<< std::setw(4) << i+1 <<": ";
      this->printCEQ(os,i);
    }
    os << std::endl;
  }

  if (meqn && madof && nnod > 0 && neq > 0)
  {
    char code[7], dofType[7];
    os <<"\n\nNode --> Equations";
    for (int n = 0; n < nnod; n++)
    {
      int dof, i = 0;
      os <<'\n'<< std::setw(4) << (minex ? minex[n] : n+1) <<":";
      for (dof = madof[n]; dof < madof[n+1]; dof++, i++)
      {
	os <<' '<< std::setw(4) << meqn[dof-1];
	code[i] = msc[dof-1] > 0 ? '0' + msc[dof-1] : '0';
	dofType[i] = !dof_type.empty() ? dof_type[dof-1] : '\0';
      }
      code[i] = dofType[i] = '\0';
      os <<'\t'<< code;
      if (n < static_cast<int>(nodeType.size()))
        os <<'\t'<< nodeType[n];
      if (strlen(dofType) > 0)
        os <<'\t' << dofType;
    }
    os << std::endl;
  }
}


void SAM::printStatusCodes (std::ostream& os) const
{
  os <<"\nNode\tDOFs\t MSC";
  for (int i = 0; i < nnod; i++)
  {
    os <<"\n"<< std::setw(4) << 1+i
       <<"\t"<< madof[i] <<" - "<< madof[i+1]-1 <<"\t";
    for (int j = madof[i]; j < madof[i+1]; j++)
      os <<" "<< msc[j-1];
  }
  os << std::endl;
}


bool SAM::initSystemEquations ()
{
#ifdef SP_DEBUG
  std::cout <<"SAM::initSystemEquations()"<< std::endl;
#endif
  if (!msc && ndof > 0) return false;
  if ((!msc || !mpmceq) && nceq > 0) return false;

  // Initialize the DOF-to-equation connectivity array
  int i, j, ierr = 0;
  meqn = new int[ndof];
  memset(meqn,0,ndof*sizeof(int));
#ifdef USE_F77SAM
  syseq_(msc,mpmceq,mmceq,6,mpar,meqn,ierr);
  for (i = j = 0; i < ndof; i++)
    if (msc[i] == 0) msc[i] = --j; // reaction force indices
#else
  int ndof1  = 0;
  int ndof2  = 0;
  int npdof  = 0;
  int nddof  = 0;
  for (i = 0; i < ndof; i++)
    if (msc[i] == 0)
      msc[i] = -(++nspdof); // reaction force indices
    else if (msc[i] == 1)
      ndof1++;
    else if (msc[i] == 2)
      ndof2++;
    else
      ierr--;

  if (ierr < 0)
    std::cerr <<" *** SAM::initSystemEquations: The MSC array has "<< -ierr
              <<" invalid values, should only contain 0, 1 or 2."<< std::endl;

  for (i = 1; i <= nceq; i++)
  {
    int ip = mpmceq[i-1];
    int jp = mpmceq[i]-1;
    int idof = mmceq[ip-1];
    if (idof < 1 || idof > ndof)
    {
      ierr--;
      std::cerr <<" *** SAM::initSystemEquations: idof = "<< idof
                <<" is out of range [1,"<< ndof <<"]."<< std::endl;
      continue;
    }
    else if (msc[idof-1] > 0)
    {
      ierr--;
      std::cerr <<" *** SAM::initSystemEquations: Invalid status code "
                << msc[idof-1] <<" for slave dof "<< idof << std::endl;
    }
    else if (jp == ip)
      npdof++; // prescribed DOF
    else if (jp > ip)
    {
      nddof++; // slave DOF
      while (ip++ < jp)
      {
	int jdof = mmceq[ip-1];
	if (jdof < 1 || jdof > ndof)
	{
          ierr--;
          std::cerr <<" *** SAM::initSystemEquations: jdof = "<< jdof
                    <<" is out of range [1,"<< ndof <<"]."<< std::endl;
	}
	else if (msc[jdof-1] < 1)
	{
          ierr--;
          std::cerr <<" *** SAM::initSystemEquations: Invalid status code "
                    << msc[jdof-1] <<" for master dof "<< jdof << std::endl;
	}
      }
    }
    else
    {
      ierr--;
      std::cerr <<" *** SAM::initSystemEquations: Logic error "<< jp
                <<" < "<< ip <<" for constraint equation "<< i << std::endl;
      break;
    }
    meqn[idof-1] = -i;
  }

  mpar[3] = ndof1;
  mpar[4] = ndof2;
  mpar[7] = nspdof - nceq;
  mpar[8] = npdof;
  mpar[9] = nddof;
  neq = ndof1 + ndof2;

  i = 1;
  j = ndof1+1;
  for (int idof = 0; idof < ndof; idof++)
    if (msc[idof] == 1)
      meqn[idof] = i++;
    else if (msc[idof] == 2)
      meqn[idof] = j++;
#endif

  if (ierr == 0) return true;

  std::cerr <<" *** SAM::initSystemEquations: Failure "<< ierr << std::endl;
#ifdef SP_DEBUG
  this->printStatusCodes(std::cerr);
#endif
  return false;
}


int SAM::getNoNodes (char dofType) const
{
  if (dofType == 'A')
    return nnod;
  else if (nodeType.empty())
    return dofType == 'D' ? nnod : 0;

  int n = 0;
  for (int i = 0; i < nnod; i++)
    if (nodeType[i] == dofType) n++;

  return n;
}


void SAM::getDofCouplings (IntVec& irow, IntVec& jcol) const
{
  irow.clear();
  jcol.clear();

  // Find the set of DOF couplings for each free DOF
  std::vector<IntSet> dofc;
  this->getDofCouplings(dofc);

  irow.resize(neq+1);
  irow.front() = 0;

  // Find total number of dof couplings or non-zeroes in the system matrix
  int i;
  for (i = 0; i < neq; i++)
    irow[i+1] = irow[i] + dofc[i].size();

  jcol.reserve(irow.back());
  for (i = 0; i < neq; i++)
    jcol.insert(jcol.end(),dofc[i].begin(),dofc[i].end());
}


void SAM::getDofCouplings (std::vector<IntSet>& dofc) const
{
  dofc.resize(neq);
  for (int iel = 1; iel <= nel; iel++)
  {
    IntVec meen;
    this->getElmEqns(meen,iel);

    for (size_t j = 0; j < meen.size(); j++)
    {
      int jeq = meen[j];
      if (jeq > 0)
      {
        dofc[jeq-1].insert(jeq);
        for (size_t i = 0; i < j; i++)
        {
          int ieq = meen[i];
          if (ieq > 0)
          {
            dofc[ieq-1].insert(jeq);
            dofc[jeq-1].insert(ieq);
          }
        }
      }
      else if (jeq < 0)
      {
        int jpmceq1 = mpmceq[-jeq-1];
        int jpmceq2 = mpmceq[-jeq]-1;
        for (int jp = jpmceq1; jp < jpmceq2; jp++)
          if (mmceq[jp] > 0)
          {
            jeq = meqn[mmceq[jp]-1];
            for (size_t i = 0; i < meen.size(); i++)
            {
              int ieq = meen[i];
              if (ieq > 0)
              {
                dofc[ieq-1].insert(jeq);
                dofc[jeq-1].insert(ieq);
              }
              else if (ieq < 0)
              {
                int ipmceq1 = mpmceq[-ieq-1];
                int ipmceq2 = mpmceq[-ieq]-1;
                for (int ip = ipmceq1; ip < ipmceq2; ip++)
                  if (mmceq[ip] > 0)
                  {
                    ieq = meqn[mmceq[ip]-1];
                    dofc[ieq-1].insert(jeq);
                  }
              }
            }
          }
      }
    }
  }
}


bool SAM::assembleSystem (SystemMatrix& sysK, SystemVector& sysRHS,
                          const Matrix& eK, int iel,
                          RealArray* reactionForces) const
{
  if (eK.empty())
    return true; // silently ignore empty element matrices

  if (reactionForces)
  {
    Vector eS;
    IntVec meen;
    if (!this->getElmEqns(meen,iel,eK.rows()))
      return false;

    // Add (appropriately weighted) terms corresponding to constrained
    // (dependent and prescribed) DOFs in eK into reactionForces
    for (size_t j = 1; j <= meen.size(); j++)
    {
      int jceq = -meen[j-1];
      if (jceq < 1) continue;

      int jp = mpmceq[jceq-1];
      Real c0 = ttcc[jp-1];

      for (size_t i = 1; i <= meen.size(); i++)
	if (meen[i-1] < 1)
	{
	  eS.resize(eK.rows());
	  eS(i) = -c0*eK(i,j);
	}
    }

    if (!eS.empty())
      this->assembleReactions(*reactionForces,eS,iel);
  }

  return sysK.assemble(eK,*this,sysRHS,iel);
}


bool SAM::assembleSystem (SystemMatrix& sysM,
			  const Matrix& eM, int iel) const
{
  return eM.empty() ? true : sysM.assemble(eM,*this,iel);
}


bool SAM::assembleSystem (SystemVector& sysRHS,
                          const Matrix& eK, int iel,
                          RealArray* reactionForces) const
{
  if (eK.empty())
    return true; // silently ignore empty element matrices

  IntVec meen;
  if (!this->getElmEqns(meen,iel,eK.rows()))
    return false;

  // Add (appropriately weighted) terms corresponding to constrained
  // (dependent and prescribed) DOFs in eK into sysRHS and reactionForces
  Vector eS;
  for (size_t j = 1; j <= meen.size(); j++)
  {
    int jceq = -meen[j-1];
    if (jceq < 1) continue;

    int jp = mpmceq[jceq-1];
    Real c0 = ttcc[jp-1];

    for (size_t i = 1; i <= meen.size(); i++)
    {
      int ieq = meen[i-1];
      this->assembleRHS(sysRHS.getPtr(),-c0*eK(i,j),ieq);
      if (reactionForces && ieq < 1)
      {
        eS.resize(eK.rows());
        eS(i) = -c0*eK(i,j);
      }
    }
  }

  if (reactionForces && !eS.empty())
    this->assembleReactions(*reactionForces,eS,iel);

  return true;
}


bool SAM::assembleSystem (SystemVector& sysRHS,
                          const RealArray& eS, int iel,
                          RealArray* reactionForces) const
{
  if (eS.empty())
    return true; // silently ignore empty element vectors

  int ierr = 0;
#ifdef USE_F77SAM
  int* work = new int[eS.size()];
  addev2_(&eS.front(), ttcc, mpar, madof, meqn, mpmnpc, mmnpc, mpmceq, mmceq,
	  iel, eS.size(), 6, sysRHS.getPtr(), work, ierr);
  delete[] work;
#else
  IntVec meen;
  if (!this->getElmEqns(meen,iel,eS.size()))
    ierr = 1;
  else for (size_t i = 0; i < meen.size() && i < eS.size(); i++)
    this->assembleRHS(sysRHS.getPtr(),eS[i],meen[i]);
#endif

  if (reactionForces)
    this->assembleReactions(*reactionForces,eS,iel);

  return ierr == 0;
}


bool SAM::assembleSystem (SystemVector& sysRHS, const Real* nS, int inod,
                          RealArray* reactionForces) const
{
  IntVec mnen;
  if (!this->getNodeEqns(mnen,inod))
    return false;

  for (size_t i = 0; i < mnen.size(); i++)
    this->assembleRHS(sysRHS.getPtr(),nS[i],mnen[i]);

  if (reactionForces)
  {
    int k = 0;
    for (int j = madof[inod-1]; j < madof[inod]; j++, k++)
    {
      int ipR = -msc[j-1];
      if (ipR > 0 && ipR <= static_cast<int>(reactionForces->size()))
        (*reactionForces)[ipR-1] += nS[k];
    }
  }

  return true;
}


bool SAM::assembleSystem (SystemVector& sysRHS, Real S,
                          const std::pair<int,int>& dof) const
{
  if (dof.first < 1 || dof.first > nnod)
  {
    std::cerr <<" *** SAM::assembleSystem: Node "<< dof.first
              <<" is out of range [1,"<< nnod <<"]."<< std::endl;
    return false;
  }

  int idof = madof[dof.first-1] + dof.second-1;
  if (dof.second < 1 || idof >= madof[dof.first])
  {
    int nndof = madof[dof.first] - madof[dof.first-1];
    std::cerr <<" *** SAM::assembleSystem: Local dof "<< dof.second
              <<" is out of range [1,"<< nndof <<"]."<< std::endl;
    return false;
  }

  this->assembleRHS(sysRHS.getPtr(),S,meqn[idof-1]);

  return true;
}


void SAM::addToRHS (SystemVector& sysRHS, const RealArray& S) const
{
  if (ndof < 1 || S.empty()) return;

  int n = S.size();
  if (n > ndof) n = ndof;

  for (int i = 0; i < n; i++)
    this->assembleRHS(sysRHS.getPtr(),S[i],meqn[i]);
}


void SAM::assembleRHS (Real* RHS, Real value, int ieq) const
{
  int iceq = -ieq;
  if (ieq > 0)
    RHS[ieq-1] += value;
  else if (iceq > 0)
    for (int ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ip++)
      if (mmceq[ip] > 0)
      {
        ieq = meqn[mmceq[ip]-1];
        RHS[ieq-1] += ttcc[ip]*value;
      }
}


void SAM::assembleReactions (RealArray& rf, const RealArray& eS, int iel) const
{
  size_t k = 0;
  int ip = mpmnpc[iel-1];
  int nenod = mpmnpc[iel] - ip;
  for (int i = 0; i < nenod; i++, ip++)
  {
    int node = mmnpc[ip-1];
    if (node < 0)
      k += madof[-node] - madof[-node-1];
    else if (node > 0)
      for (int j = madof[node-1]; j < madof[node] && k < eS.size(); j++, k++)
      {
        int ipR = -msc[j-1];
        if (ipR > 0 && ipR <= static_cast<int>(rf.size()))
          rf[ipR-1] += eS[k];
      }
  }
}


bool SAM::getElmNodes (IntVec& mnpc, int iel) const
{
  mnpc.clear();
  if (iel < 1 || iel > nel)
  {
    std::cerr <<" *** SAM::getElmNodes: Element "<< iel <<" is out of range [1,"
              << nel <<"]."<< std::endl;
    return false;
  }

  int ip = mpmnpc[iel-1]-1;
  int jp = mpmnpc[iel]-1;
  if (jp <= ip) return true;

  mnpc.reserve(jp-ip);
  mnpc.insert(mnpc.end(),mmnpc+ip,mmnpc+jp);
  return true;
}


bool SAM::getElmEqns (IntVec& meen, int iel, size_t nedof) const
{
  meen.clear();
  if (iel < 1 || iel > nel)
  {
    std::cerr <<" *** SAM::getElmEqns: Element "<< iel <<" is out of range [1,"
              << nel <<"]."<< std::endl;
    return false;
  }

  int ip = mpmnpc[iel-1];
  int nenod = mpmnpc[iel] - (ip--);
  if (nenod < 1) return true;

  int neldof = 0;
  if (nedof == 0) // Calculate element size from the nodal connectivity
    neldof = std::accumulate(mmnpc+ip, mmnpc+ip+nenod, 0,
                             [this](const int& offs, const int& inod) {
                               int node = inod > 0 ? inod : -inod;
                               return offs + madof[node] - madof[node-1]; });

#ifdef USE_F77SAM
  int neslv, neprd;
  meen.resize(nedof > 0 ? nedof : neldof, 0);
  elmeq_(madof,mmnpc+ip,mpmceq,meqn,nenod,&meen.front(),neldof,neslv,neprd);
  if (neldof < static_cast<int>(meen.size())) meen.resize(neldof);
#else
  meen.reserve(nedof > 0 ? nedof : neldof);
  for (int i = 0; i < nenod; i++, ip++)
  {
    int node = mmnpc[ip];
    if (node > 0)
      meen.insert(meen.end(),meqn+madof[node-1]-1,meqn+madof[node]-1);
    else if (node < 0)
      meen.insert(meen.end(),madof[-node]-madof[-node-1],0);
  }
  neldof = meen.size();
#endif
  if (nedof == 0 || neldof == static_cast<int>(nedof)) return true;

  std::cerr <<" *** SAM::getElmEqns: Invalid element matrix dimension "
            << nedof <<" for element "<< iel
            <<" (should have been "<< neldof <<")."<< std::endl;
  return false;
}


void SAM::getUniqueEqns (IntSet& meen, int iel) const
{
  meen.clear();

  IntVec tmp;
  if (!this->getElmEqns(tmp,iel))
    return; // only if element index is out of range

  for (int jeq : tmp)
    if (jeq < 0)
    {
      int jpmceq1 = mpmceq[-jeq-1];
      int jpmceq2 = mpmceq[-jeq]-1;
      for (int jp = jpmceq1; jp < jpmceq2; ++jp)
        if (mmceq[jp] > 0)
          meen.insert(meqn[mmceq[jp]-1]);
    }
    else if (jeq > 0)
      meen.insert(jeq);
}


bool SAM::getNodeEqns (IntVec& mnen, int inod) const
{
  mnen.clear();
  if (inod < 1 || inod > nnod)
  {
    std::cerr <<" *** SAM::getNodeEqns: Node "<< inod
              <<" is out of range [1,"<< nnod <<"]."<< std::endl;
    return false;
  }

  mnen.reserve(madof[inod]-madof[inod-1]);
  mnen.insert(mnen.end(),meqn+madof[inod-1]-1,meqn+madof[inod]-1);

  return true;
}


std::pair<int,int> SAM::getNodeDOFs (int inod) const
{
  if (inod < 1 || inod > nnod) return { 0, 0 };

  return { madof[inod-1], madof[inod]-1 };
}


int SAM::getEquation (int inod, int ldof) const
{
  if (inod < 1 || inod > nnod || ldof < 1) return -1;

  int idof = madof[inod-1] + ldof-1;
  if (idof >= madof[inod]) return -2;

  int ieq = meqn[idof-1];
  return ieq > 0 ? ieq : 0;
}


bool SAM::expandSolution (const SystemVector& solVec, Vector& dofVec,
			  Real scaleSD) const
{
  if (static_cast<int>(solVec.dim()) < neq) return false;

  return this->expandVector(solVec.getRef(),dofVec,scaleSD);
}


bool SAM::expandVector (const Vector& solVec, Vector& dofVec) const
{
  if (static_cast<int>(solVec.size()) < neq) return false;

  return this->expandVector(solVec.ptr(),dofVec,Real(0));
}


bool SAM::expandVector (const Real* solVec, Vector& dofVec, Real scaleSD) const
{
  if (!meqn) return false;

  dofVec.resize(ndof,true);
#ifdef USE_F77SAM
  expand_(solVec, ttcc, mpmceq, mmceq, meqn, Real(1), scaleSD, ndof, neq,
	  dofVec.ptr());
#else
  for (int idof = 0; idof < ndof; idof++)
  {
    int ieq = meqn[idof];
    int iceq = -ieq;
    if (ieq > 0)
      dofVec[idof] += solVec[ieq-1];
    else if (iceq > 0)
    {
      int ip = mpmceq[iceq-1];
      dofVec[idof] += scaleSD*ttcc[ip-1];
      for (; ip < mpmceq[iceq]-1; ip++)
	if (mmceq[ip] > 0)
	{
	  ieq = meqn[mmceq[ip]-1];
	  if (ieq > 0 && ieq <= neq)
	    dofVec[idof] += ttcc[ip]*solVec[ieq-1];
	}
    }
  }
#endif

  return true;
}


bool SAM::getSolVec (RealArray& solVec, const RealArray& dofVec) const
{
  if (!meqn) return false;

  solVec.resize(neq,Real(0));
  for (int idof = 0; idof < ndof; idof++)
  {
    int ieq = meqn[idof];
    if (ieq > 0)
      solVec[ieq-1] = dofVec[idof];
  }

  return true;
}


bool SAM::applyDirichlet (Vector& dofVec) const
{
  if (!meqn) return false;

  for (int idof = 0; idof < ndof; idof++)
  {
    int iceq = -meqn[idof];
    if (iceq > 0)
    {
      int ip = mpmceq[iceq-1];
      dofVec[idof] = ttcc[ip-1];
    }
    else if (iceq == 0)
      dofVec[idof] = Real(0);
  }

  return true;
}


Real SAM::dot (const Vector& x, const Vector& y, char dofType) const
{
  if ((nodeType.empty() && dof_type.empty()) || dofType == 'A')
    return x.dot(y); // All DOFs are of the same type, or consider all of them

  // Consider only the dofType DOFs
  int i, j, n = x.size() < y.size() ? x.size() : y.size();
  Real retVal = Real(0);
  for (i = 0; i < nnod; i++)
    if (nodeType.empty() || nodeType[i] == dofType)
      for (j = madof[i]-1; j < madof[i+1]-1 && j < n; j++)
        if (dof_type.empty() || dof_type[j] == dofType)
          retVal += x[j]*y[j];

  return retVal;
}


Real SAM::normL2 (const Vector& x, char dofType) const
{
  if (x.empty())
    return Real(0);
  else if (nodeType.empty() && dof_type.empty())
    return x.norm2()/sqrt(x.size()); // All DOFs are of the same type

  // Consider only the dofType DOFs
  int count, i, j, n = x.size();
  Real retVal = Real(0);
  for (count = i = 0; i < nnod; i++)
    if (nodeType.empty() || nodeType[i] == dofType)
      for (j = madof[i]-1; j < madof[i+1]-1 && j < n; j++)
        if (dof_type.empty() || dof_type[j] == dofType)
        {
          retVal += x[j]*x[j];
          count ++;
        }

  return count > 0 ? sqrt(retVal/count) : retVal;
}


Real SAM::normInf (const Vector& x, size_t& comp, char dofType) const
{
  if (x.empty() || comp < 1)
    return Real(0);
  else if (nodeType.empty() && dof_type.empty())
  {
    // All DOFs are of the same type and the number of nodal DOFs is constant
    int nndof = madof[1] - madof[0];
    if (static_cast<int>(comp) <= nndof)
      return x.normInf(--comp,nndof);
    else
      return Real(0);
  }

  // Consider only the dofType DOFs
  size_t icmp = comp;
  Real retVal = Real(0);
  for (int i = 0; i < nnod; i++)
    if (!dof_type.empty())
    {
      size_t k = 0;
      for (int idof = madof[i]-1; idof+1 < madof[i+1]; idof++)
        if (dof_type[idof] == dofType && ++k == icmp)
          if (fabs(x[idof]) > retVal)
          {
            comp = i+1;
            retVal = fabs(x[idof]);
          }
    }
    else if (nodeType[i] == dofType)
    {
      int idof = madof[i] + icmp - 2;
      if (idof >= 0 && idof < madof[i+1]-1)
        if (fabs(x[idof]) > retVal)
        {
          comp = i+1;
          retVal = fabs(x[idof]);
        }
    }

  return retVal;
}


Real SAM::normReact (const RealArray& u, const RealArray& rf) const
{
  Real retVal = Real(0);

  for (int i = 0; i < ndof; i++)
    if (meqn[i] < 0 && msc[i] < 0 && -msc[i] <= static_cast<int>(rf.size()))
      if (mpmceq[-meqn[i]] - mpmceq[-meqn[i]-1] == 1) // Only prescribed DOFs
      {
        Real RF = rf[-msc[i]-1];
        retVal += u[i]*RF;
#if SP_DEBUG > 1
        std::cout <<"SAM::normReact: idof="<< i+1 <<" SC="<< msc[i]
                  <<" u="<< u[i] <<" RF="<< RF <<" --> "<< retVal << std::endl;
#endif
      }

  return retVal;
}


bool SAM::haveReaction (int dir, const IntVec* nodes) const
{
  if (dir > 0 && nspdof > 0)
    for (int i = 0; i < nnod; i++)
      if (!nodes || std::find(nodes->begin(),nodes->end(),i+1) != nodes->end())
      {
        int idof = madof[i]+dir-2;
        if (idof < madof[i+1]-1 && msc[idof] < 0 && -msc[idof] <= nspdof)
          return true;
      }

  return false;
}


Real SAM::getReaction (int dir, const RealArray& rf, const IntVec* nodes) const
{
  Real retVal = Real(0);

  if (dir > 0 && nspdof > 0)
    for (int i = 0; i < nnod; i++)
      if (!nodes || std::find(nodes->begin(),nodes->end(),i+1) != nodes->end())
      {
        int idof = madof[i]+dir-2;
        int ipr = idof < madof[i+1]-1 ? -msc[idof]-1 : -1;
        if (ipr >= 0 && ipr < static_cast<int>(rf.size()))
        {
          retVal += rf[ipr];
          if (nodes) // clear this force component to avoid it being added twice
            const_cast<RealArray&>(rf)[ipr] = 0.0;
        }
      }

  return retVal;
}


bool SAM::getNodalReactions (int inod, const RealArray& rf,
                             RealArray& nrf) const
{
  if (inod < 1 || inod > nnod)
  {
    std::cerr <<" *** SAM::getNodalReactions: Node "<< inod
              <<" is out of range [1,"<< nnod <<"]."<< std::endl;
    return false;
  }

  bool haveRF = false;
  int ip = madof[inod-1]-1;
  nrf.resize(madof[inod]-1-ip);
  for (size_t i = 0; i < nrf.size(); i++, ip++)
    if (msc[ip] < 0 && -msc[ip] <= static_cast<int>(rf.size()))
    {
      haveRF = true;
      nrf[i] = rf[-msc[ip]-1];
    }
    else
      nrf[i] = 0.0;

  return haveRF;
}


IntSet SAM::getEquations (char dType, int dof) const
{
  IntSet result;
  auto&& getEquation = [this,dType,&result](int idof)
  {
    if (dof_type.empty() || dof_type[idof] == dType)
      if (meqn[idof] > 0)
        result.insert(meqn[idof]);
  };

  for (int inod = 0; inod < nnod; inod++)
    if ((nodeType.empty() ? 'D' : nodeType[inod]) == dType || !dof_type.empty())
    {
      if (dof > 0 && dof <= madof[inod+1]-madof[inod])
        getEquation(madof[inod]+dof-2);
      else for (int idof = madof[inod]; idof < madof[inod+1]; idof++)
        getEquation(idof-1);
    }

  return result;
}
