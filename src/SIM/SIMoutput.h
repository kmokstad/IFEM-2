// $Id$
//==============================================================================
//!
//! \file SIMoutput.h
//!
//! \date Sep 05 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Sub-class with functionality for result output to %VTF and terminal.
//!
//==============================================================================

#ifndef _SIM_OUTPUT_H
#define _SIM_OUTPUT_H

#include "SIMinput.h"
#include "Vec3.h"

class ElementBlock;
class VTF;


/*!
  \brief Sub-class with additional functionality for result output.
  \details This class extends the SIMbase class with some added functionalities
  for dumping simulation results to %VTF and ASCII files, and terminal printout.
  These items are put in a separate sub-class to hide them from the SIMbase
  class, which contains the main simulation driver.
*/

class SIMoutput : public SIMinput
{
protected:
  //! \brief The constructor just forwards to the base class constructor.
  explicit SIMoutput(IntegrandBase* itg);

public:
  //! \brief The destructor frees the dynamically allocated VTF object.
  virtual ~SIMoutput();

  //! \brief Initializes the property containers of the model.
  virtual void clearProperties();

  using SIMinput::parse;
  //! \brief Parses a data section from an input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is);

  //! \brief Merges the global equation system of \a that simulator with this.
  //! \param that The simulator whose equation system is to be merged
  //! \param[in] old2new Global node number mapping
  //! \param[in] poff Global patch index offset
  virtual bool merge(SIMbase* that, const std::map<int,int>* old2new, int poff);

  //! \brief Retrieves the specified element set as a scalar field.
  //! \param[in] iset 1-based index of the element set to retrieve
  //! \param[out] name Name of the retrieved element set
  //! \param[out] elSet 1-based element index array indicating the set members
  //! \return \e false if the \a iset index is out of range, otherwise \e true
  bool getElementSet(int iset, std::string& name, RealArray& elSet) const;

protected:
  //! \brief Parses a subelement of the \a resultoutput XML-tag.
  virtual bool parseOutputTag(const tinyxml2::XMLElement* elem);

  //! \brief Tesselates the specified patch.
  virtual ElementBlock* tesselatePatch(size_t pidx) const;

public:
  //! \brief Writes current model geometry to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] inpFile File name used to construct the VTF-file name from
  //! \param[in] doClear If \e true, clear geometry block if \a inpFile is null
  //!
  //! \details The spline patches are tesselated into linear finite elements
  //! with a fixed number of elements within each knot-span of non-zero length.
  //! The solution fields are then evaluated at the nodal points of the
  //! generated FE mesh and written to the VTF-file as vector and scalar fields
  //! by the other \a writeGlv* methods.
  virtual bool writeGlvG(int& nBlock, const char* inpFile, bool doClear = true);

  //! \brief Writes current model geometry to the currently open VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] time The time from which this (new) geometry applies,
  //! in case of time-evolution. If negative, the current geometry is appended
  //! to the already existing geometry blocks.
  bool writeGlvG(int& nBlock, double time);

  //! \brief Writes additional, problem-specific, results to the VTF-file.
  virtual bool writeGlvA(int&, int, int = 1) const { return true; }

  //! \brief Writes boundary conditions as scalar fields to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param[in] iStep Load/time step identifier
  virtual bool writeGlvBC(int& nBlock, int iStep = 1) const;

  //! \brief Writes global node numbers as scalar fields to the VTF-file.
  //! \param nBlock Running result block counter
  //! \param idBlock Result block ID number
  //! \param[in] maxBlock Highest block ID number
  //! \param[in] iStep Load/time step identifier
  //!
  //! \details If you have a complex model with lots of nodal points,
  //! this might be used as a tool to identify the location of the nodes
  //! by visualizing their global number and the original node numbers.
  bool writeGlvNo(int& nBlock, int& idBlock,
                  int maxBlock = 20, int iStep = 1) const;

  //! \brief Writes boundary tractions for a given time step to the VTF-file.
  //! \param[in] iStep Load/time step identifier
  //! \param geoBlk Running geometry block counter
  //! \param nBlock Running result block counter
  bool writeGlvT(int iStep, int& geoBlk, int& nBlock) const;

  //! \brief Writes a vector field for a given load/time step to the VTF-file.
  //! \param[in] vec The vector field to output (nodal values in DOF-order)
  //! \param[in] fieldName Name identifying the vector field
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Result block ID number
  //! \param[in] ncmp Number of components in vector field
  bool writeGlvV(const RealArray& vec, const char* fieldName,
                 int iStep, int& nBlock, int idBlock = 2, int ncmp = 0) const;

  //! \brief Writes a scalar field for a given load/time step to the VTF-file.
  //! \param[in] scl The scalar field to output (nodal values)
  //! \param[in] fieldName Name identifying the scalar field
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Result block ID number
  bool writeGlvS(const Vector& scl, const char* fieldName,
                 int iStep, int& nBlock, int idBlock = 2) const;

  //! \brief Writes solution fields for a given load/time step to the VTF-file.
  //! \param[in] psol Primary solution vector
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] time Load/time step parameter
  //! \param[in] pvecName Optional name of the primary vector field solution
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] psolComps Optional number of primary solution components
  bool writeGlvS(const Vector& psol, int iStep, int& nBlock, double time = 0.0,
                 const char* pvecName = nullptr, int idBlock = 10,
                 int psolComps = 0);

  //! \brief Writes primary solution for a given load/time step to the VTF-file.
  //! \param[in] psol Primary solution vector
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] time Load/time step parameter
  //! \param[in] pvecName Optional name of the primary vector field solution
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] psolComps Optional number of primary solution components
  //! \param[in] scalarOnly If \e true, write vector as scalar components only
  virtual int writeGlvS1(const Vector& psol, int iStep, int& nBlock,
                         double time = 0.0, const char* pvecName = nullptr,
                         int idBlock = 10, int psolComps = 0,
                         bool scalarOnly = false);

  //! \brief Writes secondary solution for a load/time step to the VTF-file.
  //! \param[in] psol Primary solution vector
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] time Load/time step parameter
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] psolComps Optional number of primary solution components
  virtual bool writeGlvS2(const Vector& psol, int iStep, int& nBlock,
                          double time = 0.0, int idBlock = 20,
                          int psolComps = 0);

  //! \brief Evaluates the secondary solution for a given load/time step.
  //! \param[in] psol Primary solution vector
  //! \param[in] time Load/time step parameter
  //! \param[in] psolComps Optional number of primary solution components
  //!
  //! \details This method only evaluates the solutions fields, and does not
  //! return any data. It is used only for load/time steps that are not saved
  //! when the solution has to be evaluated at every increment in any case to
  //! ensure consistency (e.g., when material models with history variables
  //! are in use).
  bool eval2ndSolution(const Vector& psol, double time, int psolComps = 0);

  //! \brief Writes projected solutions for a given time step to the VTF-file.
  //! \param[in] ssol Secondary solution vector (control point values)
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] prefix Common prefix for the field components
  //! \param[in] maxVal Optional array of maximum values
  bool writeGlvP(const RealArray& ssol, int iStep, int& nBlock,
                 int idBlock = 100, const char* prefix = "Global projected",
                 std::vector<PointValues>* maxVal = nullptr);

  //! \brief Writes a mode shape to the VTF-file.
  //! \param[in] mode The mode shape eigenvector and associated eigenvalue
  //! \param[in] freq \e true if the eigenvalue is a frequency
  //! \param nBlock Running result block counter
  //!
  //! \details The eigenvalue is used as a label on the step state info.
  bool writeGlvM(const Mode& mode, bool freq, int& nBlock);

  //! \brief Writes element field for a given load/time step to the VTF-file.
  //! \param[in] field The element field to output
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] name Name of field
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] internalOrder If \e true, the data in \a field are assumed to
  //! be ordered w.r.t. the internal element ordering
  bool writeGlvE(const Vector& field, int iStep, int& nBlock,
                 const char* name, int idBlock = 300,
                 bool internalOrder = false) const;

  //! \brief Writes element norms for a given load/time step to the VTF-file.
  //! \param[in] norms The element norms to output
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] prefix Prefices for projected solutions
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] dualPrefix Prefix for dual solution norms
  bool writeGlvN(const Matrix& norms, int iStep, int& nBlock,
                 const std::vector<std::string>& prefix = {},
                 int idBlock = 200, const char* dualPrefix = nullptr);

  //! \brief Writes a scalar function to the VTF-file.
  //! \param[in] f The function to output
  //! \param[in] fname Name of the function
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  //! \param[in] idBlock Starting value of result block numbering
  //! \param[in] time Load/time step parameter
  bool writeGlvF(const RealFunc& f, const char* fname,
                 int iStep, int& nBlock, int idBlock = 50, double time = 0.0);

  //! \brief Writes load/time step info to the VTF-file.
  //! \param[in] iStep Load/time step identifier
  //! \param[in] value Load parameter or time of the step
  //! \param[in] itype Type identifier of the step
  bool writeGlvStep(int iStep, double value = 0.0, int itype = 0);

  //! \brief Closes the current VTF-file.
  void closeGlv();

  //! \brief Returns the current VTF-file object.
  VTF* getVTF() const { return myVtf; }
  //! \brief Defines the VTF-file for subsequent results output.
  void setVTF(VTF* vtf) { myVtf = vtf; }

  //! \brief Returns the initial geometry block index.
  int getStartGeo() const { return myGeomID; }
  //! \brief Initializes the geometry block counter.
  void setStartGeo(int gID) { myGeomID = gID; }

  //! \brief Dumps the FE model to Matlab format.
  //! \param os Output stream to write the grid data to
  //! \param[in] name Name of the Matlab function returning the grid
  //! \param[in] sets Names of nodal sets to be printed out in addition
  //! \param[in] scale Scaling factor
  bool dumpMatlabGrid(std::ostream& os, const std::string& name = "IFEM_Mesh",
                      const std::vector<std::string>& sets = {},
                      double scale = 1.0) const;

  //! \brief Dumps the (possibly refined) spline geometry in g2-format.
  //! \param os Output stream to write the geometry data to
  bool dumpGeometry(std::ostream& os) const;

  //! \brief Dumps the primary solution in ASCII format for inspection.
  //! \param[in] psol Primary solution vector
  //! \param os Output stream to write the solution data to
  //! \param[in] withID If \e true, write node ID and coordinates too
  void dumpPrimSol(const Vector& psol, utl::LogStream& os,
                   bool withID = true) const;
  //! \brief Dumps the entire solution in ASCII format.
  //! \param[in] psol Primary solution vector to derive other quantities from
  //! \param os Output stream to write the solution data to
  bool dumpSolution(const Vector& psol, utl::LogStream& os) const;
  //! \brief Dumps solution results at specified points in ASCII format.
  //! \param[in] psol Primary solution vector to derive other quantities from
  //! \param[in] time Load/time step parameter
  //! \param os Output stream to write the solution data to
  //! \param[in] formatted If \e false, write all result points on a single line
  //!            without point identifications, but with time as first column
  //! \param[in] precision Number of digits after the decimal point
  bool dumpResults(const Vector& psol, double time, utl::LogStream& os,
                   bool formatted = false, std::streamsize precision = 3) const;
  //! \brief Dumps vector solution at specified points in ASCII format.
  //! \param[in] vsol Solution vector
  //! \param[in] fname Name of vector field
  //! \param os Output stream to write the solution data to
  //! \param[in] precision Number of digits after the decimal point
  bool dumpVector(const Vector& vsol, const char* fname,
                  utl::LogStream& os, std::streamsize precision = 3) const;
  //! \brief Dumps additional problem-specific results in ASCII format.
  virtual void dumpMoreResults(double, utl::LogStream&,
                               std::streamsize = 3) const {}

  //! \brief Saves point results to output file for a given time step.
  //! \param[in] psol Primary solution vector
  //! \param[in] time Load/time step parameter
  //! \param[in] step Load/time step counter
  bool savePoints(const Vector& psol, double time, int step) const;

  //! \brief Saves result components to output files for a given time step.
  //! \param[in] psol Primary solution vectors
  //! \param[in] time Load/time step parameter
  //! \param[in] step Load/time step counter
  bool saveResults(const Vectors& psol, double time, int step) const;

  //! \brief Checks whether result points have been defined or not.
  bool hasResultPoints() const { return !myPoints.empty(); }
  //! \brief Checks whether point result files have been defined or not.
  bool hasPointResultFile() const;
  //! \brief Sets the file name for result point output.
  //! \param[in] filename The file name prefix (optionally with extension)
  //! \param[in] dumpCoord If \e true, write point coordinates to separate file
  void setPointResultFile(const std::string& filename, bool dumpCoord = false);

  //! \brief Serialization support.
  virtual bool serialize(std::map<std::string,std::string>&) const;

  //! \brief Returns the reference norm to base mesh adaptation upon.
  virtual double getReferenceNorm(const Vectors&, size_t) const = 0;
  //! \brief Returns the global effectivity index.
  virtual double getEffectivityIndex(const Vectors&, size_t, size_t) const = 0;
  //! \brief Prints integrated solution norms to the log stream.
  virtual void printNorms(const Vectors&, size_t = 36) const = 0;
  //! \brief Prints out interface force resultants to the log stream.
  virtual void printIFforces(const Vector&, RealArray&) {}
  //! \brief Prints out the nodal reaction forces to the log stream.
  virtual int printNRforces(const std::vector<int>& glbNodes = {}) const;

  //! \brief Writes out the additional functions to VTF-file.
  virtual bool writeAddFuncs(int iStep, int& nBlock, int idBlock, double time);

protected:
  //! \brief Adds an additional function for VTF-file output.
  void addAddFunc(const std::string& name, RealFunc* f);

private:
  //! \brief Private helper to initialize patch for solution evaluation.
  bool initPatchForEvaluation(int patchNo) const;

  //! \brief Private helper to extract patch-level solution vectors.
  bool extractNodeVec(const RealArray& glbVec, Vector& locVec,
                      const ASMbase* patch, int nodalCmps,
                      bool& emptyPatches) const;

  //! \brief Private helper to write out scalar fields to VTF-file.
  bool writeScalarFields(const Matrix& field, int geomID,
                         int& nBlock, std::vector< std::vector<int> >& sID,
                         size_t* nScl = nullptr,
                         ASM::ResultClass resClass = ASM::PRIMARY);

protected:
  //! \brief Struct defining a result sampling point.
  struct ResultPoint
  {
    short int    npar;  //!< Number of parameters
    unsigned int patch; //!< Patch index [1,nPatch]
    int          inod;  //!< Local node number of the closest node
    double       u[3];  //!< Parameters of the point (u,v,w)
    Vec3         X;     //!< Spatial coordinates of the point

    //! \brief Default constructor.
    ResultPoint() : npar(3), patch(1), inod(0) { u[0] = u[1] = u[2] = 0.0; }
  };

  //! \brief Result point container.
  using ResPointVec = std::vector<ResultPoint>;
  //! \brief File name to result point group mapping.
  using ResPtPair = std::pair<std::string,ResPointVec>;

  std::vector<ResPtPair> myPoints; //!< User-defined result sampling points

  //! \brief Preprocesses the result sampling points.
  virtual void preprocessResultPoints();

private:
  //! \brief Preprocesses a result sampling point group.
  //! \param ptFile Name of file that these result points are dumped to
  //! \param points Group of result points that are dumped to the given file
  void preprocessResPtGroup(std::string& ptFile, ResPointVec& points);

  //! \brief Dumps solution results at the given points in ASCII format.
  //! \param[in] psol Primary solution vector to derive other quantities from
  //! \param[in] time Load/time step parameter
  //! \param os Output stream to write the solution data to
  //! \param[in] gPoints Group of result points to write solution data for
  //! \param[in] formatted If \e false, write all result points on a single line
  //!            without point identifications, but with time as first column
  //! \param[in] precision Number of digits after the decimal point
  bool dumpResults(const Vector& psol, double time,
                   utl::LogStream& os, const ResPointVec& gPoints,
                   bool formatted, std::streamsize precision) const;

  //! \brief Evaluates solution results at specified points for a given patch.
  //! \param[in] psol Primary solution vectors to derive other quantities from
  //! \param[in] gPoints Result point definitions
  //! \param[in] patch The patch to evaluate result points for
  //! \param[out] points List of result points within this patch
  //! \param[out] Xp Coordinates of result points within this patch
  //! \param[out] sol1 Matrix of primary solution values at result points
  //! \param[out] sol2 Matrix of secondary solution values at result points
  //! \param[out] compNames Names of solution components in \a sol1 and \a sol2
  bool evalResults(const Vectors& psol, const ResPointVec& gPoints,
                   const ASMbase* patch, std::vector<int>& points, Vec3Vec& Xp,
                   Matrix& sol1, Matrix& sol2,
                   std::vector<std::string>* compNames = nullptr) const;

  std::map<std::string,RealFunc*> myAddScalars; //!< Scalar functions to output

  int    myPrec;   //!< Output precision for result sampling
  double myPtSize; //!< Size of result point visualization in VTF-file
  int    myGeomID; //!< Geometry block ID for the first patch in the VTF-file
  int    myGeofs1; //!< Geometry block ID offset for immersed geometry in the VTF-file
  int    myGeofs2; //!< Geometry block ID offset for extra geometry in the VTF-file
  VTF*   myVtf;    //!< VTF-file for result visualization
  bool   logRpMap; //!< If \e true, print out the result point mapping
  int    idxGrid;  //!< Index into \ref myPoints for grid result output
};

#endif
