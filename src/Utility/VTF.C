// $Id$
//==============================================================================
//!
//! \file VTF.C
//!
//! \date Dec 1 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Output of FE model and results to VTF file.
//!
//==============================================================================

#include "VTF.h"
#if HAS_VTFAPI == 1
#include "VTFAPI.h"
#include "VTOAPIPropertyIDs.h"
#elif HAS_VTFAPI == 2
#include "VTFXAPI.h"
#include "VTOAPIPropertyIDs.h"
#define VTFA_FAILURE VTFXA_FAILURE //!< Convenience macro
#define VTFA_SUCCESS VTFXA_SUCCESS //!< Convenience macro
#else
#define VTFA_FAILURE(x) x <= 0 //!< Convenience macro
#define VTFA_SUCCESS(x) x > 0  //!< Convenience macro
#endif
#include "ElementBlock.h"
#include "Function.h"
#include "Tensor.h"
#include <iostream>
#include <cstdio>


Real VTF::vecOffset[3] = { 0.0, 0.0, 0.0 };


/*!
  \brief Static helper printing an error message to \a std::cerr.
  \param[in] msg The message to print
  \param[in] ID1 If non-negative, the value is appended to the message
  \param[in] ID2 If non-negative, the value is appended to the message
  \return \e false (always)
*/

static bool showError (const char* msg, int ID1 = -1, int ID2 = -2)
{
  std::cerr <<"VTF: "<< msg;
  if (ID1 >= 0) std::cerr <<" "<< ID1;
  if (ID2 >= 0) std::cerr <<" ("<< ID2 <<")";
  std::cerr << std::endl;
  return false;
}


VTF::VTF (const char* filename, int type)
{
  myFile = nullptr;
  myState = nullptr;
  myGBlock = nullptr;
  myPartID = pointGeoID = 0;
  lastStep = 0;
  if (!filename) return;

#if HAS_VTFAPI == 1
  // Create the VTF file object
  myFile = new VTFAFile();
  // Enable debug info to stderr/console
  myFile->SetOutputDebugError(1);

  if (!VTFA_FAILURE(myFile->CreateVTFFile(filename,type > 0)))
    return;

  delete myFile;
  showError("Error creating VTF file");
#elif HAS_VTFAPI == 2
  myFile = new VTFXAFile();
  VTFXAFileSettings settings;
  VTFXAInitFileSettings(&settings);
  settings.bBinary = type > 0 ? VTFXA_TRUE : VTFXA_FALSE;
  settings.pszApplicationName = "IFEM";
  settings.pszVendorName = "SINTEF ICT";
  settings.iVendorID = 1001;
  settings.pszVendorCode = "Test";
  if (!VTFA_FAILURE(myFile->CreateVTFxFile(filename,&settings)))
  {
    myDatabase = new VTFXADatabase(myFile,"Single",1);
    return;
  }
  delete myFile;
  showError("Error creating VTFx file");
#else
  showError("Not available in this version");
#endif
  myFile = nullptr;
}


VTF::~VTF ()
{
  if (!myFile) return;

  for (GridBlock& block : myBlocks)
    delete block.second;

#if HAS_VTFAPI == 1
  if (myGBlock)
  {
    if (VTFA_FAILURE(myFile->WriteBlock(myGBlock)))
      showError("Error writing Geometry Block");
    delete myGBlock;
  }

  for (VTFATransformationBlock* block : myTBlock)
    if (block)
    {
      if (VTFA_FAILURE(myFile->WriteBlock(block)))
        showError("Error writing Transformation Block");
      delete block;
    }

  for (VTFADisplacementBlock* block : myDBlock)
    if (block)
    {
      if (VTFA_FAILURE(myFile->WriteBlock(block)))
        showError("Error writing Displacement Block");
      delete block;
    }

  for (VTFAVectorBlock* block : myVBlock)
    if (block)
    {
      if (VTFA_FAILURE(myFile->WriteBlock(block)))
        showError("Error writing Vector Block");
      delete block;
    }

  for (VTFAScalarBlock* block : mySBlock)
    if (block)
    {
      if (VTFA_FAILURE(myFile->WriteBlock(block)))
        showError("Error writing Scalar Block");
      delete block;
    }

  if (myState)
  {
    if (VTFA_FAILURE(myFile->WriteBlock(myState)))
      showError("Error writing state info block");
    delete myState;
  }

  if (VTFA_FAILURE(myFile->CloseFile()))
    showError("Error closing VTF file");

  delete myFile;
#elif HAS_VTFAPI == 2
  if (myGBlock)
  {
    if (VTFA_FAILURE(myDatabase->WriteBlock(myGBlock)))
      showError("Error writing Geometry Block");
    delete myGBlock;
  }

  for (VTFATransformationBlock* block : myTBlock)
    if (block)
    {
      if (VTFA_FAILURE(myDatabase->WriteBlock(block)))
        showError("Error writing Transformation Block");
      delete block;
    }

  for (VTFADisplacementBlock* block : myDBlock)
    if (block)
    {
      if (VTFA_FAILURE(myDatabase->WriteBlock(block)))
        showError("Error writing Displacement Block");
      delete myDBlock;
    }

  for (VTFAVectorBlock* block : myVBlock)
    if (block)
    {
      if (VTFA_FAILURE(myDatabase->WriteBlock(block)))
        showError("Error writing Vector Block");
      delete block;
    }

  for (VTFAScalarBlock* block : mySBlock)
    if (block)
    {
      if (VTFA_FAILURE(myDatabase->WriteBlock(block)))
        showError("Error writing Scalar Block");
      delete block;
    }

  if (myState)
  {
    if (VTFA_FAILURE(myDatabase->WriteBlock(myState)))
      showError("Error writing state info block");
    delete myState;
  }

  VTFXACase                singleCase(myFile,"Single case",1,1);
  VTFXACasePropertiesBlock frameGeneratorProps(VT_CT_FRAME_GENERATOR_SETTINGS);

  frameGeneratorProps.AddInt(VT_PI_FG_FEM_MODEL_IDS,1);
  singleCase.WritePropertiesBlock(&frameGeneratorProps);
  for (const GridBlock& block : myBlocks)
  {
    VTFXACasePropertiesBlock partAttr(VT_CT_PART_ATTRIBUTES);
    partAttr.SetPartID(block.first.first);
    partAttr.AddBool(VT_PB_PA_MESH, VTFXA_FALSE);
    partAttr.AddBool(VT_PB_PA_DISPLACEMENTS, VTFXA_FALSE);
    singleCase.WritePropertiesBlock(&partAttr);
  }
  if (VTFA_FAILURE(myFile->CloseFile()))
    showError("Error closing VTF file");

  delete myDatabase;
  delete myFile;
#endif
}


void VTF::writeGeometryBlocks (int iStep)
{
#ifdef HAS_VTFAPI
  if (myBlocks.empty())
    return;

  std::vector<int> geomID;
  geomID.reserve(myBlocks.size());
  for (const GridBlock& block : myBlocks)
    geomID.push_back(block.first.first);

#if HAS_VTFAPI == 1
  if (!myGBlock) myGBlock = new VTFAGeometryBlock();
  myGBlock->SetGeometryElementBlocks(geomID.data(),geomID.size(),iStep);
#elif HAS_VTFAPI == 2
  if (!myGBlock) myGBlock = new VTFXAGeometryBlock();
  myGBlock->SetElementBlocksForState(geomID.data(),geomID.size(),iStep);
#endif
#endif
}


void VTF::clearGeometryBlocks ()
{
  for (GridBlock& block : myBlocks)
    delete block.second;

  myBlocks.clear();
  myPartID = pointGeoID = 0;
}


const ElementBlock* VTF::getBlock (int geomID) const
{
  if (geomID > 0 && geomID <= static_cast<int>(myBlocks.size()))
    return myBlocks[geomID-1].second;

  return nullptr;
}


int VTF::getNodeBlock (int geomID) const
{
  if (geomID > 0 && geomID <= static_cast<int>(myBlocks.size()))
    return myBlocks[geomID-1].first.second;

  return 0;
}


bool VTF::writeGrid (const ElementBlock* block, const char* gName,
                     int elemID, int nodeID)
{
  if (!myFile) return true;
  if (!block) return false;

  myBlocks.emplace_back(IntPair(elemID, nodeID > 0 ? nodeID : elemID), block);

  if (nodeID <= 0 && !this->writeNodes(nodeID = elemID))
    return showError("Error writing node block",nodeID);

  if (!this->writeElements(gName,++myPartID,elemID,nodeID))
    return showError("Error writing element block",elemID);

  return true;
}


bool VTF::writeTransformation (const Vec3& X, const Tensor& T,
                               int idBlock, int geomID)
{
  if (!myFile) return true;
  if (!this->getBlock(geomID)) return false;

#if HAS_VTFAPI == 1
  // Cast to float
  float trans[12];
  size_t i, j, k = 0;
  for (j = 1; j <= 3; j++)
    for (i = 1; i <= 3; i++)
      trans[k++] = T(i,j);
  for (j = 0; j < 3; j++)
    trans[k++] = X[j];

  VTFAMatrixResultBlock tBlock(idBlock);
  if (VTFA_FAILURE(tBlock.SetMatrix(trans)))
    return showError("Error defining result block",idBlock);

  tBlock.SetMapToElementBlockID(myBlocks[geomID-1].first.first);
  if (VTFA_FAILURE(myFile->WriteBlock(&tBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  showError("Transformation not yet implemented for VTFx");
#endif

  return true;
}


bool VTF::writeVres (const std::vector<Real>& nodeResult,
                     int idBlock, int geomID, size_t nvc)
{
  if (!myFile) return true;

  const ElementBlock* grid = this->getBlock(geomID);
  if (!grid) return false;

  const size_t nnod = grid->getNoNodes();
  const size_t nres = nodeResult.size();
  const size_t ncmp = nres/(nnod > 0 ? nnod : 1);
  if (nres != ncmp*nnod)
    return showError("Invalid size of result array",nres,ncmp*nnod);
  else if (nvc < 1 || nvc > ncmp)
    nvc = ncmp;

#ifdef HAS_VTFAPI
  // Cast to float
  std::vector<float> resVec(3*nnod,0.0f);
  if (nres == 3*nnod && nvc == 3)
    for (size_t i = 0; i < nres; i++)
      resVec[i] = nodeResult[i];
  else if (nres == nnod && nvc == 1)
    for (size_t i = 0; i < nnod; i++)
      // Writing a scalar solution as Z-deflection
      resVec[3*i+2] = nodeResult[i];
  else
    for (size_t i = 0; i < nnod; i++)
      for (size_t j = 0; j < 3; j++)
        if (j < nvc) resVec[3*i+j] = nodeResult[ncmp*i+j];

#if HAS_VTFAPI == 1
  VTFAResultBlock dBlock(idBlock,VTFA_DIM_VECTOR,VTFA_RESMAP_NODE,0);

  if (VTFA_FAILURE(dBlock.SetResults3D(resVec.data(),nnod)))
    return showError("Error defining result block",idBlock);

  dBlock.SetMapToBlockID(myBlocks[geomID-1].first.second);
  if (VTFA_FAILURE(myFile->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  VTFXAResultValuesBlock dBlock(idBlock,VTFXA_DIM_VECTOR,VTFXA_FALSE);
  dBlock.SetMapToBlockID(myBlocks[geomID-1].first.second,VTFXA_NODES);
  dBlock.SetResultValues3D(resVec.data(),nnod);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#endif
#endif

  return true;
}


/*!
  If \a elementResult is empty, the external element numbers will instead
  define the scalar field to be written.
*/

bool VTF::writeEres (const std::vector<Real>& elementResult,
                     int idBlock, int geomID)
{
  if (!myFile) return true;

  const ElementBlock* grid = this->getBlock(geomID);
  if (!grid) return false;

  size_t nels = grid->getNoElms();
  size_t nres = elementResult.size();
  if (nres > nels)
    return showError("Invalid size of result array",nres,nels);
  else if (nres > 0 && nres < nels)
    showError("Warning: Fewer element results that anticipated",nres,nels);
  else // write the external element numbers
    nres = nels;

#ifdef HAS_VTFAPI
  // Cast to float
  std::vector<float> resVec(nres);
  for (size_t i = 0; i < nres; i++)
  {
    size_t j = grid->getElmIndex(i);
    resVec[j] = elementResult.empty() ? grid->getElmId(1+i) : elementResult[i];
  }

#if HAS_VTFAPI == 1
  VTFAResultBlock dBlock(idBlock,VTFA_DIM_SCALAR,VTFA_RESMAP_ELEMENT,0);

  if (VTFA_FAILURE(dBlock.SetResults1D(resVec.data(),nres)))
    return showError("Error defining result block",idBlock);

  dBlock.SetMapToBlockID(myBlocks[geomID-1].first.first);
  if (VTFA_FAILURE(myFile->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  VTFXAResultValuesBlock dBlock(idBlock,VTFXA_DIM_SCALAR,VTFXA_FALSE);
  dBlock.SetMapToBlockID(myBlocks[geomID-1].first.first,VTFXA_ELEMENTS);
  dBlock.SetResultValues1D(resVec.data(),nres);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#endif
#endif

  return true;
}


bool VTF::writeNres (const std::vector<Real>& nodalResult,
                     int idBlock, int geomID)
{
  if (!myFile) return true;

  const ElementBlock* grid = this->getBlock(geomID);
  if (!grid) return false;

  const size_t nres = nodalResult.size();
  if (nres != grid->getNoNodes())
    return showError("Invalid size of result array",nres,grid->getNoNodes());

#ifdef HAS_VTFAPI
  // Cast to float
  std::vector<float> resVec(nres);
  for (size_t i = 0; i < nres; i++)
    resVec[i] = nodalResult[i];

#if HAS_VTFAPI == 1
  VTFAResultBlock dBlock(idBlock,VTFA_DIM_SCALAR,VTFA_RESMAP_NODE,0);

  if (VTFA_FAILURE(dBlock.SetResults1D(resVec.data(),nres)))
    return showError("Error defining result block",idBlock);

  dBlock.SetMapToBlockID(myBlocks[geomID-1].first.second);
  if (VTFA_FAILURE(myFile->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  VTFXAResultValuesBlock dBlock(idBlock,VTFXA_DIM_SCALAR,VTFXA_FALSE);
  dBlock.SetMapToBlockID(myBlocks[geomID-1].first.second,VTFXA_NODES);
  dBlock.SetResultValues1D(resVec.data(),nres);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#endif
#endif

  return true;
}


bool VTF::writeNfunc (const RealFunc& f, Real time, int idBlock, int geomID)
{
  if (!myFile) return true;

  const ElementBlock* grid = this->getBlock(geomID);
  if (!grid) return false;

#ifdef HAS_VTFAPI
  const size_t nres = grid->getNoNodes();

  // Evaluate the function at the grid points
  std::vector<float> resVec(nres);
  for (size_t i = 0; i < nres; i++)
    resVec[i] = f(Vec4(grid->getCoord(i),time,grid->getParam(i)));

#if HAS_VTFAPI == 1
  VTFAResultBlock dBlock(idBlock,VTFA_DIM_SCALAR,VTFA_RESMAP_NODE,0);

  if (VTFA_FAILURE(dBlock.SetResults1D(resVec.data(),nres)))
    return showError("Error defining result block",idBlock);

  dBlock.SetMapToBlockID(myBlocks[geomID-1].first.second);
  if (VTFA_FAILURE(myFile->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  VTFXAResultValuesBlock dBlock(idBlock,VTFXA_DIM_SCALAR,VTFXA_FALSE);
  dBlock.SetMapToBlockID(myBlocks[geomID-1].first.second,VTFXA_NODES);
  dBlock.SetResultValues1D(resVec.data(),nres);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&dBlock)))
    return showError("Error writing result block",idBlock);
#endif
#endif

  return true;
}


bool VTF::writeVectors (const std::vector<Vec3Pair>& pntResult, int& gID,
                        int idBlock, const char* resultName,
                        int iStep, int iBlock)
{
#if HAS_VTFAPI == 1
  bool writePoints = false;
  if (pointGeoID == 0)
  {
    // The big assumption here is that we have only one call to writeVectors
    // per time step, and that all subsequent calls are with the same points
    pointGeoID = gID < 0 && !myBlocks.empty() ?
      myBlocks.back().first.second+1 : ++gID;
    myBlocks.emplace_back(IntPair(pointGeoID,pointGeoID),nullptr);
    writePoints = true;
  }

  VTFANodeBlock   nBlock(pointGeoID,0);
  VTFAResultBlock rBlock(idBlock,VTFA_DIM_VECTOR,VTFA_RESMAP_NODE,0);

  // Count the nodes of the point mesh to associate the vectors with.
  // It is assumed that identically zero points are associated with
  // zero knot-span elements (due to multiple knots) which will not be
  // assigned any results. Therefore we should ignore those points.
  size_t i = 0, np = 0;
  for (const Vec3Pair& pt : pntResult)
    if (!pt.first.isZero(1.0e-15)) np++;

  if (writePoints && VTFA_FAILURE(nBlock.SetNumNodes(np)))
    return showError("Error defining node block",pointGeoID);
  else if (VTFA_FAILURE(rBlock.SetNumResults(np)))
    return showError("Error defining result block",idBlock);
  else
    rBlock.SetMapToBlockID(pointGeoID);

  std::vector<int> mnpc;
  if (writePoints) mnpc.reserve(np);
  for (const Vec3Pair& pt : pntResult)
    if (!pt.first.isZero(1.0e-15)) // skip the zero-knotspan points
    {
      if (writePoints && VTFA_FAILURE(nBlock.AddNode(vecOffset[0]+pt.first.x,
                                                     vecOffset[1]+pt.first.y,
                                                     vecOffset[2]+pt.first.z)))
        return showError("Error adding node to block",pointGeoID);
      if (VTFA_FAILURE(rBlock.AddResult(pt.second.x,pt.second.y,pt.second.z)))
        return showError("Error adding result to block",idBlock);
      else if (writePoints)
        mnpc.push_back(i++);
    }

  if (writePoints)
  {
    // We must define an element block (with point elements) also,
    // otherwise GLview does not visualize the vectors
    VTFAElementBlock eBlock(pointGeoID,0,0);
    eBlock.SetPartID(++myPartID);
    eBlock.SetNodeBlockID(pointGeoID);
    if (VTFA_FAILURE(eBlock.AddElements(VTFA_POINTS,mnpc.data(),np)))
      return showError("Error defining element block",pointGeoID);

    if (VTFA_FAILURE(myFile->WriteBlock(&nBlock)))
      return showError("Error writing node block",pointGeoID);
    else if (VTFA_FAILURE(myFile->WriteBlock(&eBlock)))
      return showError("Error writing element block",pointGeoID);
  }
  if (VTFA_FAILURE(myFile->WriteBlock(&rBlock)))
    return showError("Error writing result block",idBlock);
#elif HAS_VTFAPI == 2
  showError("Vector points are not yet implemented for VTFx");
#endif

  return iStep > 0 ? this->writeVblk(idBlock,resultName,iBlock,iStep) : true;
}


bool VTF::writePoints (const Vec3Vec& points, int& gID)
{
#if HAS_VTFAPI == 1
  VTFANodeBlock nBlock(++gID,0);
  size_t i = 0, np = points.size();
  if (VTFA_FAILURE(nBlock.SetNumNodes(np)))
    return showError("Error defining node block",gID);

  std::vector<int> mnpc; mnpc.reserve(np);
  for (const Vec3& pnt : points)
    if (VTFA_FAILURE(nBlock.AddNode(vecOffset[0]+pnt.x,
                                    vecOffset[1]+pnt.y,
                                    vecOffset[2]+pnt.z)))
      return showError("Error adding node to block",gID);
    else
      mnpc.push_back(i++);

  // We must define an element block (with point elements) also,
  // otherwise GLview does not visualize the points
  VTFAElementBlock eBlock(gID,0,0);
  eBlock.SetPartID(++myPartID);
  eBlock.SetNodeBlockID(gID);
  if (VTFA_FAILURE(eBlock.AddElements(VTFA_POINTS,mnpc.data(),np)))
    return showError("Error defining element block",gID);

  if (VTFA_FAILURE(myFile->WriteBlock(&nBlock)))
    return showError("Error writing node block",gID);
  else if (VTFA_FAILURE(myFile->WriteBlock(&eBlock)))
    return showError("Error writing element block",gID);

  myBlocks.emplace_back(IntPair(gID,gID),nullptr);
#elif HAS_VTFAPI == 2
  showError("Points are not yet implemented for VTFx");
#endif

  return true;
}


bool VTF::writeTblk (const std::vector<int>& tBlockIDs, const char* resultName,
                     int idBlock, int iStep)
{
  if ((int)myTBlock.size() < idBlock) myTBlock.resize(idBlock,0);

  int status = 1;
#if HAS_VTFAPI == 1
  if (!myTBlock[--idBlock])
  {
    myTBlock[idBlock] = new VTFATransformationBlock(idBlock+1);
    if (resultName) myTBlock[idBlock]->SetName(resultName);
    status = myTBlock[idBlock]->SetResultBlocks(tBlockIDs.data(),
                                                tBlockIDs.size(),iStep);
  }
  else for (size_t i = 0; i < tBlockIDs.size() && VTFA_SUCCESS(status); i++)
    status = myTBlock[idBlock]->AddResultBlock(tBlockIDs[i],iStep);
#elif HAS_VTFAPI == 2
  showError("Transformation not yet implemented for VTFx");
#endif
  if (VTFA_FAILURE(status))
    return showError("Error defining transformation block",idBlock);

  return true;
}


bool VTF::writeDblk (const std::vector<int>& dBlockIDs, const char* resultName,
                     int idBlock, int iStep)
{
  if ((int)myDBlock.size() < idBlock) myDBlock.resize(idBlock,0);

  int status = 1;
#if HAS_VTFAPI == 1
  if (!myDBlock[--idBlock])
  {
    myDBlock[idBlock] = new VTFADisplacementBlock(idBlock+1);
    if (resultName) myDBlock[idBlock]->SetName(resultName);
    myDBlock[idBlock]->SetRelativeDisplacementResults(1);
    status = myDBlock[idBlock]->SetResultBlocks(dBlockIDs.data(),
                                                dBlockIDs.size(),iStep);
  }
  else for (size_t i = 0; i < dBlockIDs.size() && VTFA_SUCCESS(status); i++)
    status = myDBlock[idBlock]->AddResultBlock(dBlockIDs[i],iStep);
#elif HAS_VTFAPI == 2
  if (!myDBlock[--idBlock])
  {
    myDBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
                                             VTFXA_RESTYPE_DISPLACEMENT,
                                             VTFXA_RESMAP_NODE);
    if (resultName) myDBlock[idBlock]->SetName(resultName);
  }
  myDBlock[idBlock]->SetResultID(idBlock);
  status = myDBlock[idBlock]->SetResultValuesBlocks(dBlockIDs.data(),
                                                    dBlockIDs.size(),iStep);
#endif
  if (VTFA_FAILURE(status))
    return showError("Error defining displacement block",idBlock);

  return true;
}


bool VTF::writeVblk (int vBlockID, const char* resultName,
                     int idBlock, int iStep)
{
  if ((int)myVBlock.size() < idBlock) myVBlock.resize(idBlock,0);

  int status = 1;
#if HAS_VTFAPI == 1
  if (!myVBlock[--idBlock])
  {
    myVBlock[idBlock] = new VTFAVectorBlock(idBlock+1);
    if (resultName) myVBlock[idBlock]->SetName(resultName);
    status = myVBlock[idBlock]->SetResultBlocks(&vBlockID,1,iStep);
  }
  else
    status = myVBlock[idBlock]->AddResultBlock(vBlockID,iStep);
#elif HAS_VTFAPI == 2
  if (!myVBlock[--idBlock])
  {
    myVBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
                                             VTFXA_RESTYPE_VECTOR,
                                             VTFXA_RESMAP_NODE);
    if (resultName) myVBlock[idBlock]->SetName(resultName);
  }
  myVBlock[idBlock]->SetResultID(idBlock);
  status = myVBlock[idBlock]->SetResultValuesBlocks(&vBlockID,1,iStep);
#endif
  if (VTFA_FAILURE(status))
    return showError("Error defining vector block",idBlock);

  return true;
}


bool VTF::writeVblk (const std::vector<int>& vBlockIDs, const char* resultName,
                     int idBlock, int iStep)
{
  if ((int)myVBlock.size() < idBlock) myVBlock.resize(idBlock,0);

  int status = 1;
#if HAS_VTFAPI == 1
  if (!myVBlock[--idBlock])
  {
    myVBlock[idBlock] = new VTFAVectorBlock(idBlock+1);
    if (resultName) myVBlock[idBlock]->SetName(resultName);
    status = myVBlock[idBlock]->SetResultBlocks(vBlockIDs.data(),
                                                vBlockIDs.size(),iStep);
  }
  else for (size_t i = 0; i < vBlockIDs.size() && VTFA_SUCCESS(status); i++)
    status = myVBlock[idBlock]->AddResultBlock(vBlockIDs[i],iStep);
#elif HAS_VTFAPI == 2
  if (!myVBlock[--idBlock])
  {
    myVBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
                                             VTFXA_RESTYPE_VECTOR,
                                             VTFXA_RESMAP_NODE);
    if (resultName) myVBlock[idBlock]->SetName(resultName);
  }
  myVBlock[idBlock]->SetResultID(idBlock);
  status = myVBlock[idBlock]->SetResultValuesBlocks(vBlockIDs.data(),
                                                    vBlockIDs.size(),iStep);
#endif
  if (VTFA_FAILURE(status))
    return showError("Error defining vector block",idBlock);

  return true;
}


bool VTF::writeSblk (int sBlockID, const char* resultName,
                     int idBlock, int iStep, bool elementData)
{
  if ((int)mySBlock.size() < idBlock) mySBlock.resize(idBlock,0);

  int status = 1;
#if HAS_VTFAPI == 1
  if (!mySBlock[--idBlock])
  {
    mySBlock[idBlock] = new VTFAScalarBlock(idBlock+1);
    if (resultName) mySBlock[idBlock]->SetName(resultName);
    status = mySBlock[idBlock]->SetResultBlocks(&sBlockID,1,iStep);
  }
  else
    status = mySBlock[idBlock]->AddResultBlock(sBlockID,iStep);
#elif HAS_VTFAPI == 2
  if (!mySBlock[--idBlock])
  {
    int type = elementData?VTFXA_RESMAP_ELEMENT:VTFXA_RESMAP_NODE;
    mySBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
                                             VTFXA_RESTYPE_SCALAR,type);
    if (resultName) mySBlock[idBlock]->SetName(resultName);
  }
  mySBlock[idBlock]->SetResultID(idBlock);
  status = mySBlock[idBlock]->SetResultValuesBlocks(&sBlockID,1,iStep);
#endif
  if (VTFA_FAILURE(status))
    return showError("Error defining scalar block",idBlock);

  return true;
}


bool VTF::writeSblk (const std::vector<int>& sBlockIDs, const char* resultName,
                     int idBlock, int iStep, bool elementData)
{
  if ((int)mySBlock.size() < idBlock) mySBlock.resize(idBlock,0);

  int status = 1;
#if HAS_VTFAPI == 1
  if (!mySBlock[--idBlock])
  {
    mySBlock[idBlock] = new VTFAScalarBlock(idBlock+1);
    if (resultName) mySBlock[idBlock]->SetName(resultName);
    status = mySBlock[idBlock]->SetResultBlocks(sBlockIDs.data(),
                                                sBlockIDs.size(),iStep);
  }
  else for (size_t i = 0; i < sBlockIDs.size() && VTFA_SUCCESS(status); i++)
    status = mySBlock[idBlock]->AddResultBlock(sBlockIDs[i],iStep);
#elif HAS_VTFAPI == 2
  if (!mySBlock[--idBlock])
  {
    int type = elementData?VTFXA_RESMAP_ELEMENT:VTFXA_RESMAP_NODE;
    mySBlock[idBlock] = new VTFXAResultBlock(idBlock+1,
                                             VTFXA_RESTYPE_SCALAR,type);
    if (resultName) mySBlock[idBlock]->SetName(resultName);
  }
  mySBlock[idBlock]->SetResultID(idBlock);
  status = mySBlock[idBlock]->SetResultValuesBlocks(sBlockIDs.data(),
                                                    sBlockIDs.size(),iStep);
#endif
  if (VTFA_FAILURE(status))
    return showError("Error defining scalar block",idBlock);

  return true;
}


bool VTF::writeState (int iStep, const char* fmt, Real refValue, int refType)
{
#ifdef HAS_VTFAPI
  if (iStep == lastStep)
    return true;

  lastStep = iStep;
  char stepName[32];
  sprintf(stepName,fmt,refValue);
#if HAS_VTFAPI == 1
  if (!myState) myState = new VTFAStateInfoBlock();
  if (VTFA_FAILURE(myState->SetStepData(iStep,stepName,refValue,refType)))
    return showError("Error defining state info block");
#elif HAS_VTFAPI == 2
  if (!myState) myState = new VTFXAStateInfoBlock();
  if (VTFA_FAILURE(myState->AddStateInfo(iStep,stepName,refValue,
                                         VTFXA_REFVALUETYPE_TIME)))
    return showError("Error defining state info block");
#endif
#endif

  return true;
}


bool VTF::writeNodes (int iBlockID)
{
  bool ok = true;

#ifdef HAS_VTFAPI
#if HAS_VTFAPI == 1
  VTFANodeBlock nBlock(iBlockID,0);
#elif HAS_VTFAPI == 2
  VTFXANodeBlock nBlock(iBlockID,false);
#endif

  const ElementBlock* grid = myBlocks.back().second;
  if (VTFA_FAILURE(nBlock.SetNumNodes(grid->getNoNodes())))
    ok = false;

  std::vector<Vec3>::const_iterator cit;
  for (cit = grid->begin_XYZ(); cit != grid->end_XYZ() && ok; ++cit)
    if (VTFA_FAILURE(nBlock.AddNode(cit->x, cit->y, cit->z))) ok = false;

#if HAS_VTFAPI == 1
  if (VTFA_FAILURE(myFile->WriteBlock(&nBlock))) ok = false;
#elif HAS_VTFAPI == 2
  if (VTFA_FAILURE(myDatabase->WriteBlock(&nBlock))) ok = false;
#endif
#endif

  return ok;
}


bool VTF::writeElements (const char* partName, int partID,
                         int iBlockID, int iNodeBlockID)
{
#if HAS_VTFAPI == 1
  VTFAElementBlock eBlock(iBlockID,0,0);
#elif HAS_VTFAPI == 2
  VTFXAElementBlock eBlock(iBlockID,0,0);
#else
  bool eBlock = true;
#endif
  bool ok = true;

  // Lambda function writing an element block with nen nodes per element to VTF
  auto&& addElementBlock = [&eBlock](const int* mnpc, int nel, int nen)
  {
    switch (nen) {
#if HAS_VTFAPI == 1
    case 2:
      return VTFA_SUCCESS(eBlock.AddElements(VTFA_BEAMS,mnpc,nel));
    case 3:
      return VTFA_SUCCESS(eBlock.AddElements(VTFA_TRIANGLES,mnpc,nel));
    case 4:
      return VTFA_SUCCESS(eBlock.AddElements(VTFA_QUADS,mnpc,nel));
    case 8:
      return VTFA_SUCCESS(eBlock.AddElements(VTFA_HEXAHEDRONS,mnpc,nel));
#elif HAS_VTFAPI == 2
    case 2:
      return VTFA_SUCCESS(eBlock.AddElements(VTFXA_BEAMS,mnpc,nel));
    case 3:
      return VTFA_SUCCESS(eBlock.AddElements(VTFXA_TRIANGLES,mnpc,nel));
    case 4:
      return VTFA_SUCCESS(eBlock.AddElements(VTFXA_QUADS,mnpc,nel));
    case 8:
      return VTFA_SUCCESS(eBlock.AddElements(VTFXA_HEXAHEDRONS,mnpc,nel));
#else
    case 2:
    case 3:
    case 4:
    case 8:
      return eBlock;
#endif
    }
    return false;
  };

  const ElementBlock* grid = myBlocks.back().second;
  int nenod = grid->getNoElmNodes();
  if (nenod > 0)
    ok = addElementBlock(grid->getElements(),grid->getNoElms(),nenod);
  else
  {
    std::vector<int> mnpc;
    for (nenod = 2; nenod <= 8 && ok; nenod++)
      if (nenod <= 4 || nenod == 8)
        if (grid->getElements(mnpc,nenod))
          ok = addElementBlock(mnpc.data(),mnpc.size()/nenod,nenod);
  }

#if HAS_VTFAPI == 1
  eBlock.SetPartID(partID);
  eBlock.SetPartName(partName);
  eBlock.SetNodeBlockID(iNodeBlockID);
  if (VTFA_FAILURE(myFile->WriteBlock(&eBlock)))
    ok = false;
#elif HAS_VTFAPI == 2
  eBlock.SetNodeBlockID(iNodeBlockID);
  if (VTFA_FAILURE(myDatabase->WriteBlock(&eBlock)))
    ok = false;
#endif

  return ok;
}
