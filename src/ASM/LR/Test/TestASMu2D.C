//==============================================================================
//!
//! \file TestASMu2D.C
//!
//! \date Jul 14 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Tests for driver for assembly of unstructured 2D spline FE models.
//!
//==============================================================================

#include "ASMu2D.h"
#include "SIM2D.h"

#include "gtest/gtest.h"


struct EdgeTest
{
  int edge;
  int edgeIdx;
  std::array<int,2> c1;
  std::array<int,2> c2;
};


class TestASMu2D : public testing::Test,
                   public testing::WithParamInterface<EdgeTest>
{
};


class TestSIM2D : public SIM2D
{
public:
  TestSIM2D() : SIM2D(1)
  {
    opt.discretization = ASM::LRSpline;
    EXPECT_TRUE(this->read("src/ASM/LR/Test/refdata/boundary_nodes.xinp"));
    EXPECT_TRUE(this->createFEMmodel());
  }
  virtual ~TestSIM2D() {}
};


TEST_P(TestASMu2D, ConstrainEdge)
{
  TestSIM2D sim;
  ASMu2D* pch = static_cast<ASMu2D*>(sim.getPatch(1));
  ASSERT_TRUE(pch != nullptr);

  pch->constrainEdge(GetParam().edgeIdx, false, 1, 1, 1);

  std::vector<int> glbNodes;
  pch->getBoundaryNodes(GetParam().edge, glbNodes, 1);

  for (int node : glbNodes)
    EXPECT_TRUE(pch->findMPC(node,1) != nullptr);
}


TEST_P(TestASMu2D, ConstrainEdgeOpen)
{
  TestSIM2D sim;
  ASMu2D* pch = static_cast<ASMu2D*>(sim.getPatch(1));
  ASSERT_TRUE(pch != nullptr);

  pch->constrainEdge(GetParam().edgeIdx, true, 1, 1, 1);

  std::vector<int> glbNodes;
  pch->getBoundaryNodes(GetParam().edge, glbNodes, 1);

  int crn = pch->getCorner(GetParam().c1[0], GetParam().c1[1], 1);
  EXPECT_TRUE(pch->findMPC(crn,1) == nullptr);
  glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));
  crn = pch->getCorner(GetParam().c2[0], GetParam().c2[1], 1);
  EXPECT_TRUE(pch->findMPC(crn,1) == nullptr);
  glbNodes.erase(std::find(glbNodes.begin(), glbNodes.end(), crn));

  for (int node : glbNodes)
    EXPECT_TRUE(pch->findMPC(node,1) != nullptr);
}


static const std::vector<EdgeTest> edgeTestData =
        {{{1, -1, {-1, -1}, {-1 , 1}},
         {2,  1, { 1, -1}, { 1 , 1}},
         {3, -2, {-1, -1}, { 1, -1}},
         {4,  2, {-1,  1}, { 1,  1}}}};


INSTANTIATE_TEST_CASE_P(TestASMu2D, TestASMu2D,
                        testing::ValuesIn(edgeTestData));
