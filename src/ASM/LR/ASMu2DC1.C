// $Id$
//==============================================================================
//!
//! \file ASMu2DC1.C
//!
//! \date Oct 5 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of C1-continuous 2D LR-spline FE models.
//!
//==============================================================================

#include "LRSpline/LRSplineSurface.h"

#include "ASMu2DC1.h"
#include "Function.h"
#include "Vec3Oper.h"
#include "Vec3.h"
#include "MPC.h"


bool ASMu2DC1::generateFEMTopology ()
{
  if (lrspline)
    if (lrspline->order(0) > 2 && lrspline->order(1) > 2)
      return this->ASMu2D::generateFEMTopology();

  std::cerr <<" *** ASMu2DC1::generateFEMTopology:"
            <<" The polynomial order is too low.\n    "
            <<" C1-continuity requires at least quadratic order."<< std::endl;
  return false;
}


void ASMu2DC1::constrainEdge (int dir, bool open, int dof, int code, char)
{
  // Figure out what edge we are at
  DirichletEdge de(lrspline.get(), dir, dof, code);

  // Get all basis functions on this boundary edge.
  // If dof >= 100 we also need the layer of nodes just inside the boundary.
  std::vector<LR::Basisfunction*> edgeFunctions;
  de.lr->getEdgeFunctions(edgeFunctions, de.edg, dof >= 100 ? 2 : 1);

  for (LR::Basisfunction* b : edgeFunctions)
    std::cout << b->getId() <<": "<< this->getCoord(1+b->getId()) << std::endl;

  // Add constraints for all basis functions on the edge
  for (LR::Basisfunction* b : edgeFunctions)
    if (!open || !de.isCorner(1+b->getId())) // skip corners for open boundaries
    {
      if (dof%100)
        this->prescribe(1+b->getId(),dof%100,code);
      if (dof >= 100)
      {
	int node = 1; //TODO, get the inside node here
        if (dof%100 && code == 0)
          // The edge is clamped, fix the neighboring node line
          this->prescribe(node,dof/100,0);
        else
          // The edge has a prescribed rotation, add an MPC for that
          this->add2PC(MLGN[node-1],dof/100,MLGN[b->getId()-1],code);
      }
    }
}


bool ASMu2DC1::updateDirichlet (const std::map<int,RealFunc*>& func,
                                const std::map<int,VecFunc*>& vfunc,
                                double time, const std::map<int,int>* g2l)
{
  // Update the constraint equations defining Dirichlet conditions
  // on 1st-derivatives (normal rotations). Note: We here assume that
  // all MPCs with one master DOF and non-zero slave coefficient denote
  // prescribed 1st derivatives.
  std::map<int,RealFunc*>::const_iterator fit;
  std::map<int,VecFunc*>::const_iterator vfit;
  for (MPCMap::iterator cit = dCode.begin(); cit != dCode.end(); cit++)
    if (cit->first->getNoMaster() == 1)
    {
      size_t inod = this->getNodeIndex(cit->first->getSlave().node);
      size_t jnod = this->getNodeIndex(cit->first->getMaster(0).node);
      if (inod < 1 || jnod < 1) return false;

      // Evaluate the prescribed rotation
      double theta = 0.0;
      Vec4 X(this->getCoord(jnod),time);
      if ((fit = func.find(cit->second)) != func.end())
        theta = (*fit->second)(X);
      else if ((vfit = vfunc.find(cit->second)) != vfunc.end())
        theta = (*vfit->second)(X)[cit->first->getSlave().dof-1];
      else
      {
        std::cerr <<" *** ASMu2DC1::updateDirichlet: Code "<< cit->second
                  <<" is not associated with any function."<< std::endl;
        return false;
      }
      // Update the slave coefficient, s = |X_j-X_i|*tan(theta)
      cit->first->setSlaveCoeff((this->getCoord(inod) - X).length()*tan(theta));
    }

  return this->ASMu2D::updateDirichlet(func,vfunc,time,g2l);
}
