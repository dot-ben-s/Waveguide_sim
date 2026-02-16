#include "dg_block_inverse.hpp"
using namespace mfem;

DGBlockInverse::DGBlockInverse(FiniteElementSpace &fes_, double constant_val) : fes(fes_) {
  const int ne = fes.GetNE();
  const int dof = fes.GetFE(0)->GetDof();
  inv_mass_blocks.SetSize(dof, dof, ne);

  ConstantCoefficient coeff(constant_val);
  MassIntegrator mass_integ(coeff);
  DenseMatrix M_elem;
  //Iterate over each element and calculate the local inverse.
  for (int i = 0; i < ne; i++) {
    mass_integ.AssembleElementMatrix(*fes.GetFE(i), *fes.GetElementTransformation(i), M_elem);
    inv_mass_blocks(i) = M_elem;
    inv_mass_blocks(i).Invert();
  }
}
/**
 * @brief given vector x, returns M^{-1}x in vector y
 */
void DGBlockInverse::Mult(const Vector &x, Vector &y) {
  y.SetSize(x.Size());
  const int ne = fes.GetNE();
  Array<int> dofs;
  Vector x_loc, y_loc;

  for (int i = 0; i < ne; i++) {
    fes.GetElementVDofs(i, dofs);
    x.GetSubVector(dofs, x_loc);
    y_loc.SetSize(dofs.Size());
    inv_mass_blocks(i).Mult(x_loc, y_loc);
    y.SetSubVector(dofs, y_loc);
  }
}
