#ifndef MAXWELL_EVOLUTION_H_
#define MAXWELL_EVOLUTION_H_
#include "mfem.hpp"
#include "dg_block_inverse.hpp"

using namespace mfem;

/** A time dependent operator for the right hand side of the ODE.
 *  The weak form of the PDE is M du/dt = (K + F)u + phi(t), where M is the mass-, K the
 *  curl and F the flux matrix. This is rewritten as the general ODE,
 *  du/dt = M^{-1} [(K + F)u + phi(t)], and this class returns the right hand side of
 *  equation. */

class MaxwellEvolution : public TimeDependentOperator
{
private:
  BlockOperator &K;
  DGBlockInverse &M_inv_Hx, &M_inv_Hy, &M_inv_Ez;

  LinearForm *f_inlet_ez;
  LinearForm *f_inlet_hy;
  Array<int> offsets;
  double f0, epsilon0, mu0, Z0;

  mutable Vector z;
public:
  MaxwellEvolution(BlockOperator &K_,
                   DGBlockInverse &MInvHx, DGBlockInverse &MInvHy, DGBlockInverse &MInvEz,
                   LinearForm *f_inlet_ez,LinearForm *f_inlet_hy, Array<int> offsets_);

  virtual ~MaxwellEvolution() = default;

  //Explicit step: y = f(x,t)
  virtual void Mult(const Vector &x, Vector &y) const;
};

#endif // MAXWELL_EVOLUTION_H_
