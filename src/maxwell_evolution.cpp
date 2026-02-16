#include "maxwell_evolution.hpp"
using namespace mfem;

MaxwellEvolution::MaxwellEvolution(BlockOperator &K_,
                   DGBlockInverse &MInvHx, DGBlockInverse &MInvHy, DGBlockInverse &MInvEz,
                   LinearForm *f_inlet_ez,LinearForm *f_inlet_hy, Array<int> offsets_)
    : TimeDependentOperator(K_.Height(),0.0),
      K(K_),
      M_inv_Hx(MInvHx), M_inv_Hy(MInvHy), M_inv_Ez(MInvEz),
      f_inlet_ez(f_inlet_ez),
      f_inlet_hy(f_inlet_hy),
      offsets(offsets_),
      z(K_.Height())
    {
      f0 = 11.0e9;
      mu0 = 1.2566370614e-6;
      epsilon0 = 8.8541878176e-12;
      Z0 = sqrt(mu0 / epsilon0);
    }

void MaxwellEvolution::Mult(const Vector &x, Vector &y) const{

    z = 0.0;
    K.Mult(x,z);
    z.Neg();
    //z = -(K+F)u

    double sigma_t = 1.0/f0;
    double t0 = 4.0 * sigma_t;
    double signal = exp(-pow((t-t0) / sigma_t, 2)) * sin(2.0 *M_PI * f0 * t);
    double source_val =  (1.0 / Z0) * signal;
    Vector z_Ez(z.GetData() + offsets[2], offsets[3] - offsets[2]);
    Vector z_Hy(z.GetData() + offsets[1], offsets[2]- offsets[1]);

    z_Ez.Add(source_val, *f_inlet_ez);
    z_Hy.Add(-Z0 * source_val, *f_inlet_hy);
    //z = -(K+F)u + phi(t)

    int size_Hx = offsets[1]- offsets[0];
    Vector z_Hx_block(z.GetData() + offsets[0], size_Hx);
    Vector y_Hx_block(y.GetData() + offsets[0], size_Hx);
    M_inv_Hx.Mult(z_Hx_block, y_Hx_block);

    int size_Hy = offsets[2]- offsets[1];
    Vector z_Hy_block(z.GetData() + offsets[1], size_Hy);
    Vector y_Hy_block(y.GetData() + offsets[1], size_Hy);
    M_inv_Hy.Mult(z_Hy_block, y_Hy_block);

    int size_Ez = offsets[2]- offsets[1];
    Vector z_Ez_block(z.GetData() + offsets[2], size_Ez);
    Vector y_Ez_block(y.GetData() + offsets[2], size_Ez);
    M_inv_Ez.Mult(z_Ez_block, y_Ez_block);
    //z = M^{-1}[(-K-F)u + phi(t)]
  }
