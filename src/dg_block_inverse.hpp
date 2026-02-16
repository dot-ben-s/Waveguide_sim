#ifndef DG_BLOCK_INVERSE_H_
#define DG_BLOCK_INVERSE_H_
#include "mfem.hpp"
using namespace mfem;

class DGBlockInverse {
    DenseTensor inv_mass_blocks;
    const FiniteElementSpace &fes;
public:
    DGBlockInverse(FiniteElementSpace &fes_, double constant_val);
    /**
    * @brief given vector x, returns M^{-1}x in vector y
    */
    void Mult(const Vector &x, Vector &y);
};
#endif // DG_BLOCK_INVERSE_H_
