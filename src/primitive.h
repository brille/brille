/*! \file */
#ifndef _PRIMITIVE_H_
#define _PRIMITIVE_H_

#include <array>
#include "spg_database.h"


const std::array<double,9> P_I_TRANSFORM{-0.5, 0.5, 0.5,  0.5,-0.5, 0.5,  0.5, 0.5,-0.5};
const std::array<double,9> P_F_TRANSFORM{  0 , 0.5, 0.5,  0.5,  0 , 0.5,  0.5, 0.5,  0 };
const std::array<double,9> P_A_TRANSFORM{  1 ,  0 ,  0 ,   0 , 0.5,-0.5,   0 , 0.5, 0.5};
const std::array<double,9> P_B_TRANSFORM{ 0.5,  0 ,-0.5,   0 ,  1 ,  0 ,  0.5,  0 , 0.5};
const std::array<double,9> P_C_TRANSFORM{ 0.5, 0.5,  0 , -0.5, 0.5,  0 ,   0 ,  0 ,  1 };
const std::array<double,9> P_R_TRANSFORM{ 0.66666666666666666667,-0.33333333333333333333,-0.33333333333333333333,  0.33333333333333333333, 0.33333333333333333333,-0.66666666666666666667,  0.33333333333333333333, 0.33333333333333333333, 0.33333333333333333333};
const std::array<double,9> P_P_TRANSFORM{  1 ,  1 ,  1,    1 ,  1 ,  1 ,   1 ,  1 ,  1 };
const std::array<int,9> P_P_INVTRNFRM{ 1, 0, 0,  0, 1, 0,  0, 0, 1};
const std::array<int,9> I_P_TRANSFORM{ 0, 1, 1,  1, 0, 1,  1, 1, 0};
const std::array<int,9> F_P_TRANSFORM{-1, 1, 1,  1,-1, 1,  1, 1,-1};
const std::array<int,9> A_P_TRANSFORM{ 1, 0, 0,  0, 1, 1,  0,-1, 1};
const std::array<int,9> B_P_TRANSFORM{ 1, 0, 1,  0, 1, 0, -1, 0, 1};
const std::array<int,9> C_P_TRANSFORM{ 1,-1, 0,  1, 1, 0,  0, 0, 1};
const std::array<int,9> R_P_TRANSFORM{ 1, 0, 1, -1, 1, 1,  0,-1, 1};

/*! A class to hold transformation matricies and their inverse, with the matrix
stored in an object is determined by a provided Centering type.

For each centred real-space-filling lattice there is a transformation matrix P
which converts its basis vectors into those of an equivalent primitive
space-filling lattice via,
    (aₚ,bₚ,cₚ) = (aₛ,bₛ,cₛ)P
This transformation matrix has an inverse matrix P⁻¹ which performs the opposite
conversion of basis vectors, or converts the associated reciprocal space-filling
lattice basis vectors to their equivalent primitive basis vectors.

This class holds P and P⁻¹ for a given centering type.
*/
class PrimitiveTransform{
private:
  Centering cen;            //!< The Centering enum value
  std::array<double,9> c2p; //!< The matrix taking centred lattice to the primitive lattice
  std::array<int,9> p2c;    //!< The inverse matrix taking the primitive lattice to the centred lattice
public:
  PrimitiveTransform(const Centering c): cen{c} { set_matrices(); };
  void set_matrices(void) {
    switch (cen){
      case BODY:
        c2p=P_I_TRANSFORM;
        p2c=I_P_TRANSFORM;
        break;
      case FACE:
        c2p=P_F_TRANSFORM;
        p2c=F_P_TRANSFORM;
        break;
      case A_FACE:
        c2p=P_A_TRANSFORM;
        p2c=A_P_TRANSFORM;
        break;
      case B_FACE:
        c2p=P_B_TRANSFORM;
        p2c=B_P_TRANSFORM;
        break;
      case C_FACE:
        c2p=P_C_TRANSFORM;
        p2c=C_P_TRANSFORM;
        break;
      case R_CENTER:
        c2p=P_R_TRANSFORM;
        p2c=R_P_TRANSFORM;
        break;
      default:
        c2p=P_P_TRANSFORM;
        p2c=P_P_INVTRNFRM;
    }
  };
  // const double* get_to_primitive_ptr()   const { return c2p.data(); };
  // const int* get_from_primitive_ptr() const { return p2c.data(); };
  std::array<double,9> get_to_primitive() const { return c2p; };
  std::array<int,9> get_from_primitive() const { return p2c; };
  //
  void print(){
    for (int i=0; i<3; ++i){
      printf("%3s", i==1 ? "to" :"");
      for (int j=0; j<3; ++j) printf(" % 4.2f",c2p[i*3+j]);
      printf("%5s", i==1 ? "from" : "");
      for (int j=0; j<3; ++j) printf(" % 4.2f",p2c[i*3+j]);
      printf("\n");
    }
  };
  //! A check if the Matrices *should* do anything, that *is not* if the centering is Primitive
  bool does_anything() { return (cen!=PRIMITIVE);};
  //! A check if the Matrices *shouldn't* do anything, that *is* if the centering is Primitive
  bool does_nothing() { return (cen==PRIMITIVE); };
};

// create new primitive transform traits
struct PrimitiveTransformTraits{
  using to = double;
  using from = int;
};

#endif
