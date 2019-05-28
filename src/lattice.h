/*! \file */
#ifndef _LATTICE_CLASS_H_
#define _LATTICE_CLASS_H_

#include <assert.h>

#include "linear_algebra.h"
#include "primitive.h"


// forward declare the two types of lattices so that they can be mutually-referential
class Lattice;
class Direct;
class Reciprocal;

/*! \brief A class to hold information about a space-spanning lattice in three dimensions

A space-spanning lattice in N dimensions has N basis vectors which can be
described fully by their N lengths and the `sum(1:N-1)` angles between each set
of basis vectors, or `sum(1:N)` scalars in total. Storing only their lengths and
angles is therefore always more efficient than storing all N² components of the
basis vectors in an orthonormal frame.
This class stores the 3 lengths and 3 angles required to describe a
3 dimensional space-spanning lattice, plus the volume of the lattice unit cell
and an optional Hall number if a non-Primitive lattice is desired.
*/
class Lattice{
protected:
  double len[3]; //!< basis vector lengths
  double ang[3]; //!< basis vector angles ordered θ₁₂, θ₀₂, θ₀₁, in radian
  double volume; //!< volume of the unit cell formed by the basis vectors
  int hall;      //!< Hall number of the non-Primitive lattice (`hall>1`)
public:
  //! Construct the Lattice from a matrix of the basis vectors
  Lattice(const double *, const int h=1);
  //! Construct the Lattice from a vector of the basis vector lengths and a vector of the basis vector angles
  Lattice(const double *, const double *, const int h=1);
  //! Construct the Lattice from the three scalar lengths and three scalar angles
  Lattice(const double la=1.0, const double lb=1.0, const double lc=1.0, const double al=PIOVERTWO, const double bl=PIOVERTWO, const double cl=PIOVERTWO, const int h=1);
  //! Construct the lattice from vectors, specifying an International Tables symmetry name instead of a Hall number
  Lattice(const double *, const double *, const std::string);
  //! Construct the lattice from scalars, specifying an International Tables symmetry name instead of a Hall number
  Lattice(const double, const double, const double, const double, const double, const double, const std::string);
  ~Lattice() = default;
  Lattice(const Lattice& other){
    for (int i=0; i<3; ++i){
      this->len[i] = other.len[i];
      this->ang[i] = other.ang[i];
    }
    this->volume = other.volume;
    this->hall = other.hall;
  };
  //! Return the first basis vector length
  double get_a     () const {return len[0];};
  //! Return the second basis vector length
  double get_b     () const {return len[1];};
  //! Return the third basis vector length
  double get_c     () const {return len[2];};
  //! Return the angle between the second and third basis vectors in radian
  double get_alpha () const {return ang[0];};
  //! Return the angle between the first and third basis vectors in radian
  double get_beta  () const {return ang[1];};
  //! Return the angle between the first and second basis vectors in radian
  double get_gamma () const {return ang[2];};
  //! Return the volume of the parallelpiped unit cell formed by the basis vectors
  double get_volume() const {return volume;};
  //! Calculate and return the unit cell volume
  double calculatevolume();
  /*! Calculate the metric tensor of the Lattice
  @param[out] mt Pointer to memory which can store 9 doubles
  */
  void get_metric_tensor(double * mt) const ;
  /*! Calculate the covariant metric tensor of the Lattice -- this is typically referred to as **the** metric tensor
  @param[out] mt Pointer to memory which can store 9 doubles
  */
  void get_covariant_metric_tensor(double *mt) const ;
  /*! Calculate the contravariant metric tensor of the Lattice -- the inverse of **the** metric tensor
  @param[out] mt Pointer to memory which can store 9 doubles
  */
  void get_contravariant_metric_tensor(double *mt) const ;
  // some functions don't logically make sense for this base class, but
  // do for the derived classes. define them here for funsies
  //! Determine if the passed Lattice represents the same space-spanning lattice
  bool issame(const Lattice) const; // this should really have a tolerance
  //! Print the basis vector lengths and angles to the console
  virtual void print();
  //! Return a string representation of the basis vector lengths and angles
  virtual std::string string_repr();
  //! Return the Hall number of the Lattice
  int get_hall() const {return hall;};
  //! Set the symmetry of the Lattice by changing the Hall number
  int set_hall(const int h) { check_hall_number(h); return hall; };
protected:
  double unitvolume() const;
  Lattice inner_star() const;
  void set_len_pointer(const double*);
  void set_ang_pointer(const double*);
  void set_len_scalars(const double, const double, const double);
  void set_ang_scalars(const double, const double, const double);
  void check_hall_number(const int h);
  void check_IT_name(const std::string itname);
};

/*! \brief A space-spanning Lattice that exists in real space

The Direct Lattice describes a space-spanning lattice in the real three
dimensions. The class is a logical wrapper to the Lattice class and
defines new versions of some methods.
*/
class Direct: public Lattice{
public:
  template<class ...Types> Direct(Types ... args): Lattice(args...){};
  Direct(Lattice lat): Lattice(lat){};
  //! Return the inverse Reciprocal lattice
  Reciprocal star() const;
  //! Return the basis vectors expressed in *an* orthonormal frame with a* along x
  void get_xyz_transform(double *toxyz) const;
  //! Return the basis vectors expressed in *an* orthonormal frame with a along x
  void get_lattice_matrix(double*) const;
  //! Always false
  bool isstar(const Direct) const;
  //! Determine if a Reciprocal lattice is the inverse of this lattice
  bool isstar(const Reciprocal) const;
  void print() override;
  std::string string_repr() override;
  //! For non-Primitive Direct lattices, return the equivalent Primitive lattice
  Direct primitive(void) const;
};
/*! \brief A space-spanning Lattice that exists in reciprocal space

The Reciprocal Lattice describes a space-spanning lattice in the reciprocal
three dimensions of momentum. The class is a logical wrapper to the Lattice
class and defines new versions of some methods.
*/
class Reciprocal: public Lattice{
public:
  template<class ...Types> Reciprocal(Types ... args): Lattice(args...){};
  Reciprocal(Lattice lat): Lattice(lat){};
  //! Return the inverse Direct lattice
  Direct star() const;
  //! Return the Busing-Levey B matrix http://dx.doi.org/10.1107/S0365110X67000970
  void get_B_matrix(double *mt) const;
  //! Return the basis vectors expressed in *an* orthonormal frame with a* along x
  void get_xyz_transform(double *toxyz) const;
  //! Return the basis vectors expressed in *an* orthonormal frame with a along x
  void get_lattice_matrix(double*) const;
  //! Always false
  bool isstar(const Reciprocal) const;
  //! Determine if a Direct lattice is the inverse of this lattice
  bool isstar(const Direct) const;
  void print() override;
  std::string string_repr() override;
  //! For non-Primitive Reciprocal lattices, return the equivalent Primitive Reciprocal lattice
  Reciprocal primitive(void) const;
};

/*! \brief Type information for Lattice and LatVec objects

Some templated functions require internal variables or return types which
depend on *which* subtype of Lattice or LatVec are provided. This traits struct
provides the typename of an appropriate Lattice subclass and its inverse for
those cases.

The two `using` types are `type` and `star`, defined based on the templated
typename as

| template typename | type | star |
| --- | --- | --- |
| Direct | Direct | Reciprocal |
| Reciprocal | Reciprocal | Direct |
| LDVec | Direct | Reciprocal |
| LQVec | Reciprocal | Direct |

*/
template <typename T> struct LatticeTraits{
  using type = void;
  using star = void;
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<> struct LatticeTraits<Direct>{
  using type = Direct;
  using star = Reciprocal;
};
template<> struct LatticeTraits<Reciprocal>{
  using type = Reciprocal;
  using star = Direct;
};
#endif

#endif
