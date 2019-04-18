#ifndef _LATTICE_CLASS_H_
#define _LATTICE_CLASS_H_

#include <math.h>
#include <assert.h>
#include "linear_algebra.h"
#include <string>

class Lattice;
class Direct;
class Reciprocal;

class Lattice{
protected:
	double len[3];
	double ang[3];
	double volume;
public:
	Lattice(const double *, const double *);
	Lattice(double, double, double, double, double, double, double);
	Lattice(double la=1.0, double lb=1.0, double lc=1.0, double al=PIOVERTWO, double bl=PIOVERTWO, double cl=PIOVERTWO);
	~Lattice() = default;
	// Lattice(const Lattice alat);
	double get_a     () const {return len[0];};
	double get_b     () const {return len[1];};
	double get_c     () const {return len[2];};
	double get_alpha () const {return ang[0];};
	double get_beta  () const {return ang[1];};
	double get_gamma () const {return ang[2];};
	double get_volume() const {return volume;};
	double calculatevolume();
	void get_metric_tensor(double * mt) const ;
	void get_covariant_metric_tensor(double *mt) const ;
	void get_contravariant_metric_tensor(double *mt) const ;
	// some functions don't logically make sense for this base class, but
	// do for the derived classes. define them here for funsies
	bool issame(const Lattice) const; // this should really have a tolerance
	virtual void print();
	virtual std::string string_repr();

protected:
	Lattice inner_star() const;
};

// forward declare the two types of lattices so that they can be mutually-referential

class Direct: public Lattice{
public:
	Direct(const double *l, const double *a): Lattice(l,a){};
	Direct(double a, double b, double c, double d, double e, double f, double g): Lattice(a,b,c,d,e,f,g) {};
	Direct(double la=1.0, double lb=1.0, double lc=1.0, double al=PIOVERTWO, double bl=PIOVERTWO, double cl=PIOVERTWO): Lattice(la,lb,lc,al,bl,cl) {};
	Direct(Lattice lat): Lattice(lat){};
	Reciprocal star() const;
	void get_xyz_transform(double *toxyz) const;
	//bool issame(const Direct lat) const;
	//bool issame(const Reciprocal) const;
	bool isstar(const Direct) const;
	bool isstar(const Reciprocal) const;
	void print() override;
	std::string string_repr() override;
};
class Reciprocal: public Lattice{
public:
	Reciprocal(const double *l, const double *a): Lattice(l,a){};
	Reciprocal(double a, double b, double c, double d, double e, double f, double g): Lattice(a,b,c,d,e,f,g) {};
	Reciprocal(double la=1.0, double lb=1.0, double lc=1.0, double al=PIOVERTWO, double bl=PIOVERTWO, double cl=PIOVERTWO): Lattice(la,lb,lc,al,bl,cl) {};
	Reciprocal(Lattice lat): Lattice(lat){};
	Direct star() const;
	void get_B_matrix(double *mt) const;
	void get_xyz_transform(double *toxyz) const;
	//bool issame(const Reciprocal lat) const;
	//bool issame(const Direct) const;
	bool isstar(const Reciprocal) const;
	bool isstar(const Direct) const;
	void print() override;
	std::string string_repr() override;
};

template <typename T> struct LatticeTraits{
	using type = void;
	using star = void;
};
template<> struct LatticeTraits<Direct>{
	using type = Direct;
	using star = Reciprocal;
};
template<> struct LatticeTraits<Reciprocal>{
	using type = Reciprocal;
	using star = Direct;
};

#endif
