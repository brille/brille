#ifndef _LATTICE_CLASS_H_
#define _LATTICE_CLASS_H_

#include <math.h>
#include <assert.h>
#include <string>

#include "linear_algebra.h"
#include "spg_database.h"
#include "primitive.h"

class Lattice;
class Direct;
class Reciprocal;

class Lattice{
protected:
	double len[3];
	double ang[3];
	double volume;
	int hall;
public:
	Lattice(const double *, const int h=1);
	Lattice(const double *, const double *, const int h=1);
	Lattice(const double la=1.0, const double lb=1.0, const double lc=1.0, const double al=PIOVERTWO, const double bl=PIOVERTWO, const double cl=PIOVERTWO, const int h=1);
	Lattice(const double *, const double *, const std::string);
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
	int get_hall() const {return hall;};
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

// forward declare the two types of lattices so that they can be mutually-referential

class Direct: public Lattice{
public:
	template<class ...Types> Direct(Types ... args): Lattice(args...){};
	// Direct(const double *l, const double *a): Lattice(l,a){};
	// Direct(double a, double b, double c, double d, double e, double f, double g): Lattice(a,b,c,d,e,f,g) {};
	// Direct(double la=1.0, double lb=1.0, double lc=1.0, double al=PIOVERTWO, double bl=PIOVERTWO, double cl=PIOVERTWO): Lattice(la,lb,lc,al,bl,cl) {};
	Direct(Lattice lat): Lattice(lat){};
	Reciprocal star() const;
	void get_xyz_transform(double *toxyz) const;
	void get_lattice_matrix(double*) const;
	//bool issame(const Direct lat) const;
	//bool issame(const Reciprocal) const;
	bool isstar(const Direct) const;
	bool isstar(const Reciprocal) const;
	void print() override;
	std::string string_repr() override;
	Direct primitive(void) const;
};
class Reciprocal: public Lattice{
public:
	template<class ...Types> Reciprocal(Types ... args): Lattice(args...){};
	// Reciprocal(const double *l, const double *a): Lattice(l,a){};
	// Reciprocal(double a, double b, double c, double d, double e, double f, double g): Lattice(a,b,c,d,e,f,g) {};
	// Reciprocal(double la=1.0, double lb=1.0, double lc=1.0, double al=PIOVERTWO, double bl=PIOVERTWO, double cl=PIOVERTWO): Lattice(la,lb,lc,al,bl,cl) {};
	Reciprocal(Lattice lat): Lattice(lat){};
	Direct star() const;
	void get_B_matrix(double *mt) const;
	void get_xyz_transform(double *toxyz) const;
	void get_lattice_matrix(double*) const;
	//bool issame(const Reciprocal lat) const;
	//bool issame(const Direct) const;
	bool isstar(const Reciprocal) const;
	bool isstar(const Direct) const;
	void print() override;
	std::string string_repr() override;
	Reciprocal primitive(void) const;
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
