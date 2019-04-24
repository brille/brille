#include "arrayvector.h"
#include "lattice.h"
#include "latvec.h"

#ifndef _BZ_CLASS_H_
#define _BZ_CLASS_H_

class BrillouinZone {
	Reciprocal lattice;
//	ArrayVector<double> position;
//	ArrayVector<int> type;
	ArrayVector<double> vertices;
	ArrayVector<int> faces;            // the reciprocal lattice points defining the faces -- twice the *actual* face vectors
	ArrayVector<int> faces_per_vertex;
public:
	BrillouinZone(Reciprocal lat, int extent=1);
	void set_vertices(ArrayVector<double> newverts);
	void set_faces(ArrayVector<int> newfaces);
	void set_faces_per_vertex(ArrayVector<int> newfpv);
	void determine_everything(const int extent=1);
	void determine_everything_xyz(const int extent=1);
	size_t vertices_count() const { return vertices.size();};
	size_t faces_count() const { return faces.size();};
	size_t get_vertices_bare        (const size_t max, double *out) const;
	size_t get_faces_bare           (const size_t max, int *out) const ;
	size_t get_faces_per_vertex_bare(const size_t max, int *out) const;
	const Reciprocal get_lattice() const { return this->lattice;};
	LQVec<double>    get_vertices() const;
	LQVec<int>       get_faces() const ;
	ArrayVector<int> get_faces_per_vertex() const;
	void print() const;

	template<typename T> ArrayVector<bool> isinside(const LQVec<T> *p, const double tol=1e-14);
	bool moveinto(const LQVec<double> *Q, LQVec<double> *q, LQVec<int> *tau);
};


bool between_origin_and_plane_xyz(const ArrayVector<double> *p,
	                                const ArrayVector<double> *v,
															    const ArrayVector<int> *ijk, const int idx, ArrayVector<double> *inv,
															    const int store_at=0, const double tol=1e-15);
bool three_plane_intersection_xyz(const ArrayVector<double> *n,
	                                const ArrayVector<double> *p,
															    const int i, const int j, const int k,
																	ArrayVector<double> *iat, const int idx);

bool between_origin_and_plane(const LQVec<double> *p,
	                            const LQVec<double> *v,
															const ArrayVector<int> *ijk,
															const int idx,
															      LQVec<double> *inv,
															const int store_at=0,
															const double tol=1e-15);
bool three_plane_intersection(const LQVec<double> *n,
	                            const LQVec<double> *p,
															const ArrayVector<double> *xyz,
															const int i, const int j, const int k,
															LQVec<double> *iat, const int idx=0);

int make_all_indices(LQVec<int> *ijk, const int extent=1);

#endif
