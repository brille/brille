/*! \file */
#ifndef __MUNKRES_H_
#define __MUNKRES_H_

#include <forward_list>
#include <limits>
#include <memory>
#include <vector>

enum marker {
  NORMAL,
  STARED,
  PRIMED,
};
/*! \brief A single-case implementation of the Munkres' Assignment Algorithm

    See, e.g., http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
*/
template<typename T> class Munkres{
private:
  size_t N;                  //!< dimensionality of the problem
  std::vector<T> cost;       //!< cost matrix, to be optimized
  std::vector<marker> mask;  //!< mask matrix, to indicate star and prime zeros
  int step;                  //!< the current step of the algorithm
  std::vector<int> rowcover; //!< indicator for if a row is covered
  std::vector<int> colcover; //!< indicator for if a column is covered
  size_t stored_row;         //!< stored row index
  size_t stored_col;         //!< stored column index
  bool finished;             //!< flag to indicate if the algorithm has succeeded
public:
  /*! The Munkres intializer
      @param n the dimensionality of the system, which must be square
      @param cmat [optional] the cost matrix

      If the cost matrix is provided the assignment algorithm is run automatically.
  */
  Munkres(size_t n, std::vector<T> cmat = std::vector<T>()) : N(n), step(1), cost(cmat){
    if (cost.size()<N*N) cost.resize(N*N);
    // mask, rowcover, and colcover have their values set when run_assignment is called
    mask.resize(N*N);
    rowcover.resize(N);
    colcover.resize(N);
    // try to run the assignment algorithm
    finished = run_assignment();
  }
  void reset(){
    // Make sure the mask and covers are set correctly for the start of a run
    for (size_t i=0; i<N*N; ++i) mask[i] = marker::NORMAL;
    for (size_t i=0; i<N; ++i){
      rowcover[i] = 0;
      colcover[i] = 0;
    }
    stored_row = 0;
    stored_col = 0;
    // We should start from the beginning:
    // on the first step
    step = 1;
    // and not yet finished
    finished = false;
  }
  /*! Start the assignment algorithm.

  Resets the mask matrix, cover vectors, stored indices, step, and finished state
  before attempting to run the assignment algorithm.
  @returns A flag to indicate if the algorithm finished successfully
  */
  bool run_assignment(){
    reset();
    // check that we *can* perform an assignment
    T sum = 0;
    for (size_t i=0; i < cost.size(); ++i) sum+=abs(cost[i]);
    // If there is cost information we're not done yet.
    bool done = (sum>0) ? false : true;
    while (!done) {
      // show();
      switch (step) {
        case 1: step_one();   break;
        case 2: step_two();   break;
        case 3: step_three(); break;
        case 4: step_four();  break;
        case 5: step_five();  break;
        case 6: step_six();   break;
        case 7: finished=true;
        default: done=true;
      }
    }
    return finished;
  }
  //! Get a mutable reference to the cost matrix, useful for filling.
  std::vector<T>& get_cost(void){ return cost;}
  //! Get an immutable reference to the cost matrix
  const std::vector<T>& get_cost(void) const { return cost; }
  /*! After the assignment algorithm has run, use get_assignment to retrieve
      the result
      @param[out] out A size_t pointer at which the N task assignments will be stored
      @returns The returned flag is true if the assignment was made.
  */
  bool get_assignment(size_t* out){
    if (!finished) run_assignment();
    if (finished){
      for (size_t r=0; r<N; ++r)
        for (size_t c=0; c<N; ++c)
          if (mask[r*N+c]==marker::STARED) out[r]=c;
    }
    return finished;
  }
  std::string to_string() const {
    std::string x = "X", star = "*", prime = "'", blank = "";
    std::string repr= "step=" + std::to_string(step) + "\t";
    for (size_t c=0; c<N; ++c) repr += (colcover[c] ? x : blank) + "\t";
    repr += "\n";
    for (size_t r=0; r<N; ++r){
      repr += (rowcover[r] ? x : blank) + "\t";
      for (size_t c=0; c<N; ++c)
        repr += std::to_string(cost[r*N+c]) + (mask[r*N+c]==PRIMED ? prime : mask[r*N+c]==STARED ? star : blank) + "\t";
      repr += "\n";
      }
    return repr;
  }
private:
  void show(){
    printf("s=%d\t",step);
    for (size_t c=0; c<N; ++c) printf("%s\t", colcover[c] ? "X" : "" );
    printf("\n");
    for (size_t r=0; r<N; ++r){
      printf("%s\t", rowcover[r] ? "X" : "");
      for (size_t c=0; c<N; ++c)
        printf("%5.2f%s\t", cost[r*N+c], mask[r*N+c]==PRIMED ? "'" : mask[r*N+c]==STARED ? "*" : "");
      printf("\n");
    }
  }
  void find_a_zero(long long &row, long long &col){
    row = -1;
    col = -1;
    bool notdone = true;
    for (size_t r=0; r<N && notdone; ++r)
      for (size_t c=0; c<N && notdone; ++c)
        if (cost[r*N+c]==0 && rowcover[r]==0 && colcover[c]==0){
          row = (long long) r;
          col = (long long) c;
          notdone = false;
        }
  }
  long long get_col_of_star_in_row(const long long &row){
    bool sir = false;
    long long col = -1;
    size_t r = (size_t) row;
    for (size_t c=0; c<N && !sir; ++c) if (mask[r*N+c]==marker::STARED){
      sir = true;
      col = (long long) c;
    }
    return col;
  }
  void step_one(){
    // for each row of the cost matrix, find the smallest element
    // and subtract it from every element in its row.
    T smallest;
    for (size_t r=0; r<N; ++r){
      smallest = cost[r*N];
      for (size_t c=1; c<N; ++c) if(cost[r*N+c]<smallest) smallest=cost[r*N+c];
      for (size_t c=0; c<N; ++c) cost[r*N+c]-=smallest;
    }
    // With step one finished, go on to step 2
    step = 2;
  }
  void step_two(){
    // Find a zero (Z) in the resulting matrix. If there is no starred zero
    // in its row or column, star Z. Repeat for each element in the matrix.
    T zero = 0;
    for (size_t r=0; r<N; ++r)
      for (size_t c=0; c<N; ++c)
        if (cost[r*N+c] == zero && rowcover[r]==0 && colcover[c]==0){
          mask[r*N+c] = marker::STARED;
          rowcover[r] = 1;
          colcover[c] = 1;
        }
    // reset the cover vectors, since we abused them in this step:
    for (size_t i=0; i<N; ++i){
      rowcover[i]=0;
      colcover[i]=0;
    }
    // Step two finished. Go on to three.
    step = 3;
  }
  void step_three(){
    // Cover each column containing a starred zero.
    // If N columns are covered, the starred zeros describe a complete set of
    // unique assignments, and we are done.
    for (size_t r=0; r<N; ++r)
      for (size_t c=0; c<N; ++c)
        if (mask[r*N+c]==marker::STARED) colcover[c]=1;
    size_t count=0;
    for (size_t c=0; c<N; ++c) if (colcover[c]==1) count++;
    // If we've found N starred zeros, go to step 7; otherwise step 4
    step = (count >= N) ? 7 : 4;
  }
  void step_four(){
    // Find a noncovered zero and prime it. If there is no starred zero in the
    // row containing this primed zero, go to step 5.
    // Otherwise, cover this row and uncover the column containing the starred
    // zero. Continue in this manner until there are no uncovered zeros left.
    // Save the smallest uncovered value and go to step 6.
    long long row = -1;
    long long col = -1;
    bool done = false;
    while (!done){
      find_a_zero(row,col);
      if (row < 0){
        done = true;
        step = 6;
      } else {
        size_t r = (size_t)row;
        size_t c = (size_t)col;
        mask[r*N+c] = marker::PRIMED;
        col = get_col_of_star_in_row(row);
        if (col>=0) {
          c = (size_t)col;
          rowcover[r] = 1;
          colcover[c] = 0;
        } else {
          done = true;
          step = 5;
          stored_row = r;
          stored_col = c;
        }
      }
    }
  }
  void step_five(){
    // Construct a series of alternating primed and starred zeros as follows:
    // Let Z0 represent the uncovered primed zero found in step 4.
    // Let Z1 denote the starred zero in the column of Z0 (if any).
    // Let Z2 denote the primed zero in the row of Z1 (always present if Z1 exists)
    // Continue until the series terminates at a primed zero that has no starred
    // zero in its column.
    // Unstar each starred zero of the series.
    // Star each primed zero of the series.
    // Erase all primes and uncover every line in the matrix.
    // Return to step three.

    // path contains pairs of row/column values where we have found
    // either a star or a prime that is part of the alternating sequence.
    std::forward_list<std::pair<size_t, size_t>> path {{stored_row, stored_col}};


    size_t rowcol[2] = {0,stored_col}; // the starting point is a PRIMED (found previously in step 4).
    const marker whichmark[2] = {marker::STARED, marker::PRIMED};
    // starting with i=0 --> we look first along a row for a STARED zero.
    for (size_t i=0; rowcol[i]<N; ++rowcol[i]){
      if (mask[rowcol[0]*N+rowcol[1]] == whichmark[i]){
        path.push_front({rowcol[0],rowcol[1]});
        i = (i+1)&1;    // switch between row and col, and between STARED and PRIMED
        rowcol[i] = -1; // exploit unsiged integer roll-over, and that
      }
    }

    size_t idx;
    for (const auto &i : path){
      idx = i.first*N + i.second;
      // unstar each starred zero.
      if (marker::STARED == mask[idx]) mask[idx] = marker::NORMAL;
      // star each primed zero.
      if (marker::PRIMED == mask[idx]) mask[idx] = marker::STARED;
    }
    // erase all primes.
    for (size_t i=0; i<N*N; ++i) if (marker::PRIMED==mask[i]) mask[i]==marker::NORMAL;
    // uncover every line
    for (size_t i=0; i<N; ++i){
      rowcover[i] = 0;
      colcover[i] = 0;
    }
    // return to step three;
    step = 3;
  }
  void step_six(){
    // Add the value found in step 4 to every element of each covered row
    // and subtract it from every element of each uncovered column.
    // Return to step 4 without altering any stars, primes or covered lines.
    auto smallest = std::numeric_limits<T>::max();
    for (size_t r=0; r<N; ++r)
      for (size_t c=0; c<N; ++c)
        if (rowcover[r]==0 && colcover[c]==0 && cost[r*N+c] < smallest) smallest=cost[r*N+c];
    for (size_t r=0; r<N; ++r)
      for (size_t c=0; c<N; ++c){
        if (rowcover[r] == 1) cost[r*N+c] += smallest;
        if (colcover[c] == 0) cost[r*N+c] -= smallest;
      }
    step=4;
  }

};

/*! \brief Type information for cost-matrix elements used by Munkres.

The Munkres' algorithm, as implemented here, can only handle real-valued
costs. Since it is desirable to make assignments of complex-valued data we
need means by which to identify the underlying real data type.

| template typename | type | max |
| T | T | std::numeric_limits<T>::max() |
| std::complex<T> | T | std::numeric_limits<T>::max() |
*/
template<class T> struct MunkresTraits{
  using type = T;
  constexpr static T max = std::numeric_limits<T>::max();
};
#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class T> struct MunkresTraits<std::complex<T>>{
  using type = T;
  constexpr static T max = std::numeric_limits<T>::max();
};
#endif


/*! \brief Use Munkres' Assignment algorithm to determine a permutation

For two arrays of data, located in memory at `centre` and `neighbour`, and each
representing `Nobj` sets of  `Nscl` scalars, `Nvec` vector elements, and `Nmat`
matrix elements, determine the sorting permutation that maps the elements at
`centre` onto the same global mapping as the sorting permutation already stored
at `permutations[neighbour_idx]`. The resultant global permutation is stored
at `permutations[centre_idx]`.

Each array of data to be compared must be formatted as:

      [ 0{(0…Nscl-1)(0…Nvec-1)(0…Nmat*Nmat-1)}
        1{(0…Nscl-1)(0…Nvec-1)(0…Nmat*Nmat-1)}
        ⋮
        Nobj-1{(0…Nscl-1)(0…Nvec-1)(0…Nmat*Nmat-1)} ]

The function constructs an `Nobj`×`Nobj` cost matrix, where each element is
given by

      Cᵢⱼ = Wscl×∑ₖ₋₀ᴺˢᶜˡ√(centre[i,k]-neighbour[j,k])
          + Wvec×𝔣ᵥ(centre_vec[i],neighbour_vec[j])
          + Wmat×𝔣ₘ(centre_mat[i],neighbour_mat[j])

where `Wscl`, `Wvec`, and `Wmat` are weight factors for adjusting the relative
cost of the scalar, vector, and matrix differences, respectively;
𝔣ᵥ is the vector cost function; and 𝔣ₘ is the matrix cost function.

The vector cost function is selected by `vec_cost_func` to be one of:
the angle between the vectors, the length of the difference between the vectors,
or one minus the inner product of the vectors.

The matrix cost function is the Frobenius norm of the difference between the
two matrices.

@param centre The values to be sorted by the determined permutation
@param neighbour The values against which the `centre` values are compared
@param Nscl The number of scalars per object
@param Neig The number of eigenvectors per object
@param Deig The eigenvector dimensionality
@param Nvec The number of vector elements per object
@param Nmat The square root of the number of matrix elements per object
@param Wscl The cost weight for scalar elements
@param Weig The cost weight for eigvenvectors
@param Wvec The cost weight for vectors
@param Wmat The cost weight for matrices
@param span The number of elements per object, `Nscl+Nvec+Nmat*Nmat`
@param Nobj The number of objects at each of `centre` and `neighbour`
@param[out] permutations Contains the global permutation sorting for the
                         neighbour and is where the centre permutation is stored
@param centre_idx The index into `permutations` where the output is stored
@param neighbour_idx The index into `permutations` to find the neighbour permutation
@param eig_cost_func Used to select the eigenvector cost function:
                      0 --> `vector_angle`
                      1 --> `vector_distance`
                      2 --> `1-vector_product`
                      3 --> `abs(sin(hermitian_angle))`
@param vec_cost_func Used to select the vector cost function:
                      0 --> `vector_angle`
                      1 --> `vector_distance`
                      2 --> `1-vector_product`
@returns `true` if the permutation was assigned successfully, otherwise `false`
*/
template<class T, class R,
          typename=typename std::enable_if<std::is_same<typename MunkresTraits<T>::type, R>::value>::type
        >
bool munkres_permutation(const T* centre, const T* neighbour,
                         const size_t Nscl,
                         const size_t Neig, const size_t Deig,
                         const size_t Nvec, const size_t Nmat,
                         const R Wscl, const R Weig, const R Wvec, const R Wmat,
                         const size_t span, const size_t Nobj,
                         ArrayVector<size_t>& permutations,
                         const size_t centre_idx, const size_t neighbour_idx,
                         const int eig_cost_func = 3,
                         const int vec_cost_func = 0
                       ){
// initialize variables
R s_cost{0}, e_cost{0}, v_cost{0}, m_cost{0};
Munkres<R> munkres(Nobj);
size_t* assignment = new size_t[Nobj]();
size_t any_evm = Deig*Neig + Nvec + Nmat*Nmat;
// TODO: change the input from (N,M,...,Nobj) to (Nobj,N,M,...)
//       so that the per-object information is contiguous in memory.
//       Then we can get rid of c_data and n_data completely.
T* c_data = new T[any_evm]();
T* n_data = new T[any_evm]();
// calculate costs and fill the Munkres cost matrix
for (size_t i=0; i<Nobj; ++i){
  for (size_t j=0; j<Nobj; ++j){
    s_cost = R(0);
    e_cost = R(0);
    if (Nscl){
      for (size_t k=0; k<Nscl; ++k)
        s_cost += magnitude(centre[i+k*Nobj] - neighbour[j+k*Nobj]);
    }
    if (any_evm){
      for (size_t k=0; k<any_evm; ++k){
        c_data[k] =    centre[i+(Nscl+k)*Nobj];
        n_data[k] = neighbour[j+(Nscl+k)*Nobj];
      }
    }
    if (Neig*Deig){
      for (size_t k=0; k<Neig; ++k){
        switch (eig_cost_func){
          case 0: e_cost +=     vector_angle(Deig,c_data+k*Deig,n_data+k*Deig); break;
          case 1: e_cost +=  vector_distance(Deig,c_data+k*Deig,n_data+k*Deig); break;
          case 2: e_cost += 1-vector_product(Deig,c_data+k*Deig,n_data+k*Deig); break;
          case 3: e_cost += std::abs(std::sin(hermitian_angle(Deig,c_data+k*Deig,n_data+k*Deig))); break;
        }
      }
    }
    if (Nvec){
      switch (vec_cost_func){
        case 0: v_cost =     vector_angle(Nvec,c_data+Neig*Deig,n_data+Neig*Deig); break;
        case 1: v_cost =  vector_distance(Nvec,c_data+Neig*Deig,n_data+Neig*Deig); break;
        case 2: v_cost = 1-vector_product(Nvec,c_data+Neig*Deig,n_data+Neig*Deig); break;
      }
    }
    if (Nmat){
      m_cost = frobenius_distance(Nmat,c_data+Neig*Deig+Nvec,n_data+Neig*Deig+Nvec);
    }
    // for each i we want to determine the cheapest j
    munkres.get_cost()[i*Nobj+j] = Wscl*s_cost + Weig*e_cost + Wvec*v_cost + Wmat*m_cost;
  }
}
// clear variables in case assignment fails and we exit early
delete[] c_data;
delete[] n_data;
// use the Munkres' algorithm to determine the optimal assignment
munkres.run_assignment();
if (!munkres.get_assignment(assignment)){
  delete[] assignment;
  return false;
}
/* use the fact that the neighbour objects have already had their global
   permutation saved into `permutations` to determine the global permuation
   for the centre objects too; storing the result into `permutations` as well.
*/
for (size_t i=0; i<Nobj; ++i)
  for (size_t j=0; j<Nobj; ++j)
    if (permutations.getvalue(neighbour_idx,i)==assignment[j])
      permutations.insert(j,centre_idx,i);
delete[] assignment;
return true;
}

#endif
