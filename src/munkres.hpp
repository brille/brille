/* This file is part of brille.

Copyright © 2019,2020 Greg Tucker <greg.tucker@stfc.ac.uk>

brille is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

brille is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with brille. If not, see <https://www.gnu.org/licenses/>.            */

#ifndef BRILLE_MUNKRES_H_
#define BRILLE_MUNKRES_H_
#include <forward_list>
#include <limits>
#include <memory>
#include <vector>
namespace brille::assignment {

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
  Munkres(size_t n, std::vector<T> cmat = std::vector<T>()) : N(n), cost(cmat), step(1){
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
    for (size_t i=0; i < cost.size(); ++i) sum+=std::abs(cost[i]);
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
        case 7: finished=true; done=true; break;
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
    for (size_t i=0; i<N*N; ++i) if (marker::PRIMED==mask[i]) mask[i]=marker::NORMAL;
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
    auto smallest = (std::numeric_limits<T>::max)();
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

} // end namespace brille::assignment
#endif
