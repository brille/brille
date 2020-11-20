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

/* An implementation of a solver for the stable matching problem */


#ifndef BRILLE_SMP_HPP_
#define BRILLE_SMP_HPP_
#include <iostream>
namespace brille::assignment {

template<typename T, typename P> void smp_pmat(const T dim, const P *restrict tp){
  for (T i=0; i<dim; ++i){
    for (T j=0; j<dim; ++j) std::cout << " " << std::to_string(tp[i*dim + j]);
    std::cout << std::endl;
  }
}
template<typename T, typename P> void smp_pvec(const T dim, const P *restrict tp){
  for (T j=0; j<dim; ++j) std::cout << " " << std::to_string(tp[j]);
  std::cout << std::endl;
}

template<typename T, typename P>
int smp(const T dim, const P *restrict prefs, T *restrict rsol, T *restrict csol, const bool verbose){
  T nsel=0;
  auto rpref = std::unique_ptr<T[]>(new T[dim*dim]);
  auto cpref = std::unique_ptr<T[]>(new T[dim*dim]);
  P min;
  auto tmp = std::unique_ptr<P[]>(new P[dim]);
  auto inv = std::unique_ptr<T[]>(new T[dim]);
  // Create the preference lists.
  for (T i=0; i<dim; ++i){
    // rpref are based on the rows of prefs
    for (T j=0; j<dim; ++j) tmp[j] = prefs[i*dim +j];
    for (T j=0; j<dim; ++j){
      min = tmp[0];
      rpref[i*dim+j] = T(0);
      for (T k=0; k<dim; ++k) if (tmp[k] < min){
        rpref[i*dim+j] = k;
        min = tmp[k];
      }
      tmp[rpref[i*dim+j]] = (std::numeric_limits<P>::max)();
    }
    // cpref are based on the columns of pref
    for (T j=0; j<dim; ++j) tmp[j] = prefs[j*dim + i];
    for (T j=0; j<dim; ++j){
      min = tmp[0];
      cpref[i*dim+j] = T(0);
      for (T k=0; k<dim; ++k) if (tmp[k] < min){
        cpref[i*dim+j] = k;
        min = tmp[k];
      }
      tmp[cpref[i*dim+j]] = (std::numeric_limits<P>::max)();
    }
    /* cpref contains the row indexes in preferrential order, but it needs
    to contain the prefferential index in row order. */
    for (T j=0; j<dim; ++j) for (T k=0; k<dim; ++k) if (j==cpref[i*dim+k]) inv[j] = k;
    for (T j=0; j<dim; ++j) cpref[i*dim+j] = inv[j];
  }
  if (verbose){
    std::cout << "For a preference matrix:" << std::endl;
    for (T i=0; i<dim; ++i){
      for (T j=0; j<dim; ++j) std::cout << " " << std::to_string(prefs[i*dim + j]);
      std::cout << std::endl;
    }
    std::cout << "We find rpref:" << std::endl;
    for (T i=0; i<dim; ++i){
      for (T j=0; j<dim; ++j) std::cout << " " << std::to_string(rpref[i*dim + j]);
      std::cout << std::endl;
    }
    std::cout << "and cpref:" << std::endl;
    for (T i=0; i<dim; ++i){
      for (T j=0; j<dim; ++j) std::cout << " " << std::to_string(cpref[i*dim + j]);
      std::cout << std::endl;
    }
  }

  // flags to indicate if a matching has been proposed
  auto rmatched = std::unique_ptr<bool[]>(new bool[dim]);
  auto cmatched = std::unique_ptr<bool[]>(new bool[dim]);
  auto rproposed = std::unique_ptr<T[]>(new T[dim]);
  // intialize: no proposed matchings, all rows and columns try for their preferred match
  for (T i=0; i<dim; ++i){
    rmatched[i] = false;
    cmatched[i] = false;
    rproposed[i] = T(0);
    rsol[i] = T(0);
    csol[i] = T(0);
  }
  T thisc;
  T count =0;
  // as long as not all matchings have been proposed
  while (nsel < dim && count < dim*dim*dim){
    ++count;
    // check each row in order
    for(T i=0; i<dim; ++i){
      // if it hasn't been matched
      if (!rmatched[i]){
        if (rproposed[i] >= dim) throw std::runtime_error("last choice was rejected?!");
        // find its preferred column of non-rejected proposals
        thisc = rpref[i*dim + rproposed[i]];
        // if the preferred column has a proposed match
        if (cmatched[thisc]){
          // and if the column prefers this row over the proposed row
          if (cpref[i*dim + i] < cpref[i*dim + csol[thisc]]){
            // unmatch the other row
            rmatched[csol[thisc]] = false;
            rsol[csol[thisc]] = T(0);
            // make sure it moves on to its next choice
            rproposed[csol[thisc]]++;
            // and match this one
            rmatched[i] = true;
            csol[thisc] = i;
            rsol[i] = thisc;
          } else {
            // othwerwise move on to the next preferred column
            rproposed[i]++;
          }
        } else {
          // the column hasn't been proposed yet, so tentatively accept
          rsol[i] = thisc;
          csol[thisc] = i;
          rmatched[i] = true;
          cmatched[thisc] = true;
          nsel++; // we have one more match than before
        }
      }
    }
  }
  bool ok=true;
  for (size_t i=0; i<dim; ++i) ok &= rmatched[i] && cmatched[i];
  if (!ok){
    std::cout << "Not all matches found?!" << std::endl;
  }
  if (verbose){
    std::cout << "Matches are:" << std::endl;
    std::cout << "rsol:";
    for (T j=0; j<dim; ++j) std::cout << " " << std::to_string(rsol[j]);
    std::cout << std::endl;
    std::cout << "csol:";
    for (T j=0; j<dim; ++j) std::cout << " " << std::to_string(csol[j]);
    std::cout << std::endl;
  }
  return 0;
}

} // end namespace brille::assignment
#endif
