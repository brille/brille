/* This file is part of brille.

Copyright © 2020 Greg Tucker <greg.tucker@stfc.ac.uk>

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
#ifndef BRILLE_PERMUTATION_TABLE_HPP_
#define BRILLE_PERMUTATION_TABLE_HPP_
/*! \file
    \author Greg Tucker
    \brief A class holding pairwise permutation vectors for a Grid
*/
#include <map>
#include <set>
#include <array>
#include <tuple>
#include <vector>
#include <algorithm>
namespace brille {

// /*
// For future futher optimisation we might want to use the upper triangular
// matrix to pre-allocate all possible permutation entries in the std::vector;
// There are N(N-1)/2 ordered pairs of vertices (i<j) for N total vertices which is
// still very large for large N [𝒪(0.5 N²)].
// To help with this, we will need to convert (i,j) to a linear index into the
// upper triangular part; which is complicated but already solved
// 	https://stackoverflow.com/a/27088560
// */
// static size_t upper_triangular_ij2k(size_t i, size_t j, size_t n) {
// 	return (n*(n-1)/2) - (n-i)*((n-i)-1)/2 +j -i -1;
// }
// static std::tuple<size_t, size_t> upper_triangular_k2ij(size_t k, size_t n){
// 	size_t i = n - 2 - static_cast<size_t>(std::floor(std::sqrt(-8*k + 4*n*(n-1)-7)/2 - 0.5));
// 	size_t j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
// 	return
// }

/*! \brief A class to hold pairwise permutation vectors for a connected graph

The data stored at the vertices of a pairwise connected graph may not be in the
optimal order, and reordering the data to be fully optimal may not be possible
in higher than one dimension.

To overcome this problem we can store every useful permutation vector **p**ᵢⱼ
which contains the optimal ordering of the data at vertex j given the current
order of the data at vertices `i` and `j`.
As the total number of vertices, `N`, may be large but the number of pairwise
connected vertices is likely to be considerably smaller than `N`², this class
uses a mapping of `i*N+j` → `a` where `a` ∈ (0,`number_of_connected_pairs`)
to keep track of and store only the relevant **p**ᵢⱼ vectors.
*/
class PermutationTable
{
// public:
// 	using ind_t = unsigned; // uint_fastN_t or uint_leastN_t?
protected:
	static const size_t offset{1u}; // first valid value in map
	size_t IndexSize;
	std::map<size_t,size_t> ijmap;
	std::vector<std::vector<ind_t>> permutations;
public:
  bool operator!=(const PermutationTable& other) const {
    if (IndexSize != other.IndexSize) return true;
    if (ijmap != other.ijmap) {
//      for (const auto& [key, value]: ijmap){
//        auto o_at = other.ijmap.find(key);
//        info_update_if(o_at == other.ijmap.end(), "The other map is missing key ", key);
//        info_update_if((o_at != other.ijmap.end() && value != other.ijmap.at(key)),
//          "The other value at key ",key,", ",other.ijmap.at(key),", does not match ",value);
//      }
//      for (const auto& [key, value]: other.ijmap){
//        auto o_at = ijmap.find(key);
//        info_update_if(o_at == ijmap.end(), "This map is missing key ", key);
//        info_update_if((o_at != ijmap.end() && value != ijmap.at(key)),
//          "The value at key ",key,", ",ijmap.at(key)," does not match ",value);
//      }
      return true;
    }
    if (permutations != other.permutations) return true;
    return false;
  }
  PermutationTable(size_t ni, std::map<size_t, size_t> map, std::vector<std::vector<ind_t>> perm)
      : IndexSize(ni), ijmap(std::move(map)), permutations(std::move(perm)) {}
  PermutationTable(size_t ni, size_t branches): IndexSize(ni) {
    this->add_zeroth(branches);
  };
  PermutationTable(size_t ni, size_t branches, const std::set<size_t>& kys): IndexSize(ni) {
    this->add_zeroth(branches);
    for (size_t k: kys) ijmap.emplace(k, 0u); // 0u ≡ not-yet-added value
  };
public:
	/*! \brief Reset the total number of index pairs and the length of each permuation vector

	\param ni the total number of index pairs
	\param br the length of each permutation vector
	\returns true if the number of index pairs has been changed and the
	         double-index map has been invalidated
	\note If the number of index pairs is changed by this method the double-index
	      map is cleared and this object no longer holds a valid 	permutation
				table. If the length of the permutation vectors changes the list of all
				permutations is cleared and all double-index mappings are reset to the
				identity permutation.
	*/
	bool refresh(const size_t ni, const size_t br){
		bool invalidated{false};
		if (ni != this->IndexSize){
			info_update("Resizing the PermutationTable is probably not what you wanted");
			this->ijmap.clear();
			this->IndexSize = ni;
			invalidated = true;
		}
		if (br != permutations[0].size()){
			for (auto& itr: ijmap) itr.second = 0u; // all mapped values reset to the invalid value
			permutations.clear();
			this->add_zeroth(br);
		}
		return invalidated;
	}
	//! Access the find method of the double-index map
	std::map<size_t,size_t>::const_iterator find(const size_t i, const size_t j) const {
		auto itr = this->ij2key(i,j);
		return ijmap.find(itr);
	}
	//! Check for inclusion of an index pair in the double-index map
	bool has(const size_t i, const size_t j) const {
		auto itr =this->find(i,j);
		return itr != ijmap.end() && itr->second >= offset;
	}
	//! Check if an index pair is not-present in the double-index map or has a invalid mapping
	bool value_needed(const size_t i, const size_t j) const {
		auto itr = this->find(i,j);
		return itr != ijmap.end() && itr->second < offset;
	}
	//! Insert a permutation for an index pair
	template<class I>
	size_t set(const size_t i, const size_t j, const std::vector<I>& vin){
		std::vector<ind_t> v(vin.begin(), vin.end());
		bool contains{false};
		size_t idx{0};
		size_t key = this->ij2key(i,j);
		auto itr = ijmap.find(key);
		// if the key is already present, emplace will not overwrite it anyway
		if (itr != ijmap.end()) return itr->second;
		// the ijmap does not contain key (i,j), so check for permutation duplication
		std::tie(contains, idx) = this->find_permutation(v);
		// add this permutation if it's not present
		if (!contains) permutations.push_back(v);
		// and store the mapping
		ijmap.emplace(key, idx+offset); // +offset for -offset in get
		return idx;
	}
	//! Replace an existing mapping in the double-index map
	size_t overwrite(const size_t i, const size_t j, const size_t idx){
		size_t key = this->ij2key(i,j);
		std::map<size_t,size_t>::iterator itr = ijmap.find(key);
		if (itr == ijmap.end()) throw std::runtime_error("Can not overwrite a non-existant key");
		itr->second = idx+offset;
		return idx;
	}
	//! Replace the permutation vector for an index pair in the double-index map
	template<class I>
	size_t overwrite(const size_t i, const size_t j, const std::vector<I>& vin){
		std::vector<ind_t> v(vin.begin(), vin.end());
		bool contains{false};
		size_t idx{0};
		size_t key = this->ij2key(i,j);
		std::map<size_t,size_t>::iterator itr = ijmap.find(key);
		if (itr == ijmap.end()) throw std::runtime_error("Can not overwrite a non-existant key");
		std::tie(contains, idx) = this->find_permutation(v);
		if (!contains) permutations.push_back(v);
		itr->second = idx+offset;
		return idx;
	}
	/*! \brief Retrieve the permutation vector for an index pair.

	\param i the first index of the ordered pair
	\param j the second index of the ordered pair
	\returns a permutation vector
	\note If the index pair is not present in the double-index map or it has an
	      invalid mapping the identity permutation vector is returned.
	*/
	std::vector<ind_t> safe_get(const size_t i, const size_t j) const {
		auto itr = this->find(i,j);
		bool actually_present = itr != ijmap.end() && itr->second >= offset;
		// return actually_present ? permutations[itr->second - offset] : std::vector<ind_t>();
		return permutations[actually_present ? itr->second - offset : 0];
	}
	//! Return all double-index map keys
	std::set<size_t> keys() const {
		std::set<size_t> k;
		for (auto ij: ijmap) k.insert(ij.first);
		return k;
	}
	//! Insert double-index map keys with an invalid mapping, skipping existing keys
	std::set<size_t> insert_keys(const std::set<size_t>& ks) {
		for (auto k: ks){
			std::map<size_t,size_t>::iterator itr = ijmap.find(k);
			if(itr == ijmap.end()) ijmap.emplace(k, 0u); // value=0u ≡ invalid-value
		}
		return this->keys();
	}
private:
	size_t ij2key(const size_t i, const size_t j) const { return i==j ? 0u : i*IndexSize+j; }
	//
	std::tuple<bool, size_t> find_permutation(const std::vector<ind_t>& v) const {
		size_t N = v.size();
		auto equal_perm = [&](const std::vector<ind_t>& p){
			if (p.size() != N) return false;
			for (size_t i=0; i<N; ++i)
				if (p[i] != v[i]) return false;
			return true;
		};
		auto itr = std::find_if(permutations.begin(), permutations.end(), equal_perm);
		// returns (true, found_index) or (false, permutations.size())
		return std::make_tuple(itr != permutations.end(), std::distance(permutations.begin(), itr));
	}
	void add_zeroth(size_t branches) {
		std::vector<ind_t> identity(branches);
		std::iota(identity.begin(), identity.end(), 0);
		if (permutations.size()<1) permutations.resize(1);
		permutations[0] = identity;
		ijmap[0u] = 0u + offset;// offset first valid value
	}

#ifdef USE_HIGHFIVE
      public:
        template<class R>
        std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, bool>
        to_hdf(R& obj, const std::string& dataset) const {
          auto group = overwrite_group(obj, dataset);
          bool ok{true};
          group.createAttribute("size", IndexSize);
          ok &= map_to_hdf(ijmap, group, "map");
          ok &= lists_to_hdf(permutations, group, "permutations");
          return ok;
        }
        [[nodiscard]] bool
        to_hdf(const std::string& filename, const std::string& dataset, const unsigned perm=HighFive::File::OpenOrCreate) const{
          HighFive::File file(filename, perm);
          return this->to_hdf(file, dataset);
        }
        template<class R>
        static std::enable_if_t<std::is_base_of_v<HighFive::Object, R>, PermutationTable>
        from_hdf(R& obj, const std::string& dataset){
          auto group = obj.getGroup(dataset);
          size_t isz;
          group.getAttribute("size").read(isz);
          auto m = map_from_hdf<size_t, size_t>(group, "map");
          auto p = lists_from_hdf<ind_t>(group, "permutations");
          return {isz, m, p};
        }
        static PermutationTable from_hdf(const std::string& filename, const std::string& dataset){
          HighFive::File file(filename, HighFive::File::ReadOnly);
          return PermutationTable::from_hdf(file, dataset);
        }
#endif
};

/*! \brief Build a list of all possible keys from an iterable container

If a container holding a list of indices `cont = [a,b,c,...,z]` is provided then
this function will find all unique ordered pairs of indices

{(a,b), (b,a), (a,c), (c,a), ..., (a,z), (z,a), (b,c), (c,b), ... }

with each key value (x,y) = x*n + y.

\param i_beg the starting iterator, e.g., `cont.begin()`
\param i_end the ending iterator, .e.g., `cont.end()`
\param n     the number-of-index-pairs for the double-index mapping, could be
             calculated from, e.g., `i_end - i_beg`, but may also be something
						 different.
*/
template<typename Itr>
std::set<size_t>
permutation_table_keys_from_indicies(Itr i_beg, Itr i_end, const size_t n){
	std::set<size_t> keys;
	for (Itr j=i_beg; j!=i_end; ++j) for (Itr k=j+1; k!=i_end; ++k) if (*j!=*k){
		keys.insert((*j)*n + (*k));
		keys.insert((*j) + (*k)*n);
	}
	return keys;
}

} // namespace brille
#endif
