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

class PermutationTable
{
public:
	using ind_t = unsigned; // uint_fastN_t or uint_leastN_t?
private:
	static const size_t offset{1u}; // first valid value in map
	size_t IndexSize;
	std::map<size_t,size_t> ijmap;
	std::vector<std::vector<ind_t>> permutations;
public:
	PermutationTable(size_t ni, size_t branches): IndexSize(ni) {
		this->add_zeroth(branches);
	};
	PermutationTable(size_t ni, size_t branches, const std::set<size_t>& kys): IndexSize(ni) {
		this->add_zeroth(branches);
		for (size_t k: kys) ijmap.emplace(k, 0u); // 0u ≡ not-yet-added value
	};
public:
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
	std::map<size_t,size_t>::const_iterator find(const size_t i, const size_t j) const {
		auto itr = this->ij2key(i,j);
		return ijmap.find(itr);
	}
	bool has(const size_t i, const size_t j) const {
		auto itr =this->find(i,j);
		return itr != ijmap.end() && itr->second >= offset;
	}
	bool value_needed(const size_t i, const size_t j) const {
		auto itr = this->find(i,j);
		return itr != ijmap.end() && itr->second < offset;
	}
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
	size_t overwrite(const size_t i, const size_t j, const size_t idx){
		size_t key = this->ij2key(i,j);
		std::map<size_t,size_t>::iterator itr = ijmap.find(key);
		if (itr == ijmap.end()) throw std::runtime_error("Can not overwrite a non-existant key");
		itr->second = idx+offset;
		return idx;
	}
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
	std::vector<ind_t> safe_get(const size_t i, const size_t j) const {
		auto itr = this->find(i,j);
		bool actually_present = itr != ijmap.end() && itr->second >= offset;
		// return actually_present ? permutations[itr->second - offset] : std::vector<ind_t>();
		return permutations[actually_present ? itr->second - offset : 0];
	}
	std::set<size_t> keys() const {
		std::set<size_t> k;
		for (auto ij: ijmap) k.insert(ij.first);
		return k;
	}
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
};

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
