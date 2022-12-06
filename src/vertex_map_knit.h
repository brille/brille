//
// Created by g on 02/12/22.
//

#ifndef BRILLE_VERTEX_MAP_KNIT_H
#define BRILLE_VERTEX_MAP_KNIT_H

#include "vertex_map_set.h"
#include "vertex_index_map.h"

namespace brille::vertex_maps {

/*
 * Each VertexMapSet contains:
 *    - pristine  : a shared copy of the original vertex set
 *    - appended  : a unique copy of extra vertices
 *    - preserved : a unique record of when each of pristine was encountered
 *
 *    The format of preserved is, for X = pristine.size(0),+1 something like
 *      [X, ..., X, 2, X, ..., 1, 0, ..., N, ... X, ...]
 *    it contains integers from (0,N) for the N+1 vertices in pristine
 *    encountered plus some X (which may be as small as pristine.size(0))
 *    for vertices *not* seen when checking for generated vertices against
 *    known vertices
 *
 *
 * Each VertexIndexMap contains:
 *    - a map connecting node index to a fixed-size inner map (a vector)
 *      which links node-vertex index to 'final' consolidated vertex index
 *
 *    Each inner map of a VertexIndexMap is of the form
 *      [44, 9, 0, ..., N, N+10, ... , M-3]
 *    where M - N is the total number of vertices added to the appended list.
 *
 *    Each inner map index corresponds to the *order* in which vertices were
 *    encountered, so care must be taken when combining two independently
 *    constructed VertexIndexMaps!
 *
 *
 * */

template <class T, template<class> class A>
std::pair<VertexMapSet<T,A>, VertexIndexMap>
combine(const std::pair<VertexMapSet<T,A>, VertexIndexMap>& one,
        const std::pair<VertexMapSet<T,A>, VertexIndexMap>& two){
  auto no = one.first.pristine_count();

  std::cout << "Combine\n" << one.first << one.second << "with\n" << two.first << two.second << "\n";

  // copy the preserved vectors to avoid modifying the inputs
  typename VertexMapSet<T,A>::preserve_t one_prev, two_prev;
  one_prev.reserve(no);
  two_prev.reserve(no);
  for (const auto & x: one.first.preserved()) one_prev.emplace_back(x);
  for (const auto & x: two.first.preserved()) two_prev.emplace_back(x);

  // Check if there are any preserved indexes that show up in both one and two
  std::vector<int> differ;
  differ.reserve(no);
  std::transform(one_prev.begin(), one_prev.end(),
                 two_prev.begin(), std::back_inserter(differ),
                 [no](ind_t a, ind_t b){return (a<no && b<no && a!=b) ? 1 : 0;});

  // Copy-construct new maps to avoid modifying input
  VertexIndexMap one_second(one.second.data()), two_second(two.second.data());
  // Move the appended vertex indices to account for the total number preserved
  // and, in the case of `two`, the vertices in one.first.appended
  auto duplicate_count = std::count(differ.begin(), differ.end(), true);
  auto total_preserved = one.first.preserved_count() + two.first.preserved_count() - duplicate_count;
  profile_update(duplicate_count, " differ pristine vertices -- ", total_preserved, " preserved in total");
  // Use a lambda to avoid (minimal) duplication
//  auto offset_appended = [no, total_preserved](auto & map, ind_t upto, ind_t appended){
//    for (ind_t i=no; i < no + upto; ++i)
//      map.replace(i, i + (total_preserved + appended - no));
//  };
  // we need to keep the indexes beyond the end of pristine for now
  auto offset_appended = [no, total_preserved](auto & map, ind_t upto, ind_t appended){
    for (ind_t i=no; i < no + upto; ++i)
      map.replace(i, i + appended);
  };
  auto one_appended_count = one.first.appended_count();
  auto two_appended_count = two.first.appended_count();
  //offset_appended(one_second, one_appended_count, no);
  offset_appended(two_second, two_appended_count, one_appended_count);
  profile_update("appended vertex counts: ", one_appended_count, " ", two_appended_count);

  // Replace seen second indices by their seen first equivalents
  if (duplicate_count){
    for (ind_t i=0; i<no; ++i) if (differ[i]) {
      two_second.replace(two_prev[i], one_prev[i]);
    }
  }

  // Handle preserved indexes that appear only in the second
  std::vector<int> only_second;
  only_second.reserve(no);
  std::transform(one_prev.begin(), one_prev.end(),
                 two_prev.begin(), std::back_inserter(only_second),
                 [no](ind_t a, ind_t b){return (a >= no && b < no) ? 1 : 0;});


  profile_update("we now know which",std::count(only_second.begin(), only_second.end(), true)," pristine vertices appear only in the second");

//  profile_update(only_second);
//  profile_update(two_prev);

  // Build up the partial sum of preceding vertex indices which are in both
  ind_t os_max{0};
  for (ind_t i=0; i<no; ++i) if (only_second[i] && two_prev[i]>os_max) os_max = two_prev[i];
  std::vector<ind_t> os_partial(os_max+2, 0);
  // If the preserved index differs *or* is the same and is not only in two
  // the vertex index should be included in the partial sum
  for (ind_t i=0; i<no; ++i)
    if (!only_second[i] && two_prev[i] <= no)
    //if (differ[i] || (!only_second[i] && one_prev[i] < no && one_prev[i] == two_prev[i]))
      os_partial[two_prev[i]] = 1;

  profile_update(os_partial);
  std::partial_sum(os_partial.begin(), os_partial.end(), os_partial.begin());
  profile_update(os_partial);

  std::vector<ind_t> os_value(os_max+1, 0); //no+1);
  for (ind_t i=0; i<no; ++i) if (only_second[i]) os_value[two_prev[i]] = two_prev[i] - os_partial[two_prev[i]];

  profile_update(os_value);

  profile_update("We now have the cumulative sum of sorted only second");
//  auto two_into_one_offset = one.first.preserved_count();
  auto two_into_one_offset = no + one.first.preserved_count();
  for (ind_t i=0; i<no; ++i) if (only_second[i]) {
    auto updated = os_value[two_prev[i]];
    if (two_prev[i] != updated) {
      // earlier vertices were removed; so update the map *and* the preserved list
      two_second.replace(two_prev[i], updated);
    }
    // plus update one_prev[i], which *is* currently no+1 but should be
    // the reduced two-preserved index plus the one-preserved count
    one_prev[i] = updated + two_into_one_offset;
    debug_update_if(one_prev[i] > no+1, i, " -> ", one_prev[i], " (from two_prev[i]=",two_prev[i]," -> ", updated ," + ", two_into_one_offset,")");
  }

  profile_update("All mapped indices are correct");

  // All mapped indexes are correct, so move on to combining the VertexMapSets
  // 0. The pristine vertices are as needed
  // 1. Combine the two preserved vectors into one, done above.
  // 2. Combine the two sets of appended vertices

  // The naive approach to combine the appended vertices is insufficient
  // since there *could* be equivalent vertices in the two lists :/
  // auto appended = cat(0, one.first.appended(), two.first.appended());

  // Reuse the VertexMapSet machinery to sort this out for us!
  // We're done with one_prev, so it's safe to move
  VertexMapSet<T,A> out_vms(one.first.pristine(), one.first.appended(), std::move(one_prev), one.first.relative(), one.first.digits());
  profile_update("Output VertexMapSet created");
  auto appended = two.first.appended();
  ind_t two_app_no = appended.size(0);
  // reserve space now to avoid doubling out_vms.appended (multiple times)
  out_vms.reserve_appended(out_vms.appended_count() + appended.size(0));

  std::vector<std::pair<ind_t, ind_t>> a_map;
  a_map.reserve(two_app_no);
  ind_t a_skipped{0};
  for (ind_t i=0; i<two_app_no; ++i){
    // what would have been the mapped index for this one?
    // [new_kept_pristine_count] + i
    auto j = out_vms.add(appended.view(i), AddVertexType::Crafted);
    auto x = i + total_preserved + one_appended_count;
    if (j != x) ++a_skipped;
    a_map.emplace_back(x, (j!=x ? j : x - a_skipped));
  }
  // out_vms now *has* all of our vertices!

  profile_update("Output VertexMapSet given all vertices");

  // The second map *might* have to be updated _again_
  for (const auto & pair: a_map) if (pair.first != pair.second) {
    two_second.replace(pair.first, pair.second);
  }

  profile_update("Second map indexes updated");

  // We can now combine the two maps, checking to ensure they have no shared keys
  // VertexIndexMap::merge returns true if the second data_ map is not emptied
  // which only happens if the two have equivalent keys
  if (one_second.merge(two_second)){
    throw std::runtime_error("Two VertexIndexMaps were not mutually exclusive!");
  }

  profile_update("Output maps merged");

  std::cout << "Resulting in\n" << out_vms << one_second << "\n";

  return std::make_pair(out_vms, one_second);
}

namespace steps {
template <class T, template <class> class A>
std::vector<std::pair<VertexMapSet<T, A>, VertexIndexMap>>
parallel_reduce(const std::vector<std::pair<VertexMapSet<T, A>, VertexIndexMap>> &input)
{
  std::vector<std::pair<VertexMapSet<T, A>, VertexIndexMap>> output;
  auto out_len = static_cast<int64_t>(input.size() >> 1);
  auto isodd = static_cast<size_t>(2 * out_len + 1) == input.size();
  output.reserve((isodd ? out_len + 1 : out_len));
#pragma omp parallel for default(none) shared(input, output, out_len)
  for (int64_t idx = 0; idx < out_len; ++idx) {
    output.push_back(combine(input[2 * idx], input[2 * idx + 1]));
  }
  if (isodd) output.push_back(input.back());
  return output;
}

template <class T, template <class> class A>
std::vector<std::pair<VertexMapSet<T, A>, VertexIndexMap>>
reduce(const std::vector<std::pair<VertexMapSet<T, A>, VertexIndexMap>> &input)
{
  std::vector<std::pair<VertexMapSet<T, A>, VertexIndexMap>> output;
  auto out_len = input.size() >> 1;
  auto isodd = (2 * out_len + 1) == input.size();
  output.reserve((isodd ? out_len + 1 : out_len));
  for (size_t idx = 0; idx < out_len; ++idx) {
    output.push_back(combine(input[2 * idx], input[2 * idx + 1]));
  }
  if (isodd) output.push_back(input.back());
  return output;
}

template<class T, template<class> class A>
std::pair<VertexMapSet<T,A>, VertexIndexMap>
consolidate_pair(std::pair<VertexMapSet<T,A>, VertexIndexMap>& pair){
  auto vim = pair.second;
  // shift mapped indexes above the highest pristine vertex index down
  // in preparation for consolidating the vertex map
  auto pri = pair.first.pristine_count();
  auto pre = pair.first.preserved_count();
  auto app = pair.first.appended_count();
  for (ind_t i=pri; i < pri + app; ++i) vim.replace(i, i + pre - pri);

  auto vms = pair.first.consolidate();

  return std::make_pair(vms, vim);
}

}; // namespace steps

template <class T, template<class> class A>
std::pair<VertexMapSet<T,A>, VertexIndexMap>
parallel_reduce(const std::vector<std::pair<VertexMapSet<T,A>, VertexIndexMap>>& input){
  debug_update_if(input.size() < 1, "Zero-length vector will cause an error!");
  if (input.size() == 1) return input[0];
  auto workspace = steps::parallel_reduce(input);
  while (workspace.size() > 1) workspace = steps::parallel_reduce(workspace);
  return steps::consolidate_pair(workspace[0]);
}

template <class T, template<class> class A>
std::pair<VertexMapSet<T,A>, VertexIndexMap>
reduce(const std::vector<std::pair<VertexMapSet<T,A>, VertexIndexMap>>& input){
  debug_update_if(input.size() < 1, "Zero-length vector will cause an error!");
  if (input.size() == 1) return input[0];
  auto workspace = steps::reduce(input);
  while (workspace.size() > 1) workspace = steps::reduce(workspace);
  return steps::consolidate_pair(workspace[0]);
}


};

#endif // BRILLE_VERTEX_MAP_KNIT_H
