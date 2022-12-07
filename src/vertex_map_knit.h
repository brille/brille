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

  // copy the preserved vectors to avoid modifying the inputs
  typename VertexMapSet<T,A>::preserve_t one_prev, two_prev;
  for (const auto & x: one.first.preserved()) one_prev.emplace(x);
  for (const auto & x: two.first.preserved()) two_prev.emplace(x);

  // Check if there are any preserved indexes that show up in both one and two
  // but which differ in their assigned index:
  std::map<ind_t, std::pair<ind_t, ind_t>> differ;
  const auto & tp{two_prev};
  for (const auto& p: one_prev) {
    if (auto s = tp.find(p.first); s != tp.end() && s->second != p.second)
      differ.emplace(p.first, std::make_pair(p.second, s->second));
  }

  // Copy-construct new maps to avoid modifying input
  VertexIndexMap one_second(one.second.data()), two_second(two.second.data());
  // Move the appended vertex indices to account for the total number preserved
  // and, in the case of `two`, the vertices in one.first.appended
//  auto different_count = std::count(differ.begin(), differ.end(), true); // required if a vector
  auto different_count = differ.size(); // the map *only* has the different vertexes
  auto total_preserved = one.first.preserved_count() + two.first.preserved_count() - different_count;
  profile_update(different_count, " differently indexed pristine vertices -- ", total_preserved, " preserved in total");

  // Collect the preserved indices that appeared only in the second view
  std::vector<std::pair<ind_t, ind_t>> os_vector;
  os_vector.reserve(std::max(one_prev.size(), two_prev.size()));
  const auto & op{one_prev};
  std::set_difference(tp.cbegin(), tp.cend(),
                      op.cbegin(), op.cend(), std::back_inserter(os_vector),
                      [](const auto & a, const auto &b){return a.first < b.first;});

  // Sort these by their mapped value
  std::sort(os_vector.begin(), os_vector.end(), [](auto & a, auto & b){return a.second < b.second;});

  // Create an updated map with their 'correct' index in it:
  std::map<ind_t, std::pair<ind_t, ind_t>> only_second;
  ind_t os_new_val{one.first.preserved_count()};
  for (const auto & x: os_vector) only_second.emplace(x.first, std::make_pair(x.second, os_new_val++));

  // we need to keep the indexes beyond the end of pristine for now
  auto offset_appended = [no, total_preserved](auto & map, ind_t upto, ind_t appended){
    for (ind_t i=no; i < no + upto; ++i)
      map.replace(i, i + appended);
  };
  auto one_appended_count = one.first.appended_count();
  auto two_appended_count = two.first.appended_count();
  offset_appended(two_second, two_appended_count, one_appended_count);
  profile_update("appended vertex counts: ", one_appended_count, " ", two_appended_count);

  // Replace second seen indices by their first seen equivalents
  for (const auto & ifs: differ) two_second.replace(ifs.second.second, ifs.second.first);

  profile_update("We now have the cumulative sum of sorted only second");
  //  auto two_into_one_offset = one.first.preserved_count();
  auto two_into_one_offset = no + one.first.preserved_count();

  for (const auto & p: only_second) {
    // if earlier vertices were removed, update the map and the preserved list
    if (p.second.first != p.second.second) two_second.replace(p.second.first, p.second.second);
    // And update on_prev since we'll use it again:
    one_prev[p.first] = p.second.second + two_into_one_offset;
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
  auto two_app_no = appended.size(0);
  auto ono = out_vms.preserved_count();
  // reserve space now to avoid doubling out_vms.appended (multiple times)
  out_vms.reserve_appended(out_vms.appended_count() + two_app_no);
  std::vector<std::pair<ind_t, ind_t>> app_from_to;
  app_from_to.reserve(two_app_no);
  for (ind_t i=0; i<two_app_no; ++i){
    // what is the index assigned for this vertex above?
    // its [length original pristine] + [length appended in one] + i
    auto from = no + one_appended_count + i;
    // Find an already appended equivalent vertex, or append this one
    // return the index in the appended array offset by the number of pristine
    // vertices, either way
    // *and convert it to the index *after* condensing the end result!*
    auto to = out_vms.add(appended.view(i), AddVertexType::Crafted) + ono - no;
    app_from_to.emplace_back(from, to);
  }
  // out_vms now *has* all of our vertices!

  profile_update("Output VertexMapSet given all vertices");

  // Update the second map to account fo condensing the output:
  for (const auto &x: app_from_to) two_second.replace(x.first, x.second);

  profile_update("Second map indexes updated");

  // We can now combine the two maps, checking to ensure they have no shared keys
  // VertexIndexMap::merge returns true if the second data_ map is not emptied
  // which only happens if the two have equivalent keys
  if (one_second.merge(two_second)){
    throw std::runtime_error("Two VertexIndexMaps were not mutually exclusive!");
  }

  profile_update("Output maps merged");
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
