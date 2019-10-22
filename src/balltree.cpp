#include "balltree.h"
#include "debug.h"

BallNode construct_balltree(const std::vector<BallLeaf>& leaves, const size_t max_depth){
  debug_update("Construct a tree from ",leaves.size()," BallLeaf objects");
  std::array<double,3> axis;
  BallNode root = construct_ballnode(leaves, axis);
  bifurcate_balltree(root, leaves, max_depth, root.centre(), axis);
  debug_update("Tree construction finished");
  return root;
}

BallNode construct_ballnode(const std::vector<BallLeaf>& leaves, std::array<double,3>& v){
  std::array<double,3> x, c{0.,0.,0.};
  // find the centroid of the leaf positions
  for (size_t i=0; i<3u; ++i) c[i] = 0.;
  for (auto leaf: leaves){
    x = leaf.centre();
    for (size_t i=0; i<3u; ++i) c[i] += x[i];
  }
  for (size_t i=0; i<3u; ++i) c[i] /= static_cast<double>(leaves.size());
  // plus the vector pointing to the farthest leaf from the centroid
  // and the radius fully encompasing all leaves;
  double d, dmax = 0., radius = 0.;
  for (auto leaf: leaves){
    x = leaf.centre();
    d = 0.;
    for (size_t i=0; i<3u; ++i) d += (x[i]-c[i])*(x[i]-c[i]);
    d = std::sqrt(d);
    if (d>dmax){
      dmax = d;
      for (size_t i=0; i<3u; ++i) v[i] = x[i] - c[i];
    }
    if (d+leaf.radius() > radius) radius = d+leaf.radius();
  }
  debug_update("BallNode encompasing ",leaves.size()," BallLeaf objects will be centred at ",c," with radius ",radius);
  BallNode node(c, radius);
  return node;
}

bool bifurcate_balltree(BallNode& root, const std::vector<BallLeaf>& leaves, const size_t max_depth, const std::array<double,3>& at, const std::array<double,3>& along){
  if (leaves.size() < 1) return true;
  if (max_depth == 0 || leaves.size() == 1){ // we've hit bottom or can't branch further
    debug_update("Giving up with ",max_depth," bifurcations to go and ",leaves.size()," BallLeaf objects");
    root.leaves(leaves);
    debug_update("This terminal node now contains ",root.count_leaves()," leaves");
  } else {
    debug_update("Bifurcating at ",at," along ",along);
    std::vector<BallLeaf> upper, lower;
    std::array<double,3> x;
    double dot;
    for (auto leaf: leaves){
      dot = 0.;
      x = leaf.centre();
      for (size_t i=0; i<3u; ++i) dot += (x[i]-at[i])*along[i];
      if (dot > leaf.radius() ) // was > 0.
        upper.push_back(leaf);
      else
        lower.push_back(leaf);
    }
    debug_update("BallLeafs split into ",upper.size()," and ",lower.size()," objects in the upper and lower branches");
    std::array<double,3> axis;
    // create new nodes if we might need to bifurcate further
    if (upper.size()>1){
      debug_update("Constructing upper branch");
      BallNode upper_branch = construct_ballnode(upper, axis);
      bifurcate_balltree(upper_branch, upper, max_depth-1, upper_branch.centre(), axis);
      debug_update("Adding an upper child at height ",max_depth);
      root.addchild(upper_branch);
    }
    if (lower.size()>1){
      debug_update("Constructing lower branch");
      BallNode lower_branch = construct_ballnode(lower, axis);
      bifurcate_balltree(lower_branch, lower, max_depth-1, lower_branch.centre(), axis);
      debug_update("Adding a lower child at height ",max_depth);
      root.addchild(lower_branch);
    }
    // single element vectors can not be bifurcated further and do not need their own nodes
    if (upper.size()==1) root.addleaf(upper[0]);
    if (lower.size()==1) root.addleaf(lower[0]);
  }
  debug_update("Finishing bifurcation at height ",max_depth);
  return false;
}
