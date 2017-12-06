This is a simple program to compute the stability polygon of a legged robot
over n contacts. It computes a set of extremum points via cone-programming
and builds an outer and inner polygon from them. This recursive projection
method can be extended to other convex projection problems.

See the paper by Bretl and Lall, "Testing static equilibrium for legged robots"
(IEEE Transaction on Robotics, August 2008) for more details on the 2D method.

See the following paper, https://hal-lirmm.ccsd.cnrs.fr/lirmm-01477362, for more
details on the 3D method.
