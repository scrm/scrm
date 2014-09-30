scrm Version History
========================

Version 1.3.0
------------------------
Released: Not yet

### Improvements
+ Improved the autotools configuration to support more compilers and operating
  systems (#10, #27).




Version 1.2.0
------------------------
Released: 2014-09-10

### New Features
+ New `oriented forest` summary statistic as suggested by 

    J. Kelleher, A.M. Etheridge, N.H. Barton (2014) Coalescent simulation in 
    continuous space: Algorithms for large neighbourhood size, 
    Theoretical Population Biology, Volume 95, August2014, Pages 13-23, 
    ISSN 0040-5809, http://dx.doi.org/10.1016/j.tpb.2014.05.001. 

### Improvements
+ Optimized the generation of newick trees (#22) and use a buffer (#23)


 

Version 1.1.0
------------------------
Released: 2014-08-04

### Improvements
+ Improved the handling & storage of contemporary nodes. This gives a huge
  performance boost if scrm is used with large sample sizes (>1000) (#20). 
+ Optimized scrm for use with many populations (#20). 
+ Added more automatic tests of the produced distribution of trees (#20).




Version 1.0.0 
------------------------
Released: 2014-07-09

### Bug fixes
+ Fixed an access to unmapped memory




Version 1.0-beta2
------------------------
Released: 2014-06-18

### Improvements
+ Fix file permissions
+ Remove clang warning suppression
+ Added a man page




Version 1.0-beta1
------------------------
Released: 2014-06-02

### Bug fixes:
+ Option '-es' is now ms-compatible (#16).
+ It's now possible to use 3 seeds (as in ms).

### Improvements:
+ Added help and version information.
+ Small performance tweaks.

### New Features
+ Added variable recombination & mutation rates




Version 0.9-1
------------------------

### Bug fixes
+ Very first node in the tree was assigned to wrong population
+ The input time of "-eI" option was not scaled
+ Fixed the scaling of growth rates




Version 0.9-0
------------------------

Algorithm passes all tests now, starting explicit versioning.

