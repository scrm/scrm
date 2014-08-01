scrm Version History
========================

Version 1.0.1 
------------------------
Released: Not yet

### Improvements
+ Improved the handling & storage of contemporary nodes. This gives a huge
  performance boost if scrm is used with large sample sizes (>1000) (#20). 
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



Version 0.9-2
------------------------

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
