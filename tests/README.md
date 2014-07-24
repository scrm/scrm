# Tests for scrm

## Automatic tests
scrm comes with a large number of tests and self checks which are automatically
executed after each commit. These tests fall into three classes:

- unit tests: Quick tests that ensure that the different building blocks of 
  scrm work correctly. Before committing, you should ensure that scrm passes the
  unit tests using `make test`.
- debug self checks: scrm has a lot of assertions and self checks which are 
  only active in debug mode. The little `test_binaries.sh` script executes a
  number of different scenarios in debug mode to ensure that everything works
  correctly. It also checks for memory leaks using valgrind. The script is
  automatically executed an travis-CI.  
- algorithm test: Here we check that the complete algorithm produces ARGs that
  are equal in distribution to those produced by ms. They are also automatically
  executed.

## Manual test
In the folder `manualtests`, there is a large number of expensive tests that we
execute from time to time.
