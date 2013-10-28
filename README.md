scrm
====

[![Build Status](https://travis-ci.org/paulstaab/scrm.png?branch=master)](https://travis-ci.org/paulstaab/scrm)

scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.

## Licence
You can freely use all code in this project unter the conditions of the GNU
GPLv3+ with the following restriction: 

You are not allowed not use scrm to create derived software (as defined in the GPL) about which you write an publication in an peer
revied journal without our explicit agreement. Please contact me in this case.

This restriction will be removed once we published a paper about scrm ourselfs.

## Installation
On Debian/Ubuntu based systems:
```bash
apt-get install build-essential libboost-all-dev autoconf autoconf-archive doxygen libcppunit-dev
CXXFLAGS="-O3" ./bootstrap
make
```

On Mac OS:
```bash
port install boost autoconf autoconf-archive doxygen cppunit 
./bootstrap
make

```

=======
