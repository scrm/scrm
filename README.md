scrm
====

[![Build Status](https://travis-ci.org/paulstaab/scrm.png?branch=master)](https://travis-ci.org/paulstaab/scrm)

scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.


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

## Usage
We designed scrm to behave similar to the famous program
[ms](http://home.uchicago.edu/~rhudson1/source/mksamples.html). 
The following command line options behave the same effect as in ms:
"-r", "-I", "-M", "-eM", "-m", "-em", "-G", "-eG", "-g", "-eg", 
"-ej", "-t", "-T", "-L"
so if your familiar with ms, then your are almost good to go.

The option "-es" is a bit different to it's ms equivalent. Scrm does not support
simulating a fixed number of segregating sites (ms' "-s") and gene conversions
("-c"). 

Additionally, scrm has an option that ms does not have. Use "-l" followed by an
non-negative integer to give the length of an "exact window". This is the
sequence part that scrm will simulate according to the full Ancestral
Recombination Graph. Lines of more sequence positions further away will be
ignored to speedup the simulation process. 

## Licence
You can freely use all code in this project under the conditions of the GNU
GPL Version 3 or later.
