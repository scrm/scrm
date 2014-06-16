scrm
====

[![Build Status](https://travis-ci.org/paulstaab/scrm.png?branch=master)](https://travis-ci.org/paulstaab/scrm)

scrm is an implementation of the Sequential-Coalescent-with-Recombination Model.


## Installation
On Debian/Ubuntu based systems:
```bash
apt-get install build-essential autoconf autoconf-archive doxygen libcppunit-dev
CXXFLAGS="-O3" ./bootstrap
make
```

On Mac OS:
```bash
port install autoconf autoconf-archive doxygen cppunit 
./bootstrap
make

```

## Usage
We designed scrm to behave compatible to the famous program `ms` from Richard R. Hudson. 
You can use it as a drop in replacement for `ms` if you avoid the options `-c` and `-s`. 
Details are available [in the wiki][1]. 

## Licence
You can freely use all code in this project under the conditions of the GNU
GPL Version 3 or later.

[1]: https://github.com/paulstaab/scrm/wiki/Command-Line-Options
