scrm
====

[![Build Status](https://travis-ci.org/paulstaab/scrm.png?branch=master)](https://travis-ci.org/paulstaab/scrm)

scrm is a coalescent simulator for biological sequences. Different to similar programs, 
it can approximate the Ancestral Recombination Graph as close as needed, but still has 
only linear runtime cost for long sequences. It allows you to rapidly simulate chromosome 
scale sequences with essentially correct genetic linkage.

## Installation
### Stable Release
You can download the lastest stable release of `scrm` from [scrm's homepage][1]. 
Instructions on building the binary from the source packages are available in the [wiki][3].

### Stable/Development Version From GitHub
You can also install `scrm` directly from the GitHub. 
Clone branch `stable` for the lastest stable or branch `master` for the lastet development 
version. You need to install `autoconf` to build `scrm`:  

On Debian/Ubuntu based systems:
```bash
apt-get install build-essential autoconf autoconf-archive doxygen libcppunit-dev
CXXFLAGS="-O3" ./bootstrap
make
```

On Mac OS:
```bash
port install automake autoconf autoconf-archive doxygen cppunit 
CXXFLAGS="-O3" ./bootstrap
make
```

## Usage
We designed scrm to be compatible to the famous program `ms` from Richard R. Hudson. 
You can use it as a drop in replacement for `ms` if you avoid the options `-c` and `-s`. 
Details are available [in the wiki][2]. 

## Licence
You can freely use all code in this project under the conditions of the GNU
GPL Version 3 or later.

[1]: https://scrm.github.io
[2]: https://github.com/paulstaab/scrm/wiki/Command-Line-Options
[3]: https://github.com/scrm/scrm/wiki/Installation
