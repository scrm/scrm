This is the jz_stable branch

scrm
====

scrm is a coalescent simulator for biological sequences. Different to similar programs, 
it can approximate the Ancestral Recombination Graph as close as needed, but still has 
only linear runtime cost for long sequences. It allows you to rapidly simulate chromosome 
scale sequences with essentially correct genetic linkage.

## Installation
### Stable Release (recommended) 
You can download the lastest stable release of `scrm` from [scrm's homepage][1]. 
Instructions on building the binary from the source packages are available in the [wiki][3].

### Stable/Development Version From GitHub

Version             | Branch | Build Status
------------------- | ------ | -----------------
Stable Version      | stable | [![Build Status](https://travis-ci.org/scrm/scrm.png?branch=stable)](https://travis-ci.org/scrm/scrm)
Development Version | master | [![Build Status](https://travis-ci.org/scrm/scrm.png?branch=master)](https://travis-ci.org/scrm/scrm)

You can also install `scrm` directly from the git repository. Here, you need to install `autoconf` first:  

On Debian/Ubuntu based systems:
```bash
apt-get install build-essential autoconf autoconf-archive doxygen libcppunit-dev
```

On Mac OS:
```bash
port install automake autoconf autoconf-archive doxygen cppunit 
```

Afterwards you can build the binary using 
```bash
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
