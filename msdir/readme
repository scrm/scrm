ms

The full documentation on how to use ms is in msdoc.pdf.  Here is 
a brief summary:

The files in this directory are for a program, ms, which generates 
random independent samples according to a simple Wright-Fisher neutral model.
If it is invoked with a minimum of options, it produces samples under 
a panmictic, equilibrium model without recombination.  By specifying 
various options on the command line the model can include recombination,
island-model type structure, gene conversion and simple population
size changes in the past.  

The essential files from this directory are archived in the ms.tar file.
So one need only download the tar file from which the files can be 
extracted by typing:

tar xvf ms.tar 

To compile type:  gcc -O3 -o ms ms.c streec.c rand1.c -lm
  or:             gcc -O3 -o ms ms.c streec.c rand2.c -lm
  or:             gcc -O3 -o ms ms.c streec.c rand3.c -lm

(depending on which random number generator you want.)

Example usage (simplest case):

ms 5  2 -t 6.0 >ms.out

which produces two samples, each of size 5, assuming theta = 6.0 .
Here is the resulting output (i.e. the contents of ms.out ) :

ms 5 2 -t 6.0 
3579 27011 59243

//
segsites: 7
positions: 0.1516 0.2276 0.4854 0.5467 0.7896 0.8501 0.8636 
0100010
0000000
1001001
0110110
0000000

//
segsites: 5
positions: 0.0227 0.0269 0.0856 0.2972 0.5252 
00100
00100
00001
01010
10000


