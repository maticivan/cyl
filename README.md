# A computer-assisted proof of Kuperberg's six-cylinder conjecture

Ivan Matic and Rados Radoicic

For the most current version of the paper, please visit
<https://imomath.com/maticivan>.

## Building and running

The program requires a C++23 compiler.

    make
    ./cyl -i cyl_inputData.txt -o cyl_outputData

This writes `cyl_outputData.tex`, a summary of the certificate, and
`cyl_outputData.txt`, the box-by-box printout of all cases.

## The parallel program is the default

The parallel program uses all available cores and finishes in under 5 minutes
on most modern computers. The non-parallel program may take up to 20 minutes.

If you are reading the source rather than running it, you may prefer the
non-parallel version, which is simpler. To use it, rename the folders before
compiling:

    mv src/cylinders src/cylinders_parallel
    mv src/cylinders_nonParallel src/cylinders
    make

## The large output file

`cyl_outputData.txt` is too large for GitHub. Generate it as above, or
download a copy:

<https://www.dropbox.com/scl/fi/29c747kertknl681f434j/cyl_outputData.txt?rlkey=56ndbd9ybaghwqvqvzxnyum5s&st=usi0x79j&dl=0>
 

## License

MIT. See `LICENSE`.
