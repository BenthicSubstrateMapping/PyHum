# shell script to compile PyHum cython code
# Daniel Buscombe, June 2014
echo "Compiling file reading module"
cython pyread.pyx
gcc -c -fPIC -I/usr/include/python2.7/ pyread.c
gcc -shared pyread.o -o pyread.so

echo "Compiling wavelet computations module"
cython cwt.pyx
gcc -c -fPIC -I/usr/include/python2.7/ cwt.c
gcc -shared cwt.o -o cwt.so

echo "Compiling nan infilling module"
cython replace_nans.pyx
gcc -c -fPIC -I/usr/include/python2.7/ replace_nans.c
gcc -shared replace_nans.o -o replace_nans.so

echo "Compiling spectral noise module"
cython spec_noise.pyx
gcc -c -fPIC -I/usr/include/python2.7/ spec_noise.c
gcc -shared spec_noise.o -o spec_noise.so

echo "Compiling the phase preserving dynamic range compression module"
cython ppdrc.pyx
gcc -c -fPIC -I/usr/include/python2.7/ ppdrc.c
gcc -shared ppdrc.o -o ppdrc.so
