# chess-data-processing

## Compile the cpython code

https://www.msys2.org/

Get msys2 and install MinGW 64 bit

___

Compile with:
```
gcc -c -march=native -fopenmp -pipe -Wall  -fPIC -O3 -lm hkl.c
gcc -shared -Wall -march=native -fopenmp -pipe -O3  -lm hkl.o -o libhkl.dll
```
___

Use conda to install python 3.6 and activate with:

```conda activate python3_6```

___

Now the compiled .dll should work and multithreading should work
