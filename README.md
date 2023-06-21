# Getting started

## Compile the cpython code

Get msys2 and install MinGW 64 bit

https://www.msys2.org/

Open MinGW 64 bit and navigate to the directory with the file "hkl.c", which should be created automatically upon
running the code:
```commandline
cd C:/Users/User/chess_hkl_scripts
```

Compile the files "hkl.o" and "libhkl.dll" with:
```
gcc -c -march=native -fopenmp -pipe -Wall  -fPIC -O3 -lm hkl.c
gcc -shared -Wall -march=native -fopenmp -pipe -O3  -lm hkl.o -o libhkl.dll
```

---

## Using the stacking code

Create a conda environment with an installation of Python 3.6:

```conda create -n stackindexenv python=3.6```

Activate it:

```conda activate stackindexenv```

Now the compiled .dll should be functional and multithreading should work.
