# Compile

```
g++ `pkg-config --cflags meep` test.cpp -o test `pkg-config --libs meep`
g++ `pkg-config --cflags meep` `pkg-config --cflags hdf5` test.cpp -o test `pkg-config --libs meep` `pkg-config --libs hdf5` -lhdf5_cpp
```

# Run

```
mpirun --allow-run-as-root -np 4 ./test
```