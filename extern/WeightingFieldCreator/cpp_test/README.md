# Compile

```
g++ `pkg-config --cflags meep` test.cpp -o test `pkg-config --libs meep`
```

# Run

```
mpirun --allow-run-as-root -np 4 ./test
```