# Build the container

```
docker build --no-cache -t meep .
```

# Run a simple test

```
docker run --rm -it -v ~/Work/Eisvogel/code/Eisvogel/:/home/eisvogel --entrypoint bash meep
mpirun --allow-run-as-root -np 4 python3 test.py
```