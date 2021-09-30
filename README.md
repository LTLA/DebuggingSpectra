# Spectra debugging efforts

Repository containing code to reproduce assorted Spectra bugs.
I'm using R to easily manage data I/O without needing to write all the boilerplate in C++.
Run the following to set up the environment for compilation of the C++ code itself.

```sh
cmake -S . -B build
ln -s build/_deps/spectra-src/include/Spectra/ Spectra
ln -s build/_deps/eigen3-src/Eigen
```
