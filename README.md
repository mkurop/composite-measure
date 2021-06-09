# composite-measure
Composite measure of speech qualilty from the book by Philipos Loizou "Speech Enhancement - Theory and Practice"

# Cloning repository
To clone repository for latter C++ or Python installation run

```
git clone --recurse-submodules https://github.com/mkurop/composite-measure.git
```

# Installation as a C++ library
Run the following in the project root directory

```
mkdir build
cd build
cmake ..
make 
sudo make install
```

In the CMakeLists.txt of your project use:

```
find_package(composite-measure REQUIRED)
```

# Installation as a Python C++ extension

First create conda environment using

```
conda create --name <env_name> --file requirements.txt
```

Then activate the environment

```
conda activate <env_name>
```

In the folder containing setup.py run

```
pip install .
```

And the package is available in your environment for import.

# Usage C++
Consult the `example/composite_example.cpp` for example usage.

# Usage Python
Consult the `example-python/composite-example.py` for example usage.


