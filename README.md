# cppgfunction
An open source toolbox for calculating and simulating ground heat exchanger thermal response g-functions 

Create a build with cmake:

```
cd /path/to/repo
mkdir build
cd build
ccmake ..   # set settings as needed, shouldn't be any of interest, configure then generate
```

On a Visual Studio generator this will create a solution you can launch into Visual Studio.  
On makefile style generators, this will create a makefile.  From the build directory, just make the project: 

```
make -j 4  # number of processors to use to build
```

Then you can run the test:

```
ctest
```

On my platform this resulted in:
```
Test project /home/edwin/Projects/cpgfunction/cmake-build-debug
    Start 1: RunTest1
1/1 Test #1: RunTest1 .........................   Passed    0.79 sec
```

This ran what used to just be the main.cpp built file.
