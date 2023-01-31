# Extra profiling output from `polystar`

Profiling the Python/C++ interface in `polystar` is challenging since neither a C++ or Python profiler can identify bottlenecks in the interface itself.
A Python profiler can not peer within the binary module and a C++ profiler can not examine the python code.
To help identify performance issues caused by, e.g., memory copies between the Python heap and C++ heap a special macro-bases profiling output is available in `polystar`.
The profiling macros also print the system time in `YYYY-MM-DD HH:MM:SS:mmm` format by default when called, this can be useful for determining durations spanning Python and C++.

## Programmer's guide to `polystar` output macros
A single logging object exists within the `polystar` module which is used to provide output to the standard output through a set of macros representing different log-levels, `info`, `debug`, `verbose`, and `profile`.

For each log level there are two macros defined, `[level]_update(...)` and `[level]_update_if(logical_expression,...)` where `[level]` should be replaced by `info`, `debug`, `verbose`, or `profile`.
The first macro takes any number of arguments and uses them as the arguments for the logging object's print method.
The second macro checks whether its first argument evaluates to `true` and forwards the remaining arguments if so.

The macros are always defined but, other than the `info` macros, their definitions are empty unless a preprocessor variable matching their level is defined, i.e., `DEBUG`, `VERBOSE`, `PROFILE`; except that `VERBOSE` implies `DEBUG` as well.
When a macro is left as an empty definition its contents are elided away at compile-time, and the compiled executable does not execute any instructions associated with the macro.
A final macro `debug_exec(...)` exists purely to elide expensive operations in order to support more-complex debugging functionality which can be left in place once an operation is verified.

The log level can be controlled by the CMake configuration variable `POLYSTAR_LOGLEVEL`, e.g., `cmake .. -DPOLYSTAR_LOGLEVEL=VERBOSE` will define the appropriate preprocessor variable at compile-time.
As profiling is potentially orthogonal to logging it is controlled by a separate CMake variable, `POLYSTAR_PROFILING`, and can be enabled via, e.g., `cmake .. -DPOLYSTAR_PROFILING=On`

## Profiling output from the binary Python module
At present creating the `_polystar` Python module with an elevated log level or profiling output requires *either* modifying `setup.py` to include the appropriate CMake variables (by adding to the `cmake_args` list) before running, e.g., `python setup.py develop`*or* configure CMake and build the module 'by hand', e.g., `cmake .. -DPOLYSTAR_PROFILING=ON; cmake --build . --target _polystar`, and then move the resulting module into the appropriate* location (*it's left as an exercise to determine what 'appropriate' means for your system).

