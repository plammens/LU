# LU

Object-oriented LU decomposition implementation in C++.

# Setup, testing and execution instructions

> Note: these instructions assume you're on a Unix-like system, but by installing `cygwin` and some additional hassle you can use these on Windows too.

## Build

For `make` to run correctly, this project needs `GNU Coreutils` (i.e. the `UNIX`-style commands such as `find`, `rm`, etc.) , and `GNU make`. So if you're running Linux or
Mac, you probably won't need to install anything. If you're running Windows, installing [`Cygwin`](https://www.cygwin.com) with the aforementioned packages should be enough.

No external dependencies need to be installed separately.

To build everything (main program and test), run the following command.

```bash
make all  # (or just `make`)
```

This should do two things: build the main executable, and build the tests.

Other useful `make` targets:

- To build the main program separately, run
```bash
make build
```

- To build the test program separately, run
```bash
make build-test
```

- To build the libraries separately, run
```bash
make libs
```

- To clean everything (build directory and output), run
```bash
make clean
```

- To clean the output (i.e. solution files), run
```bash
make clean-out
```

- To clean the build, run
```bash
make clean-build
```



## Execution of main program

To run the main program, just execute the command

```bash
make run
```

Or, if you're a do-it-yourself type (and you have already run `make build` or `make all`), you
can execute the program directly like so (from the project root):

```bash
./build/bin/main.x
```

If you need to pass arguments to the program, the latter is probably the easiest way:

```bash
./build/bin/main.x MAT00.DAT MAT01.DAT
```

You can also do:

```bash
make run ARGS="MAT00.DAT MAT01.DAT"
```





## Testing the program

This project comes with a few bundled unit tests, made with [`doctest`](https://github.com/onqtam/doctest). To run the tests:

```bash
make test
```

(this is equivalent to

```bash
make build-test
./build/bin/test.x
```

).

If the tests succeed and you want to inspect what's actually going on behind the scenes,
run

```bash
make test ARGS="--success"
```

so you'll see in detail every single test assertion. (In fact, you can pass any arguments you
like to the test program by giving the `ARGS` variable a value &mdash; e.g. `ARGS="--foo --bar"`;
refer to [doctest's docs](https://github.com/onqtam/doctest/blob/master/doc/markdown/commandline.md) for all the supported command-line arguments).