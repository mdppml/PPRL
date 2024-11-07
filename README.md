# Privacy Preserving Record Linkage


## Links to the Related Papers


## Installation

No installation is required.

Make sure to clone the repository using the "--recurse-submodules" or "--recurse" flag to initialise the submodules as well.

If you already have a local version of this repository without submodules, use the command "git submodule update --init --recursive" to initialise the submodules.

## Building

After cloning the repo into directory `PPRL`, you can build the library `PPRL` by executing the following commands.

```bash
mkdir build
cd build
```

```bash
cmake -S ../ -DCMAKE_BUILD_TYPE=Release
```

```bash
make
```

After the build completes, the output binaries can be found in `PPRL/build/` directory.

## Usage

```bash
./record_linkage Role <port of helper> <ip of helper> <port of proxy1> <port of proxy1>
./record_linkage Role <port of helper> <ip of helper> <port of proxy1> <port of proxy1>
./record_linkage Role <port of helper> <ip of helper> <port of proxy1> <port of proxy1>
```

- input = #input parties,#samples of the first input party,#samples of the second input party,...,#samples of the last input party
- delta = delta is a number that specifies how many selections are made after shuffling

```bash
./record_linkage helper 7777 127.0.0.1 8888 127.0.0.1 100 1
./record_linkage proxy1 7777 127.0.0.1 8888 127.0.0.1 1 1
./record_linkage proxy2 7777 127.0.0.1 8888 127.0.0.1 100 1
 ```

## License

[MIT](https://choosealicense.com/licenses/mit/)
