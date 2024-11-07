# Privacy Preserving Record Linkage

This is the implementation of the Privacy Preserving Record Linkage method presented in [Accelerating Privacy-Preserving Medical Record Linkage: A Three-Party MPC Approach](https://arxiv.org/abs/2410.21605)

## Installation

No installation is required.

Make sure to **clone** the repository using the "--recurse-submodules" or "--recurse" flag to initialise the submodules as well.

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
./record_linkage Role <helper_port> <helper_ip> <proxy1_port> <proxy1_ip> <record_count> <parallelism>
./record_linkage Role <helper_port> <helper_ip> <proxy1_port> <proxy1_ip> <record_count> <parallelism>
./record_linkage Role <helper_port> <helper_ip> <proxy1_port> <proxy1_ip> <record_count> <parallelism>
```

`Role`: Specifies the role of the entity running the program:
* `proxy1` — First proxy server.
* `proxy2` — Second proxy server.
* `helper` — Helper server.

`helper_port`:
Port number for the helper server.

`helper_ip`:
IP address of the helper server.

`proxy1_port`:
Port number for the first proxy server.

`proxy1_ip`:
IP address of the first proxy server.

`record_count`:
Number of records in the second database that will be linked or compared against.

`parallelism`:
Controls how many records from X (first dataset) are processed at a time to ensure efficient memory use, as storing all comparison data for large datasets in memory simultaneously may be impractical. 
### Local Test Example

For local testing, you can simulate record linkage between the datasets given:
```bash
# Start the helper
./record_linkage helper 7777 127.0.0.1 8888 127.0.0.1 100 1

# Start proxy1 (with 1 record in first database)
./record_linkage proxy1 7777 127.0.0.1 8888 127.0.0.1 1 1

# Start proxy2 (with 100 records in second database)
./record_linkage proxy2 7777 127.0.0.1 8888 127.0.0.1 100 1
 ```

This setup demonstrates how to run the protocol locally, linking a single record from one dataset against 100 records from another. You can adjust the number of records.

### Program Output

The program outputs matching results for records in the first dataset as follows:

```
Index   isMatch?   Matching
0       1          0
Total Time 0.0637471
```

- **Index**: The position of the record in the first dataset.
- **isMatch?**: Indicates if a match was found:
    - `1` if the similarity score exceeds the threshold (currently 0.8).
    - `0` if no match exceeds the threshold.
- **Matching**: The index of the matched record in the second dataset (\(r_2\)). If no match is found, this field may be empty or set to a default value.

### Viewing Exact Similarity Scores

To view the exact similarity scores for matches, modify the `GetMatches` function in `record_linkage.h`. Look for the commented lines in the function that guide where to enable score output.


## License

[MIT](https://choosealicense.com/licenses/mit/)
