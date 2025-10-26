## Black Holes on Graphs

This repository provides the implementation of algorithms for **black hole detection queries on directed graphs**.
 It includes our proposed methods (`DetectBH` and `DetectBH-P`) as well as a baseline approach (`BaseLine`) for comparison.

## Code Structure

The main implementations are located in the `src` directory, with each algorithm in its own files:

- `Graph.h/cpp`  - Index and algorithm implementation
- `tools.hpp` - Utility functions
- `type.h` - Classes and types assisting code implementation
- `main.cpp` - Main program

## Compile the codes

When you already download the codes, run the following commands to compile our codes.

```
cd src
g++ -std=c++17 -O3 -fopenmp main.cpp Graph.cpp -o main
```

After running the codes, there will be executable files called `main` in `./src` directory.

## Datasets

- **Default datasets** (e.g., `Wiki-Vote`) are included in the package.
- **Additional datasets** can be downloaded from [SNAP](https://snap.stanford.edu/data/index.html).
- Conversion scripts:
  - The `readFileTrans` function in `tools.cpp` can be used to convert datasets into the format required by this index, and the resulting files are stored in the `./dataset/transform` directory.

## Run the procedure

### Example Commands

**Perform a query for black holes of size k :**

```
./main -d Wiki-Vote -l 2 -r 10 -m 100
```

### Command Line Options

| Option | Required | Description                                                  |
| ------ | -------- | ------------------------------------------------------------ |
| `-d`   | Yes      | Dataset name (e.g., `Wiki-Vote`)                             |
| `-l`   | Yes      | Query the left interval of black holes of size k             |
| `-r`   | Yes      | Query the right interval of black holes of size k            |
| `-m`   | No       | The index retains the largest minimum black hole size (greater than or equal to r) |

