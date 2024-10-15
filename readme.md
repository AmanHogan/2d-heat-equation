# 2D Heat Equation Solver

## Summary

This project contains multiple implementations of solvers for the 2D heat equation using different computational techniques, including parallelization with OpenMP. The project includes three C programs, which are compiled and executed sequentially via a provided Bash script. It is apart of the UTA Parallel Proccesing Course

### Solvers:
1. **heat_solver.c**: Basic solver implementation.
2. **heat_solver_omp.c**: Solver using OpenMP for parallel computation.
3. **heat_solver_omp_sch.c**: Enhanced solver with OpenMP and a scheduling approach.

## Prerequisites

- **GCC** (GNU Compiler Collection) with OpenMP support.
- **Linux/Unix environment** or any system that supports Bash scripting and GCC.
  
## Project Structure

```
├── heat_solver.c               # Basic 2D heat equation solver
├── heat_solver_omp.c           # OpenMP parallelized solver
├── heat_solver_omp_sch.c       # OpenMP solver with scheduling optimizations
├── run_heat_solvers.sh         # Bash script to compile and run all solvers
└── solver_output.txt           # Output file generated after running the solvers
```

## Getting Started

### Step 1: Clone the Repository

To get started, clone the repository to your local machine:
```bash
git clone <repository-url>
cd <repository-directory>
```

### Step 2: Make the Bash Script Executable

If the script is not executable by default, you need to give it execute permissions:
```bash
chmod +x run_heat_solvers.sh
```

### Step 3: Run the Script

Run the Bash script to compile and execute all the solvers:
```bash
./run_heat_solvers.sh
```

This script will:
1. Compile each solver (`heat_solver.c`, `heat_solver_omp.c`, `heat_solver_omp_sch.c`).
2. Run the compiled executables in sequence.
3. Output the results of all three runs into a file called `solver_output.txt`.

### Step 4: View the Output

After the script has finished running, you can check the results in the `solver_output.txt` file:
```bash
cat solver_output.txt
```

## Example Output

The generated `solver_output.txt` will contain the outputs from each solver:
```
Output from heat_solver:
[Results from heat_solver...]

Output from heat_solver_omp:
[Results from heat_solver_omp...]

Output from heat_solver_omp_sch:
[Results from heat_solver_omp_sch...]
```

## Contributions

- **Aman Hogan**
- University of Texas at Arlington
---