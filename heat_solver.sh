#!/bin/bash

# Author: Aman Hogan
# Function: Runs three scripts solve 2D heat equation using
# various techniques

SOLVER1="heat_solver.c"
SOLVER2="heat_solver_omp.c"
SOLVER3="heat_solver_omp_sch.c"

EXE1="heat_solver"
EXE2="heat_solver_omp"
EXE3="heat_solver_omp_sch"

OUTPUT_FILE="solver_output.txt"

# Clean up previous logs and binaries
rm -f $EXE1 $EXE2 $EXE3 $OUTPUT_FILE

echo "Compiling $SOLVER1..."
gcc -fopenmp -o $EXE1 $SOLVER1 -lm
if [ $? -ne 0 ]; then
    echo "Compilation of $SOLVER1 failed."
    exit 1
fi

echo "Compiling $SOLVER2..."
gcc -fopenmp -o $EXE2 $SOLVER2 -lm
if [ $? -ne 0 ]; then
    echo "Compilation of $SOLVER2 failed."
    exit 1
fi

echo "Compiling $SOLVER3..."
gcc -fopenmp -o $EXE3 $SOLVER3 -lm
if [ $? -ne 0 ]; then
    echo "Compilation of $SOLVER3 failed."
    exit 1
fi

echo "Running $EXE1..."
./$EXE1 > $OUTPUT_FILE
echo "$EXE1 finished."

echo "Running $EXE2..."
./$EXE2 >> $OUTPUT_FILE
echo "$EXE2 finished."

echo "Running $EXE3..."
./$EXE3 >> $OUTPUT_FILE
echo "$EXE3 finished."

echo "All programs have been compiled and run successfully."
echo "Check the $OUTPUT_FILE for output."
