# This is README file

This project aims to build a simulator to test various database scheduling strategies by providing
framework that handles underlying DAG initial processing.

Below is the example of how to run the program:

- First, generate `dat` files for the simulator
ruby file_formatter.rb -e data/default/edges.csv -p data/default/prednodes.csv -s data/default/sequence.csv
ruby trace_formatter.rb data/default/

- Enter the following directory
cd lbl_scheduler

- Copy the `dat` files
cp ../*.dat .

- Compile program
make sbu_sched

- To run leveling scheduler with 8 simulated cores and enabled reading from traces:
./main.x -p 8 -t
