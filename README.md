This course is a VC that is assessed via coursework (aka mini projects) and student participation at the Alpen-Adria-University of Klagenfurt. The repository is organized as follows:

## Making AES
This lab gives the chance to “play a bit” with AES encryption rounds and the respective functions.

## Structural Attacks
The lab gives a taster of how differential cryptanalysis works, and a hint as to how structural attacks inspired cryptanalysts to use additional information in order to create
strong, practically feasible attack vectors.

## Task1
It is based on a demo script that performs a DPA style attack on a set of power traces. The demo is based on a differential attack, using a Hamming weight leakage model, and 
correlation as a distinguisher. The job was to investigate several questions:

1. Rewrite the script to attack bits instead of bytes, using a single bit leakage model, and using both correlation as well as a difference of means based distinguisher.
2. Rewrite the script to attack the input of the SubBytes operation, and attack this intermediate with any distinguisher. 
3. Establish the number of needed traces (for the best attack strategy on these traces) to determine the AES key with reasonable enumeration effort.

## Task2
The objective of this request is to apply either a profiling method or leakage detection on a dataset different from that used in Task 1. It has been implemented a profiling method 
using the fixed-key dataset ”ATMEGA_AES_V1”. This dataset consists of synchronized traces where all acquisitions share the same fixed key.
