# Energy-Efficient-Touch-Localization

XXX links

### Overview
This repo includes code for:
1. simulating Spiking Neural Network (SNN);
2. training SNNs with the learning algorithm presented in the paper;
3. reproducing the results of the experiments described in the paper.

The code for the simulation and training of SNN has been developed in C++ and is structured in such a way that general networks of supported point neurons can be implemented and trained with the adopted training algorithm.
In combination with this code, a Python 3 module has been implemented in order to handle the training process and import and preliminarily analyze the results of the simulations. The related documentation and tutorials are published on readthedocs.org (XXX).

To reproduce the experimental results from the paper, please refer to the src/build/run_indentation_example.ipynb Notebook. **Note:** Running this notebook requires compiling the C++ library for SNN simulation and training. For compilation instructions, see the setup section in the documentation (XXX).
