.. Energy-Efficient-Touch-Localization documentation master file, created by
   sphinx-quickstart on Sun Aug  4 11:10:09 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Energy-Efficient-Touch-Localization's documentation!
===============================================================

This is the documentation for the code published on GitHub at `this <https://github.com/Neuro-Robotic-Touch-Laboratory/Energy-Efficient-Touch-Localization>`_ page.
The code allows for the simulation and training of Spiking Neural Networks. It has been developed in C++ and is structured in such a way that general networks of supported point neurons can be implemented.
Supported neuron types are (name-conventions of `NEST-simulator <https://nest-simulator.readthedocs.io/en/v3.3/contents.html>`_ are adopted, when possible):

- `iaf_cond_alpha  <https://nest-simulator.readthedocs.io/en/v3.3/models/iaf_cond_alpha.html?highlight=iaf_cond_alpha>`_
- iaf_curr_alpha (current based version of iaf_cond_alpha)
- `aeif_cond_exp <https://nest-simulator.readthedocs.io/en/v3.3/models/aeif_cond_exp.html?highlight=aeif_cond_exp>`_
- aeif_curr_exp (current based version of aeif_cond_exp)
- aqif_cond_exp (adaptive quadratic integrate and fire)
- aqif_curr_exp (current based version of aqif_cond_exp)
- `iaf_cond_exp <https://nest-simulator.readthedocs.io/en/v3.3/models/iaf_cond_exp.html?highlight=iaf_cond_exp>`_
- iaf_curr_exp (current based version of iaf_cond_exp)
- parrot_neuron

In combination with this code, a python3 module has been implemented in order to handle the training process and import and preliminarily analyze the results of the simulations (see :ref:`examples`).

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Table of Contents
-----------------
.. toctree::
   :maxdepth: 2

   setup
   examples
   api/index
