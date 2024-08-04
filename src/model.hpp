#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <climits>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <list>
#include <random>
#include <errno.h>
#include <math.h>
#include <float.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <chrono>

#include <yaml-cpp/yaml.h>

#include "pcg-cpp/include/pcg_random.hpp"

using namespace std;

// SOME USEFULL GLOBAL CONSTANT (if more than one consider using namespace)
/// External input with frequency below this threshold (expressed in kHz) will be ignored
const double   EPSILON = 0.0001;
/// Parameter setting the precision for the precision in the bisection procedure followed by Network::find_sol_bisection (used in case of oscillatory external poissonian input)
const double   TOLL = 0.0001;
/// Parameter setting the number of warmup epochs which are neglected if optimizer is set to ADAM or SGD_MEAN
const unsigned DO_NOT_TRAIN_IN_FIRST_BATCHES = 5;


// SOME UTILITY FUNCTIONS:
/// Templated function printing a vector
template <class T>
void print_vector(const vector<T> & v);

/// Function managing config file names
string join_and_correct_config(string conf, string dest);


/// Function computing the CDF of the poissonian distribution (with some approximations)
double poissonian_CDF(unsigned x, double mu);

/// Function counting the spike number of spikes of v_t in time windows of half-width half_w with "stride" equal tostep_res.
void count_spikes(vector<unsigned> & v_out, vector<double> & v_t, unsigned lenght, double step_res, double half_w, double dur, double tt0 );


/// Class useful to generate a (pseudo)random double between 0 and 1 using the <b>pcg</b> generator
class RandomGenerator {
private:
    pcg32                                              gen;
    uniform_real_distribution<double>                  dis;
public:
    RandomGenerator(double seed=DBL_MAX) {
        if (seed==DBL_MAX) {
            pcg32                       _gen(pcg_extras::seed_seq_from<random_device>{});
            gen = _gen;
        }
        else {
            pcg32                       _gen(seed);
            gen = _gen;
        }
        uniform_real_distribution<double>              _dis(0,1);
        dis = _dis;
    }
    /// Function returning a (pseudo)random double between 0 and 1
    double getRandomUniform() {
        return  dis(gen);
    }
};


/// Class handling the neuron structure
class Neuron {
public:
    /// State vector
    vector<double>                  x;
    /// Dimenssion of the state vector #x
    unsigned                        dim;
    /// time of last spike (initialized to -1 - #SubNetwork::t_ref) [ms]
    double                          t_last_spike;
    /// First outwards neighbors.
    /// A dictionary is implemetnted with format: {target SubNetwork : list of neurons}
    map<unsigned, vector<unsigned>> neighbors;
    /// Excitatory inputs [ms].
    /// A dictionary is implemetnted with format: {source SubNetwork : list of input times}
    map<unsigned, vector<double>>   input_t_ex;
    /// Inhibitory inputs [ms].
    /// A dictionary is implemetnted with format: {source SubNetwork : list of input times}
    map<unsigned, vector<double>>   input_t_in;

    /// First outwards neighbors weight.
    /// A dictionary is implemetnted with format: {target SubNetwork index : list of weights (corresponding to the neuron-index stored in neighbors)}
    map<unsigned, vector<double>>   neighbors_out_weights;
    /// Auxiliary variables for momentum and ADAM optimezers.
    /// A dictionary is implemetnted with format: {target SubNetwork index : list of weights (corresponding to the neuron-index stored in neighbors)}
    map<unsigned, vector<double>>   weight_aux1; // neighbors_out_delta_weights;
    map<unsigned, vector<double>>   weight_aux2;

    /// Cumulative update weight (used only in case of training with quantized weights: alpha_quantized_weights>0)
    /// A dictionary is implemetnted with format: {target SubNetwork index : list of cumulative weights (corresponding to the neuron-index stored in neighbors)}
    map<unsigned, vector<double>>   weight_quant_cumul;

    /// Excitatory weights [a.u.].
    /// A dictionary is implemetnted with format: {source SubNetwork : list of input weights}
    map<unsigned, vector<double>>   input_w_ex;
    /// Inhibitory weights [a.u.].
    /// A dictionary is implemetnted with format: {source SubNetwork : list of input weights}
    map<unsigned, vector<double>>   input_w_in;

    /// Dictionary containing the index of the next relevant spike in the vector #input_t_ex for each excitatory Subnetwork
    map<unsigned, unsigned>         next_sp_ex_index;
    /// Dictionary containing the index of the next relevant spike in the vector #input_t_in for each inhibitory Subnetwork
    map<unsigned, unsigned>         next_sp_in_index;

    /// External input weight:
    /// for each neuron in the subnetwork, it is given by #SubNetwork::weights_ex['ext'] plus a random value uniformly extracted in [-#SubNetwork::dev_ext_weight, #SubNetwork::dev_ext_weight] [nS]
    double                          ext_weight;

    /// Vector of the neuron spike times [ms]
    vector<double>                  t_spikes;

    // Specific external current bool
    bool                            specific_I_e_bool;
    // Specific external current vector
    vector<double>                  specific_I_e;

    // Offset in external current (acts in both specific or general I_e)
    double                          offset_I_e;
    // Auxiliary variables for momentum and ADAM optimezers.
    double                          offset_aux1;  // delta_offset_I_e
    double                          offset_aux2;

    /// Vector containing the desired number of spikes for each time interval, if provided
    vector<double>                  desidered_output;
    /// Vector containing the actual number of spikes for each time interval, during training
    vector<unsigned>                R_actual;
    /// Vector containing the generalized error, during training
    vector<double>                  error;
    /// Training hyperparameter regulating the modulation of the error in case of low desired output
    double                          optimal_c_desid_zero;
    /// Vector containing the spike times to be replicated in case the neuron is set in parrot mode
    vector<double>                  parrot_spikes;

    /// Class constructor
    Neuron(unsigned _dim);

    /// Method printing a list of the attributes and their value
    void info();


};


/// Class handling the SubNetwork structure
class SubNetwork {
private:
    /// RandomGenerator object useful to generate random numbers from uniform a uniform distribution in [0,1]
    RandomGenerator     g;
public:
    /// Name of the SubNetwork (must be unique in the Network)
    string       name;
    /*!  Model of neurons (Name-conventions of <a href="https://nest-simulator.readthedocs.io/en/v3.3/contents.html">NEST-simulator</a> adopted).
    *   Supported models:
    *     - iaf_cond_alpha   (id_model 0)
    *     - iaf_curr_alpha   (id_model 1)
    *     - aeif_cond_exp    (id_model 2)
    *     - aeif_curr_exp    (id_model 3)
    *     - aqif_cond_exp    (id_model 4)
    *     - aqif_curr_exp    (id_model 5)
    *     - iaf_cond_exp     (id_model 6)
    *     - iaf_curr_exp     (id_model 7)
    *     - aqif2_cond_exp   (id_model 8)
    *     - parrot_neuron    (id_model 9)
    */
    string       neuron_model;
    /// Identificative number for the neuron model
    unsigned     id_model;
    /// Identificative number for the subnetwork (associated to the name)
    unsigned     id_subnet;
    /// Number of neurons in the SubNetwork
    unsigned     N;
    /// Membrain capacity [pF]
    double       C_m;
    /// Resting potential [mV]
    double       E_L;
    /// Excitatory reversal potential [mV]
    double       E_ex;
    /// Inhibitory reversal potential [mV]
    double       E_in;
    /// Reset potential [mV]
    double       V_res;
    /// Threshold potential [mV]
    double       V_th;
    /// Refractory time [ms]
    double       t_ref;
    /// External injected current [pA]
    double       I_e;
    /// Amplitude of the oscillatory part of the external injected current [pA]
    double       osc_amp;
    /// Angular frequency of the oscillatory part of the external injected current [kHz]
    double       osc_omega;
    /// Amplitude of the oscillatory part of the external input rate [kHz]
    double       osc_amp_poiss;
    /// Angular frequency of the oscillatory part of the external input rate [kHz]
    double       osc_omega_poiss;
    /// Deviation of #weights_ex['ext'] from its central value [nS]
    double       dev_ext_weight;
    /// Input rate from external source [kHz] (here simulating cortical input)
    double       ext_in_rate;
    /// characteristic time of excitatory synaptic inputs [ms]
    double       tau_syn_ex;
    /// characteristic time of inhibitory synaptic inputs [ms]
    double       tau_syn_in;
    /// file containg the eventual values of the "specific" external current for each neuron; "none" is interpreted as no specific current.
    string       specific_I_e_file;

    /// file containg the eventual values of the offset external current for each neuron; "none" is interpreted as no offset current.
    string       offset_I_e_file;

    /// Weights of excitatory input connections.
    /// A dictionary is implemetnted with format: {source SubNetwork : weight [nS]}.
    map<unsigned, double>     weights_ex;
    /// Weights of inhibitory input connections.
    /// A dictionary is implemetnted with format: {source SubNetwork : weight [NS]}; all weights are positive.
    map<unsigned, double>     weights_in;

    /// Synaptic delays from the subnet to the other subnets.
    /// A dictionary is implemetnted with format: {target SubNetwork : delay [ms]}; all weights are positive.
    map<unsigned, double>     delays_out;
    /// Connection probabilities of the SubNetwork with the other SubNetworks.
    /// A dictionary is implemetnted with format: {target SubNetwork : probability}.
    map<unsigned, double>     probabilities_out;
    /// Effect of spike on target population
    /// A dictionary is implemetnted with format: {target SubNetwork index : bool}. If true the effect of a spike is on the excitatory variable; otherwise it is on the inhibitory.
    map<unsigned, bool>       excit_out;
    /// Characterization of the effect of spike on target population
    /// A dictionary is implemetnted with format: {target SubNetwork : bool}. If true the effect of a spike on the target population is reverted (e.g. an excitatory input on a iaf_cond_exp neuron decreases g_ex)
    map<unsigned, bool>       reverse_effect;

    map<unsigned, double>     alpha_quantized_weights;

    // MODEL SPECIFIC ATTRIBUTES
    /// Subthreshold adaptation [nS] (only adaptive models)
    double       a_adaptive;
    /// Spike-triggered adaptation: step_height of adaptation variable after spike emission [pA] (only adaptive models)
    double       b_adaptive;
    /// Characteristic decay time of the adaptation variable (only adaptive models)
    double       tau_w_adaptive;
    /// Spike detection threshold (only aeif_cond_exp and aqif_cond_exp models)
    double       V_peak;

    /// Membrain leakage conductance [nS] (only aeif_cond_exp and iaf_cond_alpha)
    double       g_L;

    /// Slope factor of exponential rise (only aeif_cond_exp)
    double       delta_T_aeif_cond_exp;

    /// k parameter of Izhikevich adaptive model [pA/mV<SUP>2</SUP>]
    // C_m dV/dt = k(V-V_th)(V-E_L) + input currents - w + I_e; tau_w dw/dt = a (V-E_L) - w; w -> w+b  (only aqif_cond_exp and aqif2_cond_exp)
    double       k_aqif_cond_exp;

    /// V_b parameter of Izhikevich adaptive fast spiking interneurons model (only aqif2_cond_exp)
    double       V_b_aqif2_cond_exp;

    /// if true postsynaptic currents or conductances are alpha-shaped
    bool         alpha_PSS;

    /// if true the model is condunctance based; if false the model is current based
    bool         cond_based;

    /// if true the model contain an adaptation variable in the fourth dimension of the state vextor.
    bool         adap_4;

    /// Vector of #N #Neuron-type objects
    vector<Neuron>          pop;

    /// Vector of neurons whose state you want to save
    vector<unsigned>        to_save;

    /// if true the incoming weights are specific for each neuron in the population
    /// a dictionary is implemented with format { source SubNetwork index : bool }
    map<unsigned, bool>     specific_in_weight;
    /// if true the outgoing weights are specific for each neuron in the population
    /// a dictionary is implemented with format { target SubNetwork index : bool }
    map<unsigned, bool>     specific_out_weight;
    /// for each target population with specific out weights, the dictionary contains the configuration file
    /// a dictionary is implemented with format { target SubNetwork index : file }
    map<unsigned, string>   specific_out_weight_files;
    /// if true the outgoing weights are trained
    /// a dictionary is implemented with format { target SubNetwork index : bool }
    map<unsigned, bool>     train_out_weight;

    /// path to file containg the desidered output rates
    string                  desidered_output_path;

    /// regulates the order for the backpropagation stage of the training process.
    /// The output layer has train_order==0; set to -1 if none of the input weights are to be trained.
    int                     train_order;

    /// minimum value of output trainable weights
    double                  w_min;
    /// maximum value of output trainable weights
    double                  w_max;

    /// if true the offset current parameter is trained
    bool                    train_offset_I_e;
    /// minimum value of output trainable offset current
    double                  offset_I_e_min;
    /// maximum value of output trainable offset current
    double                  offset_I_e_max;
    /// equilibrium value of the ofsfset current, if ldecay_offset>0
    double                  offset_I_e_ref;
    /// decay factor of the offset current
    double                  ldecay_offset;

    /// Training hyperparameter regulating the modulation of the weight-decay mechanism: maximum decay factor
    double                  decay_factor_max;
    /// Training hyperparameter regulating the modulation of the weight-decay mechanism: steepnes of the modulation
    double                  decay_factor_steepness;
    /// Training hyperparameter regulating the modulation of the weight-decay mechanism: location of the modulation
    double                  decay_factor_position;

    /// path to file containg the timings of the spikes to be replicated
    string                  parrot_spike_path;

    /// if True, errors have already been computed in this epoch
    bool                    errors_computed;

    /// class constructor
    SubNetwork(string _name, string _neuron_model, int _N, double _C_m, double _E_L,  \
               double _V_res, double _V_th,                \
               double _t_ref, double _I_e, double _osc_amp, double _osc_omega,        \
               double _dev_ext_weight, double _ext_in_rate,                           \
               double _osc_amp_poiss, double _osc_omega_poiss,                        \
               double _tau_syn_ex, double _tau_syn_in, RandomGenerator _g);

    /// Method which save Neuron::t_spikes of each Neuron in #pop
    /// \param out_file output file XXX scrivi come viene stampato
    void save_t_spikes(string out_file);

    /// Method printing a list of the attributes and their value
    void info(map<unsigned, string> &subnet_index_name);
};

/// Class handling the network structure
class Network {
private:
    /// End time of simulation [ms]
    double      t_end;
    /// Current time during simulation (starts at #t=0 and ends at #t_end)
    double      t = 0;
    /// Number of calls to evolve method with argument #t_end / #n_step
    unsigned    n_step = 0;
    /// Time resolution of the simulation [ms]
    double      dt;
    /// Input yaml-file with Network composition and features.
    string      subnets_config_yaml;
    /// Input yaml-file with connection weights (positive weights are excitatory, negative weights are inhibitory)
    string      weights_config_yaml;
    /// Input yaml-file with connectivity probabilities between subnetworks and corresponding delays.
    /// Note that each delay must immediately follow the related probability
    string      connections_config_yaml;
    /// Input yaml-file with list of neurons whose state you want to save at each step.
    /// You can leave this file empty if you don't want to save any nuron state.
    string      to_save_config_yaml;
    /// Output directory of the simulation
    string      out_dir;

    /// Dictionary with format {#SubNetwork::name : related index in Network::subnets (equal to SubNetwork::id_subnet)}
    map<string, unsigned>   subnet_name_index;
    /// Dictionary with format {index in Network::subnets (equal to SubNetwork::id_subnet) : #SubNetwork::name}
    map<unsigned, string>   subnet_index_name;
    /// Vector of #SubNetwork-type objects
    vector<SubNetwork>      subnets;
    /// RandomGenerator object useful to generate random numbers from uniform uniform distribution in [0,1]
    RandomGenerator         g;

    /// External input mode:
    /// - <i>0</i> (base mode): each neuron receives an indipendent poisson signal with mean frequency = #SubNetwork::ext_in_rate and possibly with the osccillatory component
    /// - <i>2</i> (with_correlation mode): implementation of method A (ask for details, not compatible with oscillatory input)
    unsigned                input_mode;

    // EXTERNAL MODE SPECIFIC ATTRIBUTES
    /// Parameter regulating the input correlation of the striatum populations (only in input_mode=2)
    double                  rho_corr;
    /// Vector containing the subnets indices corresponding to the population with correlatated inputs (only in input_mode=2)
    vector<unsigned>        corr_pops;
    /// Map containg the last time generated by the exponential distribution for the (partially correlated) external input
    map<unsigned, double>   corr_last_time;

    /// number of time steps for which the external specific current is kept constant
    unsigned                repeat_specific_I_e;

    // TRAINING MEMBERS
    /// Input yaml-file with training variables
    string          training_config_yaml;

    /// Optimization method employed to minimize the loss; Supported methods are:
    /// - SGD      (id=0): momentum_w, momentum_c
    /// - SGD_MEAN (id=1): momentum_w, momentum_c
    /// - ADAM     (id=2): momentum_w, momentum_c, momentum2_w, momentum2_c
    string          optimization_method;
    /// Optimization method id (see optimization_method)
    unsigned        optimization_method_id;

    /// Training variable: current epoch
    unsigned        current_epoch;
    /// Training variable: path to mini-batches durations [ms]
    string          dur_batches_path;
    /// Training variable: vector containing the mini-batches durations [ms]
    vector<double>  dur_batches;
    /// Training variable: temporal amplitude of each time interval considered to train the network ms]
    double          window_res;
    /// Training variable: temporal shift between consecutive time intervals considered considered to train the network [ms]
    double          step_res;
    /// Training variable: maximum absolute value for the weights [pA]
    double          w_cut;
    /// Training variable: regulation the shape of the error function
    unsigned        POT;
    /// Training variable: probability of noisy weight-updates during training
    double          prob_noise;
    /// Training variable: amplitude of noisy weight-updates during training
    double          amp_noise;
    /// Training variable: intensity of weight decay mechanism
    double          l_decay;
    /// Training variable: learning rate for weight updates
    double          l_rate0;
    /// Training variable: learning rate for current updates
    double          l_rate0_curr;
    /// Training variable: momentum hyperparameter for weight-updates (SGD_MEAN and ADAM only)
    double          momentum_w;
    /// Training variable: momentum hyperparameter for current-updates (SGD_MEAN and ADAM only)
    double          momentum_c;
    /// Training variable: momentum hyperparameter for weight-updates (ADAM only)
    double          momentum2_w;
    /// Training variable: momentum hyperparameter for current-updates (ADAM only)
    double          momentum2_c;
    /// Training variable: slope of the FI curve in the F>0 regime
    double          dF_dI_app;
    /// Training variable: leaky multiplicative factor regulating the slope of the FI curve in the F=0 regime
    double          dF_dI_no_act_fact;
    /// Training variable: modulation factor for the linear component of the error
    double          error_lin_coeff;
    /// Training variable: see code for details
    double          epsilon_zero_act;
    /// Training variable: hyperparamenter regulating the modulation of the error in case of low desired output (see code for details)
    double          a_desid_zero;
    /// Training variable: hyperparamenter regulating the modulation of the error in case of low desired output (see code for details)
    double          c_desid_zero;
    /// Training variable: hyperparameter allowing for the modulation of the learning rate in the different layers
    vector<double>  correct_l_rates;
    /// Training variable: see code for details
    vector<int>     ord_target_subnets;

    /// Training variable: hyperparameter regulating the fine tuning of quantized weights
    int             n_quant;
    /// Training variable: hyperparameter regulating the fine tuning of quantized weights
    double          quant_delta_weight_pos;
    /// Training variable: hyperparameter regulating the fine tuning of quantized weights
    double          quant_delta_weight_neg;

public:

    /// Dictionary containing the relation between the supported Neuron models and the dimension of its state vector
    map<string, unsigned>   supported_models = { {"iaf_cond_alpha", 5}, {"aeif_cond_exp", 4}, {"aqif_cond_exp", 4}, {"aqif2_cond_exp", 4}, {"iaf_cond_exp", 3}, \
                                                 {"iaf_curr_alpha", 5}, {"aeif_curr_exp", 4}, {"aqif_curr_exp", 4},                        {"iaf_curr_exp", 3}, \
                                                 {"parrot_neuron",  3}    };
    /// Dictionary with format { SubNetwork::neuron_model : SubNetwork::id_model }
    map<string, unsigned>   subnet_model_id = { {"iaf_cond_alpha", 0}, {"aeif_cond_exp", 2}, {"aqif_cond_exp", 4}, {"aqif2_cond_exp", 8}, {"iaf_cond_exp", 6}, \
                                                {"iaf_curr_alpha", 1}, {"aeif_curr_exp", 3}, {"aqif_curr_exp", 5},                        {"iaf_curr_exp", 7}, \
                                                {"parrot_neuron",  9}    };
    /// Dictionary with format { SubNetwork::neuron_model : bool } bool is true if the model includes an adaptation variable (in the 4th dimension of the state attribute)
    map<string, bool>       model_is_adaptive = { {"iaf_cond_alpha", false}, {"aeif_cond_exp", true}, {"aqif_cond_exp", true}, {"aqif2_cond_exp", true}, {"iaf_cond_exp", false}, \
                                                  {"iaf_curr_alpha", false}, {"aeif_curr_exp", true}, {"aqif_curr_exp", true},                           {"iaf_curr_exp", false}, \
                                                  {"parrot_neuron",  false}   };

    /// Class constructor  // XXX handle errors in createSubnet
    Network(double _t_end, double _dt, unsigned _input_mode, string _subnets_config_yaml, string _weights_config_yaml, \
        string _connections_config_yaml, string _to_save_config_yaml, string _training_config_yaml,  \
        string _out_dir, RandomGenerator _g, string _input_mode_config, unsigned _n_step, unsigned _repeat_specific_I_e, bool evolve_only);

    /// Method initializing the #subnets vector using configurations files:
    /// - #subnets_config_yaml for the neurons' features;
    /// - #weights_config_yaml for the synaptic weights;
    /// - #connections_config_yaml for the connections features: probabilities and delays.
    void createSubnets(bool evolve_only);
    /// Method initializing the vector SubNetwork::pop of each SubNetwork in #subnets
    /// Neurons are initialized with Neuron::x[0] = SubNetwork::E_L of the belonging SubNetwork
    void createPops();
    /// Method initializing the specific_I_e_bool and specific_I_e attribute of all Neurons in all SubNetworks
    void load_specific_currents_and_offsets();
    /// Method initializing the weights attribute of all Neurons in all SubNetworks
    void load_specific_weights();

    /// Methos initializing the variables for the training
    void set_up_training();

    /// Method training the trainable weights
    void train(double dur);

    /// Method training the trainable weights
    unsigned compute_gen_error(int t_order, unsigned lenght, double dur);

    /// Method evolving the network for time _T
    void evolve(double _T);

    /// Method updating external input according to #input_mode
    void externalInputUpdate();

    /// Function whose solutions==0 needs to be find in case of oscillatory external input rate
    double input_func(double y_, double r0_, double A_, double omega_, double t0_, double t_);

    /// Function determining next spike time in case of oscillatory external input rate
    double find_sol_bisection (double y_, double r0_, double A_, double omega_, double t0_);

    /// Function parsing "specific" exteral currents for a subnetwork
    void parse_external_current(vector<Neuron>& vect, string file, int lenght);

    /// Function parsing "specific" weights from a subnetwork to another (target_p)
    void parse_specific_weight(vector<Neuron>& vect, unsigned target_p, string file);

    /// Function saving the current "specific" weights from a subnetwork to another
    void save_specific_weights();

    /// Function saving the current external I_e "offset"
    void save_offset_currents();

    /// Method freeing memory from past spikes and performing a control operation over the state vector of each neuron
    void free_past();

    /// Method printing:
    /// - main characteristics of the Network composition
    /// - the simulation control variables
    void info();

    /// Method printing main parameters of training
    void info_training();

    /// Auxiliary Method for the fine tuning of quantized weights
    double get_quant_weight_delta(double w_quant, double w_cum);

    /// Method for parsing spike times of parrot neurons
    void parse_parrot_spikes(vector<double>& spikes, string file);
};


#endif
