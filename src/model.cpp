#include <iostream>
#include "model.hpp"

using namespace std;

// SOME UTILITY FUNCTIONS:
template <class T>
void print_vector(const vector<T> & v) {
    if (v.size() > 0) {
        for (unsigned i=0; i<v.size()-1; i++) {
            cout << v[i] << " ";
        }
        cout << v[v.size()-1] << endl;
    }
    else cout << "empty" << endl;
}

string join_and_correct_config(string conf, string dest) {

    if (conf == "none") {
        return "none";
    }
    else {
        string      to_ret;
        to_ret = conf.substr( conf.find_last_of("/"), conf.length()-conf.find_last_of("/") );
        to_ret = dest + to_ret;

        return to_ret;
    }
}

void Network::parse_parrot_spikes(vector<double>& spikes, string file) {

    ifstream                in_file(file);
    string                  line;
    double                  previous_spike=0, current_spike=0;


    if (in_file.is_open()) {
        while (getline(in_file,line)) {
            current_spike = stod(line);
            spikes.push_back( current_spike );
            if (current_spike-previous_spike < dt ) {
                cerr << "Error: minimum ISI is less than dt in file " << file << endl;
                exit(1);
            }
        }
    }
    else {
        cerr << "Error: unable to open file " << file << endl;
        exit(1);
    }
}


void Network::parse_external_current(vector<Neuron>& vect, string file, int lenght) {

    vector<double>          temp(lenght,0.);
    string                  line;
    int                     start, neuron, i, n=0;
    size_t                  end;
    const char              del_inl='\t';
    ifstream                in_file(file);

    int                     N_neur=vect.size();

    if (in_file.is_open()) {
        while (getline(in_file,line)) {
            start=0;
            // read neuron id
            end = line.find(del_inl, start);
            neuron = stoi(line.substr(start, int(end)-start));
            if (neuron<n) {
                cerr << "Error: inconsistent (too low) neuron number in " << file << endl;
                exit(1);
            }
            if (neuron>N_neur) {
                cerr << "Error: inconsistent (too high) neuron number in " << file << endl;
                exit(1);
            }
            start = int(end)+1;

            for (;n<neuron;n++) {
                vect[n].specific_I_e_bool = false;
            }
            for (i=0; end!=string::npos; i++) {
                end = line.find(del_inl, start);
                temp[i] = stod(line.substr(start, int(end)-start));    // get number from string
                start = int(end)+1;
            }
            // note: ther MUST NOT be "\t" at the and of each row
            if (i!=lenght) {
                cerr << "Error: inconsistent lenght of external-current vector for neuron " << n << " in " << file << "\t(i="<<i<<")" << endl;
                exit(1);
            }
            vect[n].specific_I_e_bool = true;
            vect[n].specific_I_e = temp;
            n++;
        }
        for (;n<N_neur;n++) {
            vect[n].specific_I_e_bool = false;
        }
        if (n>N_neur) {
            cerr << "Error: too many raws in " << file << endl;
            exit(1);
        }
        in_file.close();
    }
    else {
        cerr << "Error: unable to open file " << file << endl;
        exit(1);
    }
}


void Network::parse_specific_weight(vector<Neuron>& vect, unsigned target_p, string file) {

    int                     lenght = vect[0].neighbors[target_p].size(), N_neur=vect.size();
    ifstream                in_file(file);
    vector<double>          temp(lenght ,0.);
    vector<double>          temp_delta(lenght ,0.);
    string                  line;
    int                     start;
    size_t                  end;
    const char              del_inl='\t';
    int                     n=0, i;

    if (in_file.is_open()) {
        while (getline(in_file,line)) {
            if (n==N_neur) {
                cerr << "Error: inconsistent (too high) number of rows (" << n+1 << ") in " << file << endl;
                exit(1);
            }
            start=0;
            end = line.find(del_inl, start);
            if (n!=stoi(line.substr(start, int(end)-start))) {
                // the weights of all source neurons (n) must be initialized (1)
                cerr << "Error: inconsistent neuron number " << stoi(line.substr(start, int(end)-start)) << " (n="<< n <<") in " << file << endl;
                exit(1);
            }
            start = int(end)+1;
            for (i=0; end!=string::npos; i++) {
                end = line.find(del_inl, start);
                temp[i] = stod(line.substr(start, int(end)-start));    // get number from string
                start = int(end)+1;
            }
            // note: ther MUST NOT be "\t" at the and of each row
            if (i != lenght) {
                // for each source neuron (n), the weights to all target neurons (i) must be initialized
                cerr << "Error: inconsistent lenght of specific weights for neuron " << n << " in " << file << endl;
                exit(1);
            }
            vect[n].neighbors_out_weights[target_p] = temp;
            vect[n].weight_aux1[target_p] = temp_delta;
            vect[n].weight_aux2[target_p] = temp_delta;
            n++;
        }
        if (n!=N_neur) {
            // the weights of all source neurons (n) must be initialized (2)
            cerr << "Error: inconsistent number of rows (" << n << ") in " << file << endl;
            exit(1);
        }
        in_file.close();
    }
    else {
        cerr << "Error: unable to open file " << file << endl;
        exit(1);
    }
}


void Network::save_specific_weights() {

    string          file_name;
    ofstream        of;

    // for each pop (source)
    for (auto k=subnets.begin(); k!=subnets.end(); k++) {

        // for each target pop (with train)
        for (auto p_k=k->train_out_weight.begin(); p_k!=k->train_out_weight.end(); p_k++ ) {
            if (p_k->second) {

                file_name = out_dir + "/specific_weights/" + k->name + "_to_" + subnets[p_k->first].name + "_" + to_string(current_epoch) + ".txt";
                of.open( file_name.c_str() );

                for (unsigned i=0; i<k->N; i++) {
                    of << i;
                    for (auto w : (k->pop)[i].neighbors_out_weights[p_k->first] )
                        of << "\t" << w;
                    of << endl;
                }
                of.close();
            }
        }
    } // end for over subnets
}


void Network::save_offset_currents() {

    string          file_name;
    ofstream        of;

    // for each pop
    for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        if (k->offset_I_e_file != "none") {
            file_name = out_dir + "/specific_weights/" + k->name + "_" + to_string(current_epoch) + ".txt";
            of.open( file_name.c_str() );
            for (unsigned i=0; i<k->N; i++) {
                of << k->pop[i].offset_I_e << endl;
            }
            of.close();
        }
    } // end for over subnets
}


void Network::set_up_training() {

    try {
        string              temp;
        ifstream            f_in_training(training_config_yaml);
        YAML::Node          training = YAML::Load(f_in_training);

        if (f_in_training.fail()) {
            cerr << "ERROR: training_config_yaml (i.e. " << training_config_yaml << ") cannot be opened" << endl;
            exit(1);
        }
        current_epoch       = training["epoch"].as<unsigned>();
        dur_batches_path    = training["dur_batches_path"].as<string>();
        optimization_method = training["optimization_method"].as<string>();
        window_res          = training["window_res"].as<double>();
        step_res            = training["step_res"].as<double>();
        w_cut               = training["w_cut"].as<double>();
        POT                 = training["POT"].as<unsigned>();
        prob_noise          = training["prob_noise"].as<double>();
        amp_noise           = training["amp_noise"].as<double>();
        l_decay             = training["l_decay"].as<double>();
        l_rate0             = training["l_rate0"].as<double>();
        l_rate0_curr        = training["l_rate0_curr"].as<double>();
        correct_l_rates     = training["correct_l_rates"].as<vector<double>>();
        error_lin_coeff     = training["error_lin_coeff"].as<double>();
        epsilon_zero_act    = training["epsilon_zero_act"].as<double>();
        a_desid_zero        = training["a_desid_zero"].as<double>();

        momentum_w          = training["momentum_w"].as<double>();
        momentum_c          = training["momentum_c"].as<double>();
        if (optimization_method=="SGD") {
            optimization_method_id = 0;
        }
        else if (optimization_method=="SGD_MEAN") {
            optimization_method_id = 1;
        }
        else if (optimization_method=="ADAM") {
            optimization_method_id = 2;
            momentum2_w  = training["momentum2_w"].as<double>();
            momentum2_c  = training["momentum2_c"].as<double>();
        }
        else {
            cerr << "optimization method (" << optimization_method << ") NOT supported" << endl;
            exit(1);
        }

        // quantized training of weights
        quant_delta_weight_pos = training["quant_delta_weight_pos"].as<double>();
        quant_delta_weight_neg = training["quant_delta_weight_neg"].as<double>();
        n_quant                = training["n_quant"].as<double>();                          // the total number of weights is 2*n_quant +1 (null weight)
        for (auto k=subnets.begin(); k!=subnets.end(); k++) {
            for (auto target_p : k->train_out_weight) {
                if (target_p.second && k->alpha_quantized_weights[target_p.first]>0) {
                    vector<double>    temp_zeros(k->pop[0].neighbors[target_p.first].size(), 0.);
                    for (unsigned i=0; i<k->N; i++){
                        k->pop[i].weight_quant_cumul[target_p.first] = temp_zeros;
                    }
                }
            }
        }

        temp              = training["c_desid_zero"].as<string>();
        if (temp=="none") {
            // import optimal cs for each output neuron (basing on trainining order==0)
            c_desid_zero = -1;
            try {
                for (auto k=subnets.begin(); k!=subnets.end(); k++) {
                    if (k->train_order == 0) {
                        ifstream         f_in_optimal_cs(training["c_desid_zero_file"].as<string>());
                        double           temp_d;
                        vector<double>   temp_vector;

                        while (f_in_optimal_cs>>temp_d) {
                            temp_vector.push_back(temp_d);
                        }
                        f_in_optimal_cs.close();
                        if ( temp_vector.size()==k->N ) {
                            for (unsigned i=0; i<k->N; i++) {
                                k->pop[i].optimal_c_desid_zero = temp_vector[i];
                            }
                        }
                        else {
                            cerr << "\tERROR inconsistent lenght of file with optimal_c_desid_zero: " << temp_vector.size() << " != " << k->N << endl;
                            exit(1);
                        }
                    }
                }
            } catch (...) {
                cerr << "\tERROR occurred while parsing file with optimal_c_desid_zero (check epsilon_zero_act)" << endl;
                exit(1);
            }
        }
        else {
            c_desid_zero = stod(temp);
        }

        // NOTE: STARTING FROM THE CONNECTIONS TO THE OUTPUT LAYER
        // l_rates are the multiplicative factor with respect to l_rate0 and the default learning rate: fix [1,1...] for the default configuration]
        dF_dI_app         = training["dF_dI_app"].as<double>();
        dF_dI_no_act_fact = training["dF_dI_no_act_fact"].as<double>();

        if (POT%2==0) {
            cerr << "Error: POT is even..." << endl;
            exit(1);
        }
        f_in_training.close();
    } catch (...) {
        cerr << "\tERROR occurred while parsing file training_config_yaml (i.e. " << training_config_yaml << ")" << endl;
        exit(1);
    }

    // initialize ord_target_subnets
    for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        ord_target_subnets.push_back( -1 );
    }
    for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        if (k->train_order>=0) {
            ord_target_subnets[k->train_order]=k->id_subnet;
        }

    }

    // create errors folder
    system(("mkdir -p " + out_dir + "/errors").c_str());

    // load dur_batches
    ifstream         f_in_batches(dur_batches_path);
    double           temp_d;
    if (f_in_batches.is_open()) {

        while (f_in_batches>>temp_d) {
            dur_batches.push_back(temp_d);
        }
        f_in_batches.close();
    }
    else {
        cerr << "Error: unable to open file " << dur_batches_path << endl;
        exit(1);
    }

    // load desidered_output
    unsigned         N_neur;
    string           line;
    int              n=0,i, lenght=floor(t_end/step_res), start;
    vector<double>   temp(lenght,0.);
    size_t           end;
    const char       del_inl='\t';
    for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        if (k->desidered_output_path!="none") {
            ifstream         f_in_desidered(k->desidered_output_path);

            if (f_in_desidered.is_open()) {
                N_neur=k->N;
                while (getline(f_in_desidered,line)) {
                    if (n==N_neur) {
                        cerr << "Error: inconsistent (too high) number of rows (" << n+1 << ") in " << k->desidered_output_path << endl;
                        exit(1);
                    }
                    start=0;
                    end = line.find(del_inl, start);
                    if ( n!=stoi(line.substr(start, int(end)-start)) ) {
                        // the output of all neurons m
                        cerr << "Error: inconsistent neuron number " << stoi(line.substr(start, int(end)-start)) << " (n="<< n <<") in " << k->desidered_output_path << endl;
                        exit(1);
                    }
                    start = int(end)+1;
                    for (i=0; end!=string::npos; i++) {
                        end = line.find(del_inl, start);
                        temp[i] = stod(line.substr(start, int(end)-start));    // get number from string
                        start = int(end)+1;
                    }
                    // note: ther MUST NOT be "\t" at the and of each row
                    if (i != lenght) {
                        cerr << "Error: inconsistent lenght of desidered output for neuron " << n << " in " << k->desidered_output_path << endl;
                        exit(1);
                    }
                    (k->pop)[n].desidered_output = temp;
                    n++;
                } // end while on input file
            }
            else {
                cerr << "Error: unable to open file " << k->desidered_output_path << endl;
                exit(1);
            }
        } // end if desidered_output_path!="none"
    } // end for over subnets
}


double poissonian_CDF(unsigned x, double mu) {
    double         sum=0, temp_num=1, temp_den=1;

    if (mu < 0.05) {
        // XXX neglect sporadic single spikes when mu is small
        // if (x<2) {
        //     return 0.5;
        // }
        return (x>0) ? 1. : 0.5;
    }
    else if (x>10*mu){
        return 1.;
    }

    for(unsigned k=0; k<=x; k++) {
        sum += temp_num/temp_den;
        temp_num = temp_num*mu;
        temp_den = temp_den*(k+1);
    }
    return pow( M_E, -mu ) * sum;
}



// Neuron
Neuron::Neuron(unsigned _dim) {

    x.resize(_dim, 0);
    dim = _dim;
}

void Neuron::info() {
    cout << "\tx\t\t";
    print_vector(x);

    cout << "input_t_ex" << endl;
    for (auto i : input_t_ex) {
        cout << "\t  source_subnet_id: " << i.first << "\t  " << i.second.size();
    }
    cout << endl;
    cout << "input_t_in" << endl;
    for (auto i : input_t_in) {
        cout << "\t  source_subnet_id: " << i.first << "\t  " << i.second.size();
    }
    cout << endl;
    cout << "input_w_ex" << endl;
    for (auto i : input_w_ex) {
        cout << "\t  source_subnet_id: " << i.first << "\t  " << i.second.size();
    }
    cout << endl;
    cout << "input_w_in" << endl;
    for (auto i : input_w_in) {
        cout << "\t  source_subnet_id: " << i.first << "\t  " << i.second.size();
    }
    cout << endl;

    cout << "\tneighbours"<<endl;
    for (auto i : neighbors) {
        cout << "\t  target_subnet_id: " << i.first << "\t  ";
        print_vector(i.second);
    }
    cout << "\texternal inp\t";
    print_vector(input_t_ex[UINT_MAX]);
}


// SubNetwork
SubNetwork::SubNetwork(string _name, string _neuron_model, int _N, double _C_m, double _E_L,       \
           double _V_res, double _V_th,                                \
           double _t_ref, double _I_e, double _osc_amp, double _osc_omega, double _dev_ext_weight, \
           double _ext_in_rate, double _osc_amp_poiss, double _osc_omega_poiss,                    \
           double _tau_syn_ex, double _tau_syn_in, RandomGenerator _g) {
    name = _name;
    neuron_model = _neuron_model;
    N = _N;
    C_m = _C_m;
    E_L = _E_L;
    V_res = _V_res;
    V_th = _V_th;
    t_ref = _t_ref;
    I_e = _I_e;
    osc_amp = _osc_amp;
    osc_omega = _osc_omega;
    dev_ext_weight = _dev_ext_weight;
    ext_in_rate = _ext_in_rate;
    osc_amp_poiss = _osc_amp_poiss;
    osc_omega_poiss = _osc_omega_poiss;
    tau_syn_ex = _tau_syn_ex;
    tau_syn_in = _tau_syn_in;
    g = _g;

}


void SubNetwork::save_t_spikes (string out_file) {
    ofstream        of(out_file, ios::app);
    of << setprecision(9);
    for (unsigned i=0; i<N; i++) {
        if ( pop[i].t_spikes.size() > 0 ) {
            of << i << "  ";
            for (auto s : pop[i].t_spikes) {
                of << s << "  ";
            }
            of << endl;
        }
    }
}


void SubNetwork::info(map<unsigned, string> &subnet_index_name) {
    unsigned    temp_count=0;
    string      temp_string, temp_string_tr;
    cout << "Subnetwork " << name <<  endl;
    cout << "\tneuron_model\t\t" << neuron_model << endl;
    cout << "\tid_model\t\t" << id_model << endl;
    cout << "\tN\t\t\t" << N << endl;
    cout << "\tC_m\t\t\t" << C_m << endl;
    cout << "\tE_L\t\t\t" << E_L << endl;
    if (cond_based) {
        cout << "\tE_ex\t\t\t" << E_ex << endl;
        cout << "\tE_in\t\t\t" << E_in << endl;
    }
    cout << "\tV_res\t\t\t" << V_res << endl;
    cout << "\tV_th\t\t\t" << V_th << endl;
    cout << "\tt_ref\t\t\t" << t_ref << endl;
    cout << "\tI_e\t\t\t" << I_e << endl;
    cout << "\tspecific_I_e_file\t" << specific_I_e_file << endl;
    cout << "\text_in_rate\t\t" << ext_in_rate << endl;
    cout << "\tdev_ext_weight\t\t" << dev_ext_weight << endl;
    cout << "\ttau_syn_ex\t\t" << tau_syn_ex << endl;
    cout << "\ttau_syn_in\t\t" << tau_syn_in << endl;

    if (adap_4) {
        cout << "\ta_adaptive\t\t" << a_adaptive << endl;
        cout << "\tb_adaptive\t\t" << b_adaptive << endl;
        cout << "\ttau_w_adaptive\t\t" << tau_w_adaptive << endl;
    }

    // other
    switch(id_model) {
        case 2:
        case 3:
            cout << "\tdelta_T_aeif_cond_exp\t" << delta_T_aeif_cond_exp << endl;
            cout << "\tg_L\t\t\t" << g_L << endl;
            cout << "\tV_peak\t\t\t" << V_peak << endl;
            break;
        case 4:
        case 5:
            cout << "\tk_aqif_cond_exp\t" << k_aqif_cond_exp << endl;
            cout << "\tV_peak\t\t\t" << V_peak << endl;
            break;
        case 8:
            cout << "\tk_aqif_cond_exp\t" << k_aqif_cond_exp << endl;
            cout << "\tV_b_aqif2_cond_exp\t" << V_b_aqif2_cond_exp << endl;
            cout << "\tV_peak\t\t\t" << V_peak << endl;
            break;
        case 0:
        case 1:
            cout << "\tg_L\t\t\t" << g_L << endl;
            cout << "\tV_peak\t\t\t" << V_peak << " (same as V_th)" << endl;
            break;
        case 6:
        case 7:
            cout << "\tg_L\t\t\t" << g_L << endl;
            cout << "\tV_peak\t\t\t" << V_peak << " (same as V_th)" << endl;
            break;
        case 9:
            cout << "\tparrot_spike_path\t" << parrot_spike_path << endl;
            cout << "\tV_peak\t\t\t" << V_peak << " (same as V_th)" << endl;
            break;
        // case :
        //
        //     break;
        default:
            break;
      }

    cout << "\tnon_spec weights_ex FROM\t";
    if (weights_ex.size()==0) cout << "empty" << endl;
    else {
        cout << endl;
        for (auto i : weights_ex) {
            if (i.first==UINT_MAX) {
                cout << "\t    " << "ext" << "\t\t" << i.second << endl;
                continue;
            }
            // cout << "\t\t    " << subnet_index_name[i.first] << "   " << subnets[subnet_index_name[i.first]].excit_out[id_subnet] <<endl;
            // temp_string = subnets[subnet_index_name[i.first]].excit_out[id_subnet] ? "\t(ex)" : "\t(in)";
            temp_string = "";
            cout << "\t    " << subnet_index_name[i.first] << "\t\t" << i.second << temp_string << endl;
        }
    }
    cout << "\tnon_spec weights_in FROM\t";
    if (weights_in.size()==0) cout << "empty" << endl;
    else {
        cout << endl;
        for (auto i : weights_in) {
            if (i.first==UINT_MAX) {
                cerr << "Error: inhibitory external input is not supported" << endl;
                exit(1);
            }
            cout << "\t    " << subnet_index_name[i.first] << "\t\t" << i.second << endl;
        }
    }

    cout << "\tspec weights FROM\t\t";
    for (auto i : specific_in_weight) {
        if (i.second){
            if (i.first==UINT_MAX) {
                cerr << "Error: external input with specific weights is not supported" << endl;
                exit(1);
            }
            cout << "\n\t    " << subnet_index_name[i.first];
            temp_count = temp_count+1;
        }
    }
    if (temp_count==0) cout << "empty";
    cout << endl;

    cout << "\tdelays_out TO" << endl;
    for (auto i : delays_out) {
        cout << "\t    " << subnet_index_name[i.first] << "\t\t" << i.second << endl;
    }

    cout << "\tconn prob TO (with specificity & trainability)" << endl;
    for (auto i : probabilities_out) {
        temp_string = ( specific_out_weight[i.first] ) ? "spec    \t" : "non_spec\t";
        temp_string_tr = ( train_out_weight[i.first] ) ? "train    \t" : "non_train\t";
        cout << "\t    " << subnet_index_name[i.first] << "\t\t" << i.second << "  " << temp_string << temp_string_tr <<  endl;
        if (specific_out_weight[i.first]) {
            cout << "\t\t    out_weight_file: " << specific_out_weight_files[i.first] << endl;
        }
    }

    cout << "\tto_save\t\t\t";
    print_vector(to_save);

    // for (auto i : pop) {
    //     i.info();
    // }

}


// Network
Network::Network(double _t_end, double _dt, unsigned _input_mode, string _subnets_config_yaml, string _weights_config_yaml, \
    string _connections_config_yaml, string _to_save_config_yaml, string _training_config_yaml, string _out_dir, RandomGenerator _g, \
    string _input_mode_config, unsigned  _n_step, unsigned _repeat_specific_I_e, bool evolve_only) {
    cout << "Network constructor called" << endl;

    if ( subnet_model_id["iaf_cond_alpha"]!=0 || subnet_model_id["iaf_curr_alpha"]!=1 || \
         subnet_model_id["aeif_cond_exp"] !=2 || subnet_model_id["aeif_curr_exp"] !=3 || \
         subnet_model_id["aqif_cond_exp"] !=4 || subnet_model_id["aqif_curr_exp"] !=5 || \
         subnet_model_id["aqif2_cond_exp"]!=8 ||                                         \
         subnet_model_id["iaf_cond_exp"]  !=6 || subnet_model_id["iaf_curr_exp"]  !=7 || \
         subnet_model_id["parrot_neuron"] !=9 ||                                         \
         subnet_model_id.size()!=10) {
        cerr << "ERROR: modify here, in cases in EVOLVE(), in CREATE_SUBNETS(), in subnetwork.INFO() and in the docstring related to neuron_model." << endl;
        exit(1);
    }

    t_end = _t_end;
    n_step = _n_step;
    dt = _dt;
    if (t_end/n_step < dt) {
        cerr << "\tERROR: t_end/n_step < dt" <<endl;
        exit(1);
    }
    input_mode = _input_mode;

    subnets_config_yaml = _subnets_config_yaml;
    weights_config_yaml = _weights_config_yaml;
    connections_config_yaml = _connections_config_yaml;
    to_save_config_yaml = _to_save_config_yaml;
    training_config_yaml = _training_config_yaml;
    out_dir = _out_dir;
    repeat_specific_I_e = _repeat_specific_I_e;
    g = _g;

    subnet_name_index["ext"] = UINT_MAX;

    // cout << "before create SubNetwork" << endl; // ABC
    createSubnets(evolve_only);
    // cout << "after create SubNetwork" << endl; // ABC

    if (input_mode == 2) {
        ifstream            f_in_corr(_input_mode_config);
        YAML::Node          n_corr = YAML::Load(f_in_corr);

        system(("cp " + _input_mode_config + " ./" + out_dir ).c_str());

        try {
            rho_corr = n_corr["rho_corr"].as<double>();
            for (auto j : n_corr["corr_pops"]) {
                corr_pops.push_back( subnet_name_index[j.as<string>()] );
            }
        } catch (...) {
            cerr << "\tERROR occurred while parsing file _input_mode_config (i.e. " << _input_mode_config << ")" << endl;
            exit(1);
        }
    }

    createPops();
    load_specific_currents_and_offsets();
    load_specific_weights();

    info();

    if (evolve_only) {
        auto start = chrono::system_clock::now();
        while (t < t_end) {
            if (t + t_end/n_step > t_end) {
                cout.precision(20);
                cout << "entered if in 'while (t < t_end)' t: " << t << "\t t_end: " << t_end << endl;
                evolve(t_end -t);
            }
            else {
                evolve(t_end/n_step);
                free_past();
            }
        }
        chrono::duration<double> elapsed_seconds = chrono::system_clock::now()-start;
        cout << "evolution required " << elapsed_seconds.count() << " s" << endl;
    }
    else {
        set_up_training();
        info_training();
        /// XXXA
        // int STOP_AFTER_ind =0;
        for (auto T : dur_batches) {
            // cout << "evolving starting from " << t << " to " <<  t+T << endl;
            evolve(T);
            // cout << "training starting from " << t-T << " to " <<  t << endl;
            train(T);
            free_past();
            /// XXXA
            // STOP_AFTER_ind+=1;
            // if (STOP_AFTER_ind == 5) break;
        }
        save_specific_weights();
        save_offset_currents();

        // check desidered
        for (unsigned i=0; i<subnets[ord_target_subnets[0]].N; i++) {
            cout << i << " lenght of des:  " << subnets[ord_target_subnets[0]].pop[i].desidered_output.size() << endl;
        }
    }

    // info();
    // cout << "check for overflow in the integrators if necessary (quadratic and exponential)!!" << endl;

}


void Network::createSubnets(bool evolve_only) {

    try {
        ifstream        f_in(subnets_config_yaml);
        YAML::Node      config = YAML::Load(f_in);
        unsigned        index=0;

        if (f_in.fail()) {
            cerr << "ERROR: subnets_config_yaml (i.e. " << subnets_config_yaml << ") cannot be opened" << endl;
            exit(1);
        }

        for (auto i : config) {

            if (subnet_model_id.find(i["neuron_model"].as<string>()) == subnet_model_id.end()){
                cerr << "\tERROR: neuron model not supported: " << i["neuron_model"].as<string>() << endl;
                exit(1);
            }

            // update some useful dictionaries
            subnet_name_index[i["name"].as<string>()] = index;
            subnet_index_name[index] = i["name"].as<string>();

            // cout << i["name"].as<string>() << " before subnet creation " << endl;
            // cout << "   name " << i["name"].as<string>()  << endl;
            // cout << "   neuron_model " << i["neuron_model"].as<string>()  << endl;
            // cout << "   N "  << i["N"].as<unsigned>()  << endl;
            // cout << "   C_m " << i["C_m"].as<double>()  << endl;
            // cout << "   E_L " << i["E_L"].as<double>()  << endl;
            // cout << "   V_res " << i["V_res"].as<double>()  << endl;
            // cout << "   V_th " << i["V_th"].as<double>()  << endl;
            // cout << "   t_ref " << i["t_ref"].as<double>()  << endl;
            // cout << "   I_e " << i["I_e"].as<double>()  << endl;
            // cout << "   osc_amp " << i["osc_amp"].as<double>()  << endl;
            // cout << "   osc_omega " << i["osc_omega"].as<double>()  << endl;
            // cout << "   dev_ext_weight " << i["dev_ext_weight"].as<double>()  << endl;
            // cout << "   ext_in_rate " << i["ext_in_rate"].as<double>()  << endl;
            // cout << "   osc_amp_poiss " << i["osc_amp_poiss"].as<double>()  << endl;
            // cout << "   osc_omega_poiss " << i["osc_omega_poiss"].as<double>()  << endl;
            // cout << "   tau_syn_ex " << i["tau_syn_ex"].as<double>()  << endl;
            // cout << "   tau_syn_in " << i["tau_syn_in"].as<double>()  << endl;

            // generate subnet
            subnets.push_back( SubNetwork(i["name"].as<string>(), i["neuron_model"].as<string>(), i["N"].as<unsigned>(),    \
                     i["C_m"].as<double>(), i["E_L"].as<double>(),                                                          \
                     i["V_res"].as<double>(), i["V_th"].as<double>(), i["t_ref"].as<double>(),                              \
                     i["I_e"].as<double>(), i["osc_amp"].as<double>(), i["osc_omega"].as<double>(),                         \
                     i["dev_ext_weight"].as<double>(), i["ext_in_rate"].as<double>(),                                       \
                     i["osc_amp_poiss"].as<double>(), i["osc_omega_poiss"].as<double>(),                                    \
                     i["tau_syn_ex"].as<double>(),i["tau_syn_in"].as<double>(), g) );

            subnets[index].id_subnet = index;
            subnets[index].id_model = subnet_model_id[i["neuron_model"].as<string>()];
            subnets[index].specific_I_e_file = i["specific_I_e_file"].as<string>();
            subnets[index].offset_I_e_file = i["offset_I_e_file"].as<string>();
            subnets[index].desidered_output_path = i["desidered_output_path"].as<string>();
            subnets[index].train_order = i["train_order"].as<int>();
            subnets[index].train_offset_I_e = i["train_offset_I_e"].as<bool>();

            if (evolve_only==false) {
                subnets[index].w_max = i["w_max"].as<double>();
                subnets[index].w_min = i["w_min"].as<double>();
                subnets[index].offset_I_e_max = i["offset_I_e_max"].as<double>();
                subnets[index].offset_I_e_min = i["offset_I_e_min"].as<double>();
                subnets[index].offset_I_e_ref = i["offset_I_e_ref"].as<double>();
                subnets[index].ldecay_offset  = i["ldecay_offset"].as<double>();
                subnets[index].decay_factor_max       = i["decay_factor_max"].as<double>();
                subnets[index].decay_factor_steepness = i["decay_factor_steepness"].as<double>();
                subnets[index].decay_factor_position  = i["decay_factor_position"].as<double>();
            }

            if (subnets[index].specific_I_e_file!="none") {
                system(("mkdir -p " + out_dir + "/specific_I_es").c_str());
                system(("cp " + subnets[index].specific_I_e_file + " ./" + out_dir + "/specific_I_es/").c_str());
            }

            // particular attributes
            // current VS condactance based models
            if ((subnets[index].neuron_model).find("cond")!=string::npos) {        // if neuron_model contains cond
                subnets[index].cond_based = true;
                subnets[index].E_ex = i["E_ex"].as<double>();
                subnets[index].E_in = i["E_in"].as<double>();
            }
            else if ((subnets[index].neuron_model).find("curr")!=string::npos) {    // if neuron_model contains curr
                subnets[index].cond_based = false;
            }
            else if ((subnets[index].neuron_model).find("parrot")!=string::npos) {  // if neuron_model contains parrot
                subnets[index].cond_based = false;
            }
            else {
                cerr << "Error: the model is neither current nor conductance based, nor parrot" << endl;
                exit(1);
            }

            // alpha_PSS post synaptic shape: true if neuron_model contains alpha
            subnets[index].alpha_PSS = ((subnets[index].neuron_model).find("alpha")!=string::npos) ? true : false;

            // adaptation
            subnets[index].adap_4 = model_is_adaptive[subnets[index].neuron_model];
            if (subnets[index].adap_4) {
                subnets[index].a_adaptive = i["a_adaptive"].as<double>();
                subnets[index].b_adaptive = i["b_adaptive"].as<double>();
                subnets[index].tau_w_adaptive = i["tau_w_adaptive"].as<double>();
            }

            // other
            switch(subnets[index].id_model) {
                case 2:
                case 3:
                    subnets[index].delta_T_aeif_cond_exp = i["delta_T_aeif_cond_exp"].as<double>();
                    subnets[index].g_L = i["g_L"].as<double>();
                    subnets[index].V_peak = i["V_peak"].as<double>();
                    break;
                case 4:
                case 5:
                    subnets[index].k_aqif_cond_exp = i["k_aqif_cond_exp"].as<double>();
                    subnets[index].V_peak = i["V_peak"].as<double>();
                    break;
                case 8:
                    subnets[index].k_aqif_cond_exp = i["k_aqif_cond_exp"].as<double>();
                    subnets[index].V_b_aqif2_cond_exp = i["V_b_aqif2_cond_exp"].as<double>();
                    subnets[index].V_peak = i["V_peak"].as<double>();
                    break;
                case 0:
                case 1:
                    subnets[index].g_L = i["g_L"].as<double>();
                    subnets[index].V_peak = subnets[index].V_th;
                    break;
                case 6:
                case 7:
                    subnets[index].g_L = i["g_L"].as<double>();
                    subnets[index].V_peak = subnets[index].V_th;
                    break;
                case 9:
                    subnets[index].parrot_spike_path = i["parrot_spike_path"].as<string>();
                    subnets[index].V_peak = subnets[index].V_th;
                    system(("mkdir -p " + out_dir + "/" + subnets[index].name).c_str());
                    system(("cp -r " + subnets[index].parrot_spike_path + " ./" + out_dir + "/" + subnets[index].name).c_str());
                    break;
                // case :
                //
                //     break;
                default:
                    break;
              }

            index += 1;
        }
        f_in.close();
    } catch (...) {
        cerr << "\tERROR occurred while parsing file subnets_config_yaml (i.e. " << subnets_config_yaml << ") (w_max/w_min) " << endl;
        exit(1);
    }

    // cout << "before to_save" << endl; // ABC
    // to_save
    try {
        ifstream            f_in_save(to_save_config_yaml);
        YAML::Node          to_save = YAML::Load(f_in_save);
        unsigned            pop_id;

        if (f_in_save.fail()) {
            cerr << "ERROR: to_save_config_yaml (i.e. " << to_save_config_yaml << ") cannot be opened" << endl;
            exit(1);
        }

        for (auto i : to_save){
            if (i.Type() == YAML::NodeType::Map) {
                for (auto j : i) {
                    pop_id = subnet_name_index[j.first.as<string>()];
                    subnets[ pop_id ].to_save = j.second.as<vector<unsigned>>();

                    if (subnets[ pop_id ].to_save.size() > 0) {
                        if ( *max_element(subnets[pop_id].to_save.begin(), subnets[pop_id].to_save.end())>subnets[pop_id].N ) {
                            cerr << "Error: to_save vector of " << subnet_name_index[j.first.as<string>()] << " is not consistent" <<endl;
                            exit(1);
                        }
                    }

                }
            }
        }
        f_in_save.close();
    } catch (...) {
        cerr << "\tERROR occurred while parsing file to_save_config_yaml (i.e. " << to_save_config_yaml << ")" << endl;
        exit(1);
    }

    // cout << "before weights" << endl; // ABC
    // weights
    try {
        ifstream        f_in_w(weights_config_yaml);
        YAML::Node      config_w = YAML::Load(f_in_w);
        double          temp, correction_alpha;
        unsigned        population_id;

        if (f_in_w.fail()) {
            cerr << "ERROR: weights_config_yaml (i.e. " << weights_config_yaml << ") cannot be opened" << endl;
            exit(1);
        }

        for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        // for subnet in subnets
            for (auto i : config_w) {
            // for node in weights_config_yaml
                if (k->name == i["name"].as<string>()) {
                    // check correspondence between subnet name and node name

                    for (auto j : i) {
                        if (j.first.as<string>() == "name") continue;

                        population_id = subnet_name_index[j.first.as<string>()];        // source population
                        if (j.second.as<string>()!="none") {        // if different (!) from "none"
                            // general weights

                            /// initialize specific_in_weight and specific_out_weight
                            k->specific_in_weight[population_id] = false;
                            if (population_id!=UINT_MAX) subnets[population_id].specific_out_weight[k->id_subnet] = false;

                            temp = j.second.as<double>();
                            if (temp < 0) {
                                // inhibitory
                                correction_alpha = (k->alpha_PSS) ? M_E/(k->tau_syn_in) : 1.;
                                k->weights_in[population_id] = - correction_alpha*temp;
                                subnets[population_id].excit_out[k->id_subnet] = false;
                            }
                            else {
                                // excitatory
                                correction_alpha = (k->alpha_PSS) ? M_E/(k->tau_syn_ex) : 1.;
                                k->weights_ex[population_id] = correction_alpha*temp;
                                if (population_id!=UINT_MAX) subnets[population_id].excit_out[k->id_subnet] = true;
                            }
                        }
                        else {
                            // specific weights
                            k->specific_in_weight[population_id] = true;
                            subnets[population_id].specific_out_weight[k->id_subnet] = true;
                        }

                    }

                    break;
                }
            }
        }
        f_in_w.close();
    } catch (...) {
        cerr << "\tERROR occurred while parsing file weights_config_yaml (i.e. " << weights_config_yaml << ")" << endl;
        exit(1);
    }

    // cout << "before connections" << endl; // ABC
    // connections probabilities and delays
    try {
        ifstream        f_in_conn(connections_config_yaml);
        YAML::Node      config_conn = YAML::Load(f_in_conn);
        double          prob;
        unsigned        population_id;
        string          file;

        if (f_in_conn.fail()) {
            cerr << "ERROR: connections_config_yaml (i.e. " << connections_config_yaml << ") cannot be opened" << endl;
            exit(1);
        }

        for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        // for subnet in subnets
            for (auto i : config_conn) {
            //for node in connections_config_yaml
                if (k->name == i["name"].as<string>()) {
                // check correspondence between subnet name and node name

                    for (auto iter=next(i.begin(),1); iter!=i.end(); iter=next(iter,5)) {
                    // jumping the "name" entry read the probability, the delay, the file for specific_weights and the bool for training
                        population_id = subnet_name_index[iter->first.as<string>()];

                        // FIRST ENTRY: prob
                        prob = iter->second.as<double>();
                        if (prob < 0) {
                            k->probabilities_out[population_id] = -prob;
                            k->reverse_effect[population_id] = true;
                        }
                        else {
                            k->probabilities_out[population_id] = prob;
                            k->reverse_effect[population_id] = false;
                        }
                        // SECOND ENTRY: delay
                        k->delays_out[population_id] = (next(iter,1))->second.as<double>();
                        // THIRD ENTRY: file
                        if ((next(iter,2))->second.as<string>()=="none") {
                            // general weights
                            // check coherence with weights_config_yaml
                            if (k->specific_out_weight[population_id]!=false || subnets[population_id].specific_in_weight[k->id_subnet]!=false) {
                                cerr << "Error: initialization of specific weights is not coherent in " << weights_config_yaml << " AND " << connections_config_yaml << endl;
                                exit(1);
                            }
                        }
                        else {
                            // specific weights
                            // check coherence with weights_config_yaml
                            if (k->specific_out_weight[population_id]!=true || subnets[population_id].specific_in_weight[k->id_subnet]!=true) {
                                cerr << "Error: initialization of specific weights is not coherent in " << weights_config_yaml << " AND " << connections_config_yaml << endl;
                                exit(1);
                            }
                            file = (next(iter,2))->second.as<string>();
                            k->specific_out_weight_files[population_id] = file;
                            system(("mkdir -p " + out_dir + "/specific_weights").c_str());
                            system(("cp " + file + " ./" + out_dir + "/specific_weights/").c_str());
                        }
                        // FOURTH ENTRY: train_bool
                        if ((next(iter,3))->second.as<string>()=="true") {
                            k->train_out_weight[population_id] = true;
                            if (k->specific_out_weight[population_id]!=true) {
                                cerr << "Error: training non specific weights from " << k->name << " to " << subnets[population_id].name << endl;
                                exit(1);
                            }
                        }
                        else {
                            k->train_out_weight[population_id] = false;
                        }
                        // FIFTH ENTRY: quantized training
                        // alpha_quantized_weights>0 means that training should be quantized for this connection
                        k->alpha_quantized_weights[population_id] = (next(iter,4))->second.as<double>();
                        if (k->alpha_quantized_weights[population_id]>0) {
                            if (k->train_out_weight[population_id]!=true) {
                                cerr << "Error: quantized training of non trainable connections from " << k->name << " to " << subnets[population_id].name << endl;
                                exit(1);
                            }
                        }
                    }
                    break;
                }
            }
        }
        cout << "Eventual specific files for currents and weights have been copied" << endl;
        f_in_conn.close();
    } catch (...) {
        cerr << "\tERROR occurred while parsing file connections_config_yaml (i.e. " << connections_config_yaml << ")" << endl;
        exit(1);
    }
}


void Network::createPops() {
    unsigned                    population_id;
    double                      probability, tmp_sum;
    vector<unsigned>            neig_temp;
    map <unsigned, unsigned>    neig_count;

    // Neuron initialization
    for (auto k=subnets.begin(); k!=subnets.end(); k++) {
    // for subnet in subnets
        for (unsigned i=0; i<k->N; i++){
        // for each neuron in the population
            (k->pop).push_back( Neuron(supported_models[k->neuron_model]) );
            (k->pop)[i].x[0] = k->E_L;
            (k->pop)[i].t_last_spike = -1 - k->t_ref;
            if (g.getRandomUniform() > 0.5) (k->pop)[i].ext_weight = k->weights_ex[UINT_MAX] + k->dev_ext_weight * g.getRandomUniform();
            else (k->pop)[i].ext_weight = k->weights_ex[UINT_MAX] - k->dev_ext_weight * g.getRandomUniform();

            for (auto item : k->weights_ex) {
                (k->pop)[i].next_sp_ex_index[item.first] = 0;
            }
            for (auto item : k->weights_in) {
                (k->pop)[i].next_sp_in_index[item.first] = 0;
            }
        }

        // load parrot spikes
        if (k->id_model==9) {
            for (unsigned i=0; i<k->N; i++){
                parse_parrot_spikes( k->pop[i].parrot_spikes, (out_dir + "/" + k->name + "/" + to_string(i) + ".txt").c_str() );
            }
        }
    }


    // These two parts must be kept separated!
    for (auto k=subnets.begin(); k!=subnets.end(); k++) {
    // for subnet in subnets

        neig_count.clear();
        for (auto j : k->probabilities_out) {
            neig_count[j.first] = 0;
        }

        for (unsigned i=0; i<k->N; i++){
        // for each neuron in the subnet population

            // (outwards) neighbors initialization
            for (auto j : k->probabilities_out) {
            //  for item in dictionary probabilities out
                neig_temp.clear();
                population_id = j.first;
                probability = j.second;
                if (probability < 0) {
                    cerr << "\tERROR: negative probability detected" << endl;
                    exit(1);
                }
                for(unsigned n=0; n<subnets[population_id].N; n++){
                    if (g.getRandomUniform()<probability) {
                        neig_temp.push_back(n);
                    }
                }
                (k->pop)[i].neighbors[population_id] = neig_temp;
                neig_count[population_id] += neig_temp.size();
            }

            //  exernal input initialization for input_mode 0 (base mode)
            if (input_mode == 0 && k->ext_in_rate > EPSILON){
                tmp_sum = 0;
                while (tmp_sum < dt) {
                    tmp_sum += (- log(g.getRandomUniform()) ) / k->ext_in_rate;
                    (k->pop)[i].input_t_ex[UINT_MAX].push_back(tmp_sum);
                }
            }
        } // end for over k->pop

        cout << k->name << " connected to " << endl;
        for (auto j : neig_count) {
            cout << "\t" << subnets[j.first].name << " with " <<  j.second << " tot conn:\t" << double(j.second)/k->N  << " per source Neuron\t" << double(j.second) / subnets[j.first].N << " per target neuron"<< endl;
        }

        if (input_mode == 2 && k->ext_in_rate > EPSILON){

            // non striatum pop
            if ( find( corr_pops.begin(), corr_pops.end(), k-subnets.begin() ) == corr_pops.end() ) {
                cout << "\t\t\t" << k->name << " IS NOT in corr_pops!" << endl;

                for (unsigned i=0; i<k->N; i++) {
                // for each neuron in the population
                    tmp_sum = 0.;
                    while (tmp_sum < dt) {
                        tmp_sum += (- log(g.getRandomUniform()) ) / k->ext_in_rate;
                        (k->pop)[i].input_t_ex[UINT_MAX].push_back(tmp_sum);
                    }
                }
            }

            // striatum pop
            else {
                cout << "\t\t\t" << k->name << " IS in corr_pops!" << endl;
                tmp_sum = 0;
                while (tmp_sum < dt) {
                    tmp_sum += (- log(g.getRandomUniform()) ) / k->ext_in_rate * rho_corr;
                    for (unsigned i=0; i<k->N; i++) {
                        if (g.getRandomUniform() < rho_corr) {
                            (k->pop)[i].input_t_ex[UINT_MAX].push_back(tmp_sum);
                        }
                    }
                }
                corr_last_time[k->id_subnet] = tmp_sum;
            }

        }

    } // end for over subnets
}


void Network::load_specific_currents_and_offsets() {

    // specific external currents
    try {
        for (auto k=subnets.begin(); k!=subnets.end(); k++) {
            // if "none" is passed, no specific current is employed
            if (k->specific_I_e_file == "none") {
                for (unsigned i=0; i<k->N; i++){
                    (k->pop)[i].specific_I_e_bool = false;
                }
            }
            else {
                parse_external_current( k->pop, k->specific_I_e_file, floor(t_end/(dt*repeat_specific_I_e)) );      // XXX floor or floor+1?
            }
        }
    } catch (...) {
        cerr << "\tERROR occurred while parsing specific external current files: (check for tabs at the end of raws)" << endl;
        exit(1);
    }
    // offset external currents
    try {
        for (auto k=subnets.begin(); k!=subnets.end(); k++) {
            if (k->offset_I_e_file == "none") {
                for (unsigned i=0; i<k->N; i++){
                    (k->pop)[i].offset_I_e = 0.;
                    (k->pop)[i].offset_aux1 = 0.;
                    (k->pop)[i].offset_aux2 = 0.;

                }
            }
            else {
                ifstream                in_file(k->offset_I_e_file);
                string                  line;
                int                     i=0;
                while (getline(in_file,line)) {
                    if (i==k->N) {
                        cerr << "\tERROR: inconsistent number of lines in specific offset file: " << k->offset_I_e_file << endl;
                        exit(1);
                    }
                    (k->pop)[i].offset_I_e = stod(line);
                    (k->pop)[i].offset_aux1 = 0.;
                    (k->pop)[i].offset_aux2 = 0.;
                    i++;
                }
                if (i!=k->N) {
                    cerr << "\tERROR: inconsistent number of lines in specific offset file: " << k->offset_I_e_file << endl;
                    exit(1);
                }
            }
        }
    } catch (...) {
        cerr << "\tERROR occurred while parsing offset external current files: (check for tabs at the end of raws)" << endl;
        exit(1);
    }

}

void Network::load_specific_weights() {

    // specific out weights
    try {
        for (auto k=subnets.begin(); k!=subnets.end(); k++) {
            // if "none" is passed, no specific current is employed
            for (auto target_p=(k->specific_out_weight_files).begin(); target_p!=(k->specific_out_weight_files).end(); target_p++) {
                parse_specific_weight(k->pop, target_p->first, target_p->second);
            }
        }
    } catch (...) {
        cerr << "\tERROR occurred while parsing specific weight files: (check for tabs at the end of raws)" << endl;
        exit(1);
    }
}


double Network::input_func(double y_, double r0_, double A_, double omega_, double t0_, double t_) {
    return t_ - t0_ + A_/omega_ * (cos(omega_*t0_) - cos(omega_*t_)) - y_/r0_;
}

double Network::find_sol_bisection(double y_, double r0_, double A_, double omega_, double t0_) {
    double          tD = t0_, tU = t0_+1./r0_, val = 1e10, ftD, ftU, tbar;
    ftD = input_func(y_,r0_,A_,omega_,t0_,tD);
    ftU = input_func(y_,r0_,A_,omega_,t0_,tU);

    while ( ftD*ftU > 0 ) {
        tD = tU;
        tU += 1./r0_;
        ftD = ftU;
        ftU = input_func(y_,r0_,A_,omega_,t0_,tU);
    }
    while (abs(val) > TOLL) {
        tbar = (tD+tU)/2.;
        val = input_func(y_,r0_,A_,omega_,t0_,tbar);
        if (val * ftD > 0) {
            tD = tbar;
        }
        else {
            tU = tbar;
        }
    }
    return tbar;
}


// INTEGRATORS with EULER METHOD: XXX check overflow if necessary...

// void integrator_iaf_cond_alpha(const state_type_iaf_cond_alpha &x, state_type_iaf_cond_alpha &dxdt, const double t,          \
//         double s_C_m, double s_g_L, double s_E_L, double s_E_ex, double s_E_in, double s_I_e, double s_osc_amp, double s_osc_omega, \
//         double s_tau_syn_ex, double s_tau_syn_in) {
//     dxdt[0] = 1/s_C_m * ( - s_g_L*(x[0]-s_E_L) - x[1]*(x[0]-s_E_ex) - x[3]*(x[0]-s_E_in) + s_I_e \
//                           + s_osc_amp*sin(s_osc_omega*t) );                                             // membrain potential
//     dxdt[1] = -x[1]/s_tau_syn_ex + x[2];                                                                // excitatory syn conductance
//     dxdt[2] = -x[2]/s_tau_syn_ex;                                                                       // excitatory backup variable
//     dxdt[3] = -x[3]/s_tau_syn_in + x[4];                                                                // inhibitory syn conductance
//     dxdt[4] = -x[4]/s_tau_syn_in;                                                                       // inhibitory backup variable
// }

void euler_step_iaf_cond_alpha(vector<double> &x, double t, double dt, \
        double s_C_m, double s_g_L, double s_E_L, double s_E_ex, double s_E_in, double s_I_e, double s_osc_amp, double s_osc_omega, \
        double s_tau_syn_ex, double s_tau_syn_in) {

    x[0] += 1/s_C_m * ( - s_g_L*(x[0]-s_E_L) - x[3]*(x[0]-s_E_ex) - x[4]*(x[0]-s_E_in) + s_I_e \
                          + s_osc_amp*sin(s_osc_omega*t) ) * dt;                                            // membrain potential
    x[3] += (-x[3]/s_tau_syn_ex + x[1]) * dt;                                                               // excitatory syn conductance
    x[1] += (-x[1]/s_tau_syn_ex) * dt;                                                                      // excitatory backup variable
    x[4] += (-x[4]/s_tau_syn_in + x[2]) * dt;                                                               // inhibitory syn conductance
    x[2] += (-x[2]/s_tau_syn_in) * dt;                                                                      // inhibitory backup variable

}

void euler_step_iaf_curr_alpha(vector<double> &x, double t, double dt, \
        double s_C_m, double s_g_L, double s_E_L, double s_I_e, double s_osc_amp, double s_osc_omega, \
        double s_tau_syn_ex, double s_tau_syn_in) {

    x[0] += 1/s_C_m * ( - s_g_L*(x[0]-s_E_L) + x[3] - x[4] + s_I_e \
                          + s_osc_amp*sin(s_osc_omega*t) ) * dt;                                             // membrain potential
    x[3] += (-x[3]/s_tau_syn_ex + x[1]) * dt;                                                                // excitatory syn current
    x[1] += (-x[1]/s_tau_syn_ex) * dt;                                                                       // excitatory backup variable
    x[4] += (-x[4]/s_tau_syn_in + x[2]) * dt;                                                                // inhibitory syn current
    x[2] += (-x[2]/s_tau_syn_in) * dt;                                                                       // inhibitory backup variable

}

void euler_step_aeif_cond_exp(vector<double> &x, double t, double dt, \
        double s_C_m, double s_g_L, double s_E_L, double s_E_ex, double s_E_in, double s_I_e, double s_osc_amp, double s_osc_omega, \
        double s_tau_syn_ex, double s_tau_syn_in, double s_V_th, double s_a, double s_tau_w, double s_delta_T, double s_V_peak) {

    vector<double>          init = x;
    double                  check_ov=(x[0]-s_V_th)/s_delta_T;

    // check for overflow
    if (check_ov<10) {
        x[0] += 1./s_C_m * ( - s_g_L*(x[0]-s_E_L) + s_g_L*s_delta_T*exp(check_ov) \
                             - x[1]*(x[0]-s_E_ex) - x[2]*(x[0]-s_E_in) - x[3] + s_I_e \
                             + s_osc_amp*sin(s_osc_omega*t) ) * dt;
    }
    else {
        x[0] = s_V_peak+EPSILON;
    }                                                                                                      // membrain potential
    x[1] += (-x[1]/s_tau_syn_ex) * dt;                                                                     // excitatory syn conductance
    x[2] += (-x[2]/s_tau_syn_in) * dt;                                                                     // inhibitory syn conductance
    x[3] += (-x[3]/s_tau_w + s_a/s_tau_w * (init[0]-s_E_L)) * dt;                                          // adaptation variable
}

void euler_step_aeif_curr_exp(vector<double> &x, double t, double dt, \
        double s_C_m, double s_g_L, double s_E_L, double s_I_e, double s_osc_amp, double s_osc_omega, \
        double s_tau_syn_ex, double s_tau_syn_in, double s_V_th, double s_a, double s_tau_w, double s_delta_T) {

    vector<double>          init = x;
    // cout << "check for overflow!" << endl;
    x[0] += 1./s_C_m * ( - s_g_L*(x[0]-s_E_L) + s_g_L*s_delta_T*exp((x[0]-s_V_th)/s_delta_T) \
                         + x[1] - x[2] - x[3] + s_I_e \
                         + s_osc_amp*sin(s_osc_omega*t) ) * dt;                                            // membrain potential
    x[1] += (-x[1]/s_tau_syn_ex) * dt;                                                                     // excitatory syn current
    x[2] += (-x[2]/s_tau_syn_in) * dt;                                                                     // inhibitory syn current
    x[3] += (-x[3]/s_tau_w + s_a/s_tau_w * (init[0]-s_E_L)) * dt;                                          // adaptation variable
}

// void integrator_aqif_cond_exp(const state_type_aqif_cond_exp &x, state_type_aqif_cond_exp &dxdt, const double t,             \
//         double s_C_m, double s_k, double s_E_L, double s_E_ex, double s_E_in, double s_I_e, double s_osc_amp, double s_osc_omega,   \
//         double s_tau_syn_ex, double s_tau_syn_in, double s_V_th, double s_a, double s_tau_w) {
//     dxdt[0] = 1/s_C_m * ( s_k*(x[0]-s_E_L)*(x[0]-s_V_th) - x[1]*(x[0]-s_E_ex) - x[2]*(x[0]-s_E_in) \
//                         - x[3] + s_I_e + s_osc_amp*sin(s_osc_omega*t) );                                // membrain potential
//     dxdt[1] = -x[1]/s_tau_syn_ex;                                                                       // excitatory syn conductance
//     dxdt[2] = -x[2]/s_tau_syn_in;                                                                       // inhibitory syn conductance
//     dxdt[3] = -x[3]/s_tau_w + s_a/s_tau_w * (x[0]-s_E_L);                                               // adaptation variable
// }

void euler_step_aqif_cond_exp(vector<double> &x, double t, double dt, \
        double s_C_m, double s_k, double s_E_L, double s_E_ex, double s_E_in, double s_I_e, double s_osc_amp, double s_osc_omega,   \
        double s_tau_syn_ex, double s_tau_syn_in, double s_V_th, double s_a, double s_tau_w) {

    vector<double>          init = x;
    // cout << "check for overflow!" << endl;
    x[0] += 1/s_C_m * ( s_k*(x[0]-s_E_L)*(x[0]-s_V_th) - x[1]*(x[0]-s_E_ex) - x[2]*(x[0]-s_E_in) \
                        - x[3] + s_I_e + s_osc_amp*sin(s_osc_omega*t) ) * dt;                              // membrain potential
    x[1] += (-x[1]/s_tau_syn_ex) * dt;                                                                     // excitatory syn conductance
    x[2] += (-x[2]/s_tau_syn_in) * dt;                                                                     // inhibitory syn conductance
    x[3] += (-x[3]/s_tau_w + s_a/s_tau_w * (init[0]-s_E_L)) * dt;                                          // adaptation variable
}

void euler_step_aqif_curr_exp(vector<double> &x, double t, double dt, \
        double s_C_m, double s_k, double s_E_L, double s_I_e, double s_osc_amp, double s_osc_omega,   \
        double s_tau_syn_ex, double s_tau_syn_in, double s_V_th, double s_a, double s_tau_w) {

    vector<double>          init = x;
    // cout << "check for overflow!" << endl;
    x[0] += 1/s_C_m * ( s_k*(x[0]-s_E_L)*(x[0]-s_V_th) + x[1] - x[2] \
                        - x[3] + s_I_e + s_osc_amp*sin(s_osc_omega*t) ) * dt;                              // membrain potential
    x[1] += (-x[1]/s_tau_syn_ex) * dt;                                                                     // excitatory syn current
    x[2] += (-x[2]/s_tau_syn_in) * dt;                                                                     // inhibitory syn current
    x[3] += (-x[3]/s_tau_w + s_a/s_tau_w * (init[0]-s_E_L)) * dt;                                          // adaptation variable
}


// XXX these could be delated...
// void integrator_aqif2_cond_exp(const state_type_aqif_cond_exp &x, state_type_aqif_cond_exp &dxdt, const double t,             \
//         double s_C_m, double s_k, double s_E_L, double s_E_ex, double s_E_in, double s_I_e, double s_osc_amp, double s_osc_omega,    \
//         double s_tau_syn_ex, double s_tau_syn_in, double s_V_th, double s_a, double s_tau_w, double s_V_b) {
//     dxdt[0] = 1/s_C_m * ( s_k*(x[0]-s_E_L)*(x[0]-s_V_th) - x[1]*(x[0]-s_E_ex) \
//                           - x[2]*(x[0]-s_E_in) - x[3] + s_I_e + s_osc_amp*sin(s_osc_omega*t) );         // membrain potential
//     dxdt[1] = -x[1]/s_tau_syn_ex;                                                                       // excitatory syn conductance
//     dxdt[2] = -x[2]/s_tau_syn_in;                                                                       // inhibitory syn conductance
//     if (x[0] < s_V_b)   dxdt[3] = -x[3]/s_tau_w + s_a/s_tau_w * pow((x[0]-s_V_b),3);
//     else                dxdt[3] = -x[3]/s_tau_w;                                                        // adaptation variable
// }

void euler_step_aqif2_cond_exp(vector<double> &x, double t, double dt, \
        double s_C_m, double s_k, double s_E_L, double s_E_ex, double s_E_in, double s_I_e, double s_osc_amp, double s_osc_omega,    \
        double s_tau_syn_ex, double s_tau_syn_in, double s_V_th, double s_a, double s_tau_w, double s_V_b) {

    vector<double>          init = x;
    // cout << "check for overflow!" << endl;
    x[0] += 1/s_C_m * ( s_k*(x[0]-s_E_L)*(x[0]-s_V_th) - x[1]*(x[0]-s_E_ex) \
                          - x[2]*(x[0]-s_E_in) - x[3] + s_I_e + s_osc_amp*sin(s_osc_omega*t) ) * dt;         // membrain potential
    x[1] += (-x[1]/s_tau_syn_ex) * dt;                                                                       // excitatory syn conductance
    x[2] += (-x[2]/s_tau_syn_in) * dt;                                                                       // inhibitory syn conductance
    if (init[0] < s_V_b)   x[3] += (-x[3]/s_tau_w + s_a/s_tau_w * pow((init[0]-s_V_b),3)) * dt;
    else                   x[3] += (-x[3]/s_tau_w) * dt;                                                     // adaptation variable
}

// void integrator_iaf_cond_exp(const state_type_iaf_cond_exp &x, state_type_iaf_cond_exp &dxdt, const double t,                  \
//         double s_C_m, double s_g_L, double s_E_L, double s_E_ex, double s_E_in, double s_I_e, double s_osc_amp, double s_osc_omega,   \
//         double s_tau_syn_ex, double s_tau_syn_in) {
//     dxdt[0] = 1/s_C_m * ( -s_g_L*(x[0]-s_E_L) - x[1]*(x[0]-s_E_ex) - x[2]*(x[0]-s_E_in) + s_I_e \
//                           + s_osc_amp*sin(s_osc_omega*t) );                                             // membrain potential
//     dxdt[1] = -x[1]/s_tau_syn_ex;                                                                       // excitatory syn conductance
//     dxdt[2] = -x[2]/s_tau_syn_in;                                                                       // inhibitory syn conductance
// }

void euler_step_iaf_cond_exp(vector<double> &x, double t, double dt, \
        double s_C_m, double s_g_L, double s_E_L, double s_E_ex, double s_E_in, double s_I_e, double s_osc_amp, double s_osc_omega,   \
        double s_tau_syn_ex, double s_tau_syn_in ) {

    x[0] += 1/s_C_m * ( -s_g_L*(x[0]-s_E_L) - x[1]*(x[0]-s_E_ex) - x[2]*(x[0]-s_E_in) + s_I_e \
                          + s_osc_amp*sin(s_osc_omega*t) ) * dt ;                                        // membrain potential
    x[1] += (-x[1]/s_tau_syn_ex) * dt;                                                                   // excitatory syn conductance
    x[2] += (-x[2]/s_tau_syn_in) * dt;                                                                   // inhibitory syn conductance

}

void euler_step_iaf_curr_exp(vector<double> &x, double t, double dt, \
        double s_C_m, double s_g_L, double s_E_L, double s_I_e, double s_osc_amp, double s_osc_omega,   \
        double s_tau_syn_ex, double s_tau_syn_in ) {

    x[0] += 1/s_C_m * ( -s_g_L*(x[0]-s_E_L) + x[1] - x[2] + s_I_e \
                          + s_osc_amp*sin(s_osc_omega*t) ) * dt ;                                        // membrain potential
    x[1] += (-x[1]/s_tau_syn_ex) * dt;                                                                   // excitatory syn current
    x[2] += (-x[2]/s_tau_syn_in) * dt;                                                                   // inhibitory syn current

}

/// parrot neuron evolve
void euler_step_parrot_neuron(vector<double> &x, double t, double dt, \
        double v_thr, vector<double> &parrot_spikes ) {
    if (parrot_spikes.size()>0) {
        if (parrot_spikes[0] < t+dt) {
            // register a spike here
            x[0] = v_thr+1;                                                                             // membrain potential
            parrot_spikes.erase(parrot_spikes.begin());
            x[1] += 1;      // count the spikes...
        }
    }
    // x[1] = 0;                                                                                           // excitatory syn current  (not used)
    // x[2] = 0;                                                                                           // inhibitory syn current  (not used)
}

void Network::evolve(double _T) {

    double          V_before, w_before;
    double          t0=t, I_ext;
    unsigned        j, i_t=0, i_current=0;
    bool            save_flag = false;
    ofstream        t_save;
    unsigned        st_dim, ind0, sub_ind, neig_size;    // ind0 is useful to keep trace of an index (must be always initialized to some constant)
    map<unsigned, vector<ofstream>>    to_save_files;

    // check and handle neurons to be saved
    for (auto k : subnets) {
        if (k.to_save.size() > 0) {
            save_flag = true;
            break;
        }
    }
    if (save_flag) {

        if (t < dt) {
        // first initialization of save
            mkdir((out_dir + "/neuron_states").c_str(), 0755);
        }
        t_save.open((out_dir + "/neuron_states/t.txt").c_str(), ios::app);

        for (auto k : subnets) {
            if (k.to_save.size() > 0) {
                for (unsigned u : k.to_save) {
                    to_save_files[k.id_subnet].push_back( ofstream((out_dir+"/neuron_states/" + k.name + "_" + to_string(u)+".txt").c_str() , ios::app) );
                }
            }
        }
    }


    while (t<t0+_T){

        // output percentage
        if ( fmod( t/t_end *100, 1) < dt/t_end * 99 ) {
            // cout.precision(20);
            cout  << round( t/t_end *100 ) <<" %   t = " << t <<endl;
        }

        externalInputUpdate();

        for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        // for subnet in subnets

            for (unsigned i=0; i<k->N; i++){
                // for neuron in the subnetwork

                // integration step execution
                I_ext = ((k->pop)[i].specific_I_e_bool) ? (k->pop)[i].specific_I_e[i_current] : k->I_e;
                I_ext += (k->pop)[i].offset_I_e;

                V_before = (k->pop)[i].x[0];
                w_before = (k->pop)[i].x[3];

                switch(k->id_model) {
                    case 2:
                        euler_step_aeif_cond_exp( (k->pop)[i].x, t, dt, \
                                              k->C_m, k->g_L, k->E_L, k->E_ex, k->E_in, I_ext, k->osc_amp, k->osc_omega, k->tau_syn_ex, k->tau_syn_in, \
                                              k->V_th, k->a_adaptive, k->tau_w_adaptive, k->delta_T_aeif_cond_exp, k->V_peak );
                        break;
                    case 3:
                        euler_step_aeif_curr_exp( (k->pop)[i].x, t, dt, \
                                              k->C_m, k->g_L, k->E_L, I_ext, k->osc_amp, k->osc_omega, k->tau_syn_ex, k->tau_syn_in, \
                                              k->V_th, k->a_adaptive, k->tau_w_adaptive, k->delta_T_aeif_cond_exp );
                        break;
                    case 4:
                        euler_step_aqif_cond_exp( (k->pop)[i].x, t, dt, \
                                                  k->C_m, k->k_aqif_cond_exp, k->E_L, k->E_ex, k->E_in, k->I_e, k->osc_amp, k->osc_omega, k->tau_syn_ex, k->tau_syn_in, \
                                                  k->V_th, k->a_adaptive, k->tau_w_adaptive );
                        break;
                    case 5:
                        euler_step_aqif_curr_exp( (k->pop)[i].x, t, dt, \
                                                  k->C_m, k->k_aqif_cond_exp, k->E_L, k->I_e, k->osc_amp, k->osc_omega, k->tau_syn_ex, k->tau_syn_in, \
                                                  k->V_th, k->a_adaptive, k->tau_w_adaptive );
                        break;
                    case 8:
                        euler_step_aqif2_cond_exp( (k->pop)[i].x, t, dt, \
                                                   k->C_m, k->k_aqif_cond_exp, k->E_L, k->E_ex, k->E_in, k->I_e, k->osc_amp, k->osc_omega, k->tau_syn_ex, k->tau_syn_in, \
                                                   k->V_th, k->a_adaptive, k->tau_w_adaptive, k->V_b_aqif2_cond_exp );
                        break;
                    case 6:
                        euler_step_iaf_cond_exp( (k->pop)[i].x, t, dt, \
                                                  k->C_m, k->g_L, k->E_L, k->E_ex, k->E_in, I_ext, k->osc_amp, k->osc_omega, k->tau_syn_ex, k->tau_syn_in );
                        break;
                    case 7:
                        euler_step_iaf_curr_exp( (k->pop)[i].x, t, dt, \
                                                  k->C_m, k->g_L, k->E_L, I_ext, k->osc_amp, k->osc_omega, k->tau_syn_ex, k->tau_syn_in );
                        break;
                    case 0:
                        euler_step_iaf_cond_alpha( (k->pop)[i].x, t, dt, \
                                                   k->C_m, k->g_L, k->E_L, k->E_ex, k->E_in, k->I_e, k->osc_amp, k->osc_omega, k->tau_syn_ex, k->tau_syn_in );
                        break;
                    case 1:
                        euler_step_iaf_curr_alpha( (k->pop)[i].x, t, dt, \
                                                   k->C_m, k->g_L, k->E_L, k->I_e, k->osc_amp, k->osc_omega, k->tau_syn_ex, k->tau_syn_in );
                        break;

                    case 9:
                        euler_step_parrot_neuron( (k->pop)[i].x, t, dt, k->V_peak, (k->pop)[i].parrot_spikes );
                        break;
                    // case :
                    //
                    //     break;
                    default:
                      cerr << "Error: non compatible neuron model: " << k->id_model << endl;
                      exit(1);
                  }

                // /// TEMPORARY modification
                // for (unsigned d=0; d<(k->pop[i]).dim; d++) {
                //     (k->pop[i]).x[d] += 2 * (g.getRandomUniform()-0.5);
                // }

                /// XXX correct here and above; insert it in the euler_integration functions; consider the possibility of passing init also in order to not initialize it...
                if (isnan((k->pop)[i].x[0]) || abs(V_before - (k->pop)[i].x[0]) > 150. ) {
                    // if (k->id_model == 2 || k->id_model == 3 || k->id_model ==4) {
                    if (true) {
                        cerr << "\t Error: x[0] is nan OR single step increase of x[0] was bigger than 150 mV (Pleae decrease dt)" << endl; // zzz
                        cerr << k->name << "  " << i << "  " << t << endl;
                        cerr << V_before << "   " << (k->pop)[i].x[0] << endl;
                        exit(1);
                    }
                    (k->pop)[i].x[0] = k->V_peak+EPSILON;
                    (k->pop)[i].x[3] = w_before + dt/k->tau_w_adaptive * (k->a_adaptive * (k->V_peak-k->E_L) - w_before);
                }

                for (auto item : (k->pop)[i].input_t_ex) {
                    if ((k->pop)[i].next_sp_ex_index[item.first] == item.second.size()) continue;
                    // questo caso sopra forse inutile, ma non fa male

                    if (t > item.second[(k->pop)[i].next_sp_ex_index[item.first]]) {
                        // zzz
                        cerr << showpoint << "Error: simulation time t > next_spike_ex pulse...: " << k->name << " <-- " << subnets[item.first].name << "  " << t << " < " << item.second[(k->pop)[i].next_sp_ex_index[item.first]] << endl;
                        exit(1);
                    }

                    j=0;
                    ind0 = (k->pop)[i].next_sp_ex_index[item.first];

                    // external input case
                    if (item.first == UINT_MAX) {
                        // while ( ind0 + j < item.second.size() && t <= item.second[ ind0 + j ] && item.second[ ind0 + j ] < t+dt ) {       // t+dt substitution
                        while ( ind0 + j < item.second.size() && t <= item.second[ ind0 + j ] && item.second[ ind0 + j ] < t0+(i_t+1)*dt ) {
                            (k->pop)[i].x[1] += (k->pop)[i].ext_weight;
                            j += 1;
                        }
                    }
                    // NON external input case
                    else {
                        // while ( ind0 + j < item.second.size() && t <= item.second[ ind0 + j ] && item.second[ ind0 + j ] < t+dt ) {       // t+dt substitution
                        while ( ind0 + j < item.second.size() && t <= item.second[ ind0 + j ] && item.second[ ind0 + j ] < t0+(i_t+1)*dt ) {

                            if (k->specific_in_weight[item.first]){
                                // specific input weight
                                // eps>1
                                if ( subnets[ item.first ].reverse_effect[ k->id_subnet ] ) {
                                    cerr << "eps>1... " << endl;
                                    exit(1);
                                    (k->pop)[i].x[1] -= (k->pop)[i].input_w_ex[ item.first ][ ind0 + j ];
                                }
                                else {
                                    (k->pop)[i].x[1] += (k->pop)[i].input_w_ex[ item.first ][ ind0 + j ];
                                }
                            }
                            else{
                                // general input weight
                                // eps>1
                                if ( subnets[ item.first ].reverse_effect[ k->id_subnet ] ) {
                                    cerr << "eps>1... " << endl;
                                    exit(1);
                                    (k->pop)[i].x[1] -= k->weights_ex[item.first];
                                }
                                else {
                                    (k->pop)[i].x[1] += k->weights_ex[item.first];
                                }
                            }
                            j += 1;
                        }
                    }

                    (k->pop)[i].next_sp_ex_index[item.first] += j;

                    // adjust all cases : potresti usare j al posto di next_spike_index
                    if ((k->pop)[i].next_sp_ex_index[item.first] != 0) {
                        (k->pop)[i].input_t_ex[item.first].erase( (k->pop)[i].input_t_ex[item.first].begin(), \
                                        (k->pop)[i].input_t_ex[item.first].begin()+(k->pop)[i].next_sp_ex_index[item.first]);
                        if (k->specific_in_weight[item.first]){
                            (k->pop)[i].input_w_ex[item.first].erase( (k->pop)[i].input_w_ex[item.first].begin(), \
                                            (k->pop)[i].input_w_ex[item.first].begin()+(k->pop)[i].next_sp_ex_index[item.first]);
                        }
                        (k->pop)[i].next_sp_ex_index[item.first] = 0;
                    }
                }

                // // adjust external input case
                // if ((k->pop)[i].next_sp_ex_index[UINT_MAX] != 0) {
                //     (k->pop)[i].input_t_ex[UINT_MAX].erase( (k->pop)[i].input_t_ex[UINT_MAX].begin(), \
                //                     (k->pop)[i].input_t_ex[UINT_MAX].begin()+(k->pop)[i].next_sp_ex_index[UINT_MAX]);
                //     (k->pop)[i].next_sp_ex_index[UINT_MAX] = 0;
                // }

                for (auto item : (k->pop)[i].input_t_in) {
                    if ((k->pop)[i].next_sp_in_index[item.first] == item.second.size()) continue;

                    if (t > item.second[(k->pop)[i].next_sp_in_index[item.first]]) {
                        // zzz
                        cerr << showpoint << "Error: simulation time t > next_spike_in pulse...: " << k->name << " <-- " << subnets[item.first].name << "  " << t << " < " << item.second[(k->pop)[i].next_sp_in_index[item.first]] << endl;
                        exit(1);
                    }

                    j=0;
                    ind0 = (k->pop)[i].next_sp_in_index[item.first];
                    // while ( ind0 + j < item.second.size() && t <= item.second[ ind0 + j ] && item.second[ ind0 + j ] < t+dt ) {          // t+dt substitution
                    while ( ind0 + j < item.second.size() && t <= item.second[ ind0 + j ] && item.second[ ind0 + j ] < t0+(i_t+1)*dt ) {

                        if (k->specific_in_weight[item.first]){
                            // specific input weight
                            // eps>1
                            if ( subnets[ item.first ].reverse_effect[ k->id_subnet ] ) {
                                cerr << "eps>1... " << endl;
                                exit(1);
                                (k->pop)[i].x[2] += (k->pop)[i].input_w_in[ item.first ][ ind0 + j ];
                            }
                            else {
                                (k->pop)[i].x[2] -= (k->pop)[i].input_w_in[ item.first ][ ind0 + j ];
                            }
                        }
                        else{
                            // general input weight
                            // eps>1
                            if ( subnets[ item.first ].reverse_effect[ k->id_subnet ] ) {
                                cerr << "eps>1... " << endl;
                                exit(1);
                                (k->pop)[i].x[2] -= k->weights_in[item.first];
                            }
                            else {
                                (k->pop)[i].x[2] += k->weights_in[item.first];
                            }
                        }
                        j += 1;
                    }

                    (k->pop)[i].next_sp_in_index[item.first] += j;

                    // adjust all cases : potresti usare j al posto di next_spike_index
                    if ((k->pop)[i].next_sp_in_index[item.first] != 0) {
                        (k->pop)[i].input_t_in[item.first].erase( (k->pop)[i].input_t_in[item.first].begin(), \
                                        (k->pop)[i].input_t_in[item.first].begin()+(k->pop)[i].next_sp_in_index[item.first]);
                        if (k->specific_in_weight[item.first]){
                            (k->pop)[i].input_w_in[item.first].erase( (k->pop)[i].input_w_in[item.first].begin(), \
                                            (k->pop)[i].input_w_in[item.first].begin()+(k->pop)[i].next_sp_in_index[item.first] );
                        }
                        (k->pop)[i].next_sp_in_index[item.first] = 0;
                    }

                }

                if ( t-(k->pop)[i].t_last_spike < k->t_ref ) {
                    (k->pop)[i].x[0] = k->V_res;
                }

                // neuron spike
                if ( (k->pop)[i].x[0] > k->V_peak ) {

                    (k->pop)[i].x[0] = k->V_res;
                    (k->pop)[i].t_last_spike = t;
                    (k->pop)[i].t_spikes.push_back(t);

                    // update adaptation variable x[3], if necessary
                    if (k->adap_4) (k->pop)[i].x[3] += k->b_adaptive;

                    // add input spike to neighbors with proper delay
                    for (auto item : (k->pop)[i].neighbors) {
                        neig_size = item.second.size();

                        if (k->specific_out_weight[item.first]) {
                            // Specific out weight
                            for (unsigned i_n=0; i_n<neig_size; i_n++){
                                if ((k->pop)[i].neighbors_out_weights[item.first][i_n] > 0) {
                                    // Excitatory
                                    subnets[ item.first ].pop[ item.second[i_n] ].input_t_ex[ k->id_subnet ].push_back(t + k->delays_out[item.first] + EPSILON*dt);
                                    subnets[ item.first ].pop[ item.second[i_n] ].input_w_ex[ k->id_subnet ].push_back( (k->pop)[i].neighbors_out_weights[item.first][i_n] );
                                }
                                else {
                                    // Inhibitory
                                    subnets[ item.first ].pop[ item.second[i_n] ].input_t_in[ k->id_subnet ].push_back(t + k->delays_out[item.first] + EPSILON*dt);
                                    subnets[ item.first ].pop[ item.second[i_n] ].input_w_in[ k->id_subnet ].push_back( (k->pop)[i].neighbors_out_weights[item.first][i_n] );
                                }
                            }
                        }
                        else {
                            // General out weight
                            if (k->excit_out[item.first]) {
                                // Excitatory
                                for (unsigned i_n=0; i_n<neig_size; i_n++){
                                    subnets[ item.first ].pop[ item.second[i_n] ].input_t_ex[ k->id_subnet ].push_back(t + k->delays_out[item.first] + EPSILON*dt);
                                }
                            }
                            else {
                                // Inhibitory
                                for (unsigned i_n=0; i_n<neig_size; i_n++){
                                    subnets[ item.first ].pop[ item.second[i_n] ].input_t_in[ k->id_subnet ].push_back(t + k->delays_out[item.first] + EPSILON*dt);
                                }
                            }
                        }
                    }
                } // end if V>V_peak

            } // end for over pop
        } // end for on subnets

        // update outfiles with current state
        // double     phi = M_PI;
        // if ( fmod(t*subnets[0].osc_omega, 2*M_PI)<phi+M_PI/48 && fmod(t*subnets[0].osc_omega, 2*M_PI)>phi-M_PI/48 ) {
        if (true) {
            if (save_flag) {
                t_save << t << endl;
                for (auto item_it = to_save_files.begin(); item_it!=to_save_files.end(); item_it++) {
                    sub_ind = item_it->first;
                    st_dim = subnets[ sub_ind ].pop[0].dim;
                    ind0 = 0;
                    for (auto f_it = item_it->second.begin(); f_it != item_it->second.end(); f_it++) {
                        for (unsigned s=0; s<st_dim; s++) {
                            (*f_it) << subnets[sub_ind].pop[ subnets[sub_ind].to_save[ind0] ].x[s] << "  ";
                        }
                        (*f_it) << endl;
                        ind0 += 1;
                    }
                }
            } // end if save_flag

        }     // end save only in phase phi

        i_t += 1;
        // t += dt;         // t+dt substitution
        t = t0 + i_t*dt;
        if (i_t % repeat_specific_I_e == 0) {   // note: you have already increased i_t !
            i_current += 1;
        }
    } // end while over t

    if (save_flag) {
        for (auto item_it = to_save_files.begin(); item_it!=to_save_files.end(); item_it++) {
            for (auto f_it = item_it->second.begin(); f_it != item_it->second.end(); f_it++) {
                (*f_it).close();
            }
        }
    }

    for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        for (unsigned i=0; i<k->N; i++) {
            if ((k->pop)[i].specific_I_e_bool) {
                (k->pop)[i].specific_I_e.erase( (k->pop)[i].specific_I_e.begin(), (k->pop)[i].specific_I_e.begin()+i_current );
            }
        }
    }

    for(auto k : subnets) {
        k.save_t_spikes(out_dir+"/"+k.name+"_spikes.txt");
    }
}



void Network::externalInputUpdate() {

    double           tmp_sum, tmp_end, tbar, tmp_y;

    if (input_mode == 0) {        // base mode
        for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        // for subnet in subnets
            if (k->ext_in_rate<EPSILON) continue;

            for (unsigned i=0; i<k->N; i++) {
            // for neuron in pop
                tmp_end = *((k->pop)[i].input_t_ex[UINT_MAX].end()-1);
                if ( tmp_end < t+1.00001*dt ) {
                    tmp_sum = 0;
                    if (k->osc_amp_poiss < EPSILON) {
                        while (tmp_sum < 1.00001*dt) {
                            tmp_sum += (- log(g.getRandomUniform()) ) / k->ext_in_rate;
                            (k->pop)[i].input_t_ex[UINT_MAX].push_back(tmp_end + tmp_sum);
                        }
                    }
                    else {
                        // ofstream        of("prova_osc_inp.txt", ios::app);
                        tbar = tmp_end;
                        while (tmp_sum < 1.00001*dt) {
                            tmp_y = - log(g.getRandomUniform());
                            tbar = find_sol_bisection( tmp_y, k->ext_in_rate, k->osc_amp_poiss, k->osc_omega_poiss, tbar );
                            (k->pop)[i].input_t_ex[UINT_MAX].push_back(tbar);
                            tmp_sum = tbar - tmp_end;
                            // of << tbar << endl;
                        }
                        // of.close();
                    }
                }
            } // end for over neurons
        } // end for over subnets
    } // end of input_mode == 0 (base mode)

    if (input_mode == 2) {        // with_correlation mode A;  oscillatory ext_in_rate not supported
        for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        // for subnet in subnets

            if (k->ext_in_rate < EPSILON) continue;

            // non striatum pop
            if ( find( corr_pops.begin(), corr_pops.end(), k-subnets.begin() ) == corr_pops.end() ) {
                for (unsigned i=0; i<k->N; i++) {
                // for each neuron in the population
                    tmp_end = *((k->pop)[i].input_t_ex[UINT_MAX].end()-1);
                    if ( tmp_end < t+1.00001*dt ) {
                        tmp_sum = 0.;
                        while (tmp_sum < 1.00001*dt) {
                            tmp_sum += (- log(g.getRandomUniform()) ) / k->ext_in_rate;
                            (k->pop)[i].input_t_ex[UINT_MAX].push_back(tmp_end + tmp_sum);
                        }
                    }
                }
            }

            // striatum pop
            else {
                if (corr_last_time[k->id_subnet] < t+1.00001*dt){
                    tmp_sum = 0;
                    tmp_end = corr_last_time[k->id_subnet];
                    while (tmp_sum < 1.00001*dt) {
                        tmp_sum += (- log(g.getRandomUniform()) ) / k->ext_in_rate * rho_corr;
                        for (unsigned i=0; i<k->N; i++) {
                            if (g.getRandomUniform() < rho_corr) {
                                (k->pop)[i].input_t_ex[UINT_MAX].push_back(tmp_end + tmp_sum);
                            }
                        }
                    }
                    corr_last_time[k->id_subnet] = tmp_end + tmp_sum;
                }
            }
        } // end for over subnets
    } // end of input_mode == 2 (with_correlation mode A)
}


void Network::free_past() {

    int         next_spike_ind=0;

    for (auto k=subnets.begin(); k!=subnets.end(); k++) {
    // for subnet in subnets
        for (unsigned i=0; i<k->N; i++){
        // for each neuron in the subnet population

            // perform a control operation on the values of the neuron neuron_states
            for (unsigned s=0; s<(k->pop)[i].dim; s++) {
                switch (fpclassify((k->pop)[i].x[s])) {
                    case FP_INFINITE:   {
                        cerr << "\tERROR: infinite number met at t=" << t << "("<< k->name << ", neuron " << i<< ", st " << s <<  ")" << endl;
                        exit(1);
                    }
                    case FP_NAN:        {
                        cerr << "\tERROR: NaN met at t=" << t << "("<< k->name << ", neuron " << i<< ", st " << s <<  ")" << endl;
                        exit(1);
                    }
                    case FP_SUBNORMAL:  {
                        (k->pop)[i].x[s] = 0.;
                    }
                }
            }

            // free t_spikes
            (k->pop)[i].t_spikes.clear();

            // free_past of input_t_ex vectors
            for (auto item_it=(k->pop)[i].input_t_ex.begin(); item_it!=(k->pop)[i].input_t_ex.end(); item_it++) {
                next_spike_ind = (k->pop)[i].next_sp_ex_index[item_it->first];
                if (next_spike_ind != 0) {
                    // zzz cout << "\t\t" << item_it->first << " pass from: " << item_it->second.size() << " to: " ;
                    (item_it->second).erase( (item_it->second).begin(), (item_it->second).begin() + next_spike_ind );
                    if (k->specific_in_weight[ item_it->first ]) {
                        (k->pop)[i].input_w_ex[ item_it->first ].erase( (k->pop)[i].input_w_ex[item_it->first].begin(), (k->pop)[i].input_w_ex[item_it->first].begin()+next_spike_ind );
                    }
                    // note: if next_sp_ex_index[item_it->first] == input_t_ex[item_it->first], this correctly returns the empty vector
                    (k->pop)[i].next_sp_ex_index[item_it->first] = 0;
                    // zzz cout << item_it->second.size() << endl;
                }
            }
            // free_past of input_t_in vectors
            for (auto item_it=(k->pop)[i].input_t_in.begin(); item_it!=(k->pop)[i].input_t_in.end(); item_it++) {
                next_spike_ind = (k->pop)[i].next_sp_in_index[item_it->first];
                if (next_spike_ind != 0) {
                    (item_it->second).erase( (item_it->second).begin(), (item_it->second).begin() + next_spike_ind );
                    // note: if next_sp_ex_index[item_it->first] == input_t_ex[item_it->first], this correctly returns the empty vector
                    if (k->specific_in_weight[ item_it->first ]) {
                        (k->pop)[i].input_w_in[ item_it->first ].erase( (k->pop)[i].input_w_in[item_it->first].begin(), (k->pop)[i].input_w_in[item_it->first].begin()+next_spike_ind );
                    }
                    (k->pop)[i].next_sp_in_index[item_it->first] = 0;
                }
            }
        } // end for over pop
    } // end for over subnets
}

void count_spikes(vector<unsigned> & v_out, vector<double> & v_t, unsigned lenght, double step_res, double half_w, double dur, double tt0 ) {

    unsigned        idx_temp=0, temp=0;
    double          tt;

    for (unsigned j=0; j<lenght; j++) {
        tt = tt0+ step_res*j;
        temp=0;

        for (auto it=v_t.begin()+idx_temp; it!=v_t.end(); it++) {

            if ( (step_res*j)-half_w<0 or (step_res*j)+half_w>dur ) {
                v_out[j] = 0;
            }
            else if ( *it < tt0+step_res*(j-1)-half_w ) {         // XXX this is just a check
                cerr << "Error during training: this should never happen..." << endl;
                exit(1);
            }
            else if ( *it < tt-half_w ) {     // spike before the j-th window
                idx_temp += 1;
            }
            else if ( *it < tt+half_w ) {     // spike in the j-th window
                temp += 1;
            }
            else {                             // spike after the j-th window
                break;
            }
        }
        v_out[j] = temp;
    }
}


unsigned Network::compute_gen_error( int t_order, unsigned lenght, double dur ) {
    unsigned                    computed_errors=0;
    vector<double>::iterator    iterator;
    vector<double>              zero_error(lenght, 0);
    double                      sum, sum_2;
    double                      dF_dI;
    // // XXXA
    // cout << "\tcompute_gen_error " << t_order << endl;
    for (unsigned k=0; k<subnets.size(); k++) {

        if ( subnets[k].train_order == t_order ) {
            // // XXXA
            // cout << "\t\t " << subnets[k].name << "  " << subnets[k].train_order << endl;

            // compute the error of the OUTPUT LAYER
            if ( k==ord_target_subnets[0] ) {
                computed_errors += 1;
                subnets[k].errors_computed = true;

                ofstream                    error_file( (out_dir + "/errors/" + "e_"+to_string(current_epoch)+".txt").c_str(), ios::app );
                double                      norm_factor=pow(2,POT-1), factor, c_des;
                unsigned                    output_pop_id;
                output_pop_id = ord_target_subnets[0];

                for (unsigned i=0; i<subnets[output_pop_id].N; i++) {
                    // for each target neuron i
                    iterator = (subnets[output_pop_id].pop)[i].desidered_output.begin();
                    if ((subnets[output_pop_id].pop)[i].desidered_output.size() < lenght) {
                        cerr << "Error: desidered_output of neuron " << i << "is less than lenght..." << endl;
                        exit(1);
                    }
                    (subnets[output_pop_id].pop)[i].error = zero_error;

                    sum = 0;
                    sum_2 = 0;
                    // compute E[desidered, actual]      (actual=(subnets[output_pop_id].pop)[i].R_actual)
                    // /// XXXA
                    // ofstream                    time_error_file( (out_dir + "/errors/" + "time_e_"+to_string(current_epoch)+".txt").c_str(), ios::app );

                    c_des = (c_desid_zero<0) ? (subnets[output_pop_id].pop)[i].optimal_c_desid_zero : c_desid_zero;
                    for (unsigned j=0; j<lenght; j++) {
                        // there was a mistake!
                        if (abs( floor( *(iterator+j)+1./3-1./(50*(*(iterator+j))) ) - (subnets[output_pop_id].pop)[i].R_actual[j] ) < 0.5) {
                            (subnets[output_pop_id].pop)[i].error[j] = 0;
                        }
                        else {
                            /// MODIFICATION: modulation of external input
                            // factor = (1-c)/( 1+ pow(M_E, -a*( (*(iterator+j))-b )) ) + c;
                            factor = 1. - (1.-c_des)*pow( M_E, -a_desid_zero* (*(iterator+j)) );
                            (subnets[output_pop_id].pop)[i].error[j] = - norm_factor*pow(0.5 - poissonian_CDF( (subnets[output_pop_id].pop)[i].R_actual[j], *(iterator+j) ), POT) * ( factor )  + error_lin_coeff * ( (subnets[output_pop_id].pop)[i].R_actual[j] - *(iterator+j)  ) ;

                            // (subnets[output_pop_id].pop)[i].error[j] = - norm_factor*pow(0.5 - poissonian_CDF( (subnets[output_pop_id].pop)[i].R_actual[j], *(iterator+j) ), POT) + error_lin_coeff * ( (subnets[output_pop_id].pop)[i].R_actual[j] - *(iterator+j)  ) ;
                        }
                        sum += abs( (subnets[output_pop_id].pop)[i].error[j] );
                        sum_2+= abs( (subnets[output_pop_id].pop)[i].R_actual[j] - (*(iterator+j)) );
                        // /// XXXA
                        // time_error_file << sum_2  << " " << (subnets[output_pop_id].pop[i]).R_actual[j] << " "<< *(iterator+j) << " " << (subnets[output_pop_id].pop)[i].error[j] << endl;
                    }
                    error_file << i << " " << sum/dur*1000 << " " << sum_2 << endl;

                    (subnets[k].pop)[i].desidered_output.erase( iterator, iterator+lenght );
                }
                error_file.close();
            } // end if output pop

            // compute the "generalized" error for hidden layers  (no output layer)
            // generalized_error_i = sum_(postsyn_neur_k) dF_dI_app*tau_ki w_ki E_k
            // note that multiple postsyn populations can be admitted BUT the weights to all of them MUST be trainable!
            // you are using the all to all assumption here
            else {
                computed_errors += 1;
                subnets[k].errors_computed = true;

                unsigned     lp1_pop_id;
                double       w_dest_source, tau;

                for (unsigned i=0; i<subnets[k].N; i++) {
                    (subnets[k].pop)[i].error = zero_error;
                }

                for (auto item : subnets[k].train_out_weight){
                    if (item.second==false) {
                        cerr << "computing generalized error with non-trainable connections to post_syn population...\n\tthis is not necessarily an error but shouldn't happen with a feedforward architecture.";
                        exit(1);
                        // continue;
                    }
                    lp1_pop_id = item.first;
                    for (unsigned i=0; i<subnets[k].N; i++) {
                        for (unsigned i_lp1=0; i_lp1<subnets[lp1_pop_id].N; i_lp1++) {
                            w_dest_source = (subnets[k].pop)[i].neighbors_out_weights[lp1_pop_id][i_lp1];
                            tau = (w_dest_source<0) ? subnets[lp1_pop_id].tau_syn_in : subnets[lp1_pop_id].tau_syn_ex;
                            for (unsigned j=0; j<lenght; j++) {
                                dF_dI = (subnets[lp1_pop_id].pop[i_lp1].R_actual[j]>0) ? dF_dI_app : dF_dI_no_act_fact*dF_dI_app;
                                //                              m    * tau * UPDATED w_ki  *            E_k
                                (subnets[k].pop)[i].error[j] += dF_dI* tau * w_dest_source * (subnets[lp1_pop_id].pop)[i_lp1].error[j];
                            }
                        }
                    }
                }
            }
            // // compute the "generalized" error of the LAYER L-1
            // // generalized_error = sum_o w_oh E_o : assuming that each presyn neuron is connected to all post_syn neurons
            // else if ( subnets[k].id_subnet==one_pop_id ) {
            //     for (unsigned i=0; i<subnets[one_pop_id].N; i++) {
            //         (subnets[one_pop_id].pop)[i].error = zero_error;
            //         for (unsigned j=0; j<lenght; j++) {
            //             for (unsigned i_out=0; i_out<subnets[output_pop_id].N; i_out++) {
            //                 //                                                          UPDATED w_oh                            *   E_o
            //                 (subnets[one_pop_id].pop)[i].error[j] +=  (subnets[one_pop_id].pop)[i].neighbors_out_weights[output_pop_id][i_out] * (subnets[output_pop_id].pop)[i_out].error[j];
            //             }
            //         }
            //     }
            // } // end if pop_id==one_pop_id
        }
    } // end for over subnets
    return computed_errors;
}


double Network::get_quant_weight_delta( double w_quant, double w_cum ) {
    double quant_weight_delta=0;
    if (w_quant<EPSILON) {
        quant_weight_delta = (w_cum>0) ? quant_delta_weight_pos : quant_delta_weight_neg;
    }
    else {
        quant_weight_delta = (w_quant>0) ? quant_delta_weight_pos : quant_delta_weight_neg;
    }
    return quant_weight_delta;
}


/// using original weights to compute the errors of hidden layers (correct gradient descent)
void Network::train( double dur ) {         // dur:     batch duration

    unsigned                    lenght=floor(dur/step_res);
    unsigned                    output_pop_id, current_t_order, computed_errors;
    vector<unsigned>            zero_actual(lenght, 0);
    double                      half_w = window_res/2;
    double                      tt0=t-dur, sum;  // tt0: batch initial time
    bool                        compute_R_actual;
    double                      w_dest_source, tau, dF_dI, temp_epsilon;
    double                      neuron_mean_activity, decay_factor;
    double                      quant_weight_delta, quant_sign;

    output_pop_id = ord_target_subnets[0];

    for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        k->errors_computed = false;
    }

    // compute the actual activities of each population, if necessary (output pop and pops with outgoing weights to be trained)
    for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        compute_R_actual = false;

        if (k->id_subnet == output_pop_id){
            // output layer
            compute_R_actual=true;
        }
        else {
            for (auto i : k->train_out_weight) {
                if (i.second==true) {
                    compute_R_actual=true;
                    break;
                }
            }
        }

        if (compute_R_actual) {
            for (unsigned i=0; i<k->N; i++) {
                (k->pop)[i].R_actual=zero_actual;
                count_spikes( (k->pop)[i].R_actual, (k->pop)[i].t_spikes, lenght, step_res, half_w, dur, tt0 );
            }
        }
    }

    // compute errors
    current_t_order = 0;
    for (auto k : ord_target_subnets) {
        // cout << current_t_order << "  " << k << endl;
        if (k==-1){
            break;
        }
        else {
            // compute error for the current train order
            computed_errors = compute_gen_error( current_t_order, lenght, dur );
            if (computed_errors != 1) {
                cerr << "Error: computed errors = " <<computed_errors<< " !=1 during current_t_order = " << current_t_order <<  endl;
                exit(1);
            }
            current_t_order += 1;
        }
    }


    // update weights
    // for each population having an incoming trainable connection
    current_t_order = 0;
    for (auto k : ord_target_subnets) {
        // cout << current_t_order << "  " << k << endl;
        if (k==-1){
            break;
        }
        else {

            // for each presynaptic population with trainable weights
            for (auto id_it=subnets[k].specific_in_weight.begin(); id_it!=subnets[k].specific_in_weight.end(); id_it++) {
                if (id_it->first==UINT_MAX) {
                    continue;
                }
                if (subnets[id_it->first].train_out_weight[subnets[k].id_subnet]) {

                    vector<double>    TEMP_dw(subnets[id_it->first].N,0);
                    ofstream          weight_file( (out_dir + "/specific_weights/" + \
                                                    "w_"+to_string(current_epoch)+"_"+subnets[id_it->first].name+"_to_"+subnets[k].name+".txt").c_str(), ios::app );

                    // for each target neuron
                    for (unsigned i=0; i<subnets[k].N; i++) {

                        if ( subnets[k].decay_factor_max > EPSILON ) {
                            neuron_mean_activity = ((subnets[k].pop)[i].t_spikes).size()/dur;       // in kilo hertz...
                            decay_factor = 1. + subnets[k].decay_factor_max / ( 1+pow(M_E, -(neuron_mean_activity-subnets[k].decay_factor_position)/subnets[k].decay_factor_steepness) );
                            // cout << subnets[k].name << "  " << neuron_mean_activity << "   " << decay_factor << endl;
                        }
                        else {
                            decay_factor = 1.;
                        }

                        // for each source neuron
                        for (unsigned i_pre=0; i_pre<subnets[id_it->first].N; i_pre++) {

                            if (current_epoch==0) {
                                // do not update weight
                                TEMP_dw[i_pre] = subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i];
                            }
                            else {
                                // update the i_pre -> i weight using the input_activity and generalized error terms: assuming that each presyn neuron is connected to all post_syn neurons
                                sum = 0;
                                w_dest_source = subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i];
                                tau = (w_dest_source<0) ? subnets[k].tau_syn_in : subnets[k].tau_syn_ex;
                                for (unsigned j=0; j<lenght; j++) {
                                    dF_dI = (subnets[k].pop[i].R_actual[j]>0) ? dF_dI_app : dF_dI_no_act_fact*dF_dI_app;
                                    // m*tau *      activity of presynaptic neuron      * generalized error of postsynaptic neuron
                                    temp_epsilon = 0. ? subnets[id_it->first].pop[i_pre].R_actual[j]>0 : epsilon_zero_act;
                                    sum += dF_dI * (subnets[id_it->first].pop[i_pre].R_actual[j]+temp_epsilon) * subnets[k].pop[i].error[j];
                                }
                                sum = sum*tau;

                                switch(optimization_method_id) {
                                    case 0:     // SGD   unchanged
                                        if (momentum_w>0.) {
                                            if (t>dur) {
                                                // in case we dont save the momentum in a specific file, at batch zero we use sum directly
                                                sum = ( sum + momentum_w*subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i] );
                                            }
                                            subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i] = sum;
                                        }
                                        break;
                                    case 1:     // SGD_MEAN
                                        if (t-EPSILON<dur) {
                                            // in case we dont save the momentum in a specific file, at batch zero we use sum directly
                                            subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i] = sum;
                                        }
                                        else {
                                            subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i] = (1-momentum_w)*sum + momentum_w*subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i];
                                        }
                                        if (t-EPSILON<dur*DO_NOT_TRAIN_IN_FIRST_BATCHES) {
                                            sum = 0;
                                        }
                                        else {
                                            sum = subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i];
                                        }
                                        break;
                                    case 2:     // ADAM
                                        // if (t-EPSILON<dur) {
                                        //     subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i] = sum;
                                        //     subnets[id_it->first].pop[i_pre].weight_aux2[subnets[k].id_subnet][i] = sum*sum;
                                        // }
                                        // else {
                                        //     subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i] = (1-momentum_w) *sum     + momentum_w *subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i];
                                        //     subnets[id_it->first].pop[i_pre].weight_aux2[subnets[k].id_subnet][i] = (1-momentum2_w)*sum*sum + momentum2_w*subnets[id_it->first].pop[i_pre].weight_aux2[subnets[k].id_subnet][i];
                                        // }
                                        subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i] = (1-momentum_w) *sum     + momentum_w *subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i];
                                        subnets[id_it->first].pop[i_pre].weight_aux2[subnets[k].id_subnet][i] = (1-momentum2_w)*sum*sum + momentum2_w*subnets[id_it->first].pop[i_pre].weight_aux2[subnets[k].id_subnet][i];
                                        sum = subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i] / ( sqrt(subnets[id_it->first].pop[i_pre].weight_aux2[subnets[k].id_subnet][i]) + EPSILON );
                                        // if (t-EPSILON<dur*DO_NOT_TRAIN_IN_FIRST_BATCHES) {
                                        //     sum = 0;
                                        // }
                                        // else {
                                        //     sum = subnets[id_it->first].pop[i_pre].weight_aux1[subnets[k].id_subnet][i] / ( sqrt(subnets[id_it->first].pop[i_pre].weight_aux2[subnets[k].id_subnet][i]) + EPSILON );
                                        // }
                                        break;
                                    default:
                                        cerr << "optimization method id not supported" << endl;
                                        exit(1);
                                 }


                                // momentum

                                // if (momentum_w>0.) {
                                //     if (t>dur) {
                                //         // in case we dont save the momentum in a specific file, at batch zero we use sum directly
                                //         sum = ( sum + momentum_w*subnets[id_it->first].pop[i_pre].neighbors_out_delta_weights[subnets[k].id_subnet][i] );
                                //     }
                                //     subnets[id_it->first].pop[i_pre].neighbors_out_delta_weights[subnets[k].id_subnet][i] = sum;
                                // }

                                sum = -sum*l_rate0*correct_l_rates[current_t_order];

                                // saturation
                                if (sum*w_dest_source>0) {
                                    TEMP_dw[i_pre] = sum \
                                                     * 1./( 1+pow(M_E, abs(subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i])-w_cut) ) \
                                                     - decay_factor*l_decay*subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i];
                                }
                                // NOT saturation
                                else {
                                    TEMP_dw[i_pre] = sum \
                                                     - decay_factor*l_decay*subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i];
                                }
                                /// apply weight update
                                if ((subnets[id_it->first]).alpha_quantized_weights[subnets[k].id_subnet]>0) {
                                    // quantized training
                                    (subnets[id_it->first]).pop[i_pre].weight_quant_cumul[subnets[k].id_subnet][i] += TEMP_dw[i_pre];
                                    quant_weight_delta = get_quant_weight_delta( (subnets[id_it->first]).pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i], \
                                                                                 (subnets[id_it->first]).pop[i_pre].weight_quant_cumul[subnets[k].id_subnet][i] );
                                    if ( abs( (subnets[id_it->first]).pop[i_pre].weight_quant_cumul[subnets[k].id_subnet][i] ) > (subnets[id_it->first]).alpha_quantized_weights[subnets[k].id_subnet]*quant_weight_delta ) {
                                        quant_sign = ( (subnets[id_it->first]).pop[i_pre].weight_quant_cumul[subnets[k].id_subnet][i]>0 ) ? 1 : -1;
                                        subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i] += quant_sign*quant_weight_delta;
                                        subnets[id_it->first].pop[i_pre].weight_quant_cumul[subnets[k].id_subnet][i] = 0;
                                        // XXX to be removed
                                        cout << "change in weight: " << subnets[id_it->first].name << " (" << i_pre <<  ") to " << subnets[k].name << " (" << i <<  ")" << endl;

                                        if (subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i] < -n_quant*quant_delta_weight_neg) {
                                            subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i] = -n_quant*quant_delta_weight_neg;
                                        }
                                        if (subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i] > n_quant*quant_delta_weight_pos) {
                                            subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i] = n_quant*quant_delta_weight_pos;
                                        }

                                    }
                                }
                                else {
                                    // NOT quantized training
                                    subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i] += TEMP_dw[i_pre];
                                }

                                // Noise
                                if (g.getRandomUniform()<prob_noise) {
                                    sum = 1-amp_noise + 2*amp_noise*g.getRandomUniform();
                                    // TEMP_dw[i_pre] = (sum-1)*subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i]+TEMP_dw[i_pre];
                                    subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i] *= sum;
                                }

                                /// MODIFICATION w_min/w_max
                                if (subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i] < subnets[id_it->first].w_min) {
                                    subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i] = subnets[id_it->first].w_min;
                                }
                                if (subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i] > subnets[id_it->first].w_max) {
                                    subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i] = subnets[id_it->first].w_max;
                                }
                                TEMP_dw[i_pre] = subnets[id_it->first].pop[i_pre].neighbors_out_weights[subnets[k].id_subnet][i];
                            }
                        } // end for over presyn neurons

                        // XXX
                        weight_file << i << " ";
                        for (const auto w : TEMP_dw) weight_file << " " << w;
                        weight_file << "\n";

                    } // end for over target neurons
                    weight_file.close();
                }   // end if pre-->post is to be trained
            }  // end for over presyn populations
            current_t_order += 1;
        }
    } // end for over ord_target_subnets

    // update external input currents
    for (auto k=subnets.begin(); k!=subnets.end(); k++) {
        if ( k->train_offset_I_e ) {
            for (unsigned i=0; i<k->N; i++) {
                sum = 0;
                if (current_epoch==0) {
                    sum = 0;
                }
                else {
                    if ( k->errors_computed==false ) {
                        computed_errors = compute_gen_error( -1, lenght, dur );
                        if (computed_errors != 1) {
                            cerr << "Error: computed errors = " <<computed_errors<< " !=1 during current_t_order = -1 (current traiing of input layer)" <<  endl;
                            exit(1);
                        }
                    }
                    for (unsigned j=0; j<lenght; j++) {
                        dF_dI = (k->pop[i].R_actual[j]>0) ? dF_dI_app : dF_dI_no_act_fact*dF_dI_app;
                        sum += dF_dI * k->pop[i].error[j];
                    }

                    switch(optimization_method_id) {
                        case 0:     // SGD   unchanged
                            if (momentum_c>0.) {
                                if (t>dur) {
                                    // in case we dont save the momentum in a specific file, at batch zero we use sum directly
                                    sum = ( sum + momentum_c* (k->pop[i]).offset_aux1 );
                                }
                                (k->pop[i]).offset_aux1 = sum;
                            }
                            break;
                        case 1:     // SGD_MEAN
                            if (t-EPSILON<dur) {
                                // in case we dont save the momentum in a specific file, at batch zero we use sum directly
                                (k->pop[i]).offset_aux1 = sum;
                            }
                            else {
                                (k->pop[i]).offset_aux1 = (1-momentum_c)*sum + momentum_c*(k->pop[i]).offset_aux1;
                            }
                            if (t-EPSILON<dur*DO_NOT_TRAIN_IN_FIRST_BATCHES) {
                                sum = 0;
                            }
                            else {
                                sum = (k->pop[i]).offset_aux1;
                            }
                            break;
                        case 2:     // ADAM
                            // if (t-EPSILON<dur) {
                            //     (k->pop[i]).offset_aux1 = sum;
                            //     (k->pop[i]).offset_aux2 = sum*sum;
                            // }
                            // else {
                            //     (k->pop[i]).offset_aux1 = (1-momentum_c) *sum     + momentum_c *(k->pop[i]).offset_aux1;
                            //     (k->pop[i]).offset_aux2 = (1-momentum2_c)*sum*sum + momentum2_c*(k->pop[i]).offset_aux2;
                            // }
                            (k->pop[i]).offset_aux1 = (1-momentum_c) *sum     + momentum_c *(k->pop[i]).offset_aux1;
                            (k->pop[i]).offset_aux2 = (1-momentum2_c)*sum*sum + momentum2_c*(k->pop[i]).offset_aux2;
                            sum = (k->pop[i]).offset_aux1 / ( sqrt((k->pop[i]).offset_aux2) + EPSILON );
                            // if (t-EPSILON<dur*DO_NOT_TRAIN_IN_FIRST_BATCHES) {
                            //     sum = 0;
                            // }
                            // else {
                            //     sum = (k->pop[i]).offset_aux1 / ( sqrt((k->pop[i]).offset_aux2) + EPSILON );
                            // }
                            break;
                        default:
                            cerr << "optimization method id not supported" << endl;
                            exit(1);
                      }





                    // momentum
                    // if (momentum_c>0.) {
                    //     if (t>dur) {
                    //         // in case we dont save the momentum in a specific file, at batch zero we use sum directly
                    //         sum = ( sum + momentum_c* (k->pop[i]).delta_offset_I_e );
                    //     }
                    //     (k->pop[i]).delta_offset_I_e = sum;
                    // }

                    sum = -sum*l_rate0_curr;
                }
                (k->pop[i]).offset_I_e += sum - ((k->pop[i]).offset_I_e - k->offset_I_e_ref)*(k->ldecay_offset);

                if ( (k->pop[i]).offset_I_e< k->offset_I_e_min ) {
                    (k->pop[i]).offset_I_e = k->offset_I_e_min;
                }
                if ( (k->pop[i]).offset_I_e> k->offset_I_e_max ) {
                    (k->pop[i]).offset_I_e = k->offset_I_e_max;
                }
            }

            ofstream          curr_file( (out_dir + "/specific_weights/" + \
                                            "c_"+to_string(current_epoch)+"_"+k->name+".txt").c_str(), ios::app );
            for (unsigned i=0; i<k->N-1; i++) curr_file << (k->pop[i]).offset_I_e << " ";
            curr_file << (k->pop[k->N-1]).offset_I_e << "\n";
            curr_file.close();
        }
    }

}



void Network::info() {
    cout << "Network info:" << endl;

    cout << "\tt_end \t\t\t" << t_end << endl;
    cout << "\tn_step \t\t\t" << n_step << endl;
    cout << "\tdt (resolut)\t\t" << dt << endl;
    cout << "\tinput mode\t\t" << input_mode << endl;
    if (input_mode == 2) {
        cout << "\t\trho_corr\t" << rho_corr << endl;
        cout << "\t\tcorr_pops\t";
        for (auto i : corr_pops) {
            cout << subnets[i].name << "  ";
        }
        cout << endl;
    }
    cout << "\tsubnets_config_yaml     " << subnets_config_yaml << endl;
    cout << "\tweights_config_yaml     " << weights_config_yaml << endl;
    cout << "\tconnections_config_yaml " << connections_config_yaml << endl;
    cout << "\tto_save_config_yaml     " << to_save_config_yaml << endl;
    cout << "\tout_dir                 " << out_dir << endl;

    cout << "\tsupported neuron models\t";
    for (auto i : subnet_model_id) cout << i.first << " [" << i.second << "]   ";
    cout << endl;

    cout << "\tsubnets ordering \n";
    for (auto i : subnet_name_index) cout << "\t\t[" << i.second << "] " << i.first << endl;

    cout << "subnet info " << endl;
    for (auto k : subnets) k.info(subnet_index_name);
}

void Network::info_training() {

    cout << "Training info:" << endl;
    cout << "\ttraining_config_yaml \t" << training_config_yaml << endl;
    cout << "\tcurrent_epoch \t\t" << current_epoch << endl;
    cout << "\tdur_batches_path \t" << dur_batches_path << endl;
    cout << "\twindow_res   \t\t" << window_res << endl;
    cout << "\tstep_res \t\t" << step_res << endl;
    cout << "\tw_cut \t\t\t" << w_cut << endl;
    cout << "\tPOT \t\t\t" << POT << endl;
    cout << "\tprob_noise \t\t" << prob_noise << endl;
    cout << "\tamp_noise \t\t" << amp_noise << endl;
    cout << "\tl_decay \t\t" << l_decay <<endl;
    cout << "\tl_rate0 \t\t" << l_rate0 <<endl;
    cout << "\tl_rate0_curr \t" << l_rate0_curr <<endl;
    cout << "\tdF_dI_app \t\t" << dF_dI_app <<endl;
    cout << "\tdF_dI_no_act_fact\t" << dF_dI_no_act_fact <<endl;
    cout << "\terror_lin_coeff  \t" << error_lin_coeff <<endl;

    cout << "\toptimization_method \t" << optimization_method << " (" << optimization_method_id << ")" <<endl;
    cout << "\t\tmomentum_w \t" << momentum_w <<endl;
    cout << "\t\tmomentum_c \t" << momentum_c <<endl;
    if (optimization_method_id==2) {
        cout << "\t\tmomentum2_w \t" << momentum2_w <<endl;
        cout << "\t\tmomentum2_c \t" << momentum2_c <<endl;
    }

    cout << "\tord_target_subnets\t";
    for (auto j : ord_target_subnets) {
        cout << j << " ";
    }
    cout << endl;
    cout << "\tcorrect_l_rates\t\t";
    for (auto j : correct_l_rates) {
        cout << j << " ";
    }
    cout << endl;

    cout << "\tto train connections:" << endl;
    for (auto k : subnets) {
        for (auto t : k.train_out_weight) {
            if (t.second) {
                cout << "\t\t" << k.name << " --> " << subnets[t.first].name << "\tw_min " << k.w_min << "\tw_max " << k.w_max << "\ttrain_order " << subnets[t.first].train_order <<"\tdes: " << subnets[t.first].desidered_output_path << endl;
            }
        }
    }

    cout << "\tto train offset current:" << endl;
    for (auto k : subnets) {
        if (k.train_offset_I_e) {
            cout << "\t\t" << k.name << "\ttrue" << "\toffset_I_e_min " << k.offset_I_e_min << "\toffset_I_e_max " <<  k.offset_I_e_max << "\n\t\t\t\toffset_I_e_ref " << k.offset_I_e_ref << "\toffset_I_e_decay " << k.ldecay_offset << "\tfrom file: " << k.offset_I_e_file << endl;
        }
        else {
            cout << "\t\t" << k.name << "\tfalse" << endl;
        }
    }
}
