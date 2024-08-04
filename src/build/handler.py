import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from subprocess import run
from tqdm import tqdm

sys.path.append('../')
import python_utils as utils

REPS = 10
out_dir = 'example/Train'

''' The out_dir should already exists with the following files:
        - config1.yaml
        - config_connections1.yaml
        - config_weights1
        - config_to_save1.yaml
        - batches_dur.txt

    These files should be initialized as indicated in the tutorials.
'''

total_datapoints = 160
single_input_dur = 10       # seconds   change also in batches.txt
n_hidden = 8               #           change also in config1.yaml

high_inp = 50
low_inp  = 30
high_out_fr = 15            # hertz
low_out_fr  = 5
noise_sigma = 0.1

time_shift  = 0.5           # seconds
time_window = 1.0           # seconds

for dir in ['specific_w', 'specific_Ie', 'desired_out']:
    if not os.path.exists(f'{out_dir}/{dir}'):
        os.makedirs(f'{out_dir}/{dir}')

for rep in range(REPS):

    res_dir = f'res{rep}'

    if os.path.exists(f'{out_dir}/{res_dir}'):
        print(f'{out_dir}/{res_dir} already exists')
        exit()
    else:
        os.makedirs(f'{out_dir}/{res_dir}')
        os.makedirs(f'{out_dir}/{res_dir}/log')

    def str_line(ind, v):
        str_i = str(ind)
        str_v = "\t".join([ f'{j:.2f}' for j in v ])
        return "\t".join([ str_i, str_v ])

    hist_acc = []

    pbar = tqdm(range(total_datapoints))
    for data_point in pbar:
        input_data = np.random.randint(0,4)
        # print(f'starting datapoint {data_point}  ({input_data})')

        # setup data points
        outT = np.zeros(int(single_input_dur/time_shift)) + low_out_fr*time_window          # out True
        outF = np.zeros(int(single_input_dur/time_shift)) + low_out_fr*time_window          # out False
        if input_data==0:
            inp0 = low_inp + np.random.normal(0, noise_sigma, size=single_input_dur*1000)
            inp1 = low_inp + np.random.normal(0, noise_sigma, size=single_input_dur*1000)
            outF = outF + high_out_fr*time_window
            des_out = 1
        elif input_data==1:
            inp0 = high_inp + np.random.normal(0, noise_sigma, size=single_input_dur*1000)
            inp1 = low_inp + np.random.normal(0, noise_sigma, size=single_input_dur*1000)
            outT = outT + high_out_fr*time_window
            des_out = 0
        elif input_data==2:
            inp0 = low_inp + np.random.normal(0, noise_sigma, size=single_input_dur*1000)
            inp1 = high_inp + np.random.normal(0, noise_sigma, size=single_input_dur*1000)
            outT = outT + high_out_fr*time_window
            des_out = 0
        elif input_data==3:
            inp0 = high_inp + np.random.normal(0, noise_sigma, size=single_input_dur*1000)
            inp1 = high_inp + np.random.normal(0, noise_sigma, size=single_input_dur*1000)
            outF = outF + high_out_fr*time_window
            des_out = 1

        # setting up input
        input_string = ''
        for i_inp, inp in enumerate([inp0, inp1]):
            input_string += str_line(i_inp, inp) + '\n'
        with open(f'{out_dir}/specific_Ie/L1.txt', 'w') as f:
            f.write(input_string)

        # setting up output
        out_string = ''
        for i_out, out in enumerate([outT, outF]):
            out_string += str_line(i_out, out) + '\n'
        with open(f'{out_dir}/desired_out/L3.txt', 'w') as f:
            f.write(out_string)

        # setting up weights
        if data_point==0:
            # L1_to_L2
            temp_str = ''
            for i_neur in range(2):
                v_zero = np.abs(np.random.normal(0,15+0.010,n_hidden))
                temp_str += str_line(i_neur,v_zero) + '\n'
            with open(f'{out_dir}/specific_w/L1_to_L2.txt', 'w') as f:
                f.write(temp_str)

            np.savetxt(f'{out_dir}/specific_w/L2.txt', np.zeros(n_hidden))

            # L2_to_L3
            temp_str = ''
            for i_neur in range(n_hidden):
                v_zero = np.abs(np.random.normal(0,10+0.010,2))
                temp_str += str_line(i_neur,v_zero) + '\n'
            with open(f'{out_dir}/specific_w/L2_to_L3.txt', 'w') as f:
                f.write(temp_str)

            np.savetxt(f'{out_dir}/specific_w/L3.txt', np.zeros(2))

        else:
            os.system(f'cp {out_dir}/{res_dir}/{data_point-1}/specific_weights/L1_to_L2_{data_point-1}.txt {out_dir}/specific_w/L1_to_L2.txt')
            os.system(f'cp {out_dir}/{res_dir}/{data_point-1}/specific_weights/L2_to_L3_{data_point-1}.txt {out_dir}/specific_w/L2_to_L3.txt')
            os.system(f'cp {out_dir}/{res_dir}/{data_point-1}/specific_weights/L2_{data_point-1}.txt {out_dir}/specific_w/L2.txt')


        training_string = f'''t_end:                          {single_input_dur*1000}      # ms
dt:                             1.   # ms
n_step:                         1    # unused
out_dir:                        {out_dir}/{res_dir}/{data_point}
training_config_yaml:           {out_dir}/config_training1.yaml
subnets_config_yaml:            {out_dir}/config1.yaml
weights_config_yaml:            {out_dir}/config_weights1.yaml
connections_config_yaml:        {out_dir}/config_connections1.yaml
to_save_config_yaml:            {out_dir}/config_to_save1.yaml
input_mode:                     0   # unused
input_mode_config:              unused
train:                          true
repeat_specific_I_e:            1
'''
        with open(f'{out_dir}/training.yaml', 'w') as f:
            f.write(training_string)

        config_training_string = f'''epoch:                   {data_point}
dur_batches_path:        {out_dir}/batches_dur.txt
window_res:              {time_window*1000}
step_res:                {time_shift*1000}
w_cut:                   100.
POT:                     1
error_lin_coeff:         0.
prob_noise:              0.
amp_noise:               0.
l_decay:                 0.
l_rate0:                 {0.5*50}
l_rate0_curr:            {50}  # {0.5*5}
optimization_method:     SGD
momentum_w:              {0.}
momentum_c:              {0.}
momentum2_w:             {0.95}     # used only in ADAM
momentum2_c:             {0.95}     # used only in ADAM
correct_l_rates:         [ 1., 1. ]
dF_dI_app:               {0.002}
dF_dI_no_act_fact:       {0.1}
epsilon_zero_act:        0.1
a_desid_zero:            0.
c_desid_zero:            1.
n_quant:                 -1.
quant_delta_weight_neg:  0.
quant_delta_weight_pos:  0.
'''
        with open(f'{out_dir}/config_training1.yaml', 'w') as f:
            f.write(config_training_string)

        with open(f'{out_dir}/{res_dir}/log/std_out_{data_point}.txt', 'w') as f_out:
            with open(f'{out_dir}/{res_dir}/log/std_err_{data_point}.txt', 'w') as f_err:
                run(f'./main {out_dir}/training.yaml', shell=True, stdout=f_out, stderr=f_err, text=True)

        sim = utils.SpikeSim( f'{out_dir}/{res_dir}/{data_point}', f'training.yaml', 0, neglect_t_end=-1, config_fname='config1.yaml' )

        if len(sim.data['L3'][0])>len(sim.data['L3'][1]):
            most_active = 0
        else:
            most_active = 1
        if most_active==des_out:
            hist_acc.append(1)
        else:
            hist_acc.append(0)

        pbar.set_description(f'rep {rep+1}/{REPS} -- accuracy (mean 20): {np.mean(hist_acc[-20:]):.3f}')

    np.savetxt(f'{out_dir}/{res_dir}/hist_acc.txt', hist_acc)
    mean_hist = 20
    hist_acc = np.convolve(hist_acc, np.ones(mean_hist)/mean_hist, 'valid')
    plt.plot(hist_acc)

plt.xlabel('presented data points')
plt.ylabel(f'accuracy (average across {mean_hist} data points)')
plt.show()
