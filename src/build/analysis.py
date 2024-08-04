import sys
sys.path.append('..')
import python_utils as utils


s = utils.SpikeSim(".example/NoTrain/sim", 'example.yaml', 500, 10000., 'config.yaml')

s.info()

s.histogram('all', res = 10.)

for (k,v) in (s.MeanActivity()).items():
    print(f'Mean activity of {k} \t {v[0]:.2f} kHz\t N_neurons = {v[1]} \t Activity per Neuron = {v[2]*1000:.4f} Hz')
