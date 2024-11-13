import numpy as np
from scipy import stats
from scipy import signal
from scipy.optimize import curve_fit
import glob
import yaml
import re
import matplotlib.pyplot as plt
import pickle

def save_pkl(obj, path):
    '''Function saving an object as a pickle file.

    :param obj: python object (list, dictionary...) to be saved
    :type obj: generic

    :param path: path to the object to be saved
    :type path: string
    '''
    with open(path, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_pkl(path):
    '''Function loading an object from a pickle file.

    :param path: path to the object to be loaded
    :type path: string
    '''

    with open(path, 'rb') as f:
        return pickle.load(f)

def readSpikes(file):
    '''Function reading spike times in the format produced by the simulation

    :param file: path to the file with spike times
    :type file: string
    '''
    l = []
    with open(file) as in_file:
        for line in in_file.readlines():
            n_list = [float(i) for i in line.split()]

            if len(l) <= round(n_list[0]):
                l.extend( [ [] for i in range( (round(n_list[0]) - len(l) + 1) ) ] )
            l[round(n_list[0])].extend(n_list[1:])

    for i in range(len(l)):
        l[i] = np.array(l[i])

    return l

def newReadSpikes(file, n):
    '''Function reading spike times in the format produced by the simulation

    :param file: path to the file with spike times
    :type file: string
    '''
    l = [ [] for i in range(n) ]
    with open(file) as in_file:
        for line in in_file.readlines():
            n_list = [float(i) for i in line.split()]

            if n <= round(n_list[0]):
                print(f'ERROR: non coherent n: {n_list[0]}>={n}')
                exit()
            # print(round(n_list[0]), n, l)
            l[round(n_list[0])].extend(n_list[1:])

    for i in range(len(l)):
        l[i] = np.array(l[i])

    return l


class SpikeSim:
    '''Class loading and parsing files given by a simulation.
    The main attributes are the simulation parameters and results:

    * end_t: end time of simulation
    * dt: time resolution of the simulation
    * input_mode: external input mode:
        - 0 (base mode): each neuron receives an indipendent poisson signal with mean frequency = SubNetwork::ext_in_rate
        - 1 (not implemented)
    * data: dictionary with spike times corresponding to each population; data['pop'] is a list of np.arrays each containing the activity of a neuron
    * subnets: a list of the SubNetworks in the simulation
    '''

    def __init__(self, path, sim_fname, neglect_t, neglect_t_end=-1, config_fname=''):
        '''Class constructor:

        :param path: path to the directory with output simulation files and configuration files
        :type path: string

        :param sim_fname: name of the simulation configuration file (inside the directory matched by path)
        :type sim_fname: string

        :param neglect_t: time (in ms) to be neglected at the beginning of the simulation
        :type sim_fname: float

        :param neglect_t_end: time (in ms) to be neglected at the end of the simulation
        :type sim_fname: float

        :param config_fname: name of the subnets_config_yaml configuration file (inside the directory matched by path)
        :type sim_fname: string
        '''
        self.input_dir = path
        self.sim_filename = sim_fname

        self.t_start = neglect_t
        self.t_end = neglect_t_end
        self.dt = 0
        self.input_mode = 0

        self.data = dict()
        self.subnets = []
        self.omegas = dict()
        self.Ns = dict()

        self.getParameterValues()

        if config_fname!='':
            with open(self.input_dir + '/' +config_fname) as file:
                in_dict = yaml.load(file, Loader=yaml.FullLoader)
            for d in in_dict:
                self.omegas[d['name']] = d['osc_omega']
                self.Ns[d['name']] = d['N']

        self.loadData()
        self.necglectTime()

    def getParameterValues(self):
        '''Method initializing the simulation parameters'''

        with open(self.input_dir + '/' +self.sim_filename) as file:
            in_dict = yaml.load(file, Loader=yaml.FullLoader)

        dict_t_end = in_dict['t_end']
        if self.t_end < 0.:
            self.t_end = dict_t_end
        elif self.t_end > dict_t_end:
            print(f'ERROR: t_end is too big: max = {dict_t_end}, passed = {self.t_end}')
            exit()
        self.dt = in_dict['dt']
        self.input_mode = in_dict['input_mode']

        if self.input_mode != 0:
            paper_configfile = in_dict['input_mode_config']
            if ( len(paper_configfile.split('/')) !=1 ):
                paper_configfile = paper_configfile.split('/')[-1]

            with open(self.input_dir + '/'+paper_configfile) as file_paper:
                paper_dict = yaml.load(file_paper, Loader=yaml.FullLoader)

    def loadData(self):
        '''Method loading spike times for each SubNetwork'''

        subnets_files = glob.glob(self.input_dir + '/*_spikes.txt')
        for f in subnets_files:
            pop = re.split('/|_', f)[-2]
            self.data[pop] = newReadSpikes(f, self.Ns[pop])
            self.subnets.append(pop)

        self.subnets = sorted(self.subnets)


    def necglectTime(self):
        '''Method removing the spikes occurring before t_start and after t_end (if > 0)'''

        for i in self.data:
            for j in range( len(self.data[i]) ):
                if self.t_end < 0.:
                    self.data[i][j] = self.data[i][j][self.data[i][j]>self.t_start] - self.t_start
                else:
                    self.data[i][j] = self.data[i][j][ np.logical_and( self.data[i][j]>self.t_start, self.data[i][j]<self.t_end ) ] - self.t_start

    def info(self):
        '''Method printing the simulation parameters'''

        print('Simulation data from: ' + self.input_dir)
        print('\t simulation config file: ' + self.sim_filename)
        print('\t subnets in the network: ', self.subnets)
        print(f'\t t_start = {self.t_start} ms')
        print(f'\t t_end = {self.t_end} ms')
        print(f'\t dt = {self.dt} ms')
        print(f'\t input_mode = {self.input_mode}')

    def saveData(self, path):
        save_pkl(self.data, path)

    def histogram(self, pop = '', res=1., save_img=''):
        '''Method showing or saving the spiking activity of a given subnet

        :param pop: desidered population; if 'all' is passed all population are showed.
        :type pop: string

        :param res: time width of each bin in the histogram
        :type pop: float

        :param save_img: path and name of the file to be saved
        :type save_img: string
        '''

        pop_passed = True
        if pop == '':
            pop_passed = False


        if pop.lower() == 'all':
            plt.figure()
            for i,p in enumerate(sorted(self.subnets)):
                print(p)
                l = len(self.subnets)
                cols = 1 if l==1 else 2
                rows = round(l/cols) if (l%2==0 or l==1) else int(l/cols)+1
                try:
                    plt.subplot(rows,cols, i+1)
                    plt.ylabel(p)
                    plt.xlabel('t [ms]')
                    plt.hist( np.concatenate(self.data[p]), bins=int((self.t_end - self.t_start)/res) )
                except Exception as e: print(e)
            # plt.tight_layout()
            # plt.savefig(self.input_dir+'/activity.png', dpi=500)
            plt.show()
        else:
            while True:

                if pop == '':
                    pop = input('histogram: enter subnetwork: ')
                    if pop.lower() == 'stop':
                        break
                    if not pop in self.subnets:
                        print(f'No subnet with name "{pop}", try again...')
                        pop = ''
                        continue

                plt.hist(np.concatenate(self.data[pop]), bins=int((self.t_end - self.t_start)/res))

                if save_img == '':
                    plt.show()
                else:
                    plt.savefig(save_img, dpi = 500)
                    plt.close()

                if pop_passed: break
                else: pop = ''

    def MeanActivity(self):
        '''Method computing the mean spiking activity of the subnets'''
        MAct = []
        Ns = []
        MActPerN = []

        for s in self.subnets:
            if len(self.data[s]) >0:
                counts = len( np.concatenate(self.data[s]) )
            else: counts = 0
            MAct.append( counts/(self.t_end - self.t_start) )
            n = len(self.data[s])
            Ns.append(n)
            if n > 0:
                MActPerN.append(counts/n/(self.t_end - self.t_start))
            else:
                MActPerN.append(0)

        return {self.subnets[i] : [ MAct[i], Ns[i], MActPerN[i] ] for i in range(len(self.subnets))}




def cartesian_to_cilindric(x,y,z, xv=0.,yv=-50.):
    r = np.sqrt((x-xv)**2+(y-yv)**2)

    if y>yv:                        # I/II quad (I:[0,pi/2], II:[pi/2,pi])
        if abs(x-xv) < 1e-8:
            theta = 1/2*np.pi
        else:
            theta = np.arctan((y-yv)/(x-xv))

        if theta<0:
            theta = np.pi+theta

        return r,theta,z
    else:                           # III/IV quad (III:[pi,3pi/2], IV:[-pi/2,0])
        if abs(x-xv) < 1e-8:
            theta = 3/2*np.pi
        else:
            theta = np.arctan((y-yv)/(x-xv))

        if theta>0:
            theta = np.pi+theta

        return r,theta,z