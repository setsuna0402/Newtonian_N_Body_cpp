import math
import numpy as np
import pylab
import os
import h5py
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm


def read_and_plot(data_filename, plot_name, low_bound = 0.0, high_bound = 1.0, SavePlot=False):
    # Read the data from the file
    file_cpu = h5py.File(str(data_filename), 'r')  # Open hdf5 snapshot file
    '''
    print('List of arrays in this file: ')
    for i in file_cpu.items():
        print(i)
    '''
    pos = np.array(file_cpu.get('Position'))
    '''
    acc = np.array(file_cpu.get('Acceleration'))
    vel = np.array(file_cpu.get('Velocity'))
    mass = np.array(file_cpu.get('Mass'))
    '''
    pe = np.array(file_cpu.get('PotentialEnergy'))
    ke = np.array(file_cpu.get('KineticEnergy'))
    pe_total = np.sum(pe) * 0.5 # 0.5 beacause we count each pair of particles twice
    ke_total = np.sum(ke)
    time = np.array(file_cpu.get('code_time'))
    delta_time = np.array(file_cpu.get('delta_t'))
    file_cpu.close()
    # print("Data read successfully!")
    # print("Number of particles: ", pos.shape[0])
    print("Time: ", time)
    print("Delta time: ", delta_time)
    print("Total Potential Energy: ", pe_total)
    print("Total Kinetic Energy: ", ke_total)
    print("Total Energy: ", pe_total + ke_total)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter3D(pos[:, 0], pos[:, 1], pos[:, 2], c="b",
                 s=2.0, alpha=1.0, cmap='Greens')
    # fix extension
    ax.set_xlim([low_bound, high_bound])
    ax.set_ylim([low_bound, high_bound])
    ax.set_zlim([low_bound, high_bound])
    if SavePlot == True:
        # print("Saving plot to file: ", str(plot_name))
        plt.savefig(str(plot_name), dpi=300)
        plt.close(fig)
    else:
        print("Plotting to screen...")
        plt.show()
    return pos


print("Here we go!")
# Read the data from the file
prefix = "N_body_8_thread_SIMD_N_1024_step_"
suffix = ".hdf5"
# step_list = [0, 10000, 20000, 30000, 38461]

# for i in range(len(step_list)):
for i in range(0, 38462, 1000):
    step = i
    data_filename = prefix + str(step) + suffix
    plot_name = "plot_" + str(step) + ".png"
    print("Reading data from file: ", data_filename)
    pos = read_and_plot(data_filename, plot_name, SavePlot=True)

# plt.hist(error_pos_8_thread_no_simd, bins='auto')
# plt.savefig("ABC")

print("May the force be with you!")
