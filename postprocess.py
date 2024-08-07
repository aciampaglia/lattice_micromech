from visualization import *
import numpy as np
import os
## Set up the analysis
output_directory = 'cyrve_ellips_2x2_d14_ar3_mat2'
os.chdir(output_directory)
rve_size = [1.4] # mm
defects_file = 'dmax_defects_2x2.txt'
defect_values = np.loadtxt(defects_file) #, skiprows=1, usecols=3)*1e-3
repetitions = 3
for i, d in enumerate(defect_values):
    os.chdir('defect_{}'.format(int(d*100)))
    for size in rve_size:
        for rep in range(repetitions):
            model_name='size_analysis_{}_{}_{}'.format(int(d*100), int(size*100), rep+1)
            os.chdir(model_name)
            odb = openOdb(model_name+'.odb')
            # Access the stress, strain, and volume data for each element
            stress_avg = np.zeros((len(odb.steps['Static_1'].frames),1))
            strain_avg = np.zeros((len(odb.steps['Static_1'].frames),1))
            total_volume = size**2*np.pi/4*2;
            for frame_i, frame in enumerate(odb.steps['Static_1'].frames):
                total_weighted_stress = 0.0;
                total_weighted_strain = 0.0;
                stressValues = frame.fieldOutputs['S'].bulkDataBlocks[0].data[:,1];
                volumeValues = frame.fieldOutputs['EVOL'].bulkDataBlocks[0].data[:,0];
                strainValues = frame.fieldOutputs['E'].bulkDataBlocks[0].data[:,1];
                total_weighted_stress = np.sum(np.multiply(stressValues, volumeValues));
                total_weighted_strain = np.sum(np.multiply(strainValues, volumeValues));
                stress_avg[frame_i] = np.divide(total_weighted_stress, total_volume);
                strain_avg[frame_i] = np.divide(total_weighted_strain, total_volume);
            # Combine vectors into a 2D array
            data_matrix = np.column_stack((stress_avg, strain_avg))
            # Create a header
            header = "Stress Strain"
            # Save data to a text file
            output_file_path = 'output_data.txt'
            np.savetxt(output_file_path, data_matrix, header=header, comments='', fmt='%1.5e')
            odb.close()
            os.chdir('..')
    os.chdir('..')
