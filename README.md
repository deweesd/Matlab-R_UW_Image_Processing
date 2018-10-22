# Matlab_CEUS_Features

Here are scripts generated via Matlab 2018a that consist of analyzing characteristics from high intensity ultrasound image processing using an injured spinal cord model. 

There are four scripts that correspond to various features of interest. Firstly, arrival time delay (seconds) is a parameter that is calculuated using the 'time_delay_plus_peaktime' function script. This function is based off of a 'curve_fit_dec' function that I use to extract adjusted y-fit values given x (i.e., timepoints). The adjusted y-values are used to threshold out a global value and given the adjusted y >= global thresh, give corresponding x value. This in turn gives me a raw integer representitive of arrival time for region of interest. (please reference curve fit images and heatmaps for visual representations).

The last script is compiling the raw values given a select region on the spinal cord of arrival times and generates a parameteric heatmap. Reference the attached .png files. 

