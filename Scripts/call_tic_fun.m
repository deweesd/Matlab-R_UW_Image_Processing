clc; clear all; close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script created by DMD
%   1) Run script
%   2) Choose .mat file(s) to load using the pop-up menu (Contrast & it's
%   associated timestamp file).
%   3) Draw the ROI. Once you close the rectangle, it automatically
%       shows the graphs and associated variables of interest. - Also saves
%       variables in your workspace that give raw values for each variable.
%   4) 5 saved figures will also be saved in the current directory you are
%   in (ROI figure, TIC plot with the 4 corresponding values highlights,
%   and finally the 3 output frames with frame number and time for each
%   measure)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
contrast = '';
timestamp = '';
path = '';

[contrast, path] = uigetfile;
[timestamp, path] = uigetfile;



if ~isempty(contrast) && ~isempty(timestamp) && ~isempty(path)
        uiopen([path contrast],1)
else
    return
end

if ~isempty(contrast) && ~isempty(timestamp) && ~isempty(path)
        uiopen([path timestamp],1)
else
    return
end
[fifty_percent_minus_thresh, plus_50, minus_50, max_index, minus_50_idx, mtt_index, SysPeakTime, at, at_idx, at_thresh, yplot1, S, ROI1_AUC] = TIC_features_fun(timeStamps, R1_contrast);