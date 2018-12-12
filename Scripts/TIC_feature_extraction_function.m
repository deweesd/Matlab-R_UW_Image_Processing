function [S, SysPeakTime, SysPeak, xDesired, xDesired_1, xDesired_2, idx, idx_1, idx_2,  yplot1, ROI1_AUC, washin, mean_transit_time] = time_delay_plus_peaktime(tic_lin, timeStamps)
x1 = timeStamps(1:120); % taking timestamps from 1 - 120 given smoothed y1 values 
y1 = smooth(tic_lin(1:120),'sgolay',1); % smoothing raw data


%% Fit: fit 7/8'.
[xData, yData] = prepareCurveData(x1', y1 );

% curvefit given adjust x,y datapoints
fitResults1 =  fit(xData, yData, 'smoothingspline');

% evaluate the fitted y-values
yplot1 = feval(fitResults1,x1');


%% Extract Rise, Peak, and time-to-peak points given yplot1,x1 %%
S = stepinfo(yplot1,x1', 'RiseTimeThreshold', [0.1 0.95]); % adjusting rist_time thresh (i.e., wash-in) given curve essement
SysPeakTime = S.PeakTime;
SysPeak = S.Peak;
washin = S.RiseTime;

%% Calculating max index range for PIV threshold used for AT extraction

% set max_index(peak) to 0
max_index = 0;

% for loop to extract range of index values up to peak point of curve
for i = 1:length(yplot1)
    if yplot1(i) == max(yplot1)
        max_index=i;
    end
end


%% Calculating mean transit time given threshold of y, taking adjacent value that is <= MTT_value (mean of y)

MTT_value = mean(yplot1);
mean_transit_time = 0;
mtt_index = 0;
for index_mtt = max_index:length(yplot1) % starting from the max_index point:+1 length of y-values

    if (yplot1(index_mtt) <= MTT_value) % calling all values post-max_index(i.e., peak values) onward.
        mtt_index = index_mtt;
        break % making sure it stops after it hits first value < MTT_value
    end
    
end

% setting mtt_index to nan if the threshold is never met due to output
% curve never settling (recirculation is too high) 
if(mtt_index == 0)
    mtt_index = length(yplot1);
end

mean_transit_time = x1(mtt_index); %extract x given MTT threshold
%disp(mean_transit_time);



%% ROI timepoint - Global AT thresh #1 %%

% call index value given x(timestamps) -- change global_thresh for each
% iteration 
global_thresh = mean(abs(yplot1(1:max_index)*.05)); % 5% of peak_I
xDesired = 0;
idx=0;
for  at = 1:length(yplot1) % starting from the max_index point:+1 length of y-values

    if (yplot1(at) >= global_thresh) % calling all values post-max_index(i.e., peak values) onward.
        idx = at;
        break % making sure it stops after it hits first value < MTT_value
    end
    
    
end
xDesired = x1( idx ); % Extracting AT



%% ROI timepoint - Global AT thresh #2 %%
global_thresh_1 = mean(abs(yplot1(1:max_index)*.1)); % 10% of peak_I

% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
xDesired_1 = fzero(@(x1) feval(fitResults1,x1) - global_thresh_1, x1(1));

% Find 1st index where y <= thresh
idx_1 = find( yplot1 >= global_thresh_1,1 ); 

% Get the x value at this index
xDesired_1 = x1( idx_1 ); % Extracting AT_2
 


%% ROI timepoint - Global AT thresh #3 %%
global_thresh_2 = mean(abs(yplot1(1:max_index)*.15)); % 15% of peak_I

% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
xDesired_2 = fzero(@(x1) feval(fitResults1,x1) - global_thresh_2, x1(1));

% Find 1st index where y <= thresh
idx_2 = find( yplot1 >= global_thresh_2,1 ); 

% Get the x value at this index
xDesired_2 = x1( idx_2 ); % Extracting AT_3






%% Get area under the curve from fitted_y values given x %%
ROI1_AUC = trapz(x1, yplot1);

%Z = trapz(X,Y) computes the integral of Y with respect to X using
%     the trapezoidal method.  X and Y must be vectors of the same
%     length, or X must be a column vector and Y an array whose first
%     non-singleton dimension is length(X).  trapz operates along this
%     dimension.
end