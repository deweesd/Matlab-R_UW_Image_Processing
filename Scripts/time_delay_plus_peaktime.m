function [S, SysRiseTime, SysPeakTime, SysPeak, xDesired, idx, yplot1] = time_delay_plus_peaktime(tic_lin, timeStamps)
x1 = timeStamps(5:100); % taking timestamps from 1 - 120 given smoothed y1 values 
y1 = smooth(tic_lin(5:100),'sgolay',1); % smoothing raw data



%% Fit: 'POLY fit 7/8'.
[xData, yData] = prepareCurveData(x1', y1 );

% Find coefficients 
fitResults1 =  fit(xData, yData, 'smoothingspline');

% evaluate the fitted y-values
yplot1 = feval(fitResults1,x1');

% set global threshold
global_thresh = mean(abs(yplot1(5:70)));

% interpolates to find yi, the values of the underlying function Y at the points in the vector or array xi. x must be a vector. 
Time_points = interp1(yplot1, x1', y1);

%% ROI 1 timepoint %%

% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
xDesired = fzero(@(x1) feval(fitResults1,x1) - global_thresh, x1(1));

% Find 1st index where y >= thresh
idx = find( yplot1 >= global_thresh,1 ); 

% Get the x value at this index
xDesired = x1( idx );

%% Extract Rise and Peak time points given yplot1,x1 %%

% S = stepinfo(sys)computes the step-response characteristics for a dynamic system model sys. The function returns the characteristics in a structure containing the fields:
% 
%     RiseTime ? Time it takes for the response to rise from 10% to 90% of the steady-state response.
% 
%     SettlingTime ? Time it takes for the error |y(t) - yfinal| between the response y(t) and the steady-state response yfinal to fall to within 2% of yfinal.
% 
%     SettlingMin ? Minimum value of y(t) once the response has risen.
% 
%     SettlingMax ? Maximum value of y(t) once the response has risen.
% 
%     Overshoot ? Percentage overshoot, relative to yfinal).
% 
%     Undershoot ? Percentage undershoot.
% 
%     Peak ? Peak absolute value of y(t)
% 
%     PeakTime ? Time at which the peak value occurs.


S = stepinfo(yplot1,x1); %
SysRiseTime = S.RiseTime;
SysPeakTime = S.PeakTime;
SysPeak = S.Peak;
end


