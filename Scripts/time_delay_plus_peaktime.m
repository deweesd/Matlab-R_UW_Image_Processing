function [S, SysRiseTime, SysPeakTime, SysPeak, xDesired, idx, yplot1] = time_delay_plus_peaktime(tic_lin, timeStamps)
x1 = timeStamps(1:120); % taking timestamps from 1 - 120 given smoothed y1 values 
y1 = smooth(tic_lin(1:120),'sgolay',1); % smoothing raw data



%% Fit: 'POLY fit 7/8'.
[xData, yData] = prepareCurveData(x1', y1 );

% Find coefficients for polynomial (order = 4 and 6, respectively)
fitResults1 =  fit(xData, yData, 'smoothingspline');

% evaluate the fitted y-values
yplot1 = feval(fitResults1,x1');

global_thresh = mean(abs(yplot1(1:60)));

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

S = stepinfo(yplot1,x1);
SysRiseTime = S.RiseTime;
SysPeakTime = S.PeakTime;
SysPeak = S.Peak;
end


