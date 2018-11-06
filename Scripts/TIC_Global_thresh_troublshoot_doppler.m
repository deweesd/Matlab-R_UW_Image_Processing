function [R1_global_thresh, R1_global_thresh_1, R1_global_thresh_2, R2_global_thresh, R2_global_thresh_1, R2_global_thresh_2, yplot1, yplot2, S, S_1] = TIC_Global_thresh_troublshoot_doppler(timeStamps, R1_Doppler)
%CREATEFIT(X,Y1)
%  Create a fit.
%
%  Data for 'smoothspline' fit:
%  INPUT: 
%      timeStamps file (timeStamps.mat)
%      Reconstructed file (R1_contrast.mat)
%      X Input : x1 (timepoints)
%      Y Input: y1 (unsmoothed data points - intensity values (AIU)
%  OUTPUT:
%      fitresult : a fit object representing the fit.
%      global_thresh: established threshold used for y (extracted from 10%
%      of Peak Intensity value from Curve Fit
%      yplot1: polyval of a fit object representing the fit for ROI 1
%      yplot2: polyval of a fit object representing the fit for ROI 2
%      Time_points: interpolated values of x given y-values (after fit) for
%      ROI 1
%      Time_points_2: interpolated values of x given y-values (after fit) for
%      ROI 2
%      xDesired: Desired x (timepoint) given ythreshold for ROI 1
%      xDesired_: Desired x (timepoint) given ythreshold for ROI 2

%% Begin reading in files and looping through set grid ROIs %%
frame = 58;
figure('Color', 'k');imagesc(10*log10(R1_Doppler(:,:,frame)),[110 138]); colormap('gray');
colorbar; title([num2str(frame) ' ' num2str(timeStamps(frame))],'Color','w');

%% ROI selection for curve fit %%

x = getrect; % 'getrect' establishes rectangle with four_output dimensions 
x = round(x);t=rectangle('Position',x);t.EdgeColor='g'; t.LineWidth = 2;
for k = 1:370
    tic_lin(k) = mean2(R1_Doppler(x(2):x(2)+x(4), x(1):x(1)+x(3),k));
    %tic_lin_flow(k) = mean2(R1_Doppler(x(2):x(2)+x(4), x(1):x(1)+x(3),k));
end

% Repeating getrect function for second grid overlap 
x_2 = getrect;
x_2 = round(x_2);t=rectangle('Position',x_2);t.EdgeColor='r';t.LineWidth = 2;
for k_1 = 1:370
    tic_lin_1(k_1) = mean2(R1_Doppler(x_2(2):x_2(2)+x_2(4),x_2(1):x_2(1)+x_2(3),k_1));
    %tic_lin_flow_1(k_1) = mean2(R1_Doppler(x_2(2):x_2(2)+x_2(4),x_2(1):x_2(1)+x_2(3),k_1));
end


savefig(['Contrast_ROI_final' '.fig'])

%% Overlaying time-intensity plots given ROI selections made over ultrasound stilled_image

% Create figure
figure('Color','k');
perf1 = plot(timeStamps(1:233),smooth(tic_lin(1:233),'sgolay',1));
perf1.Color = [0 1 0];              % set line color to red  
hold on;
perf2 = plot(timeStamps(1:233),smooth(tic_lin_1(1:233),'sgolay',1));
perf2.Color = [1 0 0];              % set line color to red
hold on;             % set line color to red
grid on;
set(findall(gca, 'Type', 'Line'),'LineWidth',1.0); 
h=xlabel('time (sec)'); %or h=get(gca,'xlabel')
set(h, 'FontSize', 16); 
set(h,'FontWeight','bold'); %bold font
set(h, 'Color', 'w'); 
y=ylabel('linear power');
set(y, 'FontSize', 16); 
set(y,'FontWeight','bold'); %bold font
set(y, 'Color', 'w');
set(gca, 'Color', 'k');
set(gca,'FontSize',20, 'FontWeight', 'bold', 'XColor',[1 1 1],'YColor',[1 1 1],...
    'ZColor',[1 1 1]);
savefig(['Contrast_bolus_final'  '.fig'])
%title(num2str(x_1),'Color','w');



%% setting x1/y1 variables %%

x1 = timeStamps(1:60); % taking timestamps from 1 - 120 given smoothed y1 values 
y1 = smooth(tic_lin(1:60),'sgolay',1);
%thresh = mean(y1(1:50));

x2 = timeStamps(1:100); % taking timestamps from 1 - 120 given smoothed y1 values 
y2 = smooth(tic_lin_1(1:100),'sgolay',1);
%thresh_2 = mean(y2(1:95));


%% Fit: 'POLY fit 7/8'.
[xData, yData] = prepareCurveData(x1', y1 );
[xData_2, yData_2] = prepareCurveData(x2', y2 );

% Find coefficients for polynomial 
fitResults1 =  fit(xData, yData, 'smoothingspline');
fitResults2 =  fit(xData_2, yData_2, 'smoothingspline');

% evaluate the fitted y-values
yplot1 = feval(fitResults1,x1');
yplot2 = feval(fitResults2,x2');


%% Extract Rise and Peak time points given yplot1,x1 %%
S = stepinfo(yplot1,x1);
SysRiseTime = S.RiseTime;
SysPeakTime = S.PeakTime;
SysPeak = S.Peak;

% take 5, 7.5, & % of peakI threshold for global_thresh for first ROI %
R1_global_thresh = abs(SysPeak*.05);
R1_global_thresh_1 = abs(SysPeak*.075);
R1_global_thresh_2 = abs(SysPeak*.1);

%% Extract Rise and Peak time points given yplot2,x2 %%
S_1 = stepinfo(yplot2,x2);
SysRiseTime_1 = S_1.RiseTime;
SysPeakTime_1 = S_1.PeakTime;
SysPeak_1 = S_1.Peak;

% take 5, 7.5, & % of peakI threshold for global_thresh for second ROI %
R2_global_thresh = abs(SysPeak_1*.05);
R2_global_thresh_1 = abs(SysPeak_1*.075);
R2_global_thresh_2 = abs(SysPeak_1*.1);



% interpolates to find yi, the values of the underlying function Y at the points in the vector or array xi. x must be a vector. 
Time_points = interp1(yplot1, x1, y1);
Time_points_2 = interp1(yplot2, x2, y2);


%% ROI 1 timepoint - first global threshold %%
% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
R1_xDesired = fzero(@(x1) feval(fitResults1,x1) - R1_global_thresh, x1(1));

% Find 1st index where y >= thresh ~ x-number of points until yplot1 >= to
% threshold
idx = find( yplot1 >= R1_global_thresh,1 ); 

% Get the x value at this index
R1_xDesired = x1( idx );

%% ROI 1 timepoint - second global threshold %%
% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
R1_xDesired_1 = fzero(@(x1) feval(fitResults1,x1) - R1_global_thresh_1, x1(1));

% Find 1st index where y >= thresh ~ x-number of points until yplot1 >= to
% threshold
R1_idx_1 = find( yplot1 >= R1_global_thresh_1,1 ); 

% Get the x value at this index
R1_xDesired_1 = x1( R1_idx_1 );

%% ROI 1 timepoint - third global threshold %%
% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
R1_xDesired_2 = fzero(@(x1) feval(fitResults1,x1) - R1_global_thresh_2, x1(1));

% Find 1st index where y >= thresh ~ x-number of points until yplot1 >= to
% threshold
R1_idx_2 = find( yplot1 >= R1_global_thresh_2,1 ); 

% Get the x value at this index
R1_xDesired_2 = x1( R1_idx_2 );


%% ROI 2 timepoint - first global threshold %%
% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
R2_xDesired = fzero(@(x2) feval(fitResults2,x2) - R1_global_thresh, x2(1));

% Find 1st index where y >= thresh
R2_idx = find( yplot2 >= R1_global_thresh,1 ); 

% Get the x value at this index
R2_xDesired = x2( R2_idx );



%% ROI 2 timepoint - Second global threshold %%
% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
R2_xDesired_1 = fzero(@(x2) feval(fitResults2,x2) - R1_global_thresh_1, x2(1));

% Find 1st index where y >= thresh
R2_idx_1 = find( yplot2 >= R1_global_thresh_1,1 ); 

% Get the x value at this index
R2_xDesired_1 = x2( R2_idx_1 );


%% ROI 2 timepoint - Third global threshold %%
% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
R2_xDesired_2 = fzero(@(x2) feval(fitResults2,x2) - R1_global_thresh_2, x2(1));

% Find 1st index where y >= thresh
R2_idx_2 = find( yplot2 >= R1_global_thresh_2,1 ); 

% Get the x value at this index
R2_xDesired_2 = x2( R2_idx_2 );



%% Plot fit with data while also plotting arrival delay point %%
figure( 'Name', 'Curvefit1_poly' );
h = plot(x1', y1);%smoothed-points
hold on;
plot(x1', yplot1);%polyfit points
hold on;
plot(x1', yplot1, '*g');%interpolated points of x given y
hold on;

plot(x2', y2);
hold on;
plot(x2', yplot2);
hold on;
plot(x2', yplot2, '*r');
hold on;
max_y = max(yplot1);
max_y_2 = max(yplot2);


R1_x_arrival = R1_xDesired;
R1_x_arrival_1 = R1_xDesired_1;
R1_x_arrival_2 = R1_xDesired_2;

R2_x_arrival = R2_xDesired;
R2_x_arrival_1 = R2_xDesired_1;
R2_x_arrival_2 = R2_xDesired_2;











R1_diff_in_delay = R2_xDesired - R1_xDesired;
R2_diff_in_delay = R2_xDesired_1 - R1_xDesired_1;
R3_diff_in_delay = R2_xDesired_2 - R1_xDesired_2;

hold on;
%% ROI 1 line troubleshoot %%
line([R1_x_arrival R1_x_arrival], [0 R1_global_thresh ], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5 ); % plot arrival point in dash vertical line 
plot([0 R1_x_arrival], [R1_global_thresh  R1_global_thresh ],  'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.25 ); % plot corresponding max(yi) horizontal line
hold on;
line([R1_x_arrival_1 R1_x_arrival_1], [0 R1_global_thresh_1], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5 ); % plot arrival point in dash vertical line 
plot([0 R1_x_arrival_1], [R1_global_thresh_1 R1_global_thresh_1],  'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.25 ); % plot corresponding max(yi) horizontal line
hold on;
line([R1_x_arrival_2 R1_x_arrival_2], [0 R1_global_thresh_2], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5 ); % plot arrival point in dash vertical line 
plot([0 R1_x_arrival_2], [R1_global_thresh_2 R1_global_thresh_2],  'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.25 ); % plot corresponding max(yi) horizontal line
%% ROI 2 line troubleshoot %%
line([R2_x_arrival R2_x_arrival], [0 R2_global_thresh ], 'Color', 'black', 'LineStyle', '-', 'LineWidth', 1.5 ); % plot arrival point in dash vertical line 
plot([R1_x_arrival R2_x_arrival], [R2_global_thresh  R2_global_thresh ],  'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.25 ); % plot corresponding max(yi) horizontal line
hold on;
line([R2_x_arrival_1 R2_x_arrival_1], [0 R2_global_thresh_1], 'Color', 'black', 'LineStyle', '-', 'LineWidth', 1.5 ); % plot arrival point in dash vertical line 
plot([R1_x_arrival_1 R2_x_arrival_1], [R2_global_thresh_1 R2_global_thresh_1],  'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.25 ); % plot corresponding max(yi) horizontal line
hold on;
line([R2_x_arrival_2 R2_x_arrival_2], [0 R2_global_thresh_2], 'Color', 'black', 'LineStyle', '-', 'LineWidth', 1.5 ); % plot arrival point in dash vertical line 
plot([R1_x_arrival_2 R2_x_arrival_2], [R2_global_thresh_2 R2_global_thresh_2],  'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.25 ); % plot corresponding max(yi) horizontal line
grid on;
legend( h, 'smoothsplinefit vs. timepoints', 'Location', 'NorthEast' );
% Label axes %
xlabel Time
ylabel Intensity
title(['Delta_T' num2str(R1_diff_in_delay) ' ' num2str(R2_diff_in_delay) ],'Color','k');
grid on
savefig(['Curvefit_ROI_Acute_AT_final' '.fig']);
end



