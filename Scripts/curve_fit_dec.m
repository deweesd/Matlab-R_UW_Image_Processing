function [global_thresh, global_thresh_1, xDesired, xDesired_2, idx, idx_2, yplot1, yplot2, S , S_1, ROI1_AUC, ROI2_AUC] = curve_fit_decomp(timeStamps, R1_contrast)
%CREATEFIT(X,Y1)
%  Create a fit.
%
%  Data for 'polynomial' fit:
%  INPUT: 
%      timeStamps file (timeStamps.mat)
%      Reconstructed file (R1_contrast.mat)
%      X Input : x1 (timepoints)
%      Y Input: y1 (unsmoothed data points - intensity values (AIU)
%  OUTPUT:
%      fitresult : a fit object representing the fit.
%      yplot1: polyval of a fit object representing the fit for ROI 1
%      yplot2: polyval of a fit object representing the fit for ROI 2
%      Time_points: interpolated values of x given y-values (after fit) for
%      ROI 1
%      Time_points_2: interpolated values of x given y-values (after fit) for
%      ROI 2
%      xDesired: Desired x (timepoint) given ythreshold for ROI 1
%      xDesired_: Desired x (timepoint) given ythreshold for ROI 2

%% Begin reading in files and looping through set grid ROIs %%
frame = 50;
figure('Color', 'k');imagesc((10*log10(R1_contrast(:,:,frame))),[110 148]); colormap('gray');
colorbar; title([num2str(frame) ' ' num2str(timeStamps(frame))],'Color','w');
% figure('Color', 'k');imagesc(fliplr(R1_contrast(:,:,idx)),[115 150]); colormap('gray');
% colorbar; title([num2str(idx) ' ' num2str(timeStamps(idx))],'Color','w');
% ax = gca;
%     set(gca,'XTickLabel',{'-4', '-3', '-2', '-1', '0', '1', '2', '3' '4' '5'} );% get the current axis
%     set(gca,'YTickLabel',1:10);
%     %set ( gca, 'xdir', 'reverse' ); % decending values 
%     ttl = ax.Title;             % get the title text object
%     ttl.FontWeight = 'bold';    % set the font to bold
%     ttl.FontSize = 28;
%     h=xlabel('Distance(mm)'); %or h=get(gca,'xlabel')
%     set(h, 'FontSize', 18); 
%     set(h,'FontWeight','bold');%bold font
%     z=ylabel('Width(mm)');
%     set(z, 'FontSize', 18);
%     set(z, 'FontWeight', 'bold');
%     box(ax,'off'); 
%     set(gca, 'XColor', 'white'); % set x-axis color to white
%     set(gca, 'YColor', 'white');
%     %set(gca,'FontSize',12, 'FontWeight', 'bold');
%     %%%%%% Reduce Whitespace around figure %%%%%%
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset; 
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3);
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];

%% ROI selection for curve fit %%

%%% 'getrect' establishes rectangle with four_output dimensions %%%     
x = getrect; 
x = round(x);t=rectangle('Position',x);t.EdgeColor='g'; t.LineWidth = 2;
for k = 1:100
    tic_lin(k) = mean2(R1_contrast(x(2):x(2)+x(4), x(1):x(1)+x(3),k));
end

% %%% second ROI selection %%%
x_1 = getrect; 
x_1 = round(x_1);t=rectangle('Position',x_1);t.EdgeColor='r'; t.LineWidth = 2;
for k_1 = 1:100
    tic_lin_1(k_1) = mean2(R1_contrast(x_1(2):x_1(2)+x_1(4), x_1(1):x_1(1)+x_1(3),k_1));
end



% %%% second ROI selection %%%
x_2 = getrect; 
x_2 = round(x_2);t=rectangle('Position',x_2);t.EdgeColor='y'; t.LineWidth = 2;
for k_1 = 1:100
    tic_lin_1(k_1) = mean2(R1_contrast(x_2(2):x_2(2)+x_2(4), x_2(1):x_2(1)+x_2(3),k_1));
end
savefig(['Frame2_ROI_Acute_DM_final' '.fig']);

%% Overlaying time-intensity plots given ROI selections made over ultrasound stilled_image

% Create figure
figure('Color', 'k');
perf1 = plot(timeStamps(1:100),smooth(tic_lin(1:100),'sgolay',1));
perf1.Color = [0 1 0];              % set line color to red  
hold on;
perf2 = plot(timeStamps(1:100),smooth(tic_lin_1(1:100),'sgolay',1));
perf2.Color = [1 0 0];              % set line color to red
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
savefig(['Contrast_bolus_acute_final'  '.fig']);
title(num2str(x_1),'Color','w');



%% setting x1/y1 variables %%

x1 = timeStamps(1:100); % taking timestamps from 1 - 120 given smoothed y1 values 
y1 = smooth(tic_lin(1:100),'sgolay',1);
%thresh = mean(y1(1:60));

x2 = timeStamps(1:100); % taking timestamps from 1 - 120 given smoothed y1 values 
y2 = smooth(tic_lin_1(1:100),'sgolay',1);
%thresh_2 = mean(y2(1:60));

%% Fit: 'POLY fit 7/8'.
[xData, yData] = prepareCurveData(x1', y1 );
[xData_2, yData_2] = prepareCurveData(x2', y2 );

% Find coefficients for polynomial 
fitResults1 =  fit(xData, yData, 'smoothingspline');
fitResults2 =  fit(xData_2, yData_2, 'smoothingspline');

% evaluate the fitted y-values
yplot1 = feval(fitResults1,x1');
yplot2 = feval(fitResults2,x2');

% Set threshold for fitted_y values
global_thresh = mean(abs(yplot1(1:75)));
global_thresh_1 = mean(abs(yplot2(1:90)));


% interpolates to find yi, the values of the underlying function Y at the points in the vector or array xi. x must be a vector. 
Time_points = interp1(yplot1, x1, y1);
Time_points_2 = interp1(yplot2, x2, y2);


%% ROI 1 timepoint %%
% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
xDesired = fzero(@(x1) feval(fitResults1,x1) - global_thresh, x1(1));

% Find 1st index where y >= thresh ~ x-number of points until yplot1 >= to
% threshold
idx = find( yplot1 >= global_thresh,1 ); 

% Get the x value at this index
xDesired = x1( idx );

% Get the rise time and peak time as well as y-peak value 
S = stepinfo(yplot1,x1);
SysRiseTime = S.RiseTime;
SysPeakTime = S.PeakTime;
SysPeak = S.Peak;
xDesired_peak = SysPeakTime;


% Get area under the curve from fitted_y values given x 
ROI1_AUC = trapz(x1, log10(yplot1));

%% ROI 2 timepoint %%
% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
xDesired_2 = fzero(@(x2) feval(fitResults2,x2) - global_thresh_1, x2(1));

% Find 1st index where y >= thresh
idx_2 = find( yplot2 >= global_thresh_1,1 ); 

% Get the x value at this index
xDesired_2 = x2( idx_2 );


% Get the rise time and peak time as well as y-peak value 
S_1 = stepinfo(yplot2,x2);
SysRiseTime_1 = S_1.RiseTime;
SysPeakTime_1 = S_1.PeakTime;
SysPeak_1 = S_1.Peak;
xDesired_peak_1 = SysPeakTime_1;

% Get area under the curve from fitted_y values given x 
ROI2_AUC = trapz(x2, log10(yplot2));

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
max_y = max(yplot1);
max_y_2 = max(yplot2);
x_arrival = xDesired;
x_arrival_2 = xDesired_2;
peak_1 = xDesired_peak;
peak_2 = xDesired_peak_1;
diff_in_peak = 10*log10(max(yplot1) - max(yplot2));
hold on;
line([peak_1 peak_1], [min(yplot1) max_y], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1.5 ); % plot arrival point in dash vertical line 
plot([0 peak_1], [SysPeak SysPeak],  'Color', 'g', 'LineStyle', '--', 'LineWidth', 1.25 ); % plot corresponding max(yi) horizontal line
hold on;
line([peak_2 peak_2], [min(yplot2) max_y_2], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1.5 ); % plot arrival point in dash vertical line 
plot([0 peak_2], [SysPeak_1 SysPeak_1],  'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.25 ); % plot corresponding max(yi) horizontal line
grid on;
legend( h, 'polyfit vs. timepoints', 'Location', 'NorthEast' );
% Label axes
xlabel Time
ylabel Intensity
title(['Delta_P' '^{log10}' num2str(diff_in_peak) ],'Color','k');
grid on
savefig(['Curvefit_ROI_BL_Peak_final' '.fig']);

%% Convert .fig files to .png files %%
% ALWAYS CHANGE DIRECT PATH BEFORE RUNNING %%

%folder = '/Users/Deweesd/Desktop/BL_Acute_arrival_time_heatmaps/Severe/100/BL/Scaled_HM/raw_peak_intensity_data/';

% Get all .fig files in the folder
files = dir(fullfile(folder, '*.fig'));
files = fullfile(folder, {files.name});


for k = 1:numel(files)
    % Get the filename
    [~, fname] = fileparts(files{k});

    % Open and display the .fig file
    hfig = openfig(files{k});
    
    
    hfig.InvertHardcopy = 'off'; % Maintain black background
   
    
    % Save as a PNG file with the same name as the .fig file
    saveas(hfig, fullfile(folder, [fname,  '.png']));

    % Close the figure again
    close(hfig)
end
end



