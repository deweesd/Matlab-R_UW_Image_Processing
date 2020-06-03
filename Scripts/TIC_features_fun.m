function [fifty_percent_minus_thresh, plus_50, minus_50, max_index, minus_50_idx, mtt_index, SysPeakTime, at, at_idx, at_thresh, yplot1, S, ROI1_AUC] = TIC_features_fun(timeStamps, R1_contrast)
%  Create a fit.
%  Data for 'polynomial' fit:

%  INPUT: 
%      timeStamps file (timeStamps.mat)
%      Reconstructed file (R1_contrast.mat)
%      X Input : x1 (timepoints)
%      Y Input: y1 (unsmoothed data points - intensity values (AIU)

%  OUTPUT:
%      fitresult : a fit object representing the fit.;
%      yplot1: polyval of a fit object representing the fit for ROI 1;
%      yplot2: polyval of a fit object representing the fit for ROI 2;
%      Time_points: interpolated values of x given y-values (after fit); 
%      xDesired: Desired x (timepoint) given ythreshold for ROI 1;
%      Thresholds (fifty_percent_minus_thresh, at_thresh);
%      Peak values (plus 50, minus 50, peak);
%      Arrival time and it's associated idx (at, at_idx);
%      S: stepinfo summary of calcuated values for each TIC variable


%% Begin reading in files and looping through set grid ROIs %%
frame = 52;
figure('Color', 'k');imagesc(10*log10(R1_contrast(:,:,frame)),[110 145]); colormap('gray');
colorbar; title([num2str(frame) ' ' num2str(timeStamps(frame))],'Color','w');
figperf = gcf; %Save figuredata as variable to workspace
axperf = figperf.CurrentAxes;
Xlab = [str2double(axperf.XTickLabel{1})...
            str2double(axperf.XTickLabel{end})]; %Dim num min/max X
      Xlab = [-6.4 6.4]; %Hardcoded stuff
Zlab = [str2double(axperf.YTickLabel{1})...
            str2double(axperf.YTickLabel{end})]; %Dim num min/max Z
      Zlab = [0 10]; %Hardcoded stuff



%% Begin reading in files and looping through set grid ROIs %%
x = getrect;
x = round(x);t=rectangle('Position',x);t.EdgeColor='r';t.LineWidth = 2;

for k = 1:200
    tic_lin(k) = mean2(R1_contrast(x(2):x(2)+x(4), x(1):x(1)+x(3),k));
    %tic_lin_1(k) = mean2(R1_contrast(x_2(2):x_2(2)+x_2(4),x_2(1):x_2(1)+x_2(3),k));
    %tic_lin_2(k) = mean2(R1_contrast(x_3(2):x_3(2)+x_3(4),x_3(1):x_3(1)+x_3(3),k));
end



%savefig(['Contrast_ROI_final' '.fig'])

%% Overlaying time-intensity plots given ROI selections made over ultrasound stilled_image

% Create figure
figure('Color','k');
perf1 = plot(timeStamps(1:200),smooth(tic_lin(1:200),'sgolay',1));
perf1.Color = [0 1 0];              % set line color to red  
hold on;
% perf2 = plot(timeStamps(1:200),smooth(tic_lin_1(1:200),'sgolay',1));
% perf2.Color = [1 0 0];              % set line color to red
% hold on;
% perf3 = plot(timeStamps(1:200),smooth(tic_lin_2(1:200),'sgolay',1));
% perf3.Color = [0 1 1];              % set line color to red
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
savefig(['Perfusion_arrival' num2str(frame) '.fig'])  % Save figure title(num2str(x),'Color','w');



%% setting x1/y1 variables %%

x1 = timeStamps(1:200); % taking timestamps from 1 - 120 given smoothed y1 values 
y1 = smooth(tic_lin(1:200),'sgolay',1);
y_1 = tic_lin(1:200);
%thresh = mean(y1(4:2:45));



%% Fit: 'POLY fit 7/8'.
[xData, yData] = prepareCurveData(x1', y1 );
% [xData_2, yData_2] = prepareCurveData(x2', y2 );
% [xData_3, yData_3] = prepareCurveData(x3', y3 );

% Find coefficients for polynomial 
fitResults1 =  fit(xData, yData, 'smoothingspline');
% fitResults2 =  fit(xData_2, yData_2, 'smoothingspline');
% fitResults3 =  fit(xData_3, yData_3, 'smoothingspline');

% evaluate the fitted y-values
yplot1 = feval(fitResults1,x1');
% yplot2 = feval(fitResults2,x2');
% yplot3 = feval(fitResults3,x3');


%% Extract Rise and Peak time points given yplot1,x1 %%

%% Extract Rise and Peak time points given yplot1,x1 %%

S = stepinfo(yplot1,x1', 'RiseTimeThreshold', [0.05 0.95]);
SysPeakTime = S.PeakTime;
SysPeak = S.Peak;
SysRiseTime = S.RiseTime;

% set max_index(peak) to 0
max_index = 0;
MTT_value = mean(yplot1);

%% for loop to extract range of index values up to peak point of curve %%
for i = 1:length(yplot1)
    if yplot1(i) == max(yplot1)
        max_index=i;
    end
end

mtt_index = 0;
for index_mtt = max_index:length(yplot1) % starting from the max_index point:+1 length of y-values

    if (yplot1(index_mtt) <= MTT_value) % calling all values post-max_index(i.e., peak values) onward.
        mtt_index = index_mtt;
        break % making sure it stops after it hits first value < MTT_value
    end
    
    
end

if(mtt_index == 0)
    mtt_index = length(yplot1);
end



plus_50 = x1(mtt_index);

%disp(mean_transit_time);
% take 10% of peakI threshold for global_thresh %
%global_thresh = abs(SysPeak*.05);
fifty_percent_minus_thresh = abs(SysPeak*.5);
%fifty_percent_plus = abs(SysPeak*1.5);
%global_thresh_2 = abs(SysPeak*.1);


% Other step response results of interest can be found by looking in
% the stepResults structure
% Find the first index where time exceeds the settling time
% To further improve this, you could interpolate between the points.
indexRiseT = find(x1 >= SysRiseTime,1,'first');
xRiseT = x1(indexRiseT);
yRiseT = yplot1(indexRiseT);

indexPeakT = find(x1 >= SysPeakTime,1,'first');
xPeakT = x1(indexPeakT);
yPeakT = yplot1(indexPeakT);




%% ROI timepoint - Global thresh #1 %%
at_thresh = abs(SysPeak*.075);
% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
%at = fzero(@(x1) feval(fitResults1,x1) - at_thresh, x1(1));

% Find 1st index where y <= thresh
at_idx = find( yplot1 >= at_thresh,1 ); 

% Get the x value at this index
at = x1( at_idx );
 





% interpolates to find yi, the values of the underlying function Y at the points in the vector or array xi. x must be a vector. 
Time_points = interp1(yplot1, x1, y1);



%% ROI 1 timepoint %%
% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
%xDesired = fzero(@(x1) feval(fitResults1,x1) - fifty_percent_minus_thresh, x1(1));

% Find 1st index where y >= thresh ~ x-number of points until yplot1 >= to
% threshold
minus_50_idx = find( yplot1 >= fifty_percent_minus_thresh,1 ); 

% Get the x value at this index
minus_50 = x1( minus_50_idx );


% Get area under the curve from fitted_y values given x 
ROI1_AUC = trapz(x1, log10(yplot1));


max_y = max(yplot1);
% max_y_2 = max(yplot2);
% max_y_3 = max(yplot3);

%x_arrival = xDesired;
% x_arrival_2 = xDesired_2;
% x_arrival_3 = xDesired_3;

mtt = x1(mtt_index);



% Plot step response
figure( 'Name', 'Curvefit1_poly' );
plot(x1',yplot1); 
hold on;
plot(x1', y_1, 'r*');
% Plot settling time point
plot(xPeakT,yPeakT,'MarkerEdgeColor',[0 0 1],'MarkerSize',23,'Marker','o');
hold on;
plot(minus_50, fifty_percent_minus_thresh, 'MarkerEdgeColor',[1 0 1],'MarkerSize',23,'Marker','o');
hold on;
plot(plus_50, MTT_value, 'MarkerEdgeColor',[0 1 1],'MarkerSize',23,'Marker','o');
hold on;
plot(at, at_thresh, 'MarkerEdgeColor',[1 0 0],'MarkerSize',23,'Marker','o');
hold on;
% Plot chart lines to settling point
plot([0 xPeakT],[yPeakT yPeakT],'r-.', 'LineWidth', 1);
plot([xPeakT xPeakT],[0 yPeakT],'black-.', 'LineWidth', 1);
hold on;
line([at at], [0 at_thresh], 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 1 ); % plot arrival point in dash vertical line 
plot([0 at], [at_thresh at_thresh],  'Color', 'g', 'LineStyle', '-.', 'LineWidth', 1 );
hold on;
line([minus_50 minus_50], [0 fifty_percent_minus_thresh], 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 1 ); % plot arrival point in dash vertical line 
plot([0 minus_50], [fifty_percent_minus_thresh fifty_percent_minus_thresh],  'Color', 'g', 'LineStyle', '-.', 'LineWidth', 1 );
hold on;
line([plus_50 plus_50], [0 MTT_value], 'Color', 'black', 'LineStyle', '-.', 'LineWidth', 1 ); % plot arrival point in dash vertical line 
plot([0 plus_50], [MTT_value MTT_value],  'Color', 'g', 'LineStyle', '-.', 'LineWidth', 1 );
hold on;
% H1=area(x,y1,'FaceColor',[1 1 1]);
% hold on
% idx=x>1&x<13;
% H=area(x(idx),y1(idx));
% set(H(1),'FaceColor',[1 0.5 0]);
grid on;
legend('sgolay','raw data', 'at', 'minus_50', 'xPeakT', 'plus_50',    'Location', 'NorthEast' );
% Label axes %
xlabel('Time (sec)');
ylabel('Linear Intensity');
title('TIC feature graph','Color','k')
hold on;
savefig(['TIC Graph'  '.fig'])

%%
idx = max_index;
figure('Color','k');
imagesc(Xlab, Zlab, 10*log10(R1_contrast(:,:,idx)),[115 145]); colormap('gray');
title([num2str(idx) ' ' num2str(timeStamps(idx))],'Color','w');
ax = gca;
h=xlabel('Lateral(mm)'); %or h=get(gca,'xlabel')
set(h, 'FontSize', 20); 
set(h,'FontWeight','bold');%bold font
z=ylabel('Depth(mm)');
set(z, 'FontSize', 20);
set(z, 'FontWeight', 'bold');
box(ax,'off'); 
set(gca, 'XColor', 'white', 'FontSize', 14); % set x-axis color to white
set(gca, 'YColor', 'white', 'FontSize', 14);
set(gca, 'FontWeight', 'bold');
hold on;
savefig(['Perfusion_arrival_peak_100%' num2str(idx) '.fig'])




idx2 = minus_50_idx;
figure('Color','k');
imagesc(Xlab, Zlab, 10*log10(R1_contrast(:,:,idx2)),[115 145]); colormap('gray');
title([num2str(idx2) ' ' num2str(timeStamps(idx2))],'Color','w');
ax1 = gca;
h=xlabel('Lateral(mm)'); %or h=get(gca,'xlabel')
set(h, 'FontSize', 20); 
set(h,'FontWeight','bold');%bold font
z=ylabel('Depth(mm)');
set(z, 'FontSize', 20);
set(z, 'FontWeight', 'bold');
box(ax1,'off'); 
set(gca, 'XColor', 'white', 'FontSize', 14); % set x-axis color to white
set(gca, 'YColor', 'white', 'FontSize', 14);
set(gca, 'FontWeight', 'bold');
hold on;
savefig(['Perfusion_peak_50%(-)' num2str(idx2) '.fig'])



idx3 = mtt_index;
figure('Color','k');
imagesc(Xlab, Zlab, 10*log10(R1_contrast(:,:,idx3)),[115 145]); colormap('gray');
title([num2str(idx3) ' ' num2str(timeStamps(idx3))],'Color','w');
ax2 = gca;
h=xlabel('Lateral(mm)'); %or h=get(gca,'xlabel')
set(h, 'FontSize', 20); 
set(h,'FontWeight','bold');%bold font
z=ylabel('Depth(mm)');
set(z, 'FontSize', 20);
set(z, 'FontWeight', 'bold');
box(ax2,'off'); 
set(gca, 'XColor', 'white', 'FontSize', 14); % set x-axis color to white
set(gca, 'YColor', 'white', 'FontSize', 14);
set(gca, 'FontWeight', 'bold');
hold on;
savefig(['Perfusion_peak_50%(+)' num2str(idx3) '.fig'])
hold off;
end


















































































































