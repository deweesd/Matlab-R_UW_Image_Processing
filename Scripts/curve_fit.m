function [xDesired, xDesired_2, idx, idx_2, yplot1, yplot2, thresh, thresh_2] = curve_fit(timeStamps, R1_contrast)
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
frame = 44;
figure('Color', 'k');imagesc((10*log10(R1_contrast(:,:,frame))),[115 150]); colormap('gray');
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
for k = 1:120
    tic_lin(k) = mean2(R1_contrast(x(2):x(2)+x(4), x(1):x(1)+x(3),k));
end

% %%% second ROI selection %%%
x_1 = getrect; 
x_1 = round(x_1);t=rectangle('Position',x_1);t.EdgeColor='r'; t.LineWidth = 2;
for k_1 = 1:120
    tic_lin_1(k_1) = mean2(R1_contrast(x_1(2):x_1(2)+x_1(4), x_1(1):x_1(1)+x_1(3),k_1));
end
savefig(['Frame2_ROI_ACUTE_DM_final' '.fig']);

%% Overlaying time-intensity plots given ROI selections made over ultrasound stilled_image

% Create figure
figure();
perf1 = plot(timeStamps(1:120),smooth(tic_lin(1:120),'sgolay',1));
perf1.Color = [0 1 0];              % set line color to red  
hold on;
perf2 = plot(timeStamps(1:120),smooth(tic_lin_1(1:120),'sgolay',1));
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
savefig(['Contrast_bolus_final'  '.fig']);
title(num2str(x_1),'Color','w');



%% setting x1/y1 variables %%

x1 = timeStamps(1:80); % taking timestamps from 1 - 120 given smoothed y1 values 
y1 = smooth(tic_lin(1:80),'sgolay',1);
thresh = mean(y1(1:60));

x2 = timeStamps(1:100); % taking timestamps from 1 - 120 given smoothed y1 values 
y2 = smooth(tic_lin_1(1:100),'sgolay',1);
thresh_2 = mean(y2(1:70));

%global_thresh = abs(thresh_2 - thresh);
%% Fit: 'POLY fit 7/8'.
[xData, yData] = prepareCurveData(x1', y1 );
[xData_2, yData_2] = prepareCurveData(x2', y2 );

% Find coefficients for polynomial (order = 4 and 6, respectively)
fitResults1 =  fit(xData, yData, 'smoothingspline');
fitResults2 =  fit(xData_2, yData_2, 'smoothingspline');

% evaluate the fitted y-values
yplot1 = feval(fitResults1,x1');
yplot2 = feval(fitResults2,x2');

% interpolates to find yi, the values of the underlying function Y at the points in the vector or array xi. x must be a vector. 
Time_points = interp1(yplot1, x1', y1);
Time_points_2 = interp1(yplot2, x2', y2);

%% ROI 1 timepoint %%

% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
xDesired = fzero(@(x1) feval(fitResults1,x1) - thresh, x1(1));

% Find 1st index where y >= thresh
idx = find( yplot1 >= thresh,1 ); 

% Get the x value at this index
xDesired = x1( idx );


%% ROI 2 timepoint %%

% Use fzero to get the root of y = a*x^n + b*x^(n-1) + ... + z when y = thresh
xDesired_2 = fzero(@(x2) feval(fitResults2,x2) - thresh_2, x2(1));

% Find 1st index where y >= thresh
idx_2 = find( yplot2 >= thresh_2,1 ); 

% Get the x value at this index
xDesired_2 = x2( idx_2 );

%% Plot fit with data while also plotting arrival delay point %%
figure( 'Name', 'Curvefit1_poly' );
h = plot(x1', y1);%smoothed-points
hold on;
plot(x1', yplot1);%polyfit points
hold on;
plot(Time_points, yplot1, '*g');%interpolated points of x given y
hold on;
plot(x2', y2);
hold on;
plot(x2', yplot2);
hold on;
plot(Time_points_2, yplot2, '*r');
max_y = max(yplot1);
max_y_2 = max(yplot2);
x_arrival = xDesired;
x_arrival_2 = xDesired_2;
diff_in_delay = xDesired_2 - xDesired;
hold on;
line([x_arrival x_arrival], [min(yplot1) max(yplot1)], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1.5 ); % plot arrival point in dash vertical line 
plot([x_arrival x_arrival], [idx idx],  'Color', 'r', 'LineStyle', '--', 'LineWidth', 1 ); % plot corresponding max(yi) horizontal line
hold on;
line([x_arrival_2 x_arrival_2], [min(yplot2) max(yplot2)], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1.5 ); % plot arrival point in dash vertical line 
plot([x_arrival_2 x_arrival_2], [idx_2 idx_2],  'Color', 'g', 'LineStyle', '--', 'LineWidth', 1 ); % plot corresponding max(yi) horizontal line
grid on;
legend( h, 'polyfit vs. timepoints', 'Location', 'NorthEast' );
% Label axes
xlabel Time
ylabel Intensity
title(['Delta_T' num2str(diff_in_delay) ],'Color','k');
grid on
savefig(['Curvefit_ROI_ACUTE_DM_final' '.fig']);


end


