% load('/Volumes/CEUS_Data_UW/Survival_group/Acute/Rat_97/97_acute_demustache/20171026T112833_thresh400_videos/R1_contrast.mat');
% load('/Volumes/CEUS_Data_UW/Survival_group/Acute/Rat_97/20171026T112833_movie_inj/20171026T112833_timestamps.mat');


%% pick frame of choice and run grid overlay %%
img = 42;
figure;imagesc(10*log10(R1_contrast(:,:,img)),[115 150]); % calling figure, log-transforming matrix-points, and setting dynamic range for image [# #]
colormap(gray); % setting colormap to gray



ix=1; iy=1; xDesired = 0; SysPeakTime=0; % global variables initiated
for x_img =  40:4:480 % x-axis range (by units of 4)
    for z_img = 165:4:300 % y-axis or depth range (again, by units of 4)
        x = [x_img z_img 4 4];t=rectangle('Position',x);t.EdgeColor='r';
        for k = 1:120 % 388 'timepoints' from timestamp file
            tic_lin(k) = mean2(R1_contrast(x(2):x(2)+x(4),x(1):x(1)+x(3),k)); %Taking the average of each RGB point relative to said frame
        end
        
        %call time_delay function to get time_point delay values
        [S, SysRiseTime, SysPeakTime, SysPeak, xDesired, idx, yplot1] = time_delay_plus_peaktime(tic_lin, timeStamps)
       
        tp(iy,ix)  = SysPeakTime;
        %mtt(iy,ix) = MTT;
        at(iy,ix)  = xDesired;
        iy = iy + 1;
    end
    disp(x);
    iy = 1;
    ix = ix + 1;
    disp(at);
end 
savefig(['Acute_grid_demustached_peak' '.fig']);



% output arrival_time parametric heatmap %%
figure;imagesc(at);colorbar;title('arrival time');
colormap('jet');caxis([2.5 5]);
savefig(['Acute_Demustached_HM_1_55_peak' '.fig']);

% output peak_time parametric heatmap %
figure;imagesc(tp);colorbar;title('arrival time');
colormap('jet');caxis([4 6]);
savefig(['Acute_Demustached_HM_1_55_peaktime' '.fig']);



%% output arrival_time parametric heatmap %%
figure;imagesc(at);colorbar;title('arrival time');
colormap('jet');caxis([2 4]);
%savefig(['Acute_1e13_draft_DRAFT_ignore' '.fig']);


%% code to blend heatmap %%
minv = 2.5;%min(min(R1_perf(:,:,29)));
maxv = 5;%max(max(R1_perf(:,:,2 t9)));
map=colormap('jet');
ncol = size(map,1);
s = round(1+(ncol-1)*(at-minv)/(maxv-minv)); % Taking arrival time values and rounding differences
rgb_at = ind2rgb(s,map);
rgb_at = imresize(rgb_at,5);
rgb_perf = ind2rgb(s,map);
rgb_perf = imresize(rgb_perf,5);
rgb_at_scale  = imresize(rgb_at,[40 480],'nearest');
%rgb_at_scale_2  = imresize(rgb_at,[170 220],'nearest');
toto          = zeros(size(rgb_at_scale));
alpha = 0.65;
rgb_blend = fliplr(alpha * rgb_at_scale + (1 - alpha) * toto);






%% code to blend heatmap %%
minv = 4;%min(min(R1_perf(:,:,29)));
maxv = 6;%max(max(R1_perf(:,:,2 t9)));
map=colormap('jet');
ncol = size(map,1);
s = round(1+(ncol-1)*(tp-minv)/(maxv-minv)); % Taking arrival time values and rounding differences
rgb_at = ind2rgb(s,map);
rgb_at = imresize(rgb_at,6);
rgb_perf = ind2rgb(s,map);
rgb_perf = imresize(rgb_perf,6);
rgb_at_scale  = imresize(rgb_at,[40 480],'nearest');
%rgb_at_scale_2  = imresize(rgb_at,[170 220],'nearest');
toto          = zeros(size(rgb_at_scale));
alpha = 0.65;
rgb_blend = fliplr(alpha * rgb_at_scale + (1 - alpha) * toto);


%% Blending heatmap with orginal_image %%
img = 42;
figure('Color', 'k');
ax1 = axes;
imagesc(10*log10(fliplr(R1_contrast(:,:,img))),[115 150]);
title([num2str(img) ' ' num2str(timeStamps(img))],'Color','w');
colormap(ax1,'gray');
ax = gca;
set(gca,'XTickLabel',{'-4', '-3', '-2', '-1', '0', '1', '2', '3' '4' '5'}, 'FontSize',20,'FontWeight','bold');% get the current axis
set(gca,'YTickLabel',1:10, 'FontSize',20,'FontWeight','bold');                   % get the current axis
h=xlabel('Distance(mm)'); %or h=get(gca,'xlabel')
set(h, 'FontSize', 22); 
set(h,'FontWeight','bold');%bold font
z=ylabel('Width(mm)');
set(z, 'FontSize', 22);
set(z, 'FontWeight', 'bold');
box(ax,'off'); 
set(gca, 'XColor', 'white'); % set x-axis color to white
set(gca, 'YColor', 'white');
set(gca, 'FontWeight', 'bold');
%ttl = ax.Title;
%set(ttl, 'FontSize', 28);
hold on;
%%blended_img_additions%%
im = image(rgb_blend, 'XData', [40 470], 'YData', [165 300]);
ax2 = axes;
im.AlphaData = max(rgb_blend,[],3);% Setting the alpha blend to it's highest transparency value '3'
colormap(ax2, 'jet');caxis([4 6]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
% Create colorbar
c = colorbar('peer',ax2,'FontWeight', 'bold','FontSize',20,...
    'TickLabels',{'4','4.5','5','5.5','6'},...
    'Color',[0.941176470588235 0.941176470588235 0.941176470588235]);
c.Label.String = 'seconds';
savefig(['test_blend_mask_peak' '.fig']);
hold off




%%Save Frame as fig file %%
img = 42;
figure('Color', 'k');
ax1 = axes;
imagesc(10*log10(fliplr(R1_contrast(:,:,img))),[115 150]);
title([num2str(img) ' ' num2str(timeStamps(img))],'Color','w');
colormap(ax1,'gray');
ax = gca;
set(gca,'XTickLabel',{'-4', '-3', '-2', '-1', '0', '1', '2', '3' '4' '5'}, 'FontSize',20,'FontWeight','bold');% get the current axis
set(gca,'YTickLabel',1:10, 'FontSize',20,'FontWeight','bold');                   % get the current axis
h=xlabel('Distance(mm)'); %or h=get(gca,'xlabel')
set(h, 'FontSize', 22); 
set(h,'FontWeight','bold');%bold font
z=ylabel('Width(mm)');
set(z, 'FontSize', 22);
set(z, 'FontWeight', 'bold');
box(ax,'off'); 
set(gca, 'XColor', 'white'); % set x-axis color to white
set(gca, 'YColor', 'white');
set(gca, 'FontWeight', 'bold');
ttl = ax.Title;
set(ttl, 'FontSize', 28);
colorbar('Color','w');
%%%%%% Reduce Whitespace around figure %%%%%%
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
savefig(['Acute_frame' '.fig']);





%'TickLabels',{'2.5','3','3.5','4','4.5','5'}
% 
% at_2 = at;
% % % thresholding values in at that are <= 1 and setting to 0
% at_2(find(at_2 <= 0.5)) = 1.75;