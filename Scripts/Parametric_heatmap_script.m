%% Scaling image from pixel distance to mm - based on input file %%

imagedim=PData(2).Size;
Z=PData(2).PDelta(3)*[0:imagedim(1)-1];
X=PData(2).PDelta(1)*[0:imagedim(2)-1]+PData(2).Origin(1);
UFXAxis=X;UFZAxis=Z; 
X_1 = X./10;
Z_1 = Z/10;
close gcf force


%% pick frame of choice and run grid overlay %%
img = 37; % setting frame number to overlay heatmap onto
figure;imagesc(10*log10(R1_contrast(:,:,img)),[110 148]); % calling figure, log-transforming matrix-points, and setting dynamic range for image [# #]
colormap(gray); % setting colormap to gray

% global variables set %
ix=1; iy=1; xDesired = 0; xDesired_1 = 0; xDesired_2 = 0; SysPeakTime=0; SysPeak = 0; ROI1_AUC=0; washin = 0; mean_transit_time =0; % global variables initiated
for x_img =  55:4:470 % x-axis range (by units of 4)
    for z_img = 140:4:305 % y-axis or depth range (again, by units of 4)
        x = [x_img z_img 4 4];t=rectangle('Position',x);t.EdgeColor='r';
        for k = 1:120 % 388 'timepoints' from timestamp file
            tic_lin(k) = mean2(R1_contrast(x(2):x(2)+x(4),x(1):x(1)+x(3),k)); %Taking the average of each RGB point relative to said frame
        end
        
        %call time_delay function to get time_point delay values
        [S, SysPeakTime, SysPeak, xDesired, xDesired_1, xDesired_2, idx, idx_1, idx_2,  yplot1, ROI1_AUC, washin, mean_transit_time] = time_delay_plus_peaktime(tic_lin, timeStamps);
       
        tp(iy,ix)  = SysPeakTime;
        peak(iy,ix) = SysPeak;
        WTT(iy,ix) = washin;
        at(iy,ix)  = xDesired; % global_thresh
        at_1(iy,ix)  = xDesired_1; % global_thresh 2
        at_2(iy,ix)  = xDesired_2; % global_thresh 3
        AUC(iy,ix) = ROI1_AUC;
        iy = iy + 1;
    end
    disp(x)
    iy = 1;
    ix = ix + 1;
    disp(at);
end 
%savefig(['' '.fig']);


AUC_log = log10(AUC); % log-transforming AIU values for scale

% TTP troubleshoot %
% TTP = tp - at;
% TTP_1 = tp - at_1;
% TTP_2 = tp - at_2;


%% Subplot of each feature of interest - unblended %%

%% AT %%
% output arrival_time parametric heatmap %%
subplot(4,1,1)
imagesc(at);
c = colorbar; c.Label.String = 'Arrival Time (sec)';
title('arrival time (at) at 5% ');
colormap('jet');caxis([1 3.5]);

%% TTP %%

subplot(4,1,2)
imagesc(tp);
c = colorbar; c.Label.String = 'Time to Peak (sec)';
title('TTP');
colormap('jet');caxis([3.5 8.5]);

%% WIT %%

subplot(4,1,3)
imagesc(WIT);
c = colorbar; c.Label.String = 'Wash-in Time (sec)';
title('WIT%');
colormap('jet');caxis([0 2.5]);

%% AUC %%

subplot(4,1,4)
imagesc(AUC_log);
c = colorbar; c.Label.String = 'log10 (AIU)';
title('AUC');
colormap('jet');caxis([12 15]);



%% output AUC parametric heatmap - unblended %%
figure;imagesc(AUC_log);colorbar;%title('Peak Intensity');
colormap('jet');caxis([11.9 14.9]);
%savefig(['BL_Demustached_HM_1_60' '.fig']);

%% mask over regions that want to be left as nan (transparent) %%
h_2 = imfreehand();
position_2 = wait(h_2);
map_1 = createMask(h_2);
AUC_log(map_1) = nan;
imagesc(AUC_log);




%% code to blend heatmap %%
minv = 1;%min(min(R1_perf(:,:,29)));
maxv = 6;%max(max(R1_perf(:,:,2 t9)));
map=colormap('jet');
ncol = size(map,1);
s = round(1+(ncol-1)*(TTP-minv)/(maxv-minv)); % Taking arrival time values and rounding differences
rgb_at = ind2rgb(s,map);
rgb_at = imresize(rgb_at, 6);

%% freehand over blend image to remove any regions with R, G, or B values that are <= 0.5 %%
imshow(rgb_at);
h_2 = imfreehand();
position_2 = wait(h_2);
map_2 = createMask(h_2);
rgb_at(map_2) = nan;
imagesc(rgb_at);



%% Resizing image as well as blending to a particular 'transparency' value (i.e., 'alpha'). 
rgb_at_scale  = imresize(rgb_at,[60 480],'nearest');
rgb_at_scale_2 = rgb_at_scale;
toto          = zeros(size(rgb_at_scale_2));
alpha = 0.64;
rgb_blend = fliplr(alpha * (rgb_at_scale) + (1 - alpha) * toto.*(0.5));
rgb_blend_final = rgb_blend;


% loop over all rows and columns to adjust RGB channel values after
% blending above - remove any leftover values <= # 
for ii=1:size(rgb_blend_final,1)
    for jj=1:size(rgb_blend_final,2)
        % get pixel value
        pixel_r=rgb_blend_final(ii,jj,1);
        pixel_g=rgb_blend_final(ii,jj,2);
        pixel_b=rgb_blend_final(ii,jj,3);
        pixel_rgb =rgb_blend_final(ii,jj);
        
          % check pixel value and assign new value
          if (pixel_r<0.3251 && pixel_r>0.3249 && pixel_g ==0 && pixel_g == 0)
              new_pixel=nan;
          
          %elseif (pixel_r<0.51 && pixel_r>0.49 && pixel_g ==0 && pixel_g == 0)
              %new_pixel=nan;
          else
              new_pixel = pixel_rgb;
          end

          % save new pixel value in thresholded image
          rgb_blend_final(ii,jj)=new_pixel;
     end
end




%% Blending heatmap with orginal_image %%
img = 42;
figure('Color', 'k');
ax1 = axes;
imagesc(X_1, Z_1, 10*log10(fliplr(R1_contrast(:,:,img))),[110 148]);
colormap(ax1,'gray');
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
colorbar;
hold on;
% blended_img_additions % 
im = image(rgb_blend_final, 'XData', [-5 5.2], 'YData', [3.87 7.78]); % adjusting coordinates of HM over original image
ax2 = axes;
im.AlphaData = max(rgb_blend_final,[],3);% Setting the alpha blend to it's highest transparency value '3'
colormap(ax2, 'jet');caxis([12 15]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
% Create colorbar
c = colorbar('peer',ax2,'FontWeight', 'bold','FontSize',18,...
    'Color',[0.941176470588235 0.941176470588235 0.941176470588235]);
c.Label.String = 'log10 (AIU)';
savefig(['AUC_heatmap_test' '.fig']);
hold off




