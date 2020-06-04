%% Scaling image from pixel distance to mm - based on input file %%
imagedim=PData(2).Size;
Z=PData(2).PDelta(3)*(0:imagedim(1)-1);
X=PData(2).PDelta(1)*(0:imagedim(2)-1)+PData(2).Origin(1);
UFXAxis=X;UFZAxis=Z; 
X_1 = X./10;
Z_1 = Z/10;
close gcf force


%% pick frame of choice and run grid overlay %%
img = 44; % setting frame number to overlay heatmap onto
figure;imagesc(10*log10(R1_contrast(:,:,img)),[110 145]); % calling figure, log-transforming matrix-points, and setting dynamic range for image [# #]
colormap(gray); % setting colormap to gray


% global variables set %
ix=1; iy=1; xDesired = 0; xDesired_1 = 0; xDesired_2 = 0; SysPeakTime=0; SysPeak = 0; ROI1_AUC=0; washin = 0; mean_transit_time =0; % global variables initiated
for x_img =  30:3:470 % x-axis range (by units of 4)
    for z_img = 175:3:330 % y-axis or depth range (again, by units of 4)
        x = [x_img z_img 3 3];t=rectangle('Position',x);t.EdgeColor='r';
        for k = 1:150 % 388 'timepoints' from timestamp file
            tic_lin(k) = mean2(R1_contrast(x(2):x(2)+x(4),x(1):x(1)+x(3),k)); %Taking the average of each RGB point relative to said frame
        end
        
        %call time_delay function to get time_point delay values 
        [S, SysPeakTime, SysPeak, xD, xD1, xD2, xD3, xDesired, xDesired_1, xDesired_2, xDesired_3, idx, idx_1, idx_2, idx_3,  yplot1, ROI1_AUC, washin, mean_transit_time] = time_delay_plus_peaktime_100_acute(tic_lin, timeStamps);
       
        tp(iy,ix)  = SysPeakTime;
        peak(iy,ix) = SysPeak;
        %WTT(iy,ix) = washin;
        at(iy,ix)  = xDesired; % global_thresh - 15% of 1:maxindex
        at_1(iy,ix)  = xDesired_1; % global_thresh 2 - mean(abs(20:95))
        at_2(iy,ix)  = xDesired_2; % global_thresh 3 - mean(abs(SysPeak*.075))
        at_3(iy,ix)  = xDesired_3; % global_thresh 4 - mean(abs(SysPeak*.125))
        %at_4(iy,ix)  = xDesired_4; %global_thresh 0.5 - mean(abs(1:120))
        %AUC(iy,ix) = ROI1_AUC;
        iy = iy + 1;
    end
    disp(x)
    iy = 1;
    ix = ix + 1;
    disp(at);
end 
%savefig(['' '.fig']);
%AUC_log = log10(AUC); % log-transforming AIU values for scale

% saving min value of AT matrix for thresholding %
min_mat_value_3 = min(at_3(:));
min_mat_value_2 = min(at_2(:));
min_mat_value_1 = min(at_1(:));
min_mat_value = min(at(:));


%% output parametric heatmap - unblended %%
figure;imagesc(at_2);colorbar;%title('Peak Intensity');
colormap('jet');caxis([1 3.5]);
at_2(at_2 <= min_mat_value) = NaN;
%at_3(at_3 <= 0.3) = NaN;
%imagesc(at_2);
imAlpha=ones(size(at_2));
imAlpha(isnan(at_2))=NaN;
imagesc(at_2,'AlphaData',imAlpha);

%%%% Viewing different at values %%%%

% AT %
% output arrival_time parametric heatmap %%
subplot(4,1,1)
imagesc(at);
c = colorbar; c.Label.String = 'Arrival Time (sec)';
title('arrival time (at) at 5% ');
colormap('jet');caxis([1 3.5]);

% AT 2 %

subplot(4,1,2)
imagesc(at_1);
c = colorbar; c.Label.String = 'Time to Peak (sec)';
title('TTP');
colormap('jet');caxis([1.1 3.6]);

% AT 3 %

subplot(4,1,3)
imagesc(at_2);
c = colorbar; c.Label.String = 'Wash-in Time (sec)';
title('WIT%');
colormap('jet');caxis([1.5 4]);

% AT 4 %

subplot(4,1,4)
imagesc(at_3);
c = colorbar; c.Label.String = 'log10 (AIU)';
title('AUC');
colormap('jet');caxis([1 3.5]);




%% code to blend heatmap %%
minv = 1;%min(min(R1_perf(:,:,29)));
maxv = 3.5;%max(max(R1_perf(:,:,2 t9)));
map=colormap('jet');
ncol = size(map,1);
s = round(1+(ncol-1)*(at_3-minv)/(maxv-minv)); % Taking arrival time values and rounding differences
rgb_at = ind2rgb(s,map);
rgb_at = imresize(rgb_at, 3.5);
figure; imagesc(rgb_at);



%% Resizing image as well as blending to a particular 'transparency' value (i.e., 'alpha'). 
%rgb_at_scale  = imresize(rgb_at,[30 460],'nearest');
rgb_at_scale_2 = rgb_at;
toto          = zeros(size(rgb_at_scale_2));
alpha = 0.6;
rgb_blend = fliplr(alpha * (rgb_at_scale_2) + (1 - alpha) * toto.*(0.55));
rgb_blend_final = rgb_blend;
rgb_blend_final (rgb_blend_final <=0.4) = NaN;
figure;imagesc(rgb_blend_final);



%% loop over all rows and columns to adjust RGB channel values after %%
% blending above - remove any leftover values <= # 
for ii=1:size(rgb_blend_final,1)
    for jj=1:size(rgb_blend_final,2)
        % get pixel value
        pixel_r=rgb_blend_final(ii,jj,1);
        pixel_g=rgb_blend_final(ii,jj,2);
        pixel_b=rgb_blend_final(ii,jj,3);
        pixel_rgb =rgb_blend_final(ii,jj);
        
          % check pixel value and assign new value 
          if (pixel_r<=0.325 && isnan(pixel_g) && isnan(pixel_b))
              new_pixel=nan;
              
          elseif (pixel_r<=0.1 && pixel_g <=0.1 && pixel_b <= 0.1)
              new_pixel=nan;
         
          elseif (pixel_r<=0.3 && isnan(pixel_g) && pixel_b <=0.3)
              new_pixel=nan;
           else
              new_pixel = pixel_rgb;
          end

          % save new pixel value in thresholded image
          rgb_blend_final(ii,jj)=new_pixel;
     end
end


figure;imagesc(rgb_blend_final);

%% Blending heatmap with orginal_image %%
img = 44;
figure('Color', 'k');
ax1 = axes;
imagesc(X_1, Z_1, 10*log10(fliplr(R1_contrast(:,:,img))),[115 145]);
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
hold on;
% blended_img_additions %
im = imagesc(rgb_blend_final, 'XData', [-5.35 4.81], 'YData', [4.51 8.42]); % adjusting coordinates of HM over original image
%set(im, 'AlphaData', ~isnan(rgb_blend_final));
ax2 = axes;
im.AlphaData = max(rgb_blend_final,[],3);% Setting the alpha blend to it's highest transparency value '3'
colormap(ax2, 'jet');caxis([0 2.5]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
% Create colorbar
c = colorbar('peer',ax2,'FontWeight', 'bold','FontSize',20,...
    'Color',[0.941176470588235 0.941176470588235 0.941176470588235]);
c.Label.String = 'Time (sec)';
%savefig(['at_heatmap_test_2' '.fig']);
hold off




