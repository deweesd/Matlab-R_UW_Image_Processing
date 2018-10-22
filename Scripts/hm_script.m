%% set up scaling grid %%
imagedim=PData(2).Size;
Z=PData(2).PDelta(3)*[0:imagedim(1)-1];
X=PData(2).PDelta(1)*[0:imagedim(2)-1]+PData(2).Origin(1);
UFXAxis=X;UFZAxis=Z; 
X_1 = X./10;
Z_1 = Z/10;
%close gcf force

%% pick frame of choice and run grid overlay %%
img = 42;
figure;imagesc(10*log10(R1_contrast(:,:,img)),[110 148]); % calling figure, log-transforming matrix-points, and setting dynamic range for image [# #]
colormap(gray); % setting colormap to gray

ix=1; iy=1; xDesired = 0; SysPeakTime=0; SysPeak = 0; ROI1_AUC=0; % global variables initiated
for x_img =  40:4:470 % x-axis range (by units of 4)
    for z_img = 140:4:290 % y-axis or depth range (again, by units of 4)
        x = [x_img z_img 4 4];t=rectangle('Position',x);t.EdgeColor='r';
        for k = 5:100 % 388 'timepoints' from timestamp file
            tic_lin(k) = mean2(R1_contrast(x(2):x(2)+x(4),x(1):x(1)+x(3),k)); %Taking the average of each RGB point relative to said frame
        end
        
        %call time_delay function to get time_point delay values
        [S, SysRiseTime, SysPeakTime, SysPeak, xDesired, idx, yplot1] = time_delay_plus_peaktime(tic_lin, timeStamps);
       
        %tp(iy,ix)  = SysPeakTime;
        peak(iy,ix) = SysPeak;
        %mtt(iy,ix) = MTT;
        at(iy,ix)  = xDesired;
        %AUC(iy,ix) = ROI1_AUC;
        iy = iy + 1;
    end
    disp(x)
    iy = 1;
    ix = ix + 1;
    disp(at);
end 
%savefig(['BL_grid_demustached_AUC' '.fig']);


% output arrival_time parametric heatmap %%
figure;imagesc(at);colorbar;title('arrival time');
colormap('jet');caxis([3 5.5]);
%savefig(['BL_Demustached_HM_1_60' '.fig']);
%% output arrival_time parametric heatmap %%
h_1 = imfreehand();
position_1 = wait(h_1);
map = createMask(h_1);
at(map) = nan;
imagesc(at);




% set up for peak intensity HM %
[maxValue, linearIndexesOfMaxes] = max(T(:));
[lowestValue, linearIndexesOfLowest] = min(peak(:));
T = 10*log10(peak);

% output arrival_time parametric heatmap %%
figure;imagesc(T);colorbar;title('Peak Intensity');
colormap('jet');caxis([123 141]);
%savefig(['BL_Demustached_HM_1_60' '.fig']);
%% output arrival_time parametric heatmap %%
h_2 = imfreehand();
position_2 = wait(h_2);
map = createMask(h_2);
raw_pos = h_2.getPosition();
T(map) = nan;
imagesc(T);









%% code to blend heatmap %%
minv = 3;%min(min(R1_perf(:,:,29)));
maxv = 5.5;%max(max(R1_perf(:,:,2 t9)));
map=colormap('jet');
ncol = size(map,1);
s = round(1+(ncol-1)*(at-minv)/(maxv-minv)); % Taking arrival time values and rounding differences
rgb_at = ind2rgb(s,map);
rgb_at = imresize(rgb_at, 5.5);
imshow(rgb_at);
h_3 = imfreehand();
position_3 = wait(h_3);
map_1 = createMask(h_3);
raw_pos_2 = h_3.getPosition();
rgb_at(map_1) = nan;
imagesc(rgb_at);

% resizing %
rgb_at_scale  = imresize(rgb_at,[40 470],'nearest');
rgb_at_scale_2 = rgb_at_scale;
toto          = zeros(size(rgb_at_scale_2));
alpha = 0.625;
rgb_blend = alpha * (rgb_at_scale) + (1 - alpha) * toto.*(0.5);
rgb_blend_final = rgb_blend;
%imshow(rgb_blend_final);
%imshow(rgb_at_scale_2);

% loop over all rows and columns
for ii=1:size(rgb_blend_final,1)
    for jj=1:size(rgb_blend_final,2)
        % get pixel value
        pixel_r=rgb_blend_final(ii,jj,1);
        pixel_g=rgb_blend_final(ii,jj,2);
        pixel_b=rgb_blend_final(ii,jj,3);
        pixel_rgb =rgb_blend_final(ii,jj);
        
          % check pixel value and assign new value
          if (pixel_r<0.31251 && pixel_r>0.31249 && pixel_g ==0 && pixel_g == 0)
              new_pixel=nan;
          
          elseif (pixel_r<0.3353 && pixel_r>0.2951 && pixel_g ==0 && pixel_g == 0)
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
imagesc(X_1,Z_1,10*log10(R1_contrast(:,:,img)),[110 148]);
title([num2str(img) ' ' num2str(timeStamps(img))],'Color','w');
colormap(ax1,'gray');
ax = gca;
%axis equal;
h=xlabel('Lateral(mm)'); %or h=get(gca,'xlabel')
set(h, 'FontSize', 22); 
set(h,'FontWeight','bold');%bold font
z=ylabel('Depth(mm)');
set(z, 'FontSize', 22);
set(z, 'FontWeight', 'bold');
box(ax,'off'); 
set(gca, 'XColor', 'white'); % set x-axis color to white
set(gca, 'YColor', 'white');
set(gca, 'FontWeight', 'bold');
hold on;
%%blended_img_additions%%
im = image(rgb_blend_final, 'XData', [-5.3 4.85], 'YData', [3.5 7.6]);
ax2 = axes;
im.AlphaData = max(rgb_blend_final,[],3);% Setting the alpha blend to it's highest transparency value '3'
colormap(ax2, 'jet');caxis([0 2.5]);
ax2.Visible = 'off';
linkprop([ax1 ax2],'Position');
% Create colorbar
c = colorbar('peer',ax2,'FontWeight', 'bold','FontSize',20,...
    'Color',[0.941176470588235 0.941176470588235 0.941176470588235]);
c.Label.String = 'seconds';
% h = imfreehand; position = wait(h); uiwait(msgbox('Locate the point'));
% [x,y] = ginput(1); hold on; % Prevent image from being blown away.
% plot(x,y,'r+', 'MarkerSize',30,'LineWidth',2);
savefig(['at_blend' '.fig']);
hold off



