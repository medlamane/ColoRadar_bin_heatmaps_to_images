% This routine reads the binary files in the cascade radar sensor dataset.

% **********************************************************************

%% Define User values

source_dir           = '../12_21_2020_ec_hallways_run3/';
device_name          = 'cascade/';
adc_data_path        = 'adc_samples/data/';
adc_data_suffix      = '.bin';
adc_time_path        = 'adc_samples/';
adc_time_suffix      = '.txt';

heatmap_path         = 'heatmaps/data/';
heatmap_suffix       = '.bin';
heatmap_time_path    = 'heatmaps/';
heatmap_time_suffix  = '.txt';

%% Define cascade radar sensor parameters (from calib.zip)

% Define the heatmap configuration (from the heatmap_cfg.txt file)
num_range_bins       = 128;
num_elevation_bins   = 32;
num_azimuth_bins     = 128;
heatmap_parm         = [num_range_bins, num_elevation_bins, num_azimuth_bins];
% for plotting:
range_bin_width      = 0.0592943951488;

%% Define frames to read
    n=1; % first frame to read
	nframe = 1; % reading from n to nframe; ( nframe>=n )
	%heatmap_angle_range = ones(num_range_bins, num_azimuth_bins,nframe-n);
% process individual frames
for frame_index    = n:nframe
      
   %% Define the full filenames
   
   adc_data_filename       = [source_dir,device_name,adc_data_path,'frame_',num2str(frame_index),adc_data_suffix];
   adc_time_filename       = [source_dir,device_name,adc_time_path,'timestamps',adc_time_suffix];
   heatmap_data_filename   = [source_dir,device_name,heatmap_path,'heatmap_',num2str(frame_index),heatmap_suffix];
   heatmap_time_filename   = [source_dir,device_name,heatmap_time_path,'timestamps',heatmap_time_suffix];
   
   %% Proceed if the files exist
   good_heatmap_data_filename = exist(heatmap_data_filename,'file');
   good_heatmap_time_filename = exist(heatmap_time_filename,'file');
   
   if((good_heatmap_data_filename == 2) && (good_heatmap_time_filename == 2))
      
       fid   = fopen(adc_data_filename,'r');
      [input_data, num_cnt] = fread(fid,Inf,'int16');
     
      
      
      %% Read the heatmap file
      % Read the binary file as a string of 16 bit integers   
      % open the file
      fid   = fopen(heatmap_data_filename,'r');
      [input_heatmap_data, num_cnt] = fread(fid,Inf,'single'); % 32 bit float
      
      % There should be two value per heatmap location: intensity & range rate
      expected_num_samples = 2*(num_range_bins * num_elevation_bins * num_azimuth_bins);
      
  %    disp(['Expected number of samples: ',num2str(expected_num_samples)]);
  %    disp(['Read number of samples:     ',num2str(num_cnt)]);
      
      %% Convert from serial format to matrix format     
      heatmap_intensity  = ones(num_range_bins, num_azimuth_bins, num_elevation_bins) .* NaN;
      heatmap_range_rate = ones(num_range_bins, num_azimuth_bins, num_elevation_bins) .* NaN;    
      for index_range = 1:num_range_bins
         % get 0-ref index
         r = index_range - 1;
         for index_azimuth = 1:num_azimuth_bins
            % get 0-ref index
            a = index_azimuth - 1;
            for index_elevation = 1:num_elevation_bins
               % get 0-base index
               e = index_elevation - 1;
               index_intensity = 2* (r + num_range_bins * (a + num_azimuth_bins*e)) + 1;
               index_range_rate = index_intensity + 1;
               % save the recorded value
               heatmap_intensity(index_range, index_azimuth, index_elevation) = input_heatmap_data(index_intensity);
               heatmap_range_rate(index_range, index_azimuth, index_elevation) = input_heatmap_data(index_range_rate);
			   %heatmap_angle_range(index_range, index_azimuth, frame_index-n+1) = heatmap_intensity(index_range, index_azimuth, num_elevation_bins/32);
                
            end % end for index_elevation
         end % end for index_azimuth
      end % end for index_range
      %% Read the ADC sample time 
   end % end if(good_filename == 2)
     A = heatmap_intensity(:, :, num_elevation_bins/32);
B = zeros(128,128);
    for i=1:128
    for k=1:128 
        azimuth = 0.4+i*1.05*2*pi/360; % 0.4 radian = 22.5° is offset to center the y axis || 1.05 is Azimuth resolution
        r = 0.117*k; % 0.117 is Range resolution
   [x,y,z] = sph2cart(azimuth,0,r);
   x=64+round(64*x/15);
   y=1+round(127*y/15);
B(y,x) = A(k,i);
    end
    end
    y=linspace(0,-15,129);
    x=linspace(-7.5,7.5,128);
    B(1,1)=0;
    B(1,2)=9*10^5;  
ha=get(gcf,'children');
    h=imagesc(x,y,B);
    colormap jet
title('X-Y Heatmap')
ylabel('Y axis')
xlabel('X axis')
filenm = strcat('imgrd/',int2str(frame_index),'.jpg');
saveas(h,filenm);
close all
end % end for frame_index loop






