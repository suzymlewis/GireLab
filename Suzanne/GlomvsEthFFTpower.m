%% IMPORT GLOM DATA AND ETH DATA

% Import glomerular data file
[fg, fdg] = uigetfile('.mat', 'CHOOSE GLOMULAR DATA - .mat');
cd(fdg);
g_data = load(fg);

trials = g_data.trials;
ODOR_ON = g_data.ODOR_ON; % Index for relevant analysis time of video 
    % Keeps baseline 10 seconds before odor presentation    
    % Dumps a few seconds at tail end when no odor present
        % This ensures all trials have the same length
C_or = g_data.C_or; % Temporal Matrix from CNMF segmentation
A_or = g_data.A_or; % Spatial Matrix from CNMF segmentation
C_dff = g_data.C_df; % Dff as defined by CNMF segmentation
session_length = g_data.y; % Number of trials to be analyzed over session
    % This number can vary if trials need to be dropped or are accidently
    % not recorded
frame_trial = g_data.tr; % What trial each frame is from

%clear g_data to save space in workspace
clear g_data;

% Extract relevant fields for glom signal
    % C_or the refined, merged ROI temporal data
    % A_or the refined, merged ROI spatial data

% Extract glomerular traces for each trial

%% LABEL TRIAL TYPES USING LOGGED TRIAL STRUCTURE

% HIGH-LOW        
turb_label(1:10) = 1;   
turb_label([11:15 21:25 31:35]) = 2;    
turb_label([16:20 26:30 36:40]) = 0;

% % LOW-HIGH
% turb_label(1:10) = 1;
% turb_label([11:15 21:25 31:35]) = 0;
% turb_label([16:20 26:30 36:40]) = 2;

trials_using = [1:40];

odor_on = [1000:1750]; %time period odor is present
%% LOAD ETH SEPERATED BY TRIAL

[fe, fde] = uigetfile('.dat', 'CHOOSE ETH DATA - .dat');
cd(fde);
fna = char(fe);

[ETT, ETT_FR] = Eth_x_trial(fe, turb_label, trials_using);


%% STFT PLUME 

[STFT, time_vector, freq_vector, Frames, TSMS] = STFT_eth (fe, turb_label, trials_using);
freq_max = 10 % Max freq for STFT in Hz;
% Find indices of STFT associated with freq range of interest
max_freq = max(find(freq_vector<=freq_max))

% Crop STFT for freq range of interest, here 1-10 Hz
STFT = STFT(1:max_freq, :,:);
%%  GET DFF GLOM CALCIUM
[tC tC_od tC_trials tC_trials_od]=DFFcal(C_or, trials ,turb_label);

%% DFF VS POWER FOR EACH GLOM

% For above STFT F(x)
    slide= 20; % 200ms slide intervals
    window=256; % 2560ms sliding window width
% Create 3D mat containing frame numbers for DFF for each trial for each
% glom to compare to frame numbers for STFT signal
for i=1:length(trials_using) % For each trial
    C_DFF_loop = C{1}; % Pull concatenated dff traces of 14 glom over all trials
    min_window = min(find(Frames(i,:))); % Find first frame associated with first window for given trial
    min_frame = Frames(i, min_window);
    max_frame = Frames(i, min_window + 349); % 350 is the total number of frames for each dff trace during trial
    max_window = max(find(Frames(i,:)<=max_frame));
    
    % Create vector of min frames associated with each window for the trial
    FR_2_WIN (i,:) = Frames(i, min_frame:max_frame);
    FR_range (i,:) = min_frame:max_frame;
    
    dff_trial(:, :, i) =C_DFF_loop(:, FR_2_WIN , i);
    %First frame of trial always ending of last trial do to looping
    %artifact so take first frame number that is not the minimum frame
    %number of the trial
%     for ii = 1:size(C_DFF_loop, 1) %for each glom, store the matrix of dff verse mean freq in cell of strucutre
%         x_coors= STFT(
%         % for iii = 1:length(win_temp) % For each window
%             % For each glom, create 2 vectors
%             % Create vector storing the DFF trace for each glom and store
%             % inside structure named F_v_freq
%             F_v_freq{ii} = (:, :, i) matrix of df and mean freq find(ETT_FR(FR_2_WIN(ii):(FR_2_WIN(ii+1)-1)
%             % Create vector storing the DFF trace for each glom and store
%             % inside structure named F_v_raw
%             %F_v_raw(iii, iii, i) = FR_2_WIN(ii):(FR_2_WIN(ii+1)-1)
%     end
    clear win_temp
end


%%

% %% FFT PLUME
% sr = 20; %estimated sample time in ms, since downsampled to match frame number, use sampling rate for CA- 20Hz
% Fs = 1/(sr/1000);
% L = length(odor_on);
% f = Fs*(0:(L/2))/L;
% 
% S_HIGH = std(ETT(turb_label==2,odor_on),[],2);
% figure, plot(S_HIGH)
% hold on
% S_LOW = std(ETT(turb_label==0,odor_on),[],2);
% plot(S_LOW)
% 
% mean(S_HIGH)
% mean(S_LOW)
% 
% title = ranksum(S_HIGH,S_LOW)
% 
% fft_L = abs(fft((ETT(turb_label==0,odor_on)-mean(ETT(turb_label==0,odor_on),2))'))';
% fft_L = fft_L(:,1:round(size(fft_L,2)/2));
% 
% fft_H = abs(fft((ETT(turb_label==2,odor_on)-mean(ETT(turb_label==2,odor_on),2))'))';
% fft_H = fft_H(:,1:round(size(fft_H,2)/2));
% 
% 
% figure('NumberTitle', 'off', 'Name', 'Avg FFT amp of "high" v. "low"')
% plot(f,mean(fft_L))
% hold on 
% plot(f,mean(fft_L)+std(fft_L)/sqrt(size(fft_L,1)),'b-')
% plot(f,mean(fft_L)-std(fft_L)/sqrt(size(fft_L,1)),'b-')
% plot(f,mean(fft_H))
% plot(f,mean(fft_H)-std(fft_H)/sqrt(size(fft_H,1)),'r-')
% plot(f,mean(fft_H)+std(fft_H)/sqrt(size(fft_H,1)),'r-')
% xlabel('Hz')
% ylabel('Power +-error bound')
% axis([0 10 0 max(mean(fft_H))+10])
% 
% for y = 1:size(fft_L,2);
%     rs(y) = ranksum(fft_L(:,y),fft_H(:,y));
% end


% Do actual multiple comparison test for rank-sums here
plot(f(rs<.001),mean(fft_H(:,rs<.001)),'r*','markersize',4)


fft_T = abs(fft((ETT(:,odor_on)-mean(ETT(:,odor_on),2))'))';
fft_T = fft_T(:,1:round(size(fft_T,2)/2));

sum_fft = sum(fft_T(:,(rs<.001)),2);
%sum_fft = sum((fft_T-mean(fft_L)).^2,2);

figure, plot(sum(sum_fft(turb_label==0),2),1+rand(1,size(fft_L,1)),'bo');
hold on
plot(sum(sum_fft(turb_label==2),2),2+rand(1,size(fft_H,1)),'ro');


% 
% % Use the following script to plot the STFT power analysis for a single trial
% eth_mxf = 3500;
% fs = 30; %Hz
% % Set number of fft points: SPECTRAL RESOLUTION
%     % Recommended to be power of 2
% nfft = 512
% % Set number of samples to slide window over, here 266.64 ms
%     % Recommended to be power of 2
% slide= 8; 
% window=64; %set length of window, here 2133.12 ms
% % To avoid spectral leaking use hamming window f(x)
%     % Because our signal is continuous, the sliding window can
%     % introduce sudden discontinuities which introduce harmonc content not
%     % present in the actual signal
% win = hamming(window, 'periodic'); 
% % Set 'window' equal to the total number of windows for looping
% period = 1/fs; % T is the sampling period (length of each window in seconds)
% signal_length= eth_mxf; % Number of timepoints sampled
%     %If signal_length is not even, drop last sample of X
%         % Signal length needs to be even valued for FFT
% time = (0:signal_length - 1)*signal_length; % Time vector of sampled timepoints in seconds
% winds = 1+fix((signal_length-window)/slide);  
% % Create stft matrix with zeros for looping
%     %Calculate number of rows and columns
% % Calculate number of freqs for row number
%     % + 1 why?
%     % Ceil() rounds towards next highest integer
% rown = ceil((1+nfft)/2);
% % Calculate number of windows for column number
%     % + 1
%     % Fix() rounds towards zero
% coln = 1+fix((signal_length-window)/slide); 
% stft = zeros(rown, coln);
% 
% 
% figure('Name',['STFT'],'NumberTitle','off')
% surf(time_vector, freq_vector(1:freq_max), STFT(1:freq_max, :, trial))
% % interp varies the color in each line segment and face by interpolating
% % the colormap index or true color value across the line or face
% shading interp
% % 'tight' : sets the axis limits to equal the range of the data so that the plot extends to the edges of the axes
% axis tight
% % 'box on' displays the box outline around the current axes
% box on
% % When using surf command and setting axis limits, view of figure can
% % shift, change the view to make sure the typical 2D plot axses are
% % maintained
% view(0, 90)
% % GCA returns the handle to the current axis
%     % It genereates one if there is no current axis
% set(gca)
% xlabel('Time, s')
% ylabel('Amplitude')
% title('Amplitude spectrogram of the signal')
% 
% handl = colorbar;
% set(handl) %, 'FontName', 'Times New Roman', 'FontSize', 14)
% ylabel(handl, 'Magnitude')
% for i = 1:length(STFT, 3)
%     % Take minimum frame number for STFT window
%        % Initialize the signal time segment index
%     temp_fr = frames(i, :)';
%     indx = 0;
%     for wind = 1:winds
%        
%         temp =  frames(indx+1:indx+window);
%         temp_frame = temp(temp ~= min(temp));
%         Frames(i,wind) = min(temp_frame);
%            % update the index
%         indx = indx + slide;
%     end
%     % For last STFT window of trial, take last frame and subtract one
%         % Last frame usually the first frame of next
%     if wind == winds
%         Frames(i, winds+1) = frames(i, end);
% end
% end






