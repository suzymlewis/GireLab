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
Y = g_data.Y; % Grab indices for ROI masks for each glomeruli
    % Y is M x N x Frames concatenated across all trials


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

% Call Eth_x_tral()
    %diff(eth) stored in 'ETT' for each frame assignment stored in 'ETT_FR'
[ETT, ETT_FR] = Eth_x_trial(fe, turb_label, trials_using);

%% AVERAGE ACROSS ETH SAMPLES WITHIN EACH FRAME
%%%%%%%%%%%%%%%% ADD 1 to all eth frame numbers to match recent changes in
%%%%%%%%%%%%%%%% CA

down_odor_on =350; % For frame matched data, need max relevant frame number during odor presentation

for i=1:length(trials_using) % LOOP TRIALS
    FRA = unique(ETT_FR(i,:)); % Create vector with all relevant frames in trial
    for ii=1:length(FRA) % LOOP UNIQUE FRAMES
        eth_data(i,ii) = mean(ETT(ETT_FR(i,:)==FRA(ii)));
    end
    clear FRA
end

eth_data = eth_data(:,1:down_odor_on); % Take eth signal only for frames relevant to odor presentation

%% THEORY REGARDING GLOMERULAR DELAY AND SNIFFING INDUCED VARIANCE
% Previous inhalation cycle drives depolarization more than 50–100 ms after the onset of inhalation
    % Concentration-dependent manner (Carey et al., 2009)
    % For modest odor concentrations:
        % MCs lag behind the odor stimulus by approximately 1/2 sniff cycle. 
    % Citation: Fukunaga, I., Berning, M., Kollo, M., Schmaltz, A., & Schaefer, A. T. (2012). Two distinct channels of olfactory bulb output. Neuron, 75(2), 320-329.
% Full sniffing Range:
    % 4 –12 Hz when sniffing
    % Citation: Oka, Y., Takai, Y., & Touhara, K. (2009). Nasal airflow rate affects the sensitivity and pattern of glomerular odorant responses in the mouse olfactory bulb. Journal of Neuroscience, 29(39), 12070-12078.
% Exploratory Sniffing Behavior:
    % Distinctive odor sampling strategy
    % High-frequency (6–10 Hz) sniff bouts lasting up to several seconds
    % Citation: Wachowiak, M. (2011). All in a sniff: olfaction as a model for active sensing. Neuron, 71(6), 962-973.
% Respiratory frequency: 
    % Increases from 2–5 Hz under normal breathing conditions
    % Citation: Oka, Y., Takai, Y., & Touhara, K. (2009). Nasal airflow rate affects the sensitivity and pattern of glomerular odorant responses in the mouse olfactory bulb. Journal of Neuroscience, 29(39), 12070-12078.
% So variance from dynamic in odor to processing in glomeruli could range ~80ms:500ms
    % If jitter consistent across a single glomeruli but differnt across glomeruli within each trial:
        % Suggest processing style of the glomeruli
        % Problem from analysis perspective is that the "trial" here is a
        % single type of change in eth dynamic, not just a change
            % So compare one feature of eth deriv across glom for a start?
    % If jitter consistent across all glomeruli in a trial:
        % Suggest dependent on sniff cycle delay
    % If different for each glom for each trial and different across glom
    % within a trial
        % For now suggest "noise"

%%  GET DFF GLOM CALCIUM

[tC tC_od tC_trials tC_trials_od]=DFFcal(C_or, trials ,turb_label);


%% LOOK AT CORRELATION BETWEEN EACH GLOM TO FIND AVG DELAY FOR EACH GLOM WITH MAX X-CORR
% Correlation to plume analysis

thr = 1; % Set threshold for std of derivation of normalized, mean centered trial data from glomeruli

tm_rn = 20; %range in frames before and after

% Set max frame from which to index
    % Use final frame analyzing in cross corr
PL_MEAN = eth_data;
PL_MEAN = PL_MEAN-repmat(mean(PL_MEAN(:,1:100),2),1,size(PL_MEAN,2));
PL_MEAN = eth_data';
%%%%%%%%%%%%%%%%%%%% FIX THIS SO VARIABLE/ VAR ASS NOT VAR/INT ASS   
mx = size(PL_MEAN,1); 

%%%%%%%%%%%%%%%%%%%%%%%%% REAL THING NEED TO DO IS GET A PL_MEAN; I HAVE
%%%%%%%%%%%%%%%%%%%%%%%%% ETT BUT THAT IS SAMPLED AT 10X RATE SO HAVE 3500
%%%%%%%%%%%%%%%%%%%%%%%%% INSTEAD OF 350 FRAMES
%%%%%%%%%%%%%%%%%%%%%% Later make this interp based on timestamp, for now
%%%%%%%%%%%%%%%%%%%%%% just average that is what david did


%%%%%%%%%%%%%%%%%%%%%%%%% DAVID REMOVES BASELINE OF ETHANOL DATA
%%%%%%%%%%%%%%%%% CONSIDER DELTA F/F ethanol as well??????????????????????


% PL_MEAN is a 350 x 40 strucutre
    % It is the raw eth over the 350 frames point for a trial (columns)
    % for each trial (rows)

xc = 0;
clear AVG gl_id AVG_gl crc tr_id

TR_MEAN = tC_trials;
TR_MEAN1= TR_MEAN;
% TR_MEAN is just C_trial(:,:,:) at this point
    % Note =  C_trial = calcium 14x350x40
   
for y = 1:size(TR_MEAN1,1); % LOOP GLOMERULI
    crr = squeeze(TR_MEAN1(y,:,:)); % % JUST RESHAPED CALCIUM
        % 350 x 40 for looping glom
        % Here please note TR_MEAN1==TR_MEAN 
        
        % Create CRR which is glomerular matrix for the current looping glomeruli    
    for y1 = 1:size(crr,2); % LOOP TRIAL WITHIN GLOMERULI
        clear f2 crr1 zs df1 f1
        crr1 = squeeze(crr(:,y1));
            % 350 x 1 signal for looping trial within looping glomeruli
        
        %%%%%%% DONT NEED TO CALC ZS BECAUSE ALREADY USING DELTA F
        %%%%%%% Extract one trial delta f for looping glom
        zs= crr1-mean(crr1(:,:));% Subtract the mean
        
%         % CALCULATE USABLE SIGNAL BY TAKING DERIV OF NORMALIZED, MEAN-CENTERED DATA
%             % NOT USING DELTA F/F
%         zs = crr1-mean(crr1(:,:)); % Subtract the mean
%         zs = zs./std(zs(:)); % Divide by the std to normalize
%         zs = diff(zs); % Take deriv of signal
        
        % FIND SIGNIFICANT CHANGES
        f1 = find(zs>thr); % Find where exceed 1 std change from last sample
            % m x 1 where m=number of sig events
        df1 = diff(f1); % Take derivative of this significant signal
        f2 = f1(df1>5); % Make sure significant points are distinct responses to temporal seperated dynamics
            % Since sampled at 30Hz
                % Events within 5 points are within a 1/6 second window
                    % So this changes 6Hz or faster
                    % Good because changes in calcium dynamics are on the range of 2Hz
            % Only use points where the change between significant points is greater than 5
         
        % CROSS CORR OF DIFF RAW ETH AND DIFF RAW CA
        %%%%%%%%%%%%%%%%%%%% As his is right now, one is across trials and
        %%%%%%%%%%%%%%%%%%%% the other across timepoints within a trial
        crc(y,:,y1) = xcorr(diff(PL_MEAN(:,y1)),diff(squeeze(TR_MEAN(y,:,y1))),'coeff');
            % diff(PL_MEAN(:, y1)) = eth for looping trial 'y1'
                % 349 x 1
                % diff(squeeze(TR_MEAN(y,:,y1)))
                % deriv ca signal for trial
                % 1 x 349
            % crc(y,:,y1) = cross correlation of eth and ca signals for trial
                % 1 x 697
                
        for u = 1:length(f2); % LOOP SIG EVENTS
            % f2 = entries of f1 with temporal seperate signifcant events
            % f1 = indices of za for sig events
                
            
            if f2(u) > tm_rn & f2(u) < mx-tm_rn; % Make sure sig event not window artifact of cross correlation
                % Make sure sig event not in first or last tm_rn (# frames suspected in range for artifact)
                % Drops any sig events in first and last 20 frames because likely a window artifact of xcorr
                
                xc = xc+1; % Index for significant event
            
                %%%%%%%%%%%%%%%%%%%%%%%%  SOME SORT OF DIMENSION MISMATCH BETWEEN THIS CODE AND DAVIDS   
%             Error using coder.internal.assert (line 17)
% When B is a vector, A must be a vector.
% 
% Error in xcorr (line 108)
%     coder.internal.assert(iscolumn(x), ...
 
                % STORE ETH DERIV CORRESPONDING TO SIG EVENT INTO 'AVG' VECTOR
            AVG(xc,:) = squeeze(diff(PL_MEAN(f2(u)-tm_rn:f2(u)+tm_rn,y1)))'; % Deriv eth of current sig event
                % Take 20 frames pre and post sign event
                % f2 is index of sig events
                % PL_MEAN is raw eth data
                % 'y1' just indexing current looping trial
                % AVG = m x 40 structure
                    % m = number of sig events not part of window artifact
                    
            % STORE DERIV OF MEAN CENTERED, NORMALIZED CA
            AVG_gl(xc,:) = zs(f2(u)-tm_rn:f2(u)+tm_rn);
                % Take 20 frames pre and post sign event    
                % zs = ca deriv, normalized and mean centered
                % m x 41
            % CREATE VECTOR INDICATING GLOMERULI AND TRIAL TAGS
            gl_id(xc) = y; % Tag glomeruli sig event comes from
            tr_id(xc) = y1; % Tag trial sig event comes from
                % Both vectors of length == number sig events
            end
            
        end
    
    end
    
end

%% PLOT GLOMERULI ROI MASKS FOR SESSION
clear ROI RROI
fnd = [1:size(A_or, 2)];
for y = 1:size(A_or, 2); % Fnd is the number of glomeruli found to be responsive to the plume
    ROI = A_or(:,y); %%%%%%%%%%%%% In future can only loop glomeruli with sig response indx(y)); % Index into A_or the spatial matrix
   % ROI = reshape(ROI,size(Y,1),size(Y,2));
   if y==1;
       RROI = ROI;
       RROI(ROI~=0) = fnd(y);
   else
    RROI(ROI~=0) = fnd(y);
   end
   
end
    RROI = reshape(RROI,size(Y,1),size(Y,2));
    
figure('name','ROI GLOMERULAR MASKS','NumberTitle','off')
imagesc(RROI)%[min(size(A_or, 2))-1 max(size(A_or, 2))+1])
%% PLOT MEAN RESPONSE OF GLOMERULI

% Avg DFF of each glomeruli across trials
resp_mean = mean(tC_trials, 3); 
    % So resp_mean is M x N
        % M = # glomeruli
        % N = # Frames within each trial

% Plot each glomeruli average response across trial 
figure, imagesc(resp_mean,[-1 1]);    

% Find max DFF value for each glomeruli across mean trial response
resp_str = max(resp_mean,[],2);
    % M x 1 column vector containing the max DFF response of each glomeruli(M)

% Sort the glomeruli by max of mean response strength   
[rs rs1] = sort(resp_str);
    % rs is the glomeruli sorted by max frame of mean activity
        % From least to most active
    % rs1 is the indices of these glomeruli in the presorted array

% Sort avg DFF response of each glomeruli across trials
    % From least to most active
resp_mean1 = resp_mean(rs1,:);
    % M x N:
        % M is glomeruli
        % N is number of frames within a trial

% Look at average glomerular response
    % NOTE: ONSET, OFFSET, and DYNAMIC RESPONSES
% Top is the average response across trials for each glomeruli sort
    % From least to most active
figure, subplot(2,1,1), 
imagesc(resp_mean1,[-0.5 1])

% Find standardized x-corr response across trials for each glomeruli
corr_resp = mean(crc(rs1,:,:),3)./std(crc(rs1,:,:),[],3);
    % Index using rs1 to sort by glomerular responsivity (less to most)
    % corr_resp = 14 x 697
corr_resp1 = corr_resp;

% Remove low level correlations
    % Remove any low value standardized x-corr
        % Set values less than 1 stdev = 0
corr_resp(abs(corr_resp)<1) = 0;

% Take the maximum cross correlation response of each glomeruli
cr1 = max(corr_resp,[],2);
    % cr1 = 14 x 1
    
% Get size of 1/2 x-corr vector and round to nearest whole int
sz = round(size(corr_resp,2)/2);

% For the Bottom half of Figure
% Plot the mean, large x-corr response of each glomeruli
subplot(2,1,2),
imagesc(corr_resp(:, sz-40:sz+40)); %imagesc(corr_resp(:,:))

%% GET MEAN DELAY FOR EACH GLOMERULI

% Find mean delay response for each glomeruli
    % Find max of x-corr mean for each glomeruli across trials
for i = 1:size(corr_resp, 1)
    if sum(corr_resp(i, sz-40:sz+40))==0 % If not large x-corr value then delay = 0
        delay(i)=0;
    else
        delay(i)=max(corr_resp(i, sz-40:sz+40)); 
    end
end
delay=delay';

%% Plot data from correlation and response analysis

corr_resp = mean(crc(rs1,:,:),3)./std(crc(rs1,:,:),[],3);
corr_resp1 = corr_resp;

corr_resp(abs(corr_resp)<1) = 0;
cr1 = max(corr_resp,[],2);

sz = round(size(corr_resp,2)./2);
subplot(2,1,2),
imagesc(corr_resp(:,sz-40:sz+40));

figure, subplot(2,1,2), imagesc(corr_resp1(cr1>0,sz-40:sz+40))
subplot(2,1,1), imagesc(resp_mean1(cr1>0,:))

LG = ([sz-40:sz+40]-mean([sz-40:sz+40]))./20;

figure, plot(LG,corr_resp1(cr1>0,sz-40:sz+40)','linewidth',2)
hold on
plot([0 0],[-2 2],'k--')

%% STFT PLUME 

[STFT, time_vector, freq_vector, Frames, TSMS] = STFT_eth (fe, turb_label, trials_using);
freq_max = 10 % Max freq for STFT in Hz;
% Find indices of STFT associated with freq range of interest
max_freq = max(find(freq_vector<=freq_max))

% Crop STFT for freq range of interest, here 1-10 Hz
STFT = STFT(1:max_freq, :,:);


%% AVERAGE ACROSS TIMEPOINTS WITHIN FRAME TO GET DERV FOR EACH FRAME


%% MATCH TO CORRESPONDING CA FRAME
%% PLOT
% % Add raw Florescence trace
% C_trial = reshape(C_or,size(C_or,1),size(C_or,2)/length(trials),length(trials));
% 
% % Need to get timepoint in trial of Calcium frame to match to ethsampe
% % Can not find exact time point for each frame of referemce
% % Will have to just reference through via first eth sample for each frame  
% % Then can finish interp1 function
% 
% % Find minimum timestamp assciated with each frame to find timestamp of
% % fram
% 
% 
% tempy = FR(find(TR==1));
% drop = find(tempy==1);
% tempy(drop)=[];
% starty = tempy(1);
% endy = starty + 349;
% rangey = starty:endy;
% for i= 1:350
%    indexy = min(find(tempy == rangey(i)));
%    stamp = TS(indexy);    
%    Frame_times(i) = stamp;
% end
% 
% % Drop any repeating time frames from timestamp interpolating
% 
% repeat = diff(Frame_times);
% dropy = find(repeat == 0);
% Frame_times(dropy) = [];
% 
% 
% 
% timestamps_forwindows= unique(TSMS(1,1:415)); 
% 
% % For point of Frame_times dropped
% % average frames for same timepoint
%  calcy = C_trial(:,:,1); % Create copy of raw calcium trace to alter for interp
%  % Take mean of neighboring frames for same timepoint and replace frame
%  % before repeat with mean
%  calcy(:,dropy-1)= mean(C_trial(:,dropy-1:dropy,1), 2);
%  %Drop repeating frame
%  calcy(:, dropy)=[];
% % NOW USE INTERP1 VALUES FOR calcy to GET FRAMES FROM COLUMN 3 THEN DONE
% % WITH DAD!!
% %Reshape dad for dropped frames
% 
% 
% 
% % Because frames do not line up perfectly with ethanol samples and therefore with ethanol STFT, interp1 CA to get STFT timepoints.
% for i =1:14
%     C_interp(i,:) = interp1(unique(Frame_times), calcy(i,:), unique(TSMS(1,1:415)), 'spline')
% end
% %NEED TO WORK OUT CALCIUM LOOP TO ADD REDUNDANT FRAMES BACK IN
% % dad= dad((1:length(calcy)*14),:,:)
% % Place in dad matrix as third column
% for glom=1:14
%     temp=undropca(glom,:);
%     indx=1;
%     for i= 1:14:(size(dad,1))
%         i=i+(glom-1);
%         dad(i, 3) = temp(indx);
%         indx=indx+1;
%     end
% end
% % undropca = C_interp(indx,
% 
% indx = 1;
% tempy = unique(TSMS(1,1:415));
% for glom = 1:14;
%     for i = 1:length(unique(TSMS(1,1:415)))
%         pos = find(TSMS(1,1:415)==tempy(i));
%         if length(pos)>1
%             for ii = 1:length(pos)
%                 undropca(glom, indx + (ii-1)) = C_interp(glom, i);
%                 indx=indx+1
%             end
%         else
%             undropca(glom, indx)=C_interp(glom ,i);
%             indx=indx+1;
%         end
%         clear pos
%     end
%     indx = 1;
% end
% % Add calcy into Dad


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
