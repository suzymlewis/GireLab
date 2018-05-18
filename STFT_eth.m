function [STFT, time_vector, freq_vector, Frames, TSMS] = STFT_eth (fn_dat, turb_label, trials_using)

%%%%%%%%%%%%%%% STFT Plume Analysis

    % INPUT
        % fn_dat = file name for ethanol .dat file from session
        % turb_label = vector containing turbulence level for each trial
          % within session
            % 0 = Low
            % 1 = Med
            % 2 = High
        % trials_using = vector containing trial numbers to be used from
        % session
        
    % OUTPUT
        % STFT: M x N x Tr matrix of short time fourier transform sliding across ethanol samples
         % from each trial
            % M = frequencies
            % N = windows across time
            % Tr = Trial number within session
        % Frames = M x N matrix where each entry corresponds to the minimum frame number from the DFF recording associated with the ethanol samples in each
         % STFT window included in the STFT matrix
            % M = window # within a trial
            % N = trial # within a session

%%%%%%%%%%%%%%% STFT Plume Analysis

%%  LOAD ETH FILE
            
fna = char(fn_dat);
fileID = fopen(fna);
data = fread(fileID,'int32','ieee-be');

%% EXTRACT ETH DATA FROM DAT FILE

% sr = 33.33; %estimated sample time in ms
% Wn = 8/50;
% [B,A] = butter(3,Wn,'low');

xcc = 0;
n_ten = find(data==-10);
TR = data(n_ten+1);
TS = data(n_ten+2);
FR = data(n_ten+3);
ETH = data(n_ten+4);
WHL = data(n_ten+5);
DST = data(n_ten+7);
SNF = data(n_ten+9);
un = unique(TR);
FR = FR-FR(1)+1;
clear fr eth_fr ts_fr dst snf unf

SNF = medfilt2(SNF,[3 1]);

%% DIVIDE ETH BY TRIAL
% Specify maximum number of ethanol time points to import
     % So all trials have equal samples and odor is not present for end of trials anyways
eth_mxf = 3500;

% % Low pass filter so no sig above 10Hz adding noise
    fc = 20; % Frequency cutoff for lowpass
    fs = 100; % Sampling rate for eth signal
    n = 4; % Fitler order
    Wn = fc/(fs/2); % Cutoff freq for filter is the frequency at which the
%     magnitude of the filter is 1/sqrt(2)
%          % Since digital filter Wn must lie between 0 and 1 where 1
%          % corresponds to Nyquist rate
%               % Here Nyquist rate = 25Hz
    [B,A] = butter(n,Wn,'low');

for i=1:length(trials_using)
    
    % Index using trial number of ethanol sample to divide ethanol by trial
    % Check trial usable by seeing if in turb_label
    % Check if length of trial is long enough
    if length(find(TR==(trials_using(i))))>450
    temp_trial = find(TR==trials_using(i));
    eth(i,:)=diff(ETH(temp_trial(1:eth_mxf+1))');
    
%   % Low pass filter out all sig above 10Hz
    eth(i,:) = filtfilt(B,A,eth(i,:));
    
    % Export frames associated with trial
    frames(i,:) = FR(temp_trial(2:eth_mxf+1))';
    TSms(i,:) = TS(temp_trial(2:eth_mxf+1))'; % Get timestamps associated with each trial
    
    % If trial in turb_label but ethanol sample is not long enough,
    % something went wrong so make trial all zeros for output and print
    % warning
    else
    eth(i,:) = zeros(size(trials_using(i), 3500));
    frames(i,:) = FR(temp_trial(2:eth_mxf+1))';
    TSms(i,:) = TS(temp_trial(2:eth_mxf+1))';
    warning (sprintf('Ethanol sample for trial %f is not of sufficient length, output zeros. Please check data.', i))
    end
    
    % clear temp_trial
end
% Create variable to indicate when odor is actually present
odor_on = [1000:1750]; %time period odor is present
eth_od = eth(:, odor_on);

%% Create STFT

%Set params for slide
% Set sampling rate of eth
fs = 100; %Hz
% Set number of fft points: SPECTRAL RESOLUTION
    % Recommended to be power of 2
nfft = 512
% Set number of samples to slide window over, here 50 ms so for each frame
    % Recommended to be power of 2
slide= 5; 
window=256; %set length of window, here 2560 ms
% To avoid spectral leaking use hamming window f(x)
    % Because our signal is continuous, the sliding window can
    % introduce sudden discontinuities which introduce harmonc content not
    % present in the actual signal
win = hamming(window, 'periodic'); 
% Set 'window' equal to the total number of windows for looping
period = 1/fs; % T is the sampling period (length of each window in seconds)
signal_length= eth_mxf; % Number of timepoints sampled
    %If signal_length is not even, drop last sample of X
        % Signal length needs to be even valued for FFT
time = (0:signal_length - 1)*signal_length; % Time vector of sampled timepoints in seconds
winds = 1+fix((signal_length-window)/slide);  
% Create stft matrix with zeros for looping
    %Calculate number of rows and columns
% Calculate number of freqs for row number
    % + 1 why?
    % Ceil() rounds towards next highest integer
rown = ceil((1+nfft)/2);
% Calculate number of windows for column number
    % + 1
    % Fix() rounds towards zero
coln = 1+fix((signal_length-window)/slide); 
stft = zeros(rown, coln);


% Rename signal to edit during SFTF code 
for i=1:size(eth,1)
    temp_sig = eth(i,:)';
    
    % Initialize the signal time segment index
    indx = 0;
    for wind = 1:winds
        
        % Section current window off from signal
        loop_win = temp_sig(indx+1:indx+window).*win;
        
        % FFT
        fft_loop = fft(loop_win, nfft);
        
        % Save looping FFT to STFT matrix
        % 'stft' representing that window's fft
        % Each row of STFT represents a window entry across all freqs(columns of stft)
        % Set all freq values of current looping window 'loop_win(1:rown)'
        % equal to column entries (freqs) of current window in stft
        stft(:, wind) = fft_loop(1:rown);
        
        % update the index
        indx = indx + slide;
    end
    STFT(:,:,i)=stft;
end

% Calculate the time and frequency vectors
time_vector= (window/2:slide:window/2+(coln-1)*slide)/fs;
freq_vector = (0:rown-1)*fs/nfft;

% Define the coherent amplification of the window
winamp = sum(hamming(window, 'periodic'))/window; % ;) lol winamp

% Take the amplitude of fft(x) and scale it, so not to be a
% function of the length of the window and its coherent amplification
for i=1:size(STFT, 3)
    temp_stft = STFT(:,:,i);
    temp_stft = abs(temp_stft/window/winamp);
    % correction of the DC & Nyquist component
    if rem(nfft, 2)                     % odd nfft excludes Nyquist point
        temp_stft(2:end, :) = temp_stft(2:end, :).*2;
    else                                % even nfft includes Nyquist point
        temp_stft(2:end-1, :) = temp_stft(2:end-1, :).*2;
    end
    STFT(:,:,i)=temp_stft;
end

%% CALCULATE MINIMUM FRAME NUMBER FOR EACH STFT WINDOW

for i = 1:size(STFT, 3)
    % Take minimum frame number for STFT window
    % Initialize the signal time segment index
    temp_fr = frames(i, :)';
    temp_TSms = TSms(i, :)';
    indx = 0;
    
    for wind = 1:winds
        loop_fr_wind = temp_fr(indx+1:indx+window);
        loop_ts_wind = temp_TSms(indx+1:indx+window);
        
         %First frame of trial always ending of last trial do to looping
         %artifact so take first frame number that is not the minimum frame
         %number of the trial
        if wind == 1
            fr_min = loop_fr_wind(loop_fr_wind~=min(loop_fr_wind));
            ts_min = loop_ts_wind(loop_ts_wind~=min(loop_ts_wind));
            fr_min= min(fr_min);
            ts_min= min(ts_min);
        else
            fr_min = min(loop_fr_wind);
            ts_min = min(loop_ts_wind);
        end
       
        Frames(i,wind) = fr_min;
        TSMS(i,wind) = ts_min;
        
        % update the index
        indx = indx + slide;
    end
end

for i = 1:size(STFT, 3)
    % Take minimum frame number for STFT window
    % Initialize the signal time segment index
    temp_fr = frames(i, :)';
    indx = 0;
    
    for wind = 1:winds
        loop_fr_wind = temp_fr(indx+1:indx+window);
        
         %First frame of trial always ending of last trial do to looping
         %artifact so take first frame number that is not the minimum frame
         %number of the trial
        if wind == 1
            fr_min = loop_fr_wind(loop_fr_wind~=min(loop_fr_wind));
            fr_min= min(fr_min);
        else
            fr_min = min(loop_fr_wind);
        end
       
        Frames(i,wind) = fr_min;
        
        % update the index
        indx = indx + slide;
    end
end
