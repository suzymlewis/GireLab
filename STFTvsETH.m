%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOOK AT FFT OF GLOMERULI SIGNAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOAD RAW DATA
[A_or C_or trials frame_trial ODOR_ON] = CalImdata();

%% LOAD COMPRESSED DATA
[raw raw_trials raw_frame_trial] = Compdata();

%% FIND ROIS FOR GLOMERULI IN TRIAL
[A ROI ROIweights] = ROImasking(A_or);

%% IMPORT ETHANOL
mxf = 350; % Set max frame for which to extract relevant ethanol signal
[Eth] = ethdata(mxf);

%% STFT OF ETH
% Resource code: https://www.mathworks.com/matlabcentral/fileexchange/45197-short-time-fourier-transformation--stft--with-matlab-implementation
    % Slide Fast Fourier Trasnform over signal for a finer scale resolution
    %spectrogram
% Rename signal to edit during SFTF code 

% LOOK AT ETH SIGNAL TO SEE HOW PREDOMINANT FREQ EVOLVE ACROSS TRIAL
X = Eth(10, 150:300)';   
fs = 20; % fs is the sample rate of the data
window=4; %set length of window, here 100 ms
% To avoid spectral leaking use hamming window f(x)
    % Because our signal is continuous, the sliding window can
    % introduce sudden discontinuities which introduce harmonc content not
    % present in the actual signal
win = hamming(window, 'periodic'); 
% Set number of fft points
    % Recommended to be power of 2
nfft = 512
% Set number of samples to slide window over, here 25 ms
    % Recommended to be power of 2
slide= 2; 
period = 1/fs; % T is the sampling period (length of each window in seconds)
signal_length= length(X); % Number of timepoints sampled
    %If signal_length is not even, drop last sample of X
        % Signal length needs to be even valued for FFT
if rem(length(X), 2)~=0
    % Drop last sample of X if it is not an even number of samples
    X(end) = [];
    % Update signal length to reflect the new even-valued signal length
    signal_length = signal_length - 1;
end
time = (0:signal_length - 1)*signal_length; % Time vector of sampled timepoints in seconds

% Take fast fourier transform of signal
Fourier = fft(X);

% Compute Single-Sided Amplitude Spectrum of 'fourier(t)'
% To get magnitude of frequencies take abs()
    % Normalize by dividing by length of signal
full_spec = abs(Fourier/signal_length); 
oneside_spec = full_spec(1:signal_length/2+1);
oneside_spec(2:end-1) = 2*oneside_spec(2:end-1);

% Define frequency range
freqrange = fs*(0:(signal_length/2))/signal_length;

% Plot Amplitude Spectrum of X(t)
figure, plot(freqrange,oneside_spec) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|oneside_spec(f)|')

% Create STFT
%Set params for slide
% Set 'window' equal to the total number of windows for looping
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


% Initialize the signal time segment index
indx = 0;
for wind = 1:winds
    
    % Section current window off from signal
    loop_win = X(indx+1:indx+window).*win;
    
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

% Calculate the time and frequency vectors
time_vector = (window/2:slide:window/2+(coln-1)*slide)/fs;
freq_vector = (0:rown-1)*fs/nfft;

% Define the coherent amplification of the window
winamp = sum(hamming(window, 'periodic'))/window; % ;) lol winamp

% Take the amplitude of fft(x) and scale it, so not to be a
% function of the length of the window and its coherent amplification
stft = abs(stft)/window/winamp;

% correction of the DC & Nyquist component
if rem(nfft, 2)                     % odd nfft excludes Nyquist point
    stft(2:end, :) = stft(2:end, :).*2;
else                                % even nfft includes Nyquist point
    stft(2:end-1, :) = stft(2:end-1, :).*2;
end

%Find indices of freq_vector to set max cutoff at 10 Hz
freq_max = min(find(freq_vector>=10));

% plot the spectrogram
figure('Name',['STFT'],'NumberTitle','off')
surf(time_vector, freq_vector(1:freq_max), stft(1:freq_max, :))
% interp varies the color in each line segment and face by interpolating
% the colormap index or true color value across the line or face
shading interp
% 'tight' : sets the axis limits to equal the range of the data so that the plot extends to the edges of the axes
axis tight
% 'box on' displays the box outline around the current axes
box on
% When using surf command and setting axis limits, view of figure can
% shift, change the view to make sure the typical 2D plot axses are
% maintained
view(0, 90)
% GCA returns the handle to the current axis
    % It genereates one if there is no current axis
set(gca)
xlabel('Time, s')
ylabel('Frequency')
title('Amplitude spectrogram of the signal')

handl = colorbar;
set(handl) %, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(handl, 'Magnitude')

%% STFT OF CA
% Resource code: https://www.mathworks.com/matlabcentral/fileexchange/45197-short-time-fourier-transformation--stft--with-matlab-implementation
    % Slide Fast Fourier Trasnform over signal for a finer scale resolution
    %spectrogram
% Rename signal to edit during SFTF code 
C_trial = reshape(C_or,size(C_or,1),size(C_or,2)/length(trials),length(trials));
% LOOK AT ETH SIGNAL TO SEE HOW PREDOMINANT FREQ EVOLVE ACROSS TRIAL
XX = squeeze(C_trial(2, 150:300 , 10))';   
fs = 20; % fs is the sample rate of the data
window=4; %set length of window, here 100 ms
% To avoid spectral leaking use hamming window f(x)
    % Because our signal is continuous, the sliding window can
    % introduce sudden discontinuities which introduce harmonc content not
    % present in the actual signal
win = hamming(window, 'periodic'); 
% Set number of fft points
    % Recommended to be power of 2
nfft = 512
% Set number of samples to slide window over, here 25 ms
    % Recommended to be power of 2
slide= 2; 
period = 1/fs; % T is the sampling period (length of each window in seconds)
signal_length= length(XX); % Number of timepoints sampled
    %If signal_length is not even, drop last sample of X
        % Signal length needs to be even valued for FFT
if rem(length(XX), 2)~=0
    % Drop last sample of X if it is not an even number of samples
    XX(end) = [];
    % Update signal length to reflect the new even-valued signal length
    signal_length = signal_length - 1;
end
time = (0:signal_length - 1)*signal_length; % Time vector of sampled timepoints in seconds

% Take fast fourier transform of signal
Fourier = fft(XX);

% Compute Single-Sided Amplitude Spectrum of 'fourier(t)'
% To get magnitude of frequencies take abs()
    % Normalize by dividing by length of signal
full_spec = abs(Fourier/signal_length); 
oneside_spec = full_spec(1:signal_length/2+1);
oneside_spec(2:end-1) = 2*oneside_spec(2:end-1);

% Define frequency range
freqrange = fs*(0:(signal_length/2))/signal_length;

% Plot Amplitude Spectrum of X(t)
figure, plot(freqrange,oneside_spec) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|oneside_spec(f)|')

% Create STFT
%Set params for slide
% Set 'window' equal to the total number of windows for looping
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


% Initialize the signal time segment index
indx = 0;
for wind = 1:winds
    
    % Section current window off from signal
    loop_win = XX(indx+1:indx+window).*win;
    
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

% Calculate the time and frequency vectors
time_vector = (window/2:slide:window/2+(coln-1)*slide)/fs;
freq_vector = (0:rown-1)*fs/nfft;

% Define the coherent amplification of the window
winamp = sum(hamming(window, 'periodic'))/window; % ;) lol winamp

% Take the amplitude of fft(x) and scale it, so not to be a
% function of the length of the window and its coherent amplification
stft = abs(stft)/window/winamp;

% correction of the DC & Nyquist component
if rem(nfft, 2)                     % odd nfft excludes Nyquist point
    stft(2:end, :) = stft(2:end, :).*2;
else                                % even nfft includes Nyquist point
    stft(2:end-1, :) = stft(2:end-1, :).*2;
end

%Find indices of freq_vector to set max cutoff at 10 Hz
freq_max = min(find(freq_vector>=10));

% plot the spectrogram
figure('Name',['STFT'],'NumberTitle','off')
surf(time_vector, freq_vector(1:freq_max), stft(1:freq_max, :))
% interp varies the color in each line segment and face by interpolating
% the colormap index or true color value across the line or face
shading interp
% 'tight' : sets the axis limits to equal the range of the data so that the plot extends to the edges of the axes
axis tight
% 'box on' displays the box outline around the current axes
box on
% When using surf command and setting axis limits, view of figure can
% shift, change the view to make sure the typical 2D plot axses are
% maintained
view(0, 90)
% GCA returns the handle to the current axis
    % It genereates one if there is no current axis
set(gca)
xlabel('Time, s')
ylabel('Frequency')
title('Amplitude spectrogram of the signal')

handl = colorbar;
set(handl) %, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(handl, 'Magnitude')

