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

%% HILBERT AMPLITUDE ENVELOPE OF DERIVE

% Bandpass butterworth filter 1-3Hz of each signal
clear ethfilt cafilt
X = diff(Eth(10, 1:349)');   
XX = diff(squeeze(C_trial(14, 1:349 , 10))'); 
% BANDPASS
Fs=20; %sampling rate
n=4; %filter order, looks at last 9 samples?
fpass=[1 3]; % Set threshold to pass for bandpass filter
[b,a]=butter(n,fpass/(Fs/2), 'bandpass'); %bandpass filter
ethfilt(:) = filtfilt(b,a,X(:));
cafilt(:)=filtfilt(b,a,XX(:));
%filter data for swrdata output
    %swrdata is now equal to the band pass filtered data, downsampled to
    %1.25kHz from the original selected channel
    
hilbert1 = hilbert(ethfilt);
hilbert2 = hilbert(cafilt);

% Plot the am[plitude envelopes of the two filtered signals
glom14_cor = figure
plot(linspace(1, length(hilbert1), length(hilbert1)), abs(hilbert1))
hold on
plot(linspace(1, length(hilbert1), length(hilbert1)), abs(hilbert2))
hold off

% Get cross corr of the two signals
[r, lags] = xcorr(abs(hilbert1), abs(hilbert2), 1250);

% Plot the cross corr
glom1_xcor = figure
plot(lags, r) % Plot cross corr of two signals
% plot(lags, r/max(r)) % Plot normalized cross corr of two signals

%% Save
% saveas(glom14_cor, 'glom14_cor.jpg')
% close all