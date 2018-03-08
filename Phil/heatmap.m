clear all

%% Load 'messages.event' file
% Load messages event file
% Paired timestamps located in messages.data
% Start times located in messages.textdata

%Load message file
[messfilename, messagelocation] = uigetfile('.events', 'PICK MESSAGE.EVENTS FILE');
% Change directory to file location to load message.events file
cd(messagelocation)
messages = importdata(messfilename);

% Get N x 2 mat for paired timestamps
    % N equals the number of timepoints for recording session
timestamp_mat = messages.data;

%Extract info on timepts for software(LabView) and processor(OE GUI)
LabView_info = char(messages.textdata(1));
LabView_info = LabView_info(24:end);
OEphy_info = char(messages.textdata(2));
OEphy_info = OEphy_info(24:end);

% Event file example uses a 2-bit word for the clock sync
% Experiment is responsible for sending regular sync signals

%% Load 'all.channels.event' file and extract 6 variables
% Event files (.events) -  Non-spike events are saved into the file "all_channels.events" as timestamps and ids.
% EVENT OUTPUT NOTES
%     DATA = EVDATA
        % 1-Dimensional array of event channels (integers)
        % 'int64' part of event file
%     TIMESTAMPS = EVTIMSTAMPS
        % 1-Dimensional array in seconds
        % 'uint64' part of event file
%     INFO = EVINFO
        % Structure with header and other information
            % info.sampleNum
                % 'int16' part of events file
            % info.nodeId 
                % 'uint8' part of events file
                % Node 102 = Signal channels with .cont files
                % Node 103 = 
                % Node 136 = 
            % info.eventType
                % 'uint8' part of events file
            % info.eventId
                % Code associated with this event, usually 1 or 0
                % 'uint8' part of event file

% Load event file
[eventfilename, eventlocation] = uigetfile('.events', 'PICK ALL_CHANNELS.EVENTS FILE')
% Change directory to file location to load message.events file
cd(eventlocation)

% Extract information with OpenEphys f(x)
[evdata, evtimestamps, evinfo] = load_open_ephys_data(eventfilename);

% info.header %run this code if need to checkdate/ samplerate/ etc.

% Index into 4 fields of 'info' structure and set each of these 4 fields equal to a variable
header = evinfo.header; %Give information like date and sampling rate
sample_number = evinfo.sampleNum;
node_id = evinfo.nodeId;
event_type = evinfo.eventType;
event_id = evinfo.eventId;
recording_number = evinfo.recordingNumber; %Always == 0


% Concatenate event information into a matrix [6 x length(timestamps)]
info_mat = vertcat(evdata', evtimestamps', sample_number', node_id', event_type', event_id', recording_number);

%% Create info_mat with only timepoints for events of Type 5

% Find indicies in row 5 (the event_type row) where event_type is 5
splice_3 = find(info_mat(5,:) == 3);
splice_5 = find(info_mat(5,:) == 5);

%Keep only the columns of info_mat that meet above properties
info_mat_type_3 = info_mat(:, splice_3);
info_mat_type_5 = info_mat(:, splice_5);

%% Read in behavioral file

[behav, behavlocation] = uigetfile('.xlsx', 'PICK XLSX BEHAVIORAL FILE');
% Change directory to file location to load message.events file
cd(behavlocation)
% Load behavioral data
behavior = xlsread(behav);

% Extract N x 3 matrix (pos_mat) with x position, y position, and timepoints
    % N equals the number of timepoints for recording session
pos_mat = behavior(:, 4:6);


%Smooth positional data by sliding window, decimate, then upsample




% Loop through pos_mat using the first and last behavioral timepoints from the paired
% timepoints matrix (unless the last timepoint is preceeded by the OE
% timepoints) to create a matrix with OE timestamps paired with behavioral
% data
paired_behavior = zeros(size(timestamp_mat, 1), 9);
paired_behavior(:, 1:2) = timestamp_mat;
x_pos = behavior(:, 4);
y_pos = behavior(:, 5);
round_x_pos = round(x_pos);
round_y_pos = round(y_pos);
round_pos_mat = horzcat(round_x_pos, round_y_pos);
pos = zeros(max(round_x_pos), max(round_y_pos));
for xp=1:size(round_pos_mat, 1)
    pos(round_pos_mat(xp,1),round_pos_mat(xp,2))= pos(round_pos_mat(xp,1),round_pos_mat(xp,2)) + 1;
end

%CREATE INDICES TO DIVIDE MATRIX INTO EQUAL SIZED BINS
pos_bin = pos(13:572, 12:431);
x_wins = [4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4];
y_wins = [3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3];
pos_binning = mat2cell(pos_bin, x_wins, y_wins);

% SUM NUMBER OF NONZERO ENTRIES IN BIN
z = cellfun(@(pos_bin) sum(pos_bin(:)), pos_binning);

% MAZE IS 5 FT x 5 FT
    % CONVERT PIXELS TO FT
x_pos_focus = x_pos-min(x_pos);
y_pos_focus = y_pos-min(y_pos);
x_convert = sqrt(25/2) / (max(x_pos_focus));
y_convert = sqrt(25/2) / (max(y_pos_focus));
convert = mean([x_convert y_convert]); %This conversion in ft per pixel
x_pos_ft = x_pos_focus.*convert;
y_pos_ft = y_pos_focus.*convert;

% PLOT TRAJECTORY
figure('Name',['Raw Trajectory'],'NumberTitle','off')
plot(x_pos_ft, y_pos_ft);
xlabel('Ft');
ylabel('Ft');
axis([0 max(x_pos_ft) 0 max(y_pos_ft)]);

% PLOT POSITIONAL HEATMAP
figure, imagesc(z) 
colormap hot

%% LOAD FILE .cont TO LOOK AT THETA, FAST GAMMA, SLOW GAMMA POWER VIA SLIDING SPECTROGRAM

[filename, contlocation] = uigetfile('.continuous', 'PICK .CONT CH FILE');
% Change directory to file location to load message.events file
cd(contlocation)

[chdata, chtimestamps, chinfo] = loadAndCorrectPhase(filename, 1); %loading data from open ephys
%.continuousXXX_CHY
    % XXX=processor and Y=ch#
    
% DOWNSAMPLE BY FACTOR OF 24 SO SAMPLING RATE 30kHz -> 1250Hz
downdata=decimate(chdata,24); %The second argument is the downsample factor, down from 30khz to 1250hz

% EXTRACT FIELDS FROM CONT FILE
header = chinfo.header; % General information from recording session
ts = chinfo.ts; % Unlike timestamps which increase by 1 each entry. ts sample @1024, so increases accordingly, but has same beginngin and end timepoints
nsamples = chinfo.nsamples; % Always 1024
recNum = chinfo.recNum; % Always == 0

% BANDPASS THETA THROUGH FAST GAMMA (4-100Hz)
Fs=1250; %sampling rate
n=4; %filter order, looks at last 9 samples?
fpass=[4 100]; % Set threshold to pass for bandpass filter
[b,a]=butter(n,fpass/(Fs/2), 'bandpass'); %bandpass filter
fildata(:) = filtfilt(b,a,downdata(:)); %filter data for swrdata output
    %swrdata is now equal to the band pass filtered data, downsampled to
    %1.25kHz from the original selected channel
%% SLIDING SPECTROGRAM
window=125; %set length of window, here 125 ms
overlap= 100; %set length of window overlap for sliding every 25ms start new window
spec = []; % create matrix to store spectrograms during loop
fs = 1250 % fs is the sample rate of the data
spec = []; % 
specfreq= [] % Cyclical frequencies
specTM = []; % time instnts at which the spectrogram is computed

% CREATE SLIDING WINDOW OF SPECTROGRAMS
[spec, specfreq, specTM] = spectrogram(fildata(1:length(fildata)), hamming(window), overlap, [], fs,'onesided');


% CREATE AVERAGED SPECTROGRAM
    % Average over all windows the overlay a single timepoint to create new
    % Averaged sliding spectrogram

% FILTER SPECTROGRAM TO ISOLATE THETA THROUGH SLOW GAMMA
% Find 1st entry of specfreq >= 4 Hz
spread(1) = min(find(specfreq>=4)); 
% Find 1st entry of specfreq >= 100 Hz
spread(2) = min(find(specfreq>=100)); 
specspread = spec(spread(1):spread(2), :); %Create a cropped matrix with only 4-100Hz
specmat = zeros(size(specspread, 1), (size(specspread, 2) - window + 1)); %Create square matrix to hold averaged spectrogram

for frequency = 1:size(specspread, 1); % For each frequency window
    for timepoint = 1:(size(specspread, 2) - window + 1); %For each timepoint sampled within frequency
            % across first dim since took transpose
            % Minus window size because can't average endpoint timepoints
            % since window cant slide past edge of timepoint ending
        dist = (specspread(frequency, timepoint:(timepoint + (window - 1))));
        specmat(frequency, timepoint)= mean(dist);
    end
end


% REMOVE NOISE ACROSS CHANNELS
% spec_mag = abs(spec);
% specmat_mag = abs(specmat);
% spec_mag(spec_mag >=5000)=5000;
% specmat_mag(specmat_mag >=5000)=5000;

% FILTER SPECTROGRAM TO ISOLATE THETA, FAST GAMMA, SLOW GAMMA
% Find 1st entry of specfreq >= 4 Hz
freq_cutoff(1) = min(find(specfreq>=4));
% Find 1st entry of specfreq >= 14 Hz
freq_cutoff(2) = min(find(specfreq>=14));
% Find 1st entry of specfreq >= 25 Hz
freq_cutoff(3) = min(find(specfreq>=25)); 
% Find 1st entry of specfreq >= 55 Hz
freq_cutoff(4) = min(find(specfreq>=55)); 
% Find 1st entry of specfreq >= 60 Hz
freq_cutoff(5) = min(find(specfreq>=60)); 
% Find 1st entry of specfreq >= 100 Hz
freq_cutoff(6) = min(find(specfreq>=100)); 

% % DEFINE THETA, SLOW GAMMA, AND FAST GAMMA CUTOFFS FOR PLOTTING SEPERATELY
% spectheta = spec(freq_cutoff(1):freq_cutoff(2), :); %Create a cropped matrix with only 4-14Hz
% specslowg = spec(freq_cutoff(3):freq_cutoff(4), :); %Create a cropped matrix with only 4-14Hz
% spectfastg = spec(freq_cutoff(5):freq_cutoff(6), :); %Create a cropped matrix with only 4-14Hz

% PLOT FULL SPEC
figure('Name',['Augmented Spectrogram'],'NumberTitle','off')
imagesc(abs(spec(spread(1):spread(2),:)),'y',specfreq(spread(1):spread(2)),'x',specTM, [200 20000])

% PLOT NON-AVERAGED VS AVERAGED SPECTROGRAM THETA
figure('Name',['Theta Non Averaged Vs Averaged Spectrogram'],'NumberTitle','off')
subplot(2,1,1)
imagesc(abs(spec(freq_cutoff(1):freq_cutoff(2),:)),'y',specfreq(freq_cutoff(1):freq_cutoff(2)),'x',specTM)
pan xon
xlabel ('Time') %x axis labeled as filename
ylabel ('4-14 Hz')
subplot(2,1,2)
imagesc(abs(specmat(freq_cutoff(1):freq_cutoff(2),:)),'y',specfreq(freq_cutoff(1):freq_cutoff(2)) ,'x',specTM, [10 500])
pan xon
ylabel ('4-14 Hz')

% PLOT NON-AVERAGED VS AVERAGED SPECTROGRAM SLOW GAMMA
figure('Name',['Slow Gamma Non Averaged Vs Averaged Spectrogram'],'NumberTitle','off')
subplot(2,1,1)
imagesc(abs(spec(freq_cutoff(3):freq_cutoff(4),:)),'y',specfreq(freq_cutoff(3):freq_cutoff(4)),'x',specTM,[200 5000])
pan xon
xlabel ('Time') %x axis labeled as filename
ylabel ('25-55 Hz')
subplot(2,1,2)
imagesc(abs(specmat(freq_cutoff(3):freq_cutoff(4),:)),'y',specfreq(freq_cutoff(3):freq_cutoff(4)),'x',specTM)
pan xon
ylabel ('25-55 Hz')

% PLOT NON-AVERAGED VS AVERAGED SPECTROGRAM FAST GAMMA
figure('Name',['Fast Gamma Non Averaged Vs Averaged Spectrogram'],'NumberTitle','off')
subplot(2,1,1)
imagesc(abs(spec(freq_cutoff(5):freq_cutoff(6),:)),'y',specfreq(freq_cutoff(5):freq_cutoff(6)),'x',specTM,[200 2500])
pan xon
xlabel ('Time') %x axis labeled as filename
ylabel ('60-100 Hz')
subplot(2,1,2)
imagesc(abs(specmat(freq_cutoff(3):freq_cutoff(4),:)), 'y' ,specfreq(freq_cutoff(5):freq_cutoff(6)),'x',specTM)
pan xon
ylabel ('60-100 Hz')






%% Short Forier for better frequency power analysis resolution
% Resource code: https://www.mathworks.com/matlabcentral/fileexchange/45197-short-time-fourier-transformation--stft--with-matlab-implementation
    % Slide Fast Fourier Trasnform over signal for a finer scale resolution
    %spectrogram
% Rename signal to edit during SFTF code 
X = fildata';   
fs = 1250; % fs is the sample rate of the data
window=128; %set length of window, here 100 ms
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
slide= 32; 
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
Fourier = fft(fildata);

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

%Find indices of freq_vector to set max cutoff at 100 Hz
freq_max = min(find(freq_vector>=100));

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
ylabel('Amplitude')
title('Amplitude spectrogram of the signal')

handl = colorbar;
set(handl) %, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(handl, 'Magnitude')
%% COMMENTED OUT PARTIAL CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% size(pos_bin)
% 
% pos_binned = zeros (28,21);
% for aloop = 29:28:560
%     for bloop = 22:21:420
%         pos_binned(((aloop-1)/28), ((bloop-1)/21)) = sum(sum(pos_bin((aloop - 28):28, (bloop - 21):bloop)));
%     end
% end
% % Smooth in 3 cm x 3cm spatial bins then upsample to match OE sample rate
% % for analysis
% % Paired time stamps do not match x,y position sampling exaclty, will have to smooth spatial data before can map trajectory 
% find(behavior(:,6) == 177411772)
% 
% for tp = 1:size(timestamp_mat, 1)
%     if timestamp_mat(1, tp) == behavior(:,6)
%     paired_behavior(timestamp_mat(1, tp)