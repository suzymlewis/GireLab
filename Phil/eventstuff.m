% load events file
clear all
eventfilename=uigetfile('.events')
% [data, timestamps, info] = loadAndCorrectPhase(filename, 1); %loading data from open ephys
[data, timestamps, info] = load_open_ephys_data(eventfilename);
eventdata=data;
eventts=timestamps;
% sampleRate=1250; %Rewrite in Hz sampling rate
% Downsample for 30kHz to 1kHz
ts=timestamps(1:1000:length(timestamps));
a=zeros(1,length(ts));
b=zeros(1,(length(ts)-1));
for i=1:length(ts);
    a(i)=ts(i)/ts(end);
    if i>=2;
        b(i)=a(i)-a(i-1);
    end
    i=i+1;
end

% figure
% plot(b)
% hold on
% plot([0,length(b)], [mean(b), mean(b)], 'k', 'Linewidth', 2)
%% LOAD SPIKES
spikefilename=uigetfile('.spikes')
[data, timestamps, info] = load_open_ephys_data(spikefilename);]
spikedata=data;
spikets=timestamps;

%% Import t file and plot spikes as a PSTH over time within a single cell
clear allb
[Fname pname] = uigetfile('*.t','multiselect','on')

samprate = 30000; %samprate is 32K

cd(pname)
S = LoadSpikes_VT(Fname);

for b = 1:length(S);
TS = S{b}.T; %this is the index for each cell

Timestamps = TS;%./samprate;

figure, hist(Timestamps,100); %plot the timestamps as a PSTH

end

