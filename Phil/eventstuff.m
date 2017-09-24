%% Set variables
samprate=30000 %sampling rate in Hz
%% load events file
clear all
eventfilename=uigetfile('.events')
% [data, timestamps, info] = loadAndCorrectPhase(filename, 1); %loading data from open ephys
[data, timestamps, info] = load_open_ephys_data(eventfilename);
eventdata=data;
eventts=timestamps;
%% Import t file and plot spikes as a PSTH over time within a single cell
clear allb
cellnumber=1;
tmatrix=zeros(500,length(eventts));
[Fname pname] = uigetfile('*.t','multiselect','on')

% samprate = 30000; %samprate is 32K

cd(pname)
S = LoadSpikes_VT(Fname);

for b = 1:length(S);
    TS = S{b}.T; %this is the index for each cell

    Timestamps = TS;%./samprate;
    Timestamps(1:10);
    size(Timestamps);
    tmatrix(cellnumber,1:length(Timestamps))=Timestamps;
% figure, hist(Timestamps,100); %plot the timestamps as a PSTH
cellnumber=cellnumber+1
end


%% Drop unused rows from tmatrix
cellnumber=cellnumber-1;
tmatrix=tmatrix(1:cellnumber, :);
findmaxrow=zeros(1,cellnumber);
for i=1:cellnumber;
    findmaxrow(i)=max(length(nonzeros(tmatrix(i, :))));
end
tmatrix=tmatrix(:,1:max(findmaxrow));
%% 










% %% CREATE LOOP THAT LOADS .t FILES into a matrix of arrays
% tmatrix=zeros(cellnumber,length(eventtimestamps))
% for k=1:cellnumber;
%     tfilename=strcat('TT',  num2str(date), '_', num2str(k))
%     if exist(tfilename, 'file')
%         datacurrent=load(tfilename);
%     else
%         fprintf('File %s does not exist./n', tfilename)
%     end
%     tmatrix(k,:)=datacurrent
%     %now that a .t file has been loaded, extract its data into a array and
%     %simultaneously concatenate vertically to .t matrix for given date of
%     %interest
%     k=k+1
% end
