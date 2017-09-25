%%EVENTS
%This includes both events that are generated by external devices (such as TTLs)
% and events that are emitted by modules within the application.
% timestamps = timestamps./info.header.sampleRate; % convert to seconds

%The messages.events file is just a plain ascii file with the messages text, 
    % not intended to be loaded by the load_open_ephys_data function. 
%The message events themselves, without the text, are stored in the all_channels.events too
    % being those of eventType = 5
%% Set variables

clear all

samprate=30000 %sampling rate in Hz, open message.events file for recording session to see start time and sampling rate
%% load events file
% Non-spike events are saved in a different format. Events generated by all channels are dumped into the file "all_channels.events" with the following data for each event:
% int64 timestamp (to align with timestamps from the continuous records)
% int16 sample position within a buffer
% uint8 event type (all the events that are saved have type TTL = 3 ; Network Event = 5)
% uint8 processor ID (the processor this event originated from)
% uint8 event ID (code associated with this event, usually 1 or 0)
% uint8 event channel (the channel this event is associated with)
% One uint16 recording number (version 0.2 and higher)
%EVENT OUTPUT
%     data = zeros(MAX_NUMBER_OF_EVENTS, 1); 'int64' part of event file
%     timestamps = zeros(MAX_NUMBER_OF_EVENTS, 1); 'uint64' part of event
%     file
%     info.sampleNum = zeros(MAX_NUMBER_OF_EVENTS, 1); 'int16' part of events file
%     info.nodeId = zeros(MAX_NUMBER_OF_EVENTS, 1); 'uint8' part of events
%       file
%     info.eventType = zeros(MAX_NUMBER_OF_EVENTS, 1); 'uint8' part of
%       events file
%     info.eventId = zeros(MAX_NUMBER_OF_EVENTS, 1); 'uint8' part of event
%       file
clear all
eventfilename=uigetfile('.events')
% [data, timestamps, info] = loadAndCorrectPhase(filename, 1); %loading data from open ephys
[data, timestamps, info] = load_open_ephys_data(eventfilename);
eventdata=data;
eventts=timestamps;
eventinfo=info;
sampleNum=getfield(eventinfo, 'sampleNum');
nodeId=getfield(eventinfo, 'nodeId');
eventType=getfield(eventinfo, 'eventType');
eventId=getfield(eventinfo, 'eventId');
%Make unique variable lists into array string or single string entry
%depending on length
eventinfovariables=[sampleNum, nodeId, eventType, eventId];
space='  '
for i=1:size(eventinfovariables, 2)
    eventinfovariable=unique(single(eventinfovariables(:, i)));
    if size(eventinfovariable,1)>1;
        infovariables{i}=num2str(eventinfovariable);
        strcatenation=' '
        for ii=1:length(infovariables{i})
            strcatenation=[strcatenation space infovariables{i}(ii, :)]
        end
        infov{i}=strcatenation
    else;
        infovariables{i}=int2str(eventinfovariable);
        infov{i}=infovariables{i}
    end
end
sprintf('The sample position within a buffer are %s', infov{1})
sprintf('The processors originating events in this session are %s', infov{2})
sprintf('The types of codes associated with this event [options 0,1] are %s', infov{3})
sprintf('The channels which have associated events are %s', infov{4})

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
cellnumber=cellnumber+1;
end


%% Drop unused rows from tmatrix
cellnumber=cellnumber-1;
tmatrix=tmatrix(1:cellnumber, :);
findmaxrow=zeros(1,cellnumber);
for i=1:cellnumber;
    findmaxrow(i)=max(length(nonzeros(tmatrix(i, :))));
end
tmatrix=tmatrix(:,1:max(findmaxrow));
%% NOW HAVE SPIKE TIMES FOR CELL
clear category
category=1;
eventts_by_cat=zeros(1, length(eventts));
%Check to make sure event file reasonable
for i=1:length(eventts);
    currentevent=[]
    if eventts(i+1)<(eventts(i)+1);
        eventts_by_cat(i)=category; %Creating a dictionary to pair with eventts in 2Xn matrix classify events with an event number
        currentevent=horzcat(currentevent,eventts(i));
        
        % disp('same event') % Use to troubleshoot code
    else
        category=category+1
        eventts_by_cat(i)=category;
        currentevent=horzcat(currentevent,eventts(i));
        % disp('new event') % Use to troubleshoot code
    end
    PSTH_curr=zeros(length(tmatrix(:,1), length(tmatrix(1,currentevent(start):currentevent(end))
    for i=1:length(tmatrix(:,1))
        for ii=1:length(currentevent)
            if ismember(currentevents(ii), tmatrix(i,:))
                PSTH_curr(i,:)=tmatrix(i,currentevent(start):currentevent(end))
            else
                PSTH_curr(i,:)=0
            end
        end
    end
end