function [Eth] = ethdata(mxf);
%% Import the ethanol data
%%%%%%%%%%%%%%% LOADS RELEVANT VARIABLES FROM SEGMENTED DATA INTO WORKSPACE
    % INPUT (OPTIONAL)
        % mxf: max frame to which corresponding ethnol data retrieved
            % Default: 400 frames, which is the majority of the trial
            % Do not take sample for all frames automatically since frame number can vary across trial
    % OUTPUT
        % Eth: Derivative of ethanol signal
            % Samples averaged across each frame
%%%%%%%%%%%%%%% 

% DETERMINE MAX FRAME FROM WHICH TO GATHER CORRESPONDING ETH SIGNAL
if nargin ==1 % Check to see if user input to designate max frame
    mxf = mxf;
else % If no input designate towards end of imaging time for each trial
    mxf = 400;
end

[fn pn] = uigetfile('*.dat', 'CHOOSE ETH DATA - .dat');
cd(pn)
fileID = fopen(fn);

data = fread(fileID,'int32','ieee-be');

n_ten = find(data==-10);
TR = data(n_ten+1);
TS = data(n_ten+2);
FR = data(n_ten+3);
ETH = data(n_ten+4);
DST = data(n_ten+7);
un = unique(TR);
FR = FR-FR(1)+1;
clear fr eth_fr ts_fr dst unf

for p = 1:length(un);
    frr = FR(TR==un(p));
    frr = frr-frr(1)+1;
    eth = ETH(TR==un(p));
    ts = TS(TR==un(p));
    ds = DST(TR==un(p));
    
    unf = unique(frr);
    for t = 1:length(unf);%min([mxf+25 max(unf)]);%length(unf);
        fr(p,t) = t;%mean(frr(frr==unf(t)));
        data2(p,t) = mean(eth(frr==unf(t)));
        ts_fr(p,t) = mean(ts(frr==unf(t)));
        dst(p,t) = mean(ds(frr==unf(t)));
        tr(p,t) = mean(ds(frr==unf(t))).*0+un(p);
    end
    data2a = data2(p,min(find(ts_fr(p,:)>1)):end);
    data2b(p,:) = data2(p,1:mxf);
    
    t2 = ts_fr(p,min(find(ts_fr(p,:)>1)):end);
    t2b(p,:) = t2(1:mxf);
end

clear data2 ts_fr
fls = [1:40];%files used for imaging; 2-20 for 7-21-2017_thy1-35; 2-24 for thy1-27; 1-10 for retro-AAV
data2 = data2b(fls,:);
dst = dst(fls,:);

data2 = data2-repmat(mean(data2(:,1:100),2),1,size(data2,2));
ts_fr = t2b;

 figure, imagesc(diff(data2,[],2));
 xlabel('Ethanol Sensor Data')
 
 Eth=diff(data2,[],2);
 trial_frame = fr;
