function [ETT, ETT_FR, TSms] = Eth_x_trial (fn_dat, turb_label, trials_using)
%% LOAD AND EXTRACT ETH BY TRIAL TYPE

% Uncomment correct trial format (HIGH-LOW or LOW-HIGH) and make necessary alterations
    % See behavioral excel file to determine session structure
    

fna = char(fn_dat);
fileID = fopen(fna);
data = fread(fileID,'int32','ieee-be');
    
n_ten = find(data==-10);
TR = data(n_ten+1); % Trial number
TS = data(n_ten+2); % Timestamps in ms
FR = data(n_ten+3); % Frame number
ETH = data(n_ten+4); % Ethanol signal
RUN = data(n_ten+5); % Running as measured by wheel position
% DST = data(n_ten+7); % Not measured - Empty Input
% SNF = data(n_ten+9);  % Not measured - Empty Input
    
un = unique(TR); % Create vector with trial numbers of file
FR = FR-FR(1)+1; % Start frame number at 1
RUN = medfilt2(RUN,[3 1]); % Take median filter across every 3 samples within each trial vector for the wheel position

% Set filter params for zero phase filtering of ethanol signal
Wn = 8/50;
[B,A] = butter(3,Wn,'low');

% Maximum number of ethanol time points to import
eth_mxf = 3500; 

%start xcc counter for looping through eth trials
xcc = 0;

for p = 1:length(trials_using); % Loop through trials
    frr = FR(TR==un(p));%
    frt = frr-frr(1)+1;
    ETT_FR(p, :) = frr(1:eth_mxf)';
    eth = ETH(TR==un(p)); 
    
    if length(eth)>1000;
        eth = filtfilt(B,A,eth);  % Perform zerp-phase digital filtering by processing the eth signal in both the forward and reverse directions 
    end
 
    eth = diff(eth); % Take derivative of filtered eth signal
    
    ts = TS(TR==un(p)); % Get timestamps associated with each trial
    
    if length(eth)>=eth_mxf & p<=length(turb_label)
        xcc = xcc+1;
        ETT(xcc,:) = eth(1:eth_mxf);
        TSms(xcc, :) = TS(1:eth_mxf);
        trial_num(xcc) = p;
        turbulence(xcc) = turb_label(p);
    end
end

% clear xcc counter
clear xcc

figure('NumberTitle', 'off', 'Name', 'ETHANOL TRACES')
imagesc(ETT(:,:))

% Only odor on indices
ETT_odor = 800:1800;

