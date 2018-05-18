%% GLOM PCA
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

%clear g_data to save space in workspace
clear g_data;

% Extract relevant fields for glom signal
    % C_or the refined, merged ROI temporal data
    % A_or the refined, merged ROI spatial data

% Extract glomerular traces for each trial

%% LOAD AND EXTRACT ETH BY TRIAL TYPE

% Uncomment correct trial format (HIGH-LOW or LOW-HIGH) and make necessary alterations
    % See behavioral excel file to determine session structure
    
% HIGH-LOW        
turb_label(1:10) = 1;   
turb_label([11:15 21:25 31:35]) = 2;    
turb_label([16:20 26:30 36:40]) = 0;

% % LOW-HIGH
% turb_label(1:10) = 1;
% turb_label([11:15 21:25 31:35]) = 0;
% turb_label([16:20 26:30 36:40]) = 2;


% Load ethanol data
[fe, fde] = uigetfile('.dat', 'CHOOSE ETH DATA - .dat');
cd(fde);

fna = char(fe);
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

for p = 1:length(un); % Loop through trials
    frr = FR(TR==un(p)); %
    frr = frr-frr(1)+1;
    
    eth = ETH(TR==un(p)); 
    
    if length(eth)>1000;
        eth = filtfilt(B,A,eth);  % Perform zerp-phase digital filtering by processing the eth signal in both the forward and reverse directions 
    end
 
    eth = diff(eth); % Take derivative of filtered eth signal
    
    ts = TS(TR==un(p)); % Get timestamps associated with each trial
    
    if length(eth)>=eth_mxf & p<=length(turb_label)
        xcc = xcc+1;
        ETT(xcc,:) = eth(1:eth_mxf);
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

%% STFT PLUMES SEPERATED BY TRIAL



%% DFF GLOMERULAR CALCIUM

% 3D mat with glom calcium signal across each timepoint across each trial
C_trial = reshape(C_or,size(C_or,1),size(C_or,2)/length(trials),length(trials));

% Create subsets by turbulence category for dataset
% High turb only
    % Includes trials with audio
C_trial_high = C_trial(:,:, find(turb_label==2));
% 
% figure
% for i = 1:size(C_trial_low, 1)
%     plot(C_trial_low(i,:,1))
%     hold on
% end

% plot
% Low turb only
    % Includes trials with audio
C_trial_med = C_trial(:,:, find(turb_label==1));
% Med turb only
C_trial_low = C_trial(:,:, find(turb_label==0));

% Only high and low flow trials combined
C_highandlow = C_trial(:,:, find(turb_label~=1));

% Create structure with all data sets
C = {C_trial, C_trial_high, C_trial_med, C_trial_low, C_highandlow};


% Set params and create guassian window
rate = 20; %Hz % should be an even number
w = gausswin(rate*2+1);
w = w/sum(w);
odor_only=130:350;
% Loop over each glomruli to create F-mean for DFF
for i=1:size(C,2)
    C_temp = C{i};
    for ii = 1:size(C_temp, 1)
        for iii = 1:size(C_temp, 3)
            % Extract glom for looping
            org_sig = C_temp(ii,:, iii)';
            
            % Invert edges of lenfth window size to avoid edge artifacts in local delta F
            ext_sig = [flipud(org_sig(1:rate)); org_sig; flipud(org_sig(end-rate+1:end))];
            
            % Delta F vector
            mean_sig = conv(ext_sig,w);
            
            % Calc delta F and assign vector to data subset
            dff_sig = (org_sig-mean_sig(rate*2+1:end-rate*2))./mean_sig(rate*2+1:end-rate*2);
            C_temp(ii, :,iii) = dff_sig';
          % Extract dff only shortly before and during odor exposure
            C_odor(ii, :, iii) = dff_sig(odor_only);
        end
    end
    % Concatenate like trials within each glomeruli for PCA analysis for
    % both odor only and entire relevant time
    C_od{i} = reshape(C_odor, size(C_odor, 1), size(C_odor, 2)*size(C_odor, 3));
    C{i} = reshape(C_temp, size(C_temp, 1), size(C_temp, 2)*size(C_temp, 3));
end
clear ext_sig org_sig mean_sig dff_sig C_temp

% Look at signal from a single glom in the high flow condition across a
% trial
% aaa = C{2};
% figure, plot(aaa(1,:))
%% PCA SEGMENTAL BY TRIAL TYPE
clear EXP
for i=2:size(C_od,2)
    C_temp = C_od{i};
    C_temp = C_temp';
    %C_temp = C_temp(:,[1 3:14]);
    for j = 1:size(C_temp,2)
        C_temp(:,j) = C_temp(:,j)-mean(C_temp(:,j));
    end
    corr = C_temp'*C_temp;
    % Loop to see if corr has any NAN glom
        % If does, drops glom, stores glom number in dropped_glom vector,
        % and recalcs corr without NAN glom
        if length(corr(~any(~isnan(corr), 2),:))>=1;
            dropglom=1;
            for ii = 1:length(corr)
                if length(find(~isnan(corr(:,ii)*zeros(1,size(corr,2)))))==0;
                    C_temp(:,ii)=[];
                    dropped_glom(dropglom)=ii;
                    dropglom=dropglom+1;
                end
            end
            corr =C_temp'*C_temp;
            clear dropglom
        end
        
    [w,d] = eig(corr);
    d = diag(d)/size(C_temp,1);
    %w2 = w(:,end-1);
    act = C_temp*w;
    dim{i} = sum(d)^2/sum(d.^2);
%     [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(C_temp);
%     size(EXPLAINED)
%     EXP{i} = EXPLAINED;
end
% plot(dim)
% %%
% 
% figure, imagesc(EXP{5}, [0 22])
% figure, imagesc(EXP{2}, [0 22])
% figure, imagesc(EXP{4}, [0 22])
% 
% figure, plot(EXP{5})
% hold on
% plot(EXP{2})
% hold on
% plot(EXP{1})
% hold on
% plot(EXP{3})
% hold on
% plot(EXP{4})
% 
% for i = 1:length(corr)
%     a(i) = find(~isnan(corr(:,i)*zeros(1,size(corr,2))));
% end
%% Quick Dimensionality analysis trial-by-trial
clear EXP EXP1 spd1
for y = 1:length(trials);
    CA = C_trial(:,1:100,y);
    
    CA = CA-mean(CA')';
    CA = CA./std(CA,[],2);
    
    
    CO = C_trial(:,180:300,y);
    
    CO = CO-mean(CO')';
    CO = CO./std(CO,[],2);
    
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(CA);
    if y>1;
    EXP(y,1:size(EXPLAINED,1)) = EXPLAINED;
    else
    EXP(y,:) = EXPLAINED;
    end
    
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(CO);
    if y>1;
    EXP1(y,1:size(EXPLAINED,1)) = EXPLAINED;
    else
    EXP1(y,:) = EXPLAINED;
    end
    
    DEV_eth(y) = std(eth_data(y,180:300));
    
    run21 = fft(run2(y,:));
    run21([1:5 end-4:end]) = 0;
    run21 = real(ifft(run21));
    
    RUN(y,:) = run21;
    
    SPD(y) = std(RUN(y,1:100));%sum(abs(spd(50:60)))./sum(abs(spd));
    
    SPD1(y) = std(RUN(y,180:300));
end


[p1 p2] = sort(SPD1);

figure, imagesc(EXP1(p2,:))
figure, imagesc(RUN(p2,180:300));

figure, plot(SPD1,sum(EXP1(:,1:20),2),'bo')

for y = 1:size(EXP1,1);
CS = cumsum(EXP1(y,:));
CS1(y) = min(find(CS>90));
end

figure, plot(SPD1,CS1,'bo')

figure, plot(mean(EXP1(p2(1:10),:))-mean(EXP1(p2(end-10:end),:)))

figure, plot(SPD1,mean(abs(corr_mat)),'bo');


%% PCA ON ENTIRE DATASET