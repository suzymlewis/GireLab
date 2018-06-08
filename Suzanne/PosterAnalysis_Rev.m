%%%%%%%%%% CODE TO CREATE IMAGES AND RUN ANALYSES FOR BIO COMP CONFERENCE POSTER IN NICE, FR -- SUMMER 2018
%% LOAD CNMF DATA
[A_or C_or trials frame_trial ODOR_ON] = CalImdata();

%% LOAD COMPRESSED DATA
[raw raw_trials raw_frame_trial] = Compdata();

%% [I] FIND ROIS FOR GLOMERULI IN TRIAL
[A ROI ROIweights] = ROImasking(A_or);

% Please note that ROI does not have overlapping pixels between glomeruli
    % ROIweights will have pixels that repeat across 2+ glomeruli
        % Generally only repeated over 2

%% EXTRACT ALL ROIS INTO DATA MATRIX
% Raw compressed
[signal glomeruli] = get_signal(ROI, raw);
% 
% % Weighted compressed acording to CNMF
% [signal_w glomeruli_w] = get_weighted_signal(ROIweights, raw);

%% CLEAN WORKSPACE
clear A A_or C_or frame_trial raw raw_frame_trial ROI ROIweights

%% IMPORT ETHANOL
mxf = 450; % Set max frame for which to extract relevant ethanol signal
[Eth] = ethdata(mxf);

%% [II] AVERAGE ACROSS PIXELS WITHIN ROI

ROI_or = reshape(signal, size(signal, 1), size(signal, 2) * size(signal,3)); % Reshape to fit into C_or format for DFFcal function
for i = 1:length(unique(glomeruli(:))) % LOOP ROI
   ROI_sig(i,:) = mean(ROI_or(find(glomeruli==i),:));
end

%% [III] GET DELTA F/F
% Label trial categories
turb_label = get_turb(); % Create turbulence label for each trial of session

% Raw compressed
Sig_or = reshape(signal, size(signal, 1), size(signal, 2) * size(signal,3)); % Reshape to fit into C_or format for DFFcal function
[Sig Sig_od Sig_trials Sig_trials_od] = DFFcal(Sig_or, trials ,turb_label); % Get localized DF/F signal for each trial for each pixel

% ROI equally weighted average compressed
[Rig Rig_od Rig_trials Rig_trials_od] = DFFcal(ROI_sig, trials ,turb_label); % Get localized DF/F signal for each trial for each pixel

%% CLEAN WORKSPACE
clear Sig_od Sig_trials_od Rig_od Rig_trials_od

%% AVERAGE AND WEIGHTED AVERAGED DELTA F TRACES
% Average over raw compressed pixels within glomeruli
[Glom_Sig resp] = ROI_sig(Sig_trials, glomeruli);

% Reshape to regain trial structure if needed
Glom_Sig_trials = reshape(Glom_Sig, size(Glom_Sig, 1), size(Sig_trials, 2), size(Sig_trials, 3));

%% PARAMS FOR STA
% SET PARAMS
    % Get Relevant time frame
    ODOR_ON=[1:350];
    % Set Eth : Here using derv
    PL_MEAN = Eth(:,1:max(ODOR_ON)-1)';
    % Set Cal : Here DFF not diff(DFF)
    TR_MEAN = Rig_trials(:, 1:max(ODOR_ON), :);
    % ROI x timepoints x trials
    TR_MEAN1 = diff(Rig_trials, [], 2);

    
%% [IV] STA FEATURE DETECTION
% Set the threshold
thr = .5;
tm_rn = 20; %range in frames before and after
mx = size(PL_MEAN,1);
xc = 0;
clear AVG gl_id AVG_gl

xc = 0;
clear AVG gl_id AVG_gl crc tr_id
for i = 1:size(TR_MEAN1,1);
    crr = squeeze(TR_MEAN1(i,:,:));
    for ii = 1:size(crr,2);
        clear f2 crr1 zs df1 f1
        crr1 = squeeze(crr(:,ii));
        
        % Normalize response of each ROI so can compare
        zs = crr1-mean(crr1(:,:));
        zs = zs./std(zs(:));
        zs = diff(zs);
        
        % Make sure exceed threshold and seperate events
        f1 = find(zs>thr);
        df1 = diff(f1);
        f2 = f1(df1>5);
        
        %Get cross correlation
        crc(i,:,ii) = xcorr(PL_MEAN(:,ii), diff(squeeze(TR_MEAN(i,:,ii))), 'coeff');
        
        for u = 1:length(f2);
            
            if f2(u) > tm_rn & f2(u) < mx-tm_rn; % Drop first and last 20 frames to remove window artifacts or xcorr
                xc = xc+1;
            AVG(xc,:) = squeeze(PL_MEAN(f2(u)-tm_rn:f2(u)+tm_rn,ii))';
            AVG_gl(xc,:) = zs(f2(u)-tm_rn:f2(u)+tm_rn);
            gl_id(xc) = i;
            tr_id(xc) = ii;
            
            end
            
        end
    
    end
    
end

% Average across trials to get mean response of ROI
resp_mean = mean(TR_MEAN, 3); %mean of each cell's response over all trials
figure, imagesc(resp_mean);

% Order ROIs based on magnitude of maximum response
resp_str = max(resp_mean(:, ODOR_ON, :), [], 2);
[rs rs1] = sort(resp_str);
resp_mean1 = resp_mean(rs1,:);


% Grab normalized cross correlation responses exceeding threshold
corr_resp = mean(crc(rs1,:,:),3)./std(crc(rs1,:,:),[],3);
corr_resp1 = corr_resp;
corr_resp(abs(corr_resp)<.8) = 0;
cr1 = max(corr_resp,[],2);

% Size of signal from which take cross-correlation
sz = round(size(corr_resp,2)./2);

% Plot the mean response of ROIs with events
    % Beneath plot the mean xorr of same ROI
figure, subplot(2,1,1), 
imagesc(resp_mean(rs1,:))
subplot(2,1,2),
imagesc(corr_resp(:,sz-40:sz+40));

figure, subplot(2,1,2), imagesc(corr_resp1(cr1>0,sz-40:sz+40))
subplot(2,1,1), imagesc(resp_mean1(cr1>0,:))

LG = ([sz-40:sz+40]-mean([sz-40:sz+40]))./20;

figure, plot(LG,corr_resp1(cr1>0,sz-40:sz+40)','linewidth',2)
hold on
plot([0 0],[-2 2],'k--')

%% GET MEAN DELAY USING CROSS CORRELATION FROM STA FEATURE ANALYSIS
wind = 40; % Indicate size of frame lag in which to look for max
    % At 20 Hz imaging this is 2 seconds pre and post
delay = delay_STA(corr_resp1, wind);

%% [VI] PLOT THE CA VERSE ETH DERIVATIVES USING THE PROPER DELAY
ODOR_ON = [130:325]
plt_rltn(Rig_trials, Eth, delay, ODOR_ON)

%% LOGISTIC REGRESSION MODEL TO FIND MOST INFORMATIVE ROIs
    % https://www.biorxiv.org/content/biorxiv/early/2018/03/13/281444.full.pdf
%% REGULARIZED LOGISTIC REGRESSION MODEL (RLRM)
    % https://www.biorxiv.org/content/biorxiv/early/2018/03/13/281444.full.pdf
% Create PALS
    % Population activity Vectors from all pixels
    % Average population activity for each trial was computed as the mean
    % response across all pixels included in the FOV
    % So 15 High PALs and 15 Low PALs
    % Data distribution for decoder
        % Train : Use 80% data
        % Optimize : 10 % To find optimal lambda (enforced sparsity factor)
        % Test : 10 % 
    % Decoders performed
        % 1- Trained and tested within trial turbulence condition 
        % 2- Trained in one condition and tested in another

    pals =     


%% GRAB ALL EVENTS EXCEEDING 2 STD FROM BASELINE
stdev = std(Glom_Sig, 1, 2); % Grab st dev (1 refers to use sample size for sample not population)
events = [];

% Add indices for beginning and end or trials to vector
    % To make sure eth signal trying to grab not events in 1st or last sec
no = [1:40 size(Glom_Sig_trials,2)-40:size(Glom_Sig_trials,2)];

for i = 1:size(Glom_Sig_trials, 1) % LOOP GLOMERULI
    et = find(abs(Glom_Sig_trials(i,:, :)) > mean(Glom_Sig(i,:)) + (2.5 * stdev(i)));
        % Even though dff close to mean zero, is not perfectly centered
    trtp = ceil(et/40);    
    drop = find(ismember(trtp, no));
    if size(drop,1)>0
        et(drop) = [];
        trtp(drop) = [];
    end
    gt = ones(length(et), 1) * i; % Store ROI relevant to event
    tr = ceil(et/size(Glom_Sig_trials, 2)); % Store trial number corresponding to events
    events = vertcat(events, horzcat(et, gt, trtp, tr));
end

%% GRAB DERIV CA SIGNAL SURRONDING EVENTS
events(:, 5:85) = 0; % Create columns for 2 seconds of eth on either side event

Glom_Sig_diff = diff(Glom_Sig_trials,[],2);

for i = 1: size(events, 1) % LOOP ROIs 
    events(i, 5:85) = Glom_Sig_diff(events(i,2), events(i,3)-40:events(i,3)+40, events(i,4));
end

%% GRAB ETH SIGNAL CORRESPONDING TO EVENTS
events(:, 86:166) = 0; % Create columns for 2 seconds of eth on either side event

for i = 1: size(events, 1) % LOOP EVENTS
    events(i, 86:166) = Eth(events(i,4), events(i,3)-40:events(i,3)+40);
end

%% CROSS CORR OF TWO DERIVS
clear x_corrs
for i = 1:size(events, 1) % LOOP EVENTS
x_corrs(i, :) = xcorr(events(i, 5:85), events(i, 86:166), 'coeff');
end

%% CROSS CORRS OF TRIAL SHUFFLED AVERAGE
clear null_x_corr
for i = 1:size(events, 1) % LOOP EVENTS
    Eth_null = find(~ismember(1:size(Glom_Sig_trials, 3), events(i,4)));
    for ii=1:length(Eth_null); % LOOP OTHER TRIALS ETH SIGNAL
    null_x_corr(ii, :, i) = xcorr(events(i, 5:85), Eth(Eth_null(ii), 86:166), 'coeff');
    end
end

%% GET CONFIDENCE INTERVAL FOR TRIAL SHUFFLED ETHANOL SIGNAL
clear CI
for i = 1:size(null_x_corr, 3) % LOOP EVENTS
    for ii = 1:size(null_x_corr,2) % LOOP TIMEPOINTS
        CI(1,ii,i) = prctile(null_x_corr(:,ii,i), 97.5); % Get upper CI threshold
            % No assumed distribution
                % Just get threshold above which lies 2.5% dist density
        CI(2,ii,i) = prctile(null_x_corr(:,ii,i), 2.5); % Get lower CI threshold
            % Same but threshold below which lies 2.5% dist density
    end
end

% Plot permutations for cross-corr CI
figure, plot(CI(1,:,1000))
hold on
plot(CI(2,:, 803))
hold on
plot(x_corrs(800,:))

%% FIND EVENTS WHERE MAX XCORR EXCEEDS CONFIDENCE INTERVAL FOR SIG CA RESPONSE
clear x_corrs_delay
for i=1:max(unique(events(:,2))) % LOOP ROIs
    roix = x_corrs(find(events(:,2)==i), :); % Grab xcorr's for looping ROI's
    [lag lag_in] = max(roix(:,61:101), [], 2); % Get value and index for max cross corr
    % ROI assignment for event max xcorr
    x_corrs_delay(find(events(:,2)==i), 1) = ones(length(lag), 1)*i;
    x_corrs_delay(find(events(:,2)==i), 2) = events(find(events(:,2)==i), 3);
    % Max xcorr value for event max xcorr
    x_corrs_delay(find(events(:,2)==i),3) = lag;
    % Delay for event max xcorr (in frames) @ 20 Hz
    x_corrs_delay(find(events(:,2)==i),4) = lag_in;
    % Corresponding delay profile for ROI (in frames) @ 20 Hz
    x_corrs_delay(find(events(:,2)==i), 5) = ones(length(lag), 1)*mean(lag_in);
end