%%%%%%%%%% CODE TO CREATE IMAGES AND RUN ANALYSES FOR BIO COMP CONFERENCE POSTER IN NICE, FR -- SUMMER 2018
%% LOAD CNMF DATA
[A_or C_or trials frame_trial ODOR_ON] = CalImdata();

%% LOAD COMPRESSED DATA
[raw raw_trials raw_frame_trial] = Compdata();

%% FIND ROIS FOR GLOMERULI IN TRIAL
[A ROI ROIweights] = ROImasking(A_or);

% Please note that ROI does not have overlapping pixels between glomeruli
    % ROIweights will have pixels that repeat across 2+ glomeruli
        % Generally only repeated over 2

%% EXTRACT ALL ROIS INTO DATA MATRIX
% Raw compressed
[signal glomeruli] = get_signal(ROI, raw);

% Weighted compressed acording to CNMF
[signal_w glomeruli_w] = get_weighted_signal(ROIweights, raw);

%% CLEAN WORKSPACE
clear A A_or C_or frame_trial

%% IMPORT ETHANOL
mxf = 450; % Set max frame for which to extract relevant ethanol signal
[Eth] = ethdata(mxf);

%% GET DELTA F/F
% Label trial categories
turb_label = get_turb(); % Create turbulence label for each trial of session

% Raw compressed
Sig_or = reshape(signal, size(signal, 1), size(signal, 2) * size(signal,3)); % Reshape to fit into C_or format for DFFcal function
[Sig Sig_od Sig_trials Sig_trials_od] = DFFcal(Sig_or, trials ,turb_label); % Get localized DF/F signal for each trial for each pixel

% Weighted compressed acording to CNMF
Wig_or = reshape(signal_w, size(signal_w, 1), size(signal_w, 2) * size(signal_w,3)); % Reshape to fit into C_or format for DFFcal function
%%%%%%%%%%%%% FIX BECAUSE CURRENTLY WIG OR == SIG OR
[Wig Wig_od Wig_trials Wig_trials_od] = DFFcal(Wig_or, trials ,turb_label); % Get localized DF/F signal for each trial for each pixel


%% AVERAGE AND WEIGHTED AVERAGED DELTA F TRACES
% Average over raw compressed pixels within glomeruli
[Glom_Sig resp] = ROI_sig(Sig_trials, glomeruli);

% Average over weighted compressed pixels within glomeruli
[Glom_Sig_W] = ROI_sig(Wig_trials, glomeruli_w);

% Reshape to regain trial structure if needed
Glom_Sig_trials = reshape(Glom_sig_df, size(Glom_sig_df, 1), size(Sig_trials, 2), size(Sig_trials, 3));
Glom_Sig_W_trials = reshape(Glom_wig_df, size(Glom_wig_df, 1), size(Wig_trials, 2), size(Wig_trials, 3));

%% NOW DO ANALYSIS DAVID WAS TALKING ABOUT
% Going to take calcium imaging and look for response >1sd from basline response
% Grab eth signal corresponding to +-2 seconds from the sig ca event
% Bootstrap 1000 trial shuffled eth signal?
% Bootstrap 1000 trial shuffled ca to get signal?
% Set 95% CI
% See If sig signal
% 
% 
% Next take xcorr of all sig calcium response to eth dynamics
%     average delay within each glomeruli
%     
% Use this delay to create eth deriv by ca deriv plots
%     plot all timepoints within a single glom
%     cubic spline fit to see if signal have relation
%     maybe try this plotting all pixels within each glom seperately

%% ALSO DO STUFF USING PCA WEIGHTS AS GLOM SIG

%% FIND X-CORR OF SIGNAL
crc = get_xcor (signal, Eth, glomeruli);

%% FIND OPTIMAL DELAY FOR EACH GLOMERULI SO CAN CORRELATE WITH ETH SIGNAL
% See working get delay function, don't start from scratch



