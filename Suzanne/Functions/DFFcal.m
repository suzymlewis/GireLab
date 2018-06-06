function [C C_od C_trials C_trials_od]=DFFcal(C_or, trials ,turb_label)
%%%%%%%%%%%%%%% DFF GLOMERULAR CALCIUM
    % INPUT
        % C_or
        % trials
         % turb_label = vector containing turbulence level for each trial
          % within session
            % 0 = Low
            % 1 = Med
            % 2 = High
    % OUTPUT
        % C : 5 x 1 cell structure
                % Cell 1 : All trials
                % Cell 2 : High turbulence trials only
                % Cell 3 : Medium turbulence trials only
                % Cell 4 : Low turbulence trials only 
                % Cell 5 : High and Low turbulence trials only
        % C_od : Same as above but only includes part of odor exposure from each trial
        % C_trials : m x n x r
            % m : number of pixels across all ROI's
            % n : number of timepoints sampled in each trial
            % r : number of trials
        % C_trials_od : same as C_trials but only includes portion of odor exposure during each trial
%%%%%%%%%%%%%%%

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
    if i == 1
        C_trials=C_temp;
        C_trials_od =  C_odor;
    end
    C_od{i} = reshape(C_odor, size(C_odor, 1), size(C_odor, 2)*size(C_odor, 3));
    C{i} = reshape(C_temp, size(C_temp, 1), size(C_temp, 2)*size(C_temp, 3));
end
clear ext_sig org_sig mean_sig dff_sig C_temp
