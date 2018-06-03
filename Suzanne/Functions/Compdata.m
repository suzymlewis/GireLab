function [raw trials frame_trial] = Compdata()

%%%%%%%%%%%%%%% LOADS RELEVANT VARIABLES FROM SEGMENTED DATA INTO WORKSPACE
    % OUTPUT
        % raw: raw signal with most blank pixels corresponding to motion correction boundary dropped
            % 116pixels x 156pixels x 450ts x 40trials
        % trials: trial numbers included in data
        % frame_trial: indicates which trial each frame is from
            % Entire session in temporal matrix horizontally concatenated
%%%%%%%%%%%%%%%

% Import compressed recording
[fc, fdc] = uigetfile('.mat', 'CHOOSE COMPRESSED DATA - .mat');
cd(fdc);
comp = load(fc);

% Extract relevant variable
Y_ca = comp.clc; % Downsampled raw data using median filter for neighboring pixels
Y_ca = Y_ca (5:120, 5:160,:); % Drop areas of frame removed from motion correction to be consistent with segmented data
frame_trial = comp.tr1a; % Vector containing trial assignment for each frame

trials = max(frame_trial); % Designate number of trials for output

for i = 1:trials % LOOP THROUGH TRIALS
        k = find(frame_trial ==i); % Get indices current trials samples
        temp = Y_ca(:,:,k); % Create array with only current looping trials information
        raw(:,:,:,i) = temp(:,:,1:450);
            % Since each trial is slightly different number of frames
                % Grab the first 450 frames and store in 4D matrix
        clear temp
end
