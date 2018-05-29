function [A_or C_or trials frame_trial ODOR_ON] = CalImdata()


%%%%%%%%%%%%%%% LOADS RELEVANT VARIABLES FROM SEGMENTED DATA INTO WORKSPACE
    % OUTPUT
        % A_or: Spatial matrix from CNMF
        % C_or: Temporal matrix from CNMF
        % trials: trials included in data
        % frame_trial: indicates which trial each frame is from
            % Entire session in temporal matrix horizontally concatenated
        % ODOR_ON: index for first 350 frames of data
[fg, fdg] = uigetfile('.mat', 'CHOOSE GLOMULAR DATA - .mat');
cd(fdg);
g_data = load(fg);
trials = g_data.trials;
ODOR_ON = g_data.ODOR_ON; % Index for relevant analysis time of video 
    % Keeps baseline 10 seconds before odor presentation    
    % Dumps a few seconds at tail end when no odor present
        % This ensures all trials have the same length
A_or = g_data.A_or; % Spatial Matrix from CNMF segmentation
    % A_or is the set of glomeruli once merged and updated
        % A is all originally identified glomerular ROIs before merging
    % A_or;
        % 18096 x N matrix
            % N = number of glomeruli
            % 18096 is the number of pixels in each 116 x 156 downsampled tiff image concatentated
C_or = g_data.C_or;
    % Temporal Matrix
frame_trial = g_data.tr; % What trial each frame is from
