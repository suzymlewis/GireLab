function [signal glomeruli] = get_signal(ROI, raw)

%%%%%%%%%%%%%% PLEASE NOTE THAT THIS FUNCTION RETRIEVES NON_OVERLAPPING ROI PIXELS ONLY
                %%%%%%%%%% ALL OVERLAPPING PIXELS ARE DROPPED IN THIS SIGNAL MATRIX

%%%%%%%%%%%%%%% LOADS RELEVANT VARIABLES FROM SEGMENTED DATA INTO WORKSPACE
    % INPUT
        % ROI : cell structure where each cell is ROI as designated by CalIm CNMF
            % Within each cell: m x 1  array
                % m (rows) : all pixel associated with cell's ROI (glomeruli) 
                  % Gives pixel index within the 116 x 152 concatenated pixel frame
    % OUTPUT
        % signal: m x n matrix
            % Ca signal for all pixels lie in ROI's designated by CalIm CNMF
            % Rows :  Pixels that fall within ROI's
            % Columns : samples for each timepoint sampled at 20 Hz
        % glomeruli: n x 1 array
            % ROI assignment for each pixel row in the signal matrix
%%%%%%%%%%%%%%% 

% Reshape raw for indexing with ROI assignments
raw_=reshape(raw, [(size(raw, 1)*size(raw, 2)), size(raw, 3), size(raw,4)]);

signal=[]; % Create empty matrix to store Raw downsampled (median filtered) Ca Imaging data
    
glomeruli = []; % Create an array to hold glomerular assigment for each row of signal matrix

for i = 1:length(ROI) % LOOP OVER GLOMERULI
    temp_indx = ROI{i}; % Grab indices for ROI pixels of looping glomeruli
    % This indices correspond to a 116 x 156 pixel tif image frame
    glom_temp = (ones(1, length(ROI{i})) * i)'; % Make an array with glomeruli number as long as number of glomeruli in ROI
    glomeruli= vertcat(glomeruli, glom_temp); % Update array so will have glomeruli assignment for each row of signal
    for ii = 1:size(temp_indx, 1) % LOOP THROUGH PIXELS
            signal((size(signal,1)+1),:, :) = raw_(temp_indx(ii),:,:);
    end
end
