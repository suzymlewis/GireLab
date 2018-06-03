function [signal glomeruli] = get_signal(ROI, raw)

%%%%%%%%%%%%%% PLEASE NOTE THAT THIS FUNCTION RETRIEVES NON_OVERLAPPING ROI PIXELS ONLY
                %%%%%%%%%% ALL OVERLAPPING PIXELS ARE DROPPED IN THIS SIGNAL MATRIX

%%%%%%%%%%%%%%% LOADS RELEVANT VARIABLES FROM SEGMENTED DATA INTO WORKSPACE
    % INPUT
        % ROI 
    % OUTPUT
        % signal: m x n matrix containing Ca signal for all pixels lie
        % within ROI's designated by CalIm CNMF
            % Rows :  Pixels that fall within ROI's
            % Columns : samples for each timepoint sampled at 20 Hz
        % glomeruli: n x 1 array designating the ROI assignment for each
        % pixel row in the signal matrix
%%%%%%%%%%%%%%% 

% Reshape raw for indexing with ROI assignments
raw_=reshape(raw, [(size(raw, 1)*size(raw, 2)), size(raw, 3)]);

signal=[]; % Create empty matrix to store Raw downsampled (median filtered) Ca Imaging data
    
glomeruli = []; % Create an array to hold glomerular assigment for each row of signal matrix
for i = 1:length(ROI)
    temp_indx = ROI{i};
    glom_temp = (ones(1, length(ROI{i})) * i)';
    glomeruli= vertcat(glomeruli, glom_temp);
    for ii = 1:size(temp_indx, 1)
        signal((size(signal,1)+1),:) = raw_(temp_indx(ii),:);
    end
end
