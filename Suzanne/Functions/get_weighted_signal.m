function [signal glomeruli] = get_weighted_signal(ROIweights, raw)

%%%%%%%%%%%%%%% LOADS RAW PIXEL SIGNAL FROM COMPRESSED FILE THAT LIES WITHIN ROIs
    % INPUT
        % ROIweights : cell structure where each cell is ROI as designated by CalIm CNMF
            % Within each cell: m x 2  array
                % m (rows) : all pixel associated with cell's ROI (glomeruli) 
                % col 1 : pixel index within the 116 x 152 concatenated pixel frame
                % col 2 : weight to give pixel when averaging glomeruli
    % OUTPUT
        % signal: m x n x r matrix containing Ca signal for all pixels lie within ROI's designated by CalIm CNMF
            % m : Pixels that fall within ROI's
            % n : samples for each timepoint sampled at 20 Hz
            % r : trials
        % glomeruli: n x 2 array 
            % Col 1 : designates the ROI (glomeruli) assignment for each pixel in signal
            % Col 2 : designates weighting to that pixel
%%%%%%%%%%%%%%% 

% Reshape raw for indexing with ROI assignments
raw_=reshape(raw, [(size(raw, 1)*size(raw, 2)), size(raw, 3), size(raw,4)]);

    
glomeruli = []; % Create an array to hold glomerular assigment for each row of signal matrix
% CURRENTLY THIS LOOP DOES NOT DROP NON ROI PIXELS SO FIX FOR THAT
for i = 1:length(ROIweights) % LOOP OVER GLOMERULI
   temp = ROIweights{i}; % Grab indices for ROI pixels of looping glomerul
    % Col 1 : Indices that correspond to a 116 x 156 pixel tif image frame
    glomeruli = vertcat(glomeruli, temp);
    tempw = repmat(temp(:,2), 1, size(raw_,2), size(raw_,3));
    if i==1
        signal(1:length(glomeruli),:,:) = raw_(temp(:,1),:,:) .* tempw ;
    else 
        signal(length(glomeruli)-length(temp)+1:length(glomeruli),:,:) = raw_(temp(:,1),:,:) .* tempw;
    end
    %    signal_(length(glomeruli_w)-length(temp):length(glomeruli_w),:,:)= vertcat(signal_w, raw_(temp(:,1),:,:))
    %             signal_w(temp,:, :) = raw_(temp,:,:);
end
