function [A ROI ROIweights] = ROImasking(A_or)

%%%%%%%%%%%%%%% FINDS ROIS FROM A_OR TO MASK RAW SIGNAL

    % INPUT
        % A_or: Spatial matrix output from CalIm CNMF
        
    % OUTPUT
        % A : reshaped full spatial matrix
            % 116pixels x 156 pixels x Mglomeruli
        % glomind : dictionary with of glomerular indices
            % Each entry N x 1 with N indexed pixels for image frame
            % Indicies correspond to 116 x 156 pixel frame, not A_or
                % Meaning dont index into matrix containing all glomeruli
            % Shared pixels are dropped for equal weights across ROI
        % ROIweights: dictionary similar to ROI
            % Each entry N x 2
                % Rows: N indexed pixels for image frame
                % Column 1 : Pixel indexing into the image frame
                % Column 2 : Weight for pixel
            % This includes overlapping pixels since weighted accordingly
            
%%%%%%%%%%%%%%%

%% GRAB COORDINATES FOR ROI AND SURRONDING IMMEDIATE REGION FOR RAW DATA
% Glom1 = imagesc(A_orrr(1:13,90:100,1))
A_or_full = full(A_or);
% Reshape A_or so that the pixels regain the original tiff structure
    % A_or horizontally concatenates pixels within each frame
        % No longer maintains morphology of glomerular units
A = reshape(A_or_full, [116 156 size(A_or, 2)]);
    % 115pixels x 156pixels x Mglomeruli

% Now loop through each glomeruli and grab pixels with nonzero weight
for i = 1:size(A, 3) % LOOP GLOMERULI
    glomind{i} = find(A(:,:,i)~=0); % Grab all pixels that are weighted for each glomerulus
end

% Check for overlapping pixels between glomeruli
indx=1
for i = 1:size(A_or, 1)
    if length(find(A_or_full(i,:)))>1
        overlap(indx,1)=i
        temp_glom = find(A_or_full(i,:)~=0)
        overlap(indx,2:3) = temp_glom
        indx=indx+1;
    end
end
clear indx

%% FOR UNWEIGHTED ROI INDICES

% Drop overlapping pixels to get non-overlapping ROI indices
    % Incides refer to 115 x 156 pixel frame concatenated
        % Therefore concatenated vector = 18096 in length
drp = []; % Create vector for dropped pixels due to overlap
    % Allow for concatenation when looping 
        % Since don't know how large vector will be
for i = 1:length(glomind)
    temp_glom = glomind{i}; % Pull out temporary indices of glomeruli ROI
        % Indices correspond to a single downsampled frmae
            % 116 x 156 pixels concatenated
                % Single 18096 vector
    for ii = 1: size(overlap,1) % Loop over every pixel in overlap
        % See if overlapping pixel member of current looping glomeruli
        drpp = find(glomind{i}==overlap(ii,1));
        % If pixel exists that is member of current looping glomeruli:
        if drpp
            drp = horzcat(drp, drpp); % Add pixel to vector of pixels to drop
        end
    end
    temp_glom(drp)=[]; % drop any indices corresponding to overlapping pixels
    ROI{i} = temp_glom; % Store updated ROI mask for output
    drp = []; % Clear drp vector since will change shape each loop
end


%% FOR WEIGHT ROI INDICES

for i = 1:length(glomind) % LOOP GLOMERULI
    temp_spat = A(:,:,i); % Extract spatial matrix for frame for glomeruli
    ROIw(:,1) = glomind{i}; % Column 1 store pixel number index
        % Reference concatnated 116 x 156 downsampled frame
    ROIw(:,2)= temp_spat(glomind{i}); % Store weight of pixel
    ROIweights{i} = ROIw; % Store N x 2 vector as dictionary entry
    clear ROIw % Clear vector since will change size with each loop
end