function [Glom_Sig resp] = ROI_sig(C_trials, glomeruli)
%%%%%%%%%%%%%%% AVERAGES RESPONSE WITHIN ROI
    % INPUT
        % Sig_trials: refers to output from DFFcal() function
            % m x n x r
                % m : number of pixels across all ROI's
                % n : number of timepoints sampled in each trial
                % r : number of trials
        % glomeruli: refers to output from get_signal() function 
            % n x 1 array OR n x 2 array
            % Col 1 : ROI assignment for each pixel row in the signal matrix
            % (optional) col 2 : from weighted glomeruli output
                % This array is dropped for the function
            
    % OUTPUT
        % Glom_Sig: m x n
            % m : number ROIs
            % n : timepoints for all trials concatenated horizontally
        % glomeruli: n x 1 array
            % ROI assignment for each pixel row in the signal matrix
%%%%%%%%%%%%%%% 

if size(glomeruli, 2) == 2
    glomeruli = glomeruli(:,1)
elseif size(glomeruli, 2) > 2
    error(' Glomeruli must have 2 columns or less, please see function description.')
end

%% Average over pixels in glomeruli
Sig_dff = reshape(C_trials, size(C_trials, 1), size(C_trials, 2)*size(C_trials, 3));
for i=1:max(unique(glomeruli)) %LOOP GLOMERULI
    df_temp = Sig_dff(find(glomeruli==i), :); % Grab all pixels correponding to looping ROI 
    Glom_Sig(i,:) = (mean(df_temp,1)); % Mean raw compressed signal of looping ROI
end

%% For glomeruli ranked by magnitude of response over session
% Find the max frame response of each glomeruli
dfresp_sig = max(Glom_Sig,[],2); 
% Sort ROI's based on their max response magnitude
[drs resp] = sort(dfresp_sig);
    % Here response gives the indices of the glomeruli ranked low-high max

