%%%%%%%%%%%%%%% ROI DELAY CALCULATED VIA STA FEATURE DETECTION
    % INPUT
        % correlation_reponse:array with average xcorr response for eah ROI
            % m x n
                % m : number of ROI's
                % n : size of cross-correlation
        % wind : scalar for frames to pad pre and post event frame
            % i. e. wind = 40  would pad for 2 sec before and after @ 20Hz
            
    % OUTPUT
        % delay : array with index and value of max mean response for xcorr
            % m x n
                % m : number of frames away from 0 lag for max response
                % n : value of max mean cross correlation for ROI
%%%%%%%%%%%%%%% 
function delay = delay_STA(correlation_reponse, wind)
% Find length of samples in cross correlation
sz = round(size(corr_resp,2)./2);

for i = 1:size(correlation_reponse, 1) % LOOP ROIs
    [lag lag_in] = max(corr_resp1(i, sz-40:sz+40));
    delay(i, 1)= lag_in-40; % Lag for max corr
    delay(i, 2) = lag; % Amount of max corr
end