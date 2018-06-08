%%%%%%%%%%%%%%% PLOT THE CA VERSE ETH DERIVATIVES USING THE PROPER DELAY
    % INPUT
        % C_trials : m x n x r
            % m : number of pixels across all ROI's
            % n : number of timepoints sampled in each trial
            % r : number of trials
        % Eth: Derivative of ethanol signal
            % Samples averaged across each frame
            % NOTE: This does not have to be a derivative
                % Just that raw ethanol signal is not relevant for analyses
        % % delay : array with index and value of max mean response for xcorr
            % m x n
                % m : number of frames away from 0 lag for max response
                % n : value of max mean cross correlation for ROI
            % NOTE: can use m x 1 array if m is number of ideal frame lag
        % ODOR_ON : indices indicatin range for plotting frame relation
            % m x 1
                % m: frame numbers relevant to analysis
            
    % OUTPUT
        % Figures: plot where x axis C_trials and y axis is Eth
            % Each dot represent one timepoint from session
            % All timepoints for all trials including in 1 ROI plot
%%%%%%%%%%%%%%% 

function plt_rltn(C_trials, Eth, delay, ODOR_ON)
for i = 1:size(Rig_trials,1) % LOOP ROIs
    y = squeeze(Rig_trials(i,min(ODOR_ON):max(ODOR_ON)-1,:)); % Grab ROI signal
    
    if delay(i, 1)>= 0
        warning( 'Warning, most likely error here rather than implied predictive coding lol')
    end
    y(1:abs(delay(i, 1)),:)=[]; % Drop beginn frames from Ca to lag signal
    x = Eth(:,min(ODOR_ON):max(ODOR_ON)-1)';
    x((length(x)-abs(delay(i, 1))+1):length(x),:)=[];
    x = reshape(x.',1,[])';
       
    
    y =  reshape(y.',1,[]);
    y = y';
    
    fitobject{i} = fit(x,y,'poly2');
    figure('Name',sprintf('Response curves of glomeruli %d to eth deriv', i) ,'NumberTitle','off')
    plot(fitobject{i},x,y);
    clear x cat y
end