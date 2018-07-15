% DataSet1.xlsx
    % Rows: 14000*f where f is the number of FOVs
        % each block of 14000 rows represents the data from 1 FOV
    % Columns:
        % 1 : Frame number within session from 1:14000
        % 2 : Date: ie 3302018 = March 30, 2018
        % 3 : Thy1 Mouse number/name
        % 4 : Trial Number
        % 5 : Frame number within trial from 1:350
        % 6 : Ethanol signal from the given calcium frame number
        % 7 : Derivative of the Ethanol signal relevant to the proceeding
        % frame number, need to shift +1 if doing direct analysis
            % Since all trials are aligned by max xcorr of ROI, this
            % doesn't usually matter
        % 8 : Turb Label
            % Low = 0
            % Med = 1
            % High = 2
       % 9 - last column : designates 'glomeruli'/ROI number
% Delay_ds.xlsx
  % Rows: Glomeruli/ ROI data
  % Columns:
        % 1 : Date: ie 3302018 = March 30, 2018
        % 2 : Thy1 Mouse number/name
        % 3 : Glomeruli/ ROI number
        % 4 : Number of frames to delay for max corr
        % 5 : Magnitude of max xcorr
%% LOAD DATASETS

ds = load(uigetfile('.mat'));
data = double(ds.st_dataset);

ds_delay = load(uigetfile('.mat'));
data_delay = double(ds_delay.st_dataset);

%% STA

% Create vector with boundaries at beginning and end of each trial
    % This includes boundary samples within 20 samples of trial beginning
    % or end
trials = vertcat(0, find(diff(data(:,4))));  % Trial indices
boundaries = [];
for i = 1:41; % Only loop first 40 trial indices since create temporary data matrix for each mouse
    boundaries = vertcat(boundaries, (trials(i)-20:trials(i)+20)');
end
boundaries(boundaries<=0 | boundaries>14000)=[];

% Field of View indicies and recording date
FOVi = vertcat(0, find(diff(data(:,2)))); %Create vector with indices corresponding to the start of each FOV
FOVs = unique(data(:, 2)); % Dates for each FOV

% Create two STA matrices, one to hold STA related data and the other the
% STAs
sta = []; % To hold all information necessary to calc STA
STA = []; % To hold actual STA

for i=1:length(FOVs) % Loop through field of views (FOVs)
    % Get ethanol signal
    e_temp = data(find(data(:,2)==FOVs(i)), 6);
    % Take deriv and repeat to easily index eth to match glom signal matrix
        % It is okay that taking the derivative of all trials concatenated
        % because ignoring STA's from within 20 samples of beginning or end
        % of a trial
    e_te = diff(reshape(e_temp, 40, 350)');
    df_e = e_te(:);
    
    % Get glomerular signal
    g_temp = data(find(data(:,2)==FOVs(i)) , 9:size(data, 2)); % Grab all rows corresponding to looping FOV 
    if sum(all(~isnan(g_temp)))>0 % Unless all columns contain signal
    g_temp = g_temp(:,all(~isnan(g_temp))); % Drop Columns with all Nan
    end
    
    % Get baseline
    g_trial = reshape(g_temp, 350, 40, size(g_temp, 2));
    bs= mean(g_trial);
    df = diff(g_trial-bs);
    % thrsh = 2*std(g_trial);
    df_g = reshape(df, size(df, 1)*size(df, 2), size(g_temp,2));
    threshold = reshape(df./std(df), size(df_g, 1), size(g_temp,2));
  
    % Find events
    events = find(abs(threshold)>=2);
    events = setdiff(events, boundaries); % Remove events associated with boundary events
    
    % Grab eth corresponding to events 
    gl_events = ceil(events/size(df_g, 1)); % Vector denoting to which glom event belong
    gl = unique(gl_events); % Vector of glomeruli with calcium events in case a glom has no events
    for ii=1:length(gl) % LOOP THROUGH GLOM
        delay = data_delay(data_delay(:, 1)==FOVs(i) & data_delay(:, 3)==ii, 4);
        if ii==1
            t=events(gl_events==ii);
        else
            t=events(gl_events==ii)-14000*(gl(ii-1));
        end
        tt(1:length(t), 1) = FOVs(i); % First column is the FOV identified by date
        tt(:, 2) = gl(ii);
        for iii= 1: length(t)-1; % LOOP THROUGH GLOM EVENTS
            if t(iii)+delay>20 && t(iii)+delay + 20<length(df_e) % Make sure still can pull 20 windows with delay
                % This will toss events before odor come on anyways since
                % they are within 1 sec of the beginning of the trial
            tt(iii, 3:43) = df_e(t(iii)-20+delay:t(iii)+20+delay);
            end
        end
        sta = vertcat(sta, tt);
        STA = vertcat(STA, mean(tt));
        clear tt
    end
    
    % Visualize distribution of significant events across glomeruli
    figure, imagesc(abs(threshold)>=2);
    
end
figure('Name', 'STA +/- 1 second')
imagesc(STA(:, 3:43));

STA_=STA(:, 3:43);
for i= 1:size(STA_, 1)
    STA_(i, :) = smooth(STA_(i, :));
end
figure('Name', 'Smoothed STA +/- 1 second')
imagesc(STA_, [-8 8]);