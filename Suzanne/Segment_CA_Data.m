%Code to segment the ROIs using CNMF

%This code can be run on any data that is in a format of frames over time.
%Current version imports a .mat file of subsampled data.

%Current settings work well for the new data sets. 
%Please refer to documentation included with the CNMF functions for more details. 


%% Set parameters and import calcium data matrix

%data import variables
mxf = 450; %maximum number of frames to import per trial
trials = [2:40]; %trials to examine. Be sure that the ethanol data and calcium data are correctly matched
ODOR_ON = [1:350]; %frames per trial to analyze 

% Import data
clear Y Yn Y1
     
   %  cd('D:\Imaging\10-23-2017_ROSA_AAV-retro-Ef1a-cre_110\Aligned');
   % Y1 = importdata('AAV-retro_110_10-23-2017_compressed.mat');

    cd('G:\ForMerav\Copressed')
    Y1 = importdata('3-30-2018-6-09 PM-Suzanne-Thy1-80.mat');
    
    %DATA = Y1.DATA;
    DATA = Y1.clc;
    trial = Y1.tr1a;
    
    Y = zeros(size(DATA,1),size(DATA,2),length(trials).*length(ODOR_ON));
    xc = 1;
    for y = trials;
        if y<7 && y>3;
            break
        else
        YT = DATA(:,:,find(trial==y));
        YT = YT(:,:,ODOR_ON);
        
        Y(:,:,xc:xc+length(ODOR_ON)-1) = YT; 
        
        tr(xc:xc+length(ODOR_ON)) = y;
        xc = xc+length(ODOR_ON);
        end
        
    end

    %Yn = Y1(:,:,ODOR_ON,trials);
    
   % Y = reshape(Yn,size(Yn,1),size(Yn,2),size(Yn,3).*size(Yn,4));
    
    Y = Y(5:end,5:end,:); %cut out the artifacts from alignment
    
    [d1,d2,T] = size(Y);                                % dimensions of dataset, assume T is last
    d = d1*d2;                                          % total number of pixels per T image

%% Set parameters for segmentation algorithm

K = 100;                                           % number of components to be found
tau = 2;                                          % std of gaussian kernel (size of neuron) 
p = 0;                                            % order of autoregressive system
                                                  % (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.95;                                  % merging threshold

options = CNMFSetParms(...                      
    'd1',d1,'d2',d2,...                         % dimensions of datasets
    'search_method','dilate','dist',3,...       % search locations when updating spatial components
    'deconv_method','constrained_foopsi',...    % activity deconvolution method
    'temporal_iter',2,...                       % number of block-coordinate descent steps 
    'fudge_factor',0.98,...                     % bias correction for AR coefficients
    'merge_thr',merge_thr,...                    % merging threshold
    'gSig',tau...
    );

%options = structmerge(defoptions,options,'ErrorIfNewField',1);

%options,

%% Data pre-processing

[P,Y] = preprocess_data(Y,p);

%% fast initialization of spatial components using greedyROI and HALS

[Ain,Cin,bin,fin,center] = initialize_components(Y,K,tau,options,P);  % initialize

% display centers of found components
Cn =  correlation_image(Y);
figure;
imagesc(Cn);
axis equal; axis tight; hold all;
scatter(center(:,2),center(:,1),'mo');
title('Center of ROIs found from initialization algorithm');
drawnow;

%% manually refine components (optional)
refine_components = true;  % flag for manual refinement
if refine_components
    [Ain,Cin,center] = manually_refine_components(Y,Ain,Cin,center,Cn,tau,options);
end

%% update spatial components
Yr = reshape(Y,d,T);
[A,b,Cin] = update_spatial_components(Yr,Cin,fin,[Ain,bin],P,options);

% update temporal components
P.p = 0;    % set AR temporarily to zero for speed
[C,f,P,S,YrA] = update_temporal_components(Yr,A,b,Cin,fin,P,options);

%% classify components
[ROIvars.rval_space,ROIvars.rval_time,ROIvars.max_pr,ROIvars.sizeA,keep] = classify_components(Y,A,C,b,f,YrA,options);
    
%% run GUI for modifying component selection (optional, close twice to save values)
run_GUI = true;
if run_GUI
    Coor = plot_contours(A,Cn,options,1); close;
    GUIout = ROI_GUI(A,options,Cn,Coor,keep,ROIvars);   
    options = GUIout{2};
    keep = GUIout{3};    
end

%% merge found components
[Am,Cm,K_m,merged_ROIs,Pm,Sm] = merge_components(Yr,A(:,keep),b,C(keep,:),f,P,S,options);

% refine estimates excluding rejected components

Pm.p = p;    % restore AR value
[A2,b2,C2] = update_spatial_components(Yr,Cm,f,[Am,b],Pm,options);
[C2,f2,P2,S2,YrA2] = update_temporal_components(Yr,A2,b2,C2,f,Pm,options);

%% do some plotting

[A_or,C_or,S_or,P_or] = order_ROIs(A2,C2,S2,P2); % order components
K_m = size(C_or,1);
[C_df,~] = extract_DF_F(Yr,A_or,C_or,P_or,options); % extract DF/F values (optional)

figure;
[Coor,json_file] = plot_contours(A_or,Cn,options,1); % contour plot of spatial footprints
%savejson('jmesh',json_file,'filename');        % optional save json file with component coordinates (requires matlab json library)

%% display components

plot_components_GUI(Yr,A_or,C_or,b2,f2,Cn,options)
fn = 'Thy100_4_2_2018_glom.mat';
%%
save(fn, '-v7.3')
%% make movie

%make_patch_video(A_or,C_or,b2,f2,Yr,Coor,options)

%% Output variables

%C_or the refined, merged ROI temporal data
%A_or the refined, merged ROI spatial data
