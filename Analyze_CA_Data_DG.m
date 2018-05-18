%% an m-file to combine CNMF-segmented imaging data with ethanol sensor data

% First, run Segment_CA_Data
%C_or is the temporal data from each ROI
%A_or is the spatial data from each ROI

[fg fdg] = uigetfile('.mat');
cd(fdg);
g_data = load(fg);

% Extract relevant Glom fields
trials = g_data.trials;
C_or = g_data.C_or;
y = g_data.y;
ODOR_ON = g_data.ODOR_ON;

%% Import the ethanol data

%Import the data

%plume data

%cd('D:\Gire Lab Data\matlab tests\Imaging')
%cd('G:\Imaging Data\7-25-2017')
%cd('G:\Imaging Data\10-23-2017_ROSA_AAV-retro-Ef1a-cre_110')

cd('G:\ForMerav\Ethanol')

[fn pn] = uigetfile('*.dat');

cd(pn)

fileID = fopen(fn);

data = fread(fileID,'int32','ieee-be');

mxf=350;

n_ten = find(data==-10);
TR = data(n_ten+1);
TS = data(n_ten+2);
FR = data(n_ten+3);
ETH = data(n_ten+4);
WHL = data(n_ten+5);
DST = data(n_ten+7);
SNF = data(n_ten+9);
%figure, plot(data1)

un = unique(TR);
FR = FR-FR(1)+1;
clear fr eth_fr ts_fr dst snf unf
SNF = medfilt2(SNF,[3 1]);
for p = 1:length(un);
    frr = FR(TR==un(p));
    frr = frr-frr(1)+1;
    eth = ETH(TR==un(p));
    ts = TS(TR==un(p));
    ds = DST(TR==un(p));
    snnf = WHL(TR==un(p));
    unf = unique(frr);
    for t = 1:length(unf);min([mxf+25 max(unf)]);length(unf);
        fr(p,t) = t;mean(frr(frr==unf(t)));
        data2(p,t) = mean(eth(frr==unf(t)));
        ts_fr(p,t) = mean(ts(frr==unf(t)));
        dst(p,t) = mean(ds(frr==unf(t)));
        snf(p,t) = mean(snnf(frr==unf(t)));
        tr(p,t) = mean(ds(frr==unf(t))).*0+un(p); 
    end
    data2a = data2(p,min(find(ts_fr(p,:)>1)):end);
    data2b(p,:) = data2(p,1:mxf);
    
        run2a = snf(p,min(find(ts_fr(p,:)>1)):end);
        run2b(p,:) = snf(p,1:mxf);
    
        t2 = ts_fr(p,min(find(ts_fr(p,:)>1)):end);
    t2b(p,:) = t2(1:mxf);
end

clear data2 ts_fr
fls = [1:40];%files used for imaging; 2-20 for 7-21-2017_thy1-35; 2-24 for thy1-27; 1-10 for retro-AAV
data2 = data2b(fls,:);
dst = dst(fls,:);

run2 = run2b(fls,:);


data2 = data2-repmat(mean(data2(:,1:100),2),1,size(data2,2));
ts_fr = t2b;

%%
t = 1:50;
B = (exp(-t./24));
clear data2n
for y = 1:size(data2,1);
[Q,R] = deconv(data2(y,:),B);
data2n(y,:) = Q;
end

figure, imagesc(data2n);

drawnow

%% Select the correct range of the ethanol data to compare to the segmented CA data



eth_data = data2(trials,ODOR_ON);

eth_data1 = reshape(eth_data',size(eth_data,1).*size(eth_data,2),1)';

eth_data = diff(eth_data,[],2);


figure, plot(eth_data1)

clear corr_mat
for y = 1:size(C_or,1);
    c1 = corrcoef(C_or(y,:),eth_data1);
    corr_mat(y) = c1(3);
end

figure, bar(corr_mat);

%% Trial-by-trial analysis of correlation between plume and activity

C_trial = reshape(C_or,size(C_or,1),size(C_or,2)/length(trials),length(trials));

C_trial = diff(C_trial,[],2);

clear corr_mat
for y1 = 1:length(trials);
for y = 1:size(C_or,1);
    c1 = corrcoef(squeeze(C_trial(y,180:300,y1)),eth_data(y1,180:300)');
    corr_mat(y,y1) = c1(3);
end

end

figure, imagesc(corr_mat);


%% Quick Dimensionality analysis trial-by-trial
clear EXP EXP1 spd1
for y = 1:length(trials);
    CA = C_trial(:,1:100,y);
    
    CA = CA-mean(CA')';
    CA = CA./std(CA,[],2);
    
    
    CO = C_trial(:,180:300,y);
    
    CO = CO-mean(CO')';
    CO = CO./std(CO,[],2);
    
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(CA);
    if y>1;
    EXP(y,1:size(EXPLAINED,1)) = EXPLAINED;
    else
    EXP(y,:) = EXPLAINED;
    end
    
    [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(CO);
    if y>1;
    EXP1(y,1:size(EXPLAINED,1)) = EXPLAINED;
    else
    EXP1(y,:) = EXPLAINED;
    end
    
    DEV_eth(y) = std(eth_data(y,180:300));
    
    run21 = fft(run2(y,:));
    run21([1:5 end-4:end]) = 0;
    run21 = real(ifft(run21));
    
    RUN(y,:) = run21;
    
    SPD(y) = std(RUN(y,1:100));%sum(abs(spd(50:60)))./sum(abs(spd));
    
    SPD1(y) = std(RUN(y,180:300));
end


[p1 p2] = sort(SPD1);

figure, imagesc(EXP1(p2,:))
figure, imagesc(RUN(p2,180:300));

figure, plot(SPD1,sum(EXP1(:,1:20),2),'bo')

for y = 1:size(EXP1,1);
CS = cumsum(EXP1(y,:));
CS1(y) = min(find(CS>90));
end

figure, plot(SPD1,CS1,'bo')

figure, plot(mean(EXP1(p2(1:10),:))-mean(EXP1(p2(end-10:end),:)))

figure, plot(SPD1,mean(abs(corr_mat)),'bo');

%%
g1 = g_data.A_or;
figure,

for y = 1:23;
    
    imagesc(reshape(g1(:,y),116,156))
    axis equal
    axis off
    drawnow
    pause(0.5)
end

