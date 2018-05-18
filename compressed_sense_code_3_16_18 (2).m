%%
%compressed sensing for glomeruli in turbulence

% First, run Segment_CA_Data
%C_or is the temporal data from each ROI
%A_or is the spatial data from each ROI

[dt1 dd] = uigetfile('*.mat');
cd(dd);
Y1 = importdata(dt1);

% C_or = Y1.DATA;
 C_or = Y1.C_or;
A_or = Y1.A_or;

trials = Y1.trials;
ODOR = [1:350]; %time range during which odor is present
    rn = [1:length(ODOR)]; %variable for internal processing
 
 MPH = 600; %minimum peak height for detected peaks in wheel data
MPP = 200; %minimum peak prominence for detected peaks in wheel data

 
%% Import the ethanol data

%Import the data

%plume data
mxf = 375;
%cd('D:\Gire Lab Data\matlab tests\Imaging')
%cd('G:\Imaging Data\7-25-2017')
%cd('D:\Imaging\10-23-2017_ROSA_AAV-retro-Ef1a-cre_110')
% mxf =350;

cd('G:\ForMerav\Ethanol')

[fn pn] = uigetfile('*.dat');

cd(pn)

fileID = fopen(fn);

data = fread(fileID,'int32','ieee-be');

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

[whl_max whl_max1] = findpeaks(WHL,'MinPeakHeight',MPH,'MinPeakProminence',MPP);

figure, plot(WHL)
hold on
plot(whl_max1,whl_max,'r.','markersize',24);

df = diff(whl_max1);
df(df<5) = 5;
df(df>5000) = 5000;

dff1 = 1./df;

WHL1 = WHL.*0;
for y = 1:length(df)-1;
    WHL1(whl_max1(y):whl_max1(y+1)) = dff1(y);
end

for p = 1:length(un);
    frr = FR(TR==un(p));
    frr = frr-frr(1)+1;
    eth = ETH(TR==un(p));
    ts = TS(TR==un(p));
    ds = DST(TR==un(p));
    whll = WHL1(TR==un(p));
    
    %interpret wheel data
    
    
    
    unf = unique(frr);
    for t = 1:length(unf);%min([mxf+25 max(unf)]);%length(unf);
        fr(p,t) = t;%mean(frr(frr==unf(t)));
        data2(p,t) = mean(eth(frr==unf(t)));
        ts_fr(p,t) = mean(ts(frr==unf(t)));
        dst(p,t) = mean(ds(frr==unf(t)));
        whl(p,t) = mean(whll(frr==unf(t)));
        tr(p,t) = mean(ds(frr==unf(t))).*0+un(p); 
    end
    data2a = data2(p,min(find(ts_fr(p,:)>1)):end);
    data2b(p,:) = data2(p,1:mxf);
    
        run2a = whl(p,min(find(ts_fr(p,:)>1)):end);
        run2b(p,:) = whl(p,1:mxf);
    
        t2 = ts_fr(p,min(find(ts_fr(p,:)>1)):end);
    t2b(p,:) = t2(1:mxf);
end

clear data2 ts_fr
fls = [1:40];%files used for imaging; 2-20 for 7-21-2017_thy1-35; 2-24 for thy1-27; 1-10 for retro-AAV
data2 = data2b(fls,:);
dst = dst(fls,:);

run2 = run2b(fls,:)./max(max(run2b)); %normalized running speed


data2 = data2-repmat(mean(data2(:,1:100),2),1,size(data2,2));
ts_fr = t2b;
figure, plot(run2(1,:))
%%
% t = 1:50;
% B = (exp(-t./24));
% clear data2n
% for y = 1:size(data2,1);
% [Q,R] = deconv(data2(y,:),B);
% data2n(y,:) = Q;
% end
% 
 figure, 
 subplot(2,1,1)
 imagesc(diff(data2,[],2));
 xlabel('Ethanol Sensor Data')
 subplot(2,1,2)
 imagesc(run2);
 xlabel('Normalized Running Speed')
% 
 drawnow
 
% TR_MEAN is glomeruli
%PL_MEAN is plume
eth_data = data2(trials,ODOR);

PL_MEAN = eth_data';%reshape(eth_data',size(eth_data,1).*size(eth_data,2),1)';

C_trial = reshape(C_or,size(C_or,1),size(C_or,2)/length(trials),length(trials));
% C_trial = reshape(C_or,size(C_or,1),size(C_or,2)/length(trials),length(trials));
TR_MEAN = C_trial(:,:,:);

%CM = Y;
 
TR_MEAN1 = TR_MEAN;
 
%% Basic response analysis

baseline = mean(TR_MEAN1(:,1:50,:),2); %baseline fluorescence on each trial before odor onset

baseline(baseline<5) = 5;

TR_MEAN_df = TR_MEAN1-baseline;

TR_MEAN_df = TR_MEAN_df./baseline;

% TR_MEAN_df(TR_MEAN_df>5) = 5;

resp_mean = mean(TR_MEAN_df,3); %mean of each cell's response over all trials

figure, imagesc(resp_mean,[-1 1]);

%order cells based on strength of response
resp_str = max(resp_mean(:, ODOR, :),[],2);
[rs rs1] = sort(resp_str);

%RESP_z_score = mean(TR_MEAN,3)./std(TR_MEAN,[],3);



%% Correlation to plume analysis

thr = 1;

tm_rn = 20; %range in frames before and after

mx = size(PL_MEAN,1);

xc = 0;
clear AVG gl_id AVG_gl crc tr_id
for y = 1:size(TR_MEAN1,1);
    crr = squeeze(TR_MEAN1(y,:,:));
    for y1 = 1:size(crr,2);
        clear f2 crr1 zs df1 f1
        crr1 = squeeze(crr(:,y1));
        
        zs = crr1-mean(crr1(:,:));
        zs = zs./std(zs(:));
        zs = diff(zs);
        
        
        f1 = find(zs>thr);
        df1 = diff(f1);
        f2 = f1(df1>5);
        
        crc(y,:,y1) = xcorr(diff(PL_MEAN(:,y1)),diff(squeeze(TR_MEAN(y,:,y1))),'coeff');
        
        for u = 1:length(f2);
            
            if f2(u) > tm_rn & f2(u) < mx-tm_rn; % Drop first and last 20 frames to remove window artifacts or xcorr
                xc = xc+1;
            AVG(xc,:) = squeeze(diff(PL_MEAN(f2(u)-tm_rn:f2(u)+tm_rn,y1)))';
            AVG_gl(xc,:) = zs(f2(u)-tm_rn:f2(u)+tm_rn);
            gl_id(xc) = y;
            tr_id(xc) = y1;
            
            end
            
        end
        
        %df = diff(crr1);
    
    end
    
end

%% Plot data from correlation and response analysis
resp_mean1 = resp_mean(rs1,:);

figure, subplot(2,1,1), 
imagesc(resp_mean(rs1,:),[-0.5 1]);

corr_resp = mean(crc(rs1,:,:),3)./std(crc(rs1,:,:),[],3);
corr_resp1 = corr_resp;

corr_resp(abs(corr_resp)<1) = 0;
cr1 = max(corr_resp,[],2);

sz = round(size(corr_resp,2)./2);
subplot(2,1,2),
imagesc(corr_resp(:,sz-40:sz+40));

figure, subplot(2,1,2), imagesc(corr_resp1(cr1>0,sz-40:sz+40))
subplot(2,1,1), imagesc(resp_mean1(cr1>0,:))

LG = ([sz-40:sz+40]-mean([sz-40:sz+40]))./20;

figure, plot(LG,corr_resp1(cr1>0,sz-40:sz+40)','linewidth',2)
hold on
plot([0 0],[-2 2],'k--')



%%
%bootstrapped average
clear AVG1 gl_id1 CRC1 tr_id1
for av = 1:10;
xc1 = 0;
clear AVG11 gl_id11 crc1
for y = 1:size(TR_MEAN1,1);
    crr = squeeze(TR_MEAN1(y,:,:));
    for y1 = 1:size(crr,2);
        clear f2 crr1 zs df1 f1
        crr1 = squeeze(crr(:,y1));
        
        zs = crr1-mean(crr1(:,:));
         zs = diff(zs);
         
        zs = zs./std(zs(:));
        
        f1 = find(zs>thr);
        df1 = diff(f1);
        f2 = f1(df1>5);
        
        rp = randperm(size(crr,2));
        rp = rp(rp~=y1);
        
        crc1(y,:,y1) = xcorr(diff(PL_MEAN(:,rp(1))),diff(squeeze(TR_MEAN(y,:,y1))),'coeff');
        
        for u = 1:length(f2);
            rp = randperm(size(crr,2));
            rp = rp(rp~=y1);
            if f2(u)>tm_rn & f2(u)<mx-tm_rn ;
                xc1 = xc1+1;
            AVG11(xc1,:) = squeeze(diff(PL_MEAN(f2(u)-tm_rn:f2(u)+tm_rn,rp(1))))';
            gl_id11(xc1) = y;
            tr_id1(xc1) = y1;
            end
        end
        %df = diff(crr1);
    end
end

AVG1(av,:,:) = AVG11;
gl_id1(:,av) = gl_id11;

CRC1(av,:,:,:) = crc1;

end

%%

dt = 0.05; %delta time between frames

clear avgg vg vg1 avgg1 p avgg_gl c_avg CVgg1
xc1 = 0;
for y = 1:size(TR_MEAN1,1);
    
   % plot(AVG(gl_id==y,:)')
    avgg(y,:) = mean(AVG(gl_id==y,:));
    
    avgg_gl(y,:) = mean(AVG_gl(gl_id==y,:));
    
    avgg1(y,:) = mean(mean(AVG1(:,gl_id1(:,1)==y,:),1));
    
    
    
    vg(y,:) = std(AVG(gl_id==y,:))/sqrt(length(find(gl_id==y)));
    
    vg1(y,:) = mean(std(AVG1(:,gl_id1(:,1)==y,:),[],1))/sqrt(200);
    
    AVGT = AVG1(:,gl_id1(:,1)==y,:);
    AVGT = reshape(AVGT,size(AVGT,1)*size(AVGT,2),size(AVGT,3));
    
    for y1 = 1:size(AVG,2);
    p(y,y1) = ranksum(AVG(gl_id==y,y1),AVGT(:,y1));
    end
    
    if min(p(y,:))<0.05/size(C_or,1) & max(abs(avgg(y,:)-avgg1(y,:)))>thr;
        xc1 = xc1+1;
        c_avg(xc1,:) = avgg(y,:)-avgg1(y,:);
        indx(xc1) = y;
    end
    
if max(abs(avgg(y,:)-avgg1(y,:)))>1.5;
figure
x = ([1:size(avgg,2)]-tm_rn).*dt;
subplot(2,1,1)
    errorbar(x,avgg(y,:),vg(y,:));
hold on
     errorbar(x,avgg1(y,:),vg1(y,:));
     subplot(2,1,2)
         hold on
    plot([min(x) max(x)],[0 0],'k--')
     errorbar(x,avgg(y,:)-avgg1(y,:),vg(y,:));
     axis([min(x) max(x) -3 20])
     
     x = ([1:size(avgg_gl,2)]-tm_rn).*dt;
     plot(x,avgg_gl(y,:))
end
end

p1 = p.*0;
p1(p<0.05/40) = 1;

figure, imagesc(1-p,[1-0.05/40 1])

clear xc
for y1 = 1:size(PL_MEAN,2);
xc(:,y1) = xcorr(diff(PL_MEAN(:,y1)),'coeff');
end

figure, imagesc(xc(rn,:)','x',([1:length(rn)]-ceil(length(rn)/2)).*dt)

figure,plot(([1:length(xc)]-length(xc)/2).*dt,xc,'linewidth',2)

axis([-2 2 -0.5 1])

figure,
x = ([1:size(c_avg,2)]-tm_rn).*dt;
hold on

plot(x,c_avg,'linewidth',2);
plot([min(x) max(x)],[0 0],'k--')

%%
clear fnd
for y = 1:size(c_avg,1);
fnd(y) = find(c_avg(y,:)==max(c_avg(y,:)));
end

[pt1 pt2] = sort(fnd);

figure, imagesc(c_avg(pt2,:))

%%
clear ROI RROI
%RROI = Y(:,:,1);
Y = Y1.Y;
for y = 1:length(fnd);
    ROI = A_or(:,indx(y));
   % ROI = reshape(ROI,size(Y,1),size(Y,2));
   if y==1;
       RROI = ROI;
       RROI(ROI~=0) = fnd(y);
   else
    RROI(ROI~=0) = fnd(y);
   end
   
end
    RROI = reshape(RROI,size(Y,1),size(Y,2));
    
    figure, imagesc(RROI,[min(fnd)-1 max(fnd)+1])

