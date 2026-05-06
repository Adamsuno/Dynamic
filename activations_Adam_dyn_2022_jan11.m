%author Netphyslab. 1. Javeria Hashmi 2020 



%% 
close all
clear all



addpath /Users/jhash1/Dropbox/macros_oa
addpath /Users/jhash1/Dropbox/macros_oa/healthy_network_analysis_gs/
addpath /Users/jhash1/Dropbox/macros_oa/System_segregation_and_Network_optimization_scripts
addpath /Users/jhash1/Dropbox/macros_oa/2019_03_03_BCT
addpath /Users/jhash1/Dropbox/scripts/BCT/
addpath /Users/jhash1/Dropbox/GAC_study/gac_timeseries/


%filebase1 = '/Users/jhash1/Dropbox/CBP_study/TIMESERIES_CBP/';
%text1 = textread(['/Users/jhash1/Dropbox/CBP_study/TIMESERIES_CBP/macros/CBPsubjectlist.txt'],'%s');
%filebase1 = '/Users/jhash1/Dropbox/CBP_study/CBPts_163/';
%text1 = textread(['/Users/jhash1/Dropbox/CBP_study/CBPts_163/subjectlist_all_2021.txt'],'%s');

%labels=textread (['/Users/jhash1/Dropbox/GAC_study/gac_timeseries/roilist_pag.txt'],'%s');
labels=textread (['/Users/jhash1/Dropbox/CBP_study/CBPts_163_suspect/labels_Hashmi163.txt'],'%s');

%for Adam's resting state data
%path to the time series: 
filebase1 = '/Users/jhash1/Dropbox/CBP_study/dyn_ts/ts_all';
%path to the text file containing names of subjects: 
text1 = textread(['/Users/jhash1/Dropbox/CBP_study/dyn_ts/subject_list_dyn.txt'],'%s');
 %load /Users/jhash1/Dropbox/CBP_study/dyn_ts/subject_list_dyn.txt
nSubj = length(text1); % number of subjects



%%%%%%%filter for bold%%%%%%%%%%%%
% Wn =[0.25];
% TR=0.95;
% fs = 1/TR;
% [b,a] = butter(4,Wn/(0.5*fs),'low');
Wn = [0.07];
[b,a] = butter(4, Wn/0.2, 'low');

%low pass = 0.1662 based on original calc. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cbp 4, 13, 16, 24, 34, 37, 41,44 are missing task scans

%%
%scan1

load /Users/jhash1/Dropbox/CBP_study/dyn_ts/timingfile_6TRdelay/Scan1_rating_ts_dyn_AS.txt
rsts1=Scan1_rating_ts_dyn_AS;

for  i=1:nSubj;
    
    i
    
    %%%%%%%%%%%%%%%%%%%
     
    b1_all = load([filebase1 '/' text1{i} '_task1_Hashmi163_ts_gsr.txt']);
    
    
    
    %  ba_all=load([filebase1 '/CBP_amyg_ts/' text1{i} '_Broad_task1.txt']);
      
    %  b1_all=horzcat(b3_all,ba_all);
      
    sb = size(b1_all); lb = sb(1,1); nb = sb(1,2); b2_all = zeros(lb,nb);
    
        
        %convert to percent change
    for nodes = 1:nb;
        b1 = b1_all(:,nodes);
        b2 = 100*(b1 - mean(b1))/mean(b1); %bf = filtfilt(b,a,b2);
        b2_all(:,nodes)=b2; %change to bf if using filter
    end
    
    %allts_scan1_pag(i,:,:)=b2_all;
   
    
% null=zeros(15,163);
% 
% b2_all=vertcat(null,b2_all);
% 
% null1=zeros(15,1);
% rsts1a=vertcat(null1,rsts1);
events_1=find(rsts1>=1); %unh


events=length(events_1)
for n=1:events; %erps2) %2:length(events_10)
    
    %all_roi_ts_pc_ep(i,1:length(startBlock(i):startBlock(i+1)-1),:)=all_roi_ts_pc(startBlock(i):startBlock(i+1)-1,:);
    % all_roi_ts_pc_ep(i,:,:)=all_roi_ts_pc(events_4(i)-8:events_4(i)+6,:);
      all_roi_ts_pc_ep1(:,:,n)=b2_all(events_1(n) -28:events_1(n)+10,:); %adjust timing based on scan it keeps flat line for missed event ignore in next step
%size(all_roi_ts_pc_ep)
end
allcuescan1(:,:,:,i)=all_roi_ts_pc_ep1;

end


%%


load /Users/jhash1/Dropbox/CBP_study/dyn_ts/timingfile_6TRdelay/Scan2_rating_ts_dyn_AS.txt
rsts2=Scan2_rating_ts_dyn_AS;


for  i=1:nSubj;
    
    i
    
    %%%%%%%%%%%%%%%%%%%
    
 %b1_all = load([filebase1 'task2/' text1{i} '_new.txt']);
   
 b1_all = load([filebase1 '/' text1{i} '_task2_Hashmi163_ts_gsr.txt']);
    
      
   sb = size(b1_all); lb = sb(1,1); nb = sb(1,2); b2_all = zeros(lb,nb);
    
    for nodes = 1:nb;
        b1 = b1_all(:,nodes);
        b2 = 100*(b1 - mean(b1))/mean(b1); % bf = filtfilt(b,a,b2);
        b2_all(:,nodes)=b2;
    end
    

%    null=zeros(15,163);
% 
% b2_all=vertcat(null,b2_all);
% 
% null1=zeros(15,1);
% rsts2a=vertcat(null1,rsts2);
events_1=find(rsts2>=1); %unh


events=length(events_1)
for n=1:events; %erps2) %2:length(events_10)
    
     all_roi_ts_pc_ep2(:,:,n)=b2_all(events_1(n) -28:events_1(n)+10,:); %adjust timing based on scan it keeps flat line for missed event ignore in next step
%size(all_roi_ts_pc_ep)
end
allcuescan2(:,:,:,i)=all_roi_ts_pc_ep2;
end

 
 
%%
 



load /Users/jhash1/Dropbox/CBP_study/dyn_ts/timingfile_6TRdelay/HiLoUnknown_rating_ts_dyn_AS.txt
rsts3=HiLoUnknown_rating_ts_dyn_AS;%rsts3(426:end)=[]
%text2=text1([1:16,18:end]);
for  i=1:nSubj;
    
    i
    b1_all = load([filebase1 '/' text1{i} '_hilo_Hashmi163_ts_gsr.txt']);
      % b1_all = load([filebase1 'task3/' text1{i} '_new.txt']);
    
    
        sb = size(b1_all); lb = sb(1,1); nb = sb(1,2); b2_all = zeros(lb,nb);
    
    for nodes = 1:nb;
        b1 = b1_all(:,nodes);
        b2 = 100*(b1 - mean(b1))/mean(b1); 
        %bf = filtfilt(b,a,b2);
        b2_all(:,nodes)=b2;
    end
 
 null=zeros(15,163);
% 
 b2_all=vertcat(b2_all,null);
% 
 null1=zeros(15,1);
 rsts3a=vertcat(rsts3,null1);
events_1=find(rsts3a>=1); %unh




events=length(events_1)
for n=1:events; %erps2) %2:length(events_10)
    
    %all_roi_ts_pc_ep(i,1:length(startBlock(i):startBlock(i+1)-1),:)=all_roi_ts_pc(startBlock(i):startBlock(i+1)-1,:);
    % all_roi_ts_pc_ep(i,:,:)=all_roi_ts_pc(events_4(i)-8:events_4(i)+6,:);
      all_roi_ts_pc_ep3(:,:,n)=b2_all(events_1(n) -28:events_1(n)+10,:); %adjust timing based on scan it keeps flat line for missed event ignore in next step

      %size(all_roi_ts_pc_ep)
end
allcuescan3(:,:,:,i)=all_roi_ts_pc_ep3;

end



%%
 


%load /Users/jhash1/Dropbox/CBP_study/CBPts_163/scan3_rsts_rating_jv.txt;
load /Users/jhash1/Dropbox/CBP_study/dyn_ts/timingfile_6TRdelay/Training_rating_ts_dyn_AS.txt
rsts4=Training_rating_ts_dyn_AS;
%text2=text1([1:16,18:end]);
for  i=1:nSubj;
    
    i
    b1_all = load([filebase1 '/' text1{i} '_training_Hashmi163_ts_gsr.txt']);
      % b1_all = load([filebase1 'task3/' text1{i} '_new.txt']);
  
  
  
        sb = size(b1_all); lb = sb(1,1); nb = sb(1,2); b2_all = zeros(lb,nb);
     
    for nodes = 1:nb;
        b1 = b1_all(:,nodes);
        b2 = 100*(b1 - mean(b1))/mean(b1); 
        %bf = filtfilt(b,a,b2);
    
        b2_all(:,nodes)=b2;
    end
    
% null=zeros(15,163);
% 
% b2_all=vertcat(null,b2_all);
% 
% null1=zeros(15,1);
% rsts4a=vertcat(null1,rsts4);
events_1=find(rsts4>=1); %unh

events=length(events_1)
for n=1:events; %erps2) %2:length(events_10)
    
    %all_roi_ts_pc_ep(i,1:length(startBlock(i):startBlock(i+1)-1),:)=all_roi_ts_pc(startBlock(i):startBlock(i+1)-1,:);
    % all_roi_ts_pc_ep(i,:,:)=all_roi_ts_pc(events_4(i)-8:events_4(i)+6,:);
      all_roi_ts_pc_ep4(:,:,n)=b2_all(events_1(n) -28:events_1(n)+10,:); %adjust timing based on scan it keeps flat line for missed event ignore in next step
%size(all_roi_ts_pc_ep)
end
allcuescan4(:,:,:,i)=all_roi_ts_pc_ep4;

end

%%
%allcuescan1 for examoples contains 39 time points, 163 brain region,
%number of epochs, and 30 subjects
%extract conditions


%load timing.mat
%for averaging erps for type of events
%UNKNOWN HI LOW

%unknown low static:1,4,10
%unknown high static :2,7,9
%unknown low dynamic: 5,8,12
%unknown high dynamic:3,6,11

%extract epochs for condition
high3stat=allcuescan3(:,:,[2,7,9],:);%sq2

%average the extracted epochs:
high3mstat=squeeze(mean(high3stat,3));

%average across subjects ONLY for plotting purposes
highstatm=squeeze(mean(high3mstat,3));

high3dyn=allcuescan3(:,:,[3,6,11],:);%sq2
high3mdyn=squeeze(mean(high3dyn,3));
highdynm=squeeze(mean(high3mdyn,3));

low3stat=allcuescan3(:,:,[1,4,10],:);%sq2
low3mstat=squeeze(mean(low3stat,3));
lowstatm=squeeze(mean(low3mstat,3));

low3dyn=allcuescan3(:,:,[5,8,12],:);%sq2
low3mdyn=squeeze(mean(low3dyn,3));
lowdynm=squeeze(mean(low3mdyn,3));

high3stat=allcuescan4(:,:,[3,6,11],:);%sq2
train_high3mstat=squeeze(mean(high3stat,3));
train_highstatm=squeeze(mean(train_high3mstat,3));

high3dyn=allcuescan4(:,:,[2,7,9],:);%sq2
train_high3mdyn=squeeze(mean(high3dyn,3));
train_highdynm=squeeze(mean(train_high3mdyn,3));

low3stat=allcuescan4(:,:,[5,8,12],:);%sq2
train_low3mstat=squeeze(mean(low3stat,3));
train_lowstatm=squeeze(mean(train_low3mstat,3));

low3dyn=allcuescan4(:,:,[4,10],:);%sq2
train_low3mdyn=squeeze(mean(low3dyn,3));
train_lowdynm=squeeze(mean(train_low3mdyn,3));

low1=allcuescan1(:,:,[2,6,13],:);%sq2[3,7,8,10,13]
low1m=squeeze(mean(low1,3));

low2=allcuescan2(:,:,[6,14],:);%sq2[2,4,8,10,13]
low2m=squeeze(mean(low2,3));
test_lowcues_highdynamic=(low1m+low2m)./2;
test_lowcues_highdynamicm=squeeze(mean(test_lowcues_highdynamic,3));

low1=allcuescan1(:,:,[4,11,16],:);%sq2[3,7,8,10,13]
low1m=squeeze(mean(low1,3));

low2=allcuescan2(:,:,[11,15],:);%sq2[2,4,8,10,13]
low2m=squeeze(mean(low2,3));


test_lowcues_highstatic=(low1m+low2m)./2;
test_lowcues_highstaticm=squeeze(mean(test_lowcues_highstatic,3));
%highthreathighstim

%high cue 47degree
%90 h dynamic
high1=allcuescan1(:,:,[3,10,14],:);%sq1 
high1m=squeeze(mean(high1,3));
high2=allcuescan2(:,:,[9,16],:);%sq2
high2m=squeeze(mean(high2,3));

test_highcues_lowdynamic=(high1m+high2m)./3;
test_highcues_lowdynamicm=squeeze(mean(test_highcues_lowdynamic,3));

high1=allcuescan1(:,:,[5,8],:);%sq1 
high1m=squeeze(mean(high1,3));
high2=allcuescan2(:,:,[5,10,13],:);%sq2
high2m=squeeze(mean(high2,3));

test_highcues_lowstatic=(high1m+high2m)./3;
test_highcues_lowstaticm=squeeze(mean(test_highcues_lowstatic,3));

%%
%low vs high<-intensity coding
node=136
figure(1)
sem_low_dyn=std(squeeze(low3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(lowdynm(:,node))',sem_low_dyn,'b');
hold on


sem_high_dyn=std(squeeze(high3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(highdynm(:,node))',sem_high_dyn,'r');


hold on
 
 %plot(timing*5)

figure(2)
%% 
sem_high_stat=std(squeeze(high3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(highstatm(:,node))',sem_high_stat,'r');
hold on

sem_low_stat=std(squeeze(low3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(lowstatm(:,node))',sem_low_stat,'b');
hold on


%plot(timing*5)

%%
%compare stats vs dynamic 
node=132
figure(3)
sem_high_dyn=std(squeeze(high3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(highdynm(:,node))',sem_high_dyn,'r');

hold on

sem_high_stat=std(squeeze(high3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(highstatm(:,node))',sem_high_stat,'b');


node=52
figure(4)
sem_high_dyn=std(squeeze(high3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(highdynm(:,node))',sem_high_dyn,'r');

hold on

sem_high_stat=std(squeeze(high3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(highstatm(:,node))',sem_high_stat,'b');


%plot different rois same stim
node1=134
node2=135

sem_high_dyn1=std(squeeze(high3mdyn(:,node1,:))')./sqrt(30);
errorbar(squeeze(highdynm(:,node1))',sem_high_dyn1,'b');
hold on
sem_high_dyn2=std(squeeze(high3mdyn(:,node2,:))')./sqrt(30);
errorbar(squeeze(highdynm(:,node2))',sem_high_dyn2,'g');



node1=132
node2=133

sem_high_stat1=std(squeeze(high3mstat(:,node1,:))')./sqrt(30);
errorbar(squeeze(highstatm(:,node1))',sem_high_stat1,'b');
hold on
sem_high_stat2=std(squeeze(high3mstat(:,node2,:))')./sqrt(30);
errorbar(squeeze(highstatm(:,node2))',sem_high_stat2,'g');

% 
% figure(4)
% sem_low_stat=std(squeeze(low3mstat(:,node,:))')./sqrt(30);
% errorbar(squeeze(lowstatm(:,node))',sem_low_stat,'b');
% 
% hold on
% 
% sem_low_dyn=std(squeeze(low3mdyn(:,node,:))')./sqrt(30);
% errorbar(squeeze(lowdynm(:,node))',sem_low_dyn,'r');



%for averaging erps for type of events
%train-matched effects of positive cueing with matched stimulus:



%low cue low static	5,8,12
%low cue low dynamci	1,4,10
%high cue high static	3,6,11
%high cue high dynamic	2,7,9


%%
%average for mean activation to run stats

UnkHiStatpre=high3mstat((15:25),:,:);
UnkHiStat=squeeze(mean(UnkHiStatpre,1));

UnkHiDynpre=high3mdyn((15:25),:,:);
UnkHiDyn=squeeze(mean(UnkHiDynpre,1));


% UnkHiStatpre=high3mstat((15:25),:,:);
% UnkHiStat=squeeze(max(UnkHiStatpre,[],1));
% 
% UnkHiDynpre=high3mdyn((15:25),:,:);
% UnkHiDyn=squeeze(max(UnkHiDynpre,[],1));

UnkLowStatpre=low3mstat((15:25),:,:);
UnkLowStat=squeeze(mean(UnkLowStatpre,1));

UnkLowDynpre=low3mdyn((15:25),:,:);
UnkLowDyn=squeeze(mean(UnkLowDynpre,1));

trainHiStatpre=train_high3mstat((15:25),:,:);
trainHiStat=squeeze(mean(trainHiStatpre,1));

trainHiDynpre=train_high3mdyn((15:25),:,:);
trainHiDyn=squeeze(mean(trainHiDynpre,1));

trainLowStatpre=train_low3mstat((15:25),:,:);
trainLowStat=squeeze(mean(trainLowStatpre,1));

trainLowDynpre=train_low3mdyn((15:25),:,:);
trainLowDyn=squeeze(mean(trainLowDynpre,1));

%Positive Expectation
node=132

LocueHiStatpre=test_lowcues_highstatic((15:25),:,:);
LocueHiStat=squeeze(mean(LocueHiStatpre,1));

LocueHiDynpre=test_lowcues_highdynamic((15:25),:,:);
LocueHiDyn=squeeze(mean(LocueHiDynpre,1));

figure(1),plot(mean(squeeze(LocueHiStatpre(:,node,:)),2))
hold on

%For testing purposes
plot(mean(squeeze(LocueHiDynpre(:,node,:)),2),'r') %For testing purposes

%Negative Expectation

HicueLoStatpre=test_highcues_lowstatic((15:25),:,:);
HicueLoStat=squeeze(mean(HicueLoStatpre,1));

HicueLoDynpre=test_highcues_lowdynamic((15:25),:,:);
HicueLoDyn=squeeze(mean(HicueLoDynpre,1));

%%

[H,P,CI,STATS]=ttest(UnkHiStat',UnkHiDyn') %lots sig all four nac
[H,P,CI,STATS]=ttest(UnkHiStat',trainHiStat')% vmpfc significant
[H,P,CI,STATS]=ttest(UnkHiDyn',trainHiDyn') %sensory accumbens, lingual, pcune, ifgpo(bu) and lateralised for caudate and pcun
[H,P,CI,STATS]=ttest(LocueHiStat',trainHiStat') %not sig
[H,P,CI,STATS]=ttest(LocueHiDyn',trainHiDyn') %not sig
[H,P,CI,STATS]=ttest(trainHiStat',trainHiDyn') %nothing sig nac  {'137_NAcC_R.nii.gz'}0.04
[H,P,CI,STATS]=ttest(LocueHiStat',LocueHiDyn') %lots sig 001 acc sig after correction
[H,P,CI,STATS]=ttest(UnkHiStat',UnkLowStat')
intcodStat=UnkHiStat'- UnkLowStat';
intcodDyn=UnkHiDyn'-UnkLowDyn';
[H,P,CI,STATS]=ttest(intcodStat,intcodDyn)

pmStat=trainHiStat'- LocueHiStat';
pmDyn=trainHiDyn' - LocueHiDyn';
[H,P,CI,STATS]=ttest(pmStat,pmDyn)


[num_rejected, fdr_vec, idx]=fdr2(P(1:131),0.05);

  l1=labels(idx)
p1=P(idx)';
s1=STATS.tstat(idx)';
st=horzcat(p1,s1)


 [num_rejected, fdr_vec, idx] = fdr2(p, 0.1)

 labels(idx)
ind=find(h>0)
labels(ind)




[t,p]=ttest(UnkHiDyn',trainHiDyn') %lots sig <---
[t,p]=ttest(UnkHiStat',trainHiStat') %not sig

[t,p]=ttest(LocueHiStat',trainHiStat') %not sig even at 0.1
[t,p]=ttest(LocueHiDyn',trainHiDyn') %not sig even at 0.1


[h, crit_p, adj_p]=fdr2(p(1:138),0.05);
ind=find(h>0)
labels(ind)

[t,p]=ttest(UnkLowStat',UnkLowDyn') %nothing sig
[t,p]=ttest(trainLowStat',trainLowDyn') %nothing sig
[t,p]=ttest(HicueLoStat',HicueLoDyn') % not sign
 
[t,p]=ttest(UnkLowDyn',trainLowDyn') %not sig
[t,p]=ttest(UnkLowStat',trainLowStat') %not sig
[t,p]=ttest(HicueLoStat',trainLowStat') %not significant
 [t,p]=ttest(HicueLoDyn',trainLowDyn') %some sig: locci r, occfg r


[h, crit_p, adj_p]=fdr2(p(1:138),0.05);
ind=find(h>0)
labels(ind)

%%
%plot train stat vs dynamic in train only
node=132
figure(1)
sem_low_train_dyn=std(squeeze(train_low3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(train_lowdynm(:,node))',sem_low_train_dyn,'r');
hold on
sem_low_train_stat=std(squeeze(train_low3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(train_lowstatm(:,node))',sem_low_train_stat,'b');

figure(2)
sem_high_train_dyn=std(squeeze(train_high3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(train_highdynm(:,node))',sem_high_train_dyn,'r');
hold on
sem_high_train_stat=std(squeeze(train_high3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(train_highstatm(:,node))',sem_high_train_stat,'b');

%%
%plot effects of positive cueing with matched stimulus:
% cues vs unknown, stat vs dyn
%in low
node=132
figure(1)
sem_low_train_dyn=std(squeeze(train_low3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(train_lowdynm(:,node))',sem_low_train_dyn,'k');
hold on
sem_low_dyn=std(squeeze(low3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(lowdynm(:,node))',sem_low_dyn,'r');

figure(2)
sem_low_train_stat=std(squeeze(train_low3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(train_lowstatm(:,node))',sem_low_train_stat,'k');
hold on
sem_low_stat=std(squeeze(low3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(lowstatm(:,node))',sem_low_stat,'b');
%%

%plot effects of negative cueing with matched stimulus:
% cues vs unknown, stat vs dyn
%in high
node=133
figure(1)
sem_high_train_dyn=std(squeeze(train_high3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(train_highdynm(:,node))',sem_high_train_dyn,'k');
hold on
sem_high_dyn=std(squeeze(high3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(highdynm(:,node))',sem_high_dyn,'r');

figure(2)
sem_high_train_stat=std(squeeze(train_high3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(train_highstatm(:,node))',sem_high_train_stat,'k');
hold on
sem_high_stat=std(squeeze(high3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(highstatm(:,node))',sem_high_stat,'b');
%%
%effects of positive cueing with mismatched stimulus:
% positive cues matched vs mismtached vs unknown, stat vs dyn
%positive expectations
%low cues high static, test scan1 = 4,11,16
%low cues high static, test scan2 = 11,15
%low cues high dynamic, test scan1 =2,6,13
%low cues high dynamic, test scan1 =6,14

%negative expectations
%high cues low static, test scan1 =5,8
%high cues low static, test scan2 = 5,10,13
%high cues low dynamic, test scan1 = 3,10,14
%high cues low dynamic, test scan2 = 9,16




%%
%this!!
%plot effects of positive cueing with mismatched stimulus:
% cues vs unknown, stat vs dyn
%in high
node=135

figure(3)
sem_high_train_dyn=std(squeeze(train_high3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(train_highdynm(:,node))',sem_high_train_dyn,'k');
hold on
sem_high_dyn=std(squeeze(high3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(highdynm(:,node))',sem_high_dyn,'r');

hold on
sem_high_dyn_lowcue=std(squeeze(test_lowcues_highdynamic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_lowcues_highdynamicm(:,node))',sem_high_dyn_lowcue,'magenta');
hold on

%plot(timing*5)

figure(4)
sem_high_train_stat=std(squeeze(train_high3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(train_highstatm(:,node))',sem_high_train_stat,'k');
hold on
sem_high_stat=std(squeeze(high3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(highstatm(:,node))',sem_high_stat,'b');
hold on

sem_high_stat_lowcue=std(squeeze(test_lowcues_highstatic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_lowcues_highstaticm(:,node))',sem_high_stat_lowcue,'magenta');
hold on

%plot(timing*5)
%%
%and this!!
%plot effects of negative cueing with mismatched stimulus:
% cues vs unknown, stat vs dyn
%in high
node=133
figure(1)
sem_low_train_dyn=std(squeeze(train_low3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(train_lowdynm(:,node))',sem_low_train_dyn,'k');
hold on
sem_low_dyn=std(squeeze(low3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(lowdynm(:,node))',sem_low_dyn,'r');

hold on
sem_low_dyn_highcue=std(squeeze(test_highcues_lowdynamic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_highcues_lowdynamicm(:,node))',sem_low_dyn_highcue,'magenta');
hold on

%plot(timing*5)

figure(2)
sem_low_train_stat=std(squeeze(train_low3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(train_lowstatm(:,node))',sem_low_train_stat,'k');
hold on
sem_low_stat=std(squeeze(low3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(lowstatm(:,node))',sem_low_stat,'b');
hold on

sem_low_stat_highcue=std(squeeze(test_highcues_lowstatic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_highcues_lowstaticm(:,node))',sem_low_stat_highcue,'magenta');
hold on

%plot(timing*5)

%%
%plot diff between stats and dyn for mismatched

figure(1)
node=135
figure(3)
sem_high_dyn_lowcue=std(squeeze(test_lowcues_highdynamic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_lowcues_highdynamicm(:,node))',sem_high_dyn_lowcue,'r');
hold on

sem_high_stat_lowcue=std(squeeze(test_lowcues_highstatic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_lowcues_highstaticm(:,node))',sem_high_stat_lowcue,'b');

figure(4)
sem_low_dyn_highcue=std(squeeze(test_highcues_lowdynamic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_highcues_lowdynamicm(:,node))',sem_low_dyn_highcue,'r');
hold on

sem_low_stat_highcue=std(squeeze(test_highcues_lowstatic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_highcues_lowstaticm(:,node))',sem_low_stat_highcue,'b');



%%


%is high vs low temperature encoded differently? for dynamic
node=51

figure(1)
sem_high_train_dyn=std(squeeze(train_high3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(train_highdynm(:,node))',sem_high_train_dyn,'r');

hold on
sem_low_train_dyn=std(squeeze(train_low3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(train_lowdynm(:,node))',sem_low_train_dyn,'b');
hold on
plot(timing*5)

figure(2)

sem_high_dyn=std(squeeze(high3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(highdynm(:,node))',sem_high_dyn,'r');

hold on
sem_low_dyn=std(squeeze(low3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(lowdynm(:,node))',sem_low_dyn,'b');
hold on
plot(timing*5)

%%
node=139
figure(1)
sem_high_train_stat=std(squeeze(train_high3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(train_highstatm(:,node))',sem_high_train_stat,'r');
hold on
sem_low_train_stat=std(squeeze(train_low3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(train_lowstatm(:,node))',sem_low_train_stat,'b');
hold on
plot(timing*5)

figure(2)
sem_high_stat=std(squeeze(high3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(highstatm(:,node))',sem_high_stat,'r');
hold on

sem_low_stat=std(squeeze(low3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(lowstatm(:,node))',sem_low_stat,'b');
hold on

plot(timing*5)

%%
sem_high_stat_lowcue=std(squeeze(test_lowcues_highstatic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_lowcues_highstatic(:,node))',sem_high_stat_lowcue,'magenta');
hold on

plot(timing*5)



figure(2)
sem_low_train_stat=std(squeeze(train_low3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(train_lowstatm(:,node))',sem_low_train_stat,'k');
hold on
sem_low_stat=std(squeeze(low3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(lowstatm(:,node))',sem_low_stat,'b');
hold on

sem_low_stat_highcue=std(squeeze(test_highcues_lowstatic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_highcues_lowstatic(:,node))',sem_low_stat_highcue,'magenta');
hold on

plot(timing*5)

%%

figure(3)

sem_high_dyn_lowcue=std(squeeze(test_lowcues_highdynamic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_lowcues_highdynamic(:,node))',sem_high_dyn_lowcue,'magenta');
hold on

sem_low_dyn_highcue=std(squeeze(test_highcues_lowdynamic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_highcues_lowdynamic(:,node))',sem_low_dyn_highcue,'magenta');
hold on

plot(timing*5)


figure(1)
sem_low_train_dyn=std(squeeze(train_low3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(train_lowdynm(:,node))',sem_low_train_dyn,'k');
hold on
sem_low_dyn=std(squeeze(low3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(lowdynm(:,node))',sem_low_dyn,'r');

hold on
sem_low_dyn_highcue=std(squeeze(test_highcues_lowdynamic(:,node,:))')./sqrt(30);
errorbar(squeeze(test_highcues_lowdynamic(:,node))',sem_low_dyn_highcue,'magenta');
hold on

%%
%high vs low wm


load /Users/jhash1/Dropbox/CBP_study/TIMESERIES_CBP/macros/CBP_high_low_WM.txt
nodes = 131;
%labels=textread (['/Users/jhash1/Dropbox/GAC_study/gac_timeseries/roilist_thisone.txt'],'%s');
labels=textread (['/Users/jhash1/Dropbox/GAC_study/gac_timeseries/roilist_pag.txt'],'%s');

%for low cues

hwm=erpl(:,CBP_high_low_WM==1);
lwm=erpl(:,CBP_high_low_WM==0);
group1 = hwm';   %group1 matrices: nsubj x nroi x nroi
group2 = lwm';
[H,P,CI,STATS] = ttest2(group1,group2);
Pl=squeeze(P);
stats=squeeze(STATS.tstat);


%corrected
[num_rejected, fdr_vec, idx] = fdr2(Pl, 0.05) %correct for multiple comparisons
labels(idx)

%uncorrected for few regions
%example
Pl(66)

ind=find(Pl<0.01)


%%

Y=vertcat(erph',erpl');

for i=1:135
    
    P=anovan(Y(:,i),group_twoway_raw,'display','off')
   p1(i,:)=P(1)
   p2(i,:)=P(2)
end


ind=find(p2<0.05)
p2(ind)
labels(ind)

[num_rejected, fdr_vec, idx] = fdr2(p1, 0.1)

%[num_rejected, fdr_vec, idx] = fdr2(p2, 0.O5)


hwm=erpl(:,CBP_high_low_WM==1);
lwm=erpl(:,CBP_high_low_WM==0);
group1 = hwm';   %group1 matrices: nsubj x nroi x nroi
group2 = lwm';
[H,P,CI,STATS] = ttest2(group1,group2);
Pl=squeeze(P);
stats=squeeze(STATS.tstat);
%%
%for high cues:
hwm=erph(idx,CBP_high_low_WM==1);
lwm=erph(idx,CBP_high_low_WM==0);
group1 = hwm';   %group1 matrices: nsubj x nroi x nroi
group2 = lwm';
[H,P,CI,STATS] = ttest2(group1,group2);
Ph=squeeze(P);
stats=squeeze(STATS.tstat);

%corrected
[num_rejected, fdr_vec, idx] = fdr2(Ph, 0.05) %correct for multiple comparisons
labels(idx)

%uncorrected for few regions
%example
Ph(66)


%%
%for unknown cues:
hwm=erpu(:,CBP_high_low_WM==1);
lwm=erpu(:,CBP_high_low_WM==0);
group1 = hwm';   %group1 matrices: nsubj x nroi x nroi
group2 = lwm';
[H,P,CI,STATS] = ttest2(group1,group2);
Pu=squeeze(P);
stats=squeeze(STATS.tstat);

%corrected
[num_rejected, fdr_vec, idx] = fdr2(Pu, 0.05) %correct for multiple comparisons
labels(idx)

%uncorrected for few regions
%example
Pu(66)


%%
%havent fixed this yet
figure(1)

node2=17;
sem_high=std(squeeze(highcues(:,node2,:))')./sqrt(38);
yhigh=sem_high;
%errorbar(squeeze(highcueallmean(:,24))',sem_high,'k');
%hold on

%plot(epoch_timing)
x=[1:39]';
y=highcueallmean';
dyhigh=y(node2,:)';
plot(yhigh,'linewidth',1,'color','k')
hold on
fill([x;flipud(x)],[yhigh-dyhigh;flipud(yhigh+dyhigh)],'r','linestyle','-','linewidth',0.01,'Facealpha', 0.3);

hold on

sem_low=std(squeeze(lowcues(:,node2,:))')./sqrt(38);
dylow=sem_low';
% errorbar(squeeze(lowcueallmean(:,24))',sem_high,'k');
% hold on
%plot(epoch_timing)
y1=lowcueallmean';
ylow=y1(node2,:)';
plot(ylow,'linewidth',1,'color','k')



fill([x;flipud(x)],[ylow-dylow;flipud(ylow+dylow)],'k','linestyle','-','linewidth',0.01,'Facealpha', 0.3);
hold on

plot(epoch_timing,'k')


%%
figure(2)
node2=23;
sem_high=std(squeeze(highcues(node2,:,:))')./sqrt(38);
%errorbar(squeeze(highcueallmean(:,24))',sem_high,'k');
%hold on
dyhigh=sem_high;
%plot(epoch_timing)
x=[1:39]';
y=highcueallmean';
yhigh=y(:,node2)';
plot(yhigh,'linewidth',1,'color','k')
hold on
fill([x;flipud(x)],[yhigh-dyhigh;flipud(yhigh+dyhigh)],'r','linestyle','-','linewidth',0.01,'Facealpha', 0.3);

hold on

sem_low=std(squeeze(lowcues(node2,:,:))')./sqrt(38);
dylow=sem_low';
% errorbar(squeeze(lowcueallmean(:,24))',sem_high,'k');
% hold on
%plot(epoch_timing)
y1=lowcueallmean';
ylow=y1(node2,:)';
plot(ylow,'linewidth',1,'color','k')



fill([x;flipud(x)],[ylow-dylow;flipud(ylow+dylow)],'k','linestyle','-','linewidth',0.01,'Facealpha', 0.3);
hold on

plot(epoch_timing,'k')
%%

load /[changepath]/CBP_high_low_suggestibility.txt %code for high vs low suggestibility
nodes = 131;
%labels=textread (['/Users/jhash1/Dropbox/GAC_study/gac_timeseries/roilist_thisone.txt'],'%s');
labels=textread (['[changepath]/roilist.txt'],'%s');

%for low cues

hsugg=ccall_noz_nothresh(CBP_high_low_suggestibility==1,:,:);
lsugg=ccall_noz_nothresh(CBP_high_low_suggestibility==0,:,:);
 group1 = hsugg;   %group1 matrices: nsubj x nroi x nroi
 group2 = lsugg;
[H,P,CI,STATS] = ttest2(group1,group2);
Pl=squeeze(P);
stats=squeeze(STATS.tstat);

%corrected
[num_rejected, fdr_vec, idx] = fdr(Pl, 0.05) %correct for multiple comparisons
labels(idx)

%uncorrected for few regions
%example
Pl(66)

ind=find(Pl<0.01)



%%
%high vs low wm


%load /Users/jenni/OneDrive/Documents/CBP/CBP_high_low_WM.txt
nodes = 134;
%labels=textread (['/Users/jhash1/Dropbox/GAC_study/gac_timeseries/roilist_thisone.txt'],'%s');
labels=textread (['/Users/jenni/OneDrive/Documents/CBP/roilist_new_Hashmi_135.txt'],'%s');

%for low cues

hwm=erpl(:,CBP_high_low_WM==1);
lwm=erpl(:,CBP_high_low_WM==0);
group1 = hwm';   %group1 matrices: nsubj x nroi x nroi
group2 = lwm';
[H,P,CI,STATS] = ttest2(group1,group2);
Pl=squeeze(P);
stats=squeeze(STATS.tstat);

%corrected
[num_rejected, fdr_vec, idx] = fdr2(Pl, 0.05) %correct for multiple comparisons
labels(idx)

%uncorrected for few regions
%example
Pl(66)
%%
bla1=erph([66,67],:);
bla2=erpl([66,67],:);
bla=vertcat(bla1,bla2); bla=bla';
vertcat
bla=vertcat(erph,erpl);
 p = anova_rm(bla);
 
 
[p,tbl,stats] = anova1(bla,hwm1_lwm0)
[p,tbl,stats] = anova1(ccall_noz_nothresh,codeforimp)

[H,P,CI,STATS] = anova1(group1,group2,group3);
P=squeeze(P);
stats=squeeze(STATS.tstat);
  
%%
%for high cues:
%load path/hwm1_lwm0.mat
hwm=erph(:,hwm1_lwm0==1);
lwm=erph(:,hwm1_lwm0==0);
group1 = hwm';   %group1 matrices: nsubj x nroi x nroi
group2 = lwm';
[H,P,CI,STATS] = ttest2(group1,group2);
Ph=squeeze(P);
stats=squeeze(STATS.tstat);

%corrected
[num_rejected, fdr_vec, idx] = fdr2(Ph, 0.05) %correct for multiple comparisons
labels(idx)

%uncorrected for few regions
%example
Ph(134)


%%
%for unknown cues:
hwm=erpu(:,CBP_high_low_WM==1);
lwm=erpu(:,CBP_high_low_WM==0);
group1 = hwm';   %group1 matrices: nsubj x nroi x nroi
group2 = lwm';
[H,P,CI,STATS] = ttest2(group1,group2);
Pu=squeeze(P);
stats=squeeze(STATS.tstat);

%corrected
[num_rejected, fdr_vec, idx] = fdr2(Pu, 0.1) %correct for multiple comparisons
labels(idx)

%uncorrected for few regions
%example
Pu(134)

%%
[pvals] = anovan(cl(1,:)',{groups});


%%
%havent fixed this yet
figure(1)

node2=17;
sem_high=std(squeeze(highcues(:,node2,:))')./sqrt(38);
yhigh=sem_high;
%errorbar(squeeze(highcueallmean(:,24))',sem_high,'k');
%hold on

%plot(epoch_timing)
x=[1:39]';
y=highcueallmean';
dyhigh=y(node2,:)';
plot(yhigh,'linewidth',1,'color','k')
hold on
fill([x;flipud(x)],[yhigh-dyhigh;flipud(yhigh+dyhigh)],'r','linestyle','-','linewidth',0.01,'Facealpha', 0.3);

hold on

sem_low=std(squeeze(lowcues(:,node2,:))')./sqrt(38);
dylow=sem_low';
% errorbar(squeeze(lowcueallmean(:,24))',sem_high,'k');
% hold on
%plot(epoch_timing)
y1=lowcueallmean';
ylow=y1(node2,:)';
plot(ylow,'linewidth',1,'color','k')



fill([x;flipud(x)],[ylow-dylow;flipud(ylow+dylow)],'k','linestyle','-','linewidth',0.01,'Facealpha', 0.3);
hold on

plot(epoch_timing,'k')


%%
figure(2)
node2=23;
sem_high=std(squeeze(highcues(node2,:,:))')./sqrt(38);
%errorbar(squeeze(highcueallmean(:,24))',sem_high,'k');
%hold on
dyhigh=sem_high;
%plot(epoch_timing)
x=[1:39]';
y=highcueallmean';
yhigh=y(:,node2)';
plot(yhigh,'linewidth',1,'color','k')
hold on
fill([x;flipud(x)],[yhigh-dyhigh;flipud(yhigh+dyhigh)],'r','linestyle','-','linewidth',0.01,'Facealpha', 0.3);

hold on

sem_low=std(squeeze(lowcues(node2,:,:))')./sqrt(38);
dylow=sem_low';
% errorbar(squeeze(lowcueallmean(:,24))',sem_high,'k');
% hold on
%plot(epoch_timing)
y1=lowcueallmean';
ylow=y1(node2,:)';
plot(ylow,'linewidth',1,'color','k')



fill([x;flipud(x)],[ylow-dylow;flipud(ylow+dylow)],'k','linestyle','-','linewidth',0.01,'Facealpha', 0.3);
hold on

plot(epoch_timing,'k')

%%
%%

sem_high_dyn=std(squeeze(high3mdyn(:,node,:))')./sqrt(30);
errorbar(squeeze(highdynm(:,node))',sem_high_dyn,'r');

figure(2)
sem_high_train_stat=std(squeeze(train_high3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(train_highstatm(:,node))',sem_high_train_stat,'k');
hold on
sem_high_stat=std(squeeze(high3mstat(:,node,:))')./sqrt(30);
errorbar(squeeze(highstatm(:,node))',sem_high_stat,'b');
%%


%FINAL FIGURE CODE unknown stat vs dynamic
%final plot

%change: 1 node number, 2. name of inputs (mean across subjects, and all subjects
%for sem, 3. name of title and 4. name of file to be saved. see if you need to change 
% xlim, ylim, explore color options to personalise

figure(3) %vlpag is 24 dlpag is 23
node2=80;
%node 80 for visual

sem_high=std(squeeze(high3mstat(:,node2,:))')./sqrt(30);
dyhigh=sem_high';

x=[1:39]';
sem_low=std(squeeze(high3mdyn(:,node2,:))')./sqrt(30);
dylow=sem_low';

y1=highstatm';
ylow=y1(node2,:)';
plot(ylow,'linewidth',5,'color','yellow')
hold on
%plot(dldeg3,'color',[0.85 0.33 0.1],'linewidth',2)
fill([x;flipud(x)],[ylow-dylow;flipud(ylow+dylow)], [0.11 0.48 0.72],'linestyle','none','linewidth',0.001,'Facealpha', 0.7);

hold on
y=highdynm';
yhigh=y(node2,:)';
plot(yhigh,'linewidth',4,'color',[0.7,0.7,0.7])
hold on
%[0.83 0.31 0.18]reddish, [0.3010,0.7450,0.9330]
%bluish,[0.8500,0.3250,0.0980]rustish
fill([x;flipud(x)],[yhigh-dyhigh;flipud(yhigh+dyhigh)],[0.4660,0.6740,0.1880],'linestyle','none','linewidth',0.001,'Facealpha', 0.7);

ylim([-40.0 45.0])
set(gca,'YTick', -40:10:45)
xlabel('Time points')
xlim([0 37.0])
set(gca,'XTick', 0:10:37)
ylabel('% fMRI activation')

hold on
bar=zeros(1,39);
bar=bar*NaN;
bar(4:9)=-35;
bar(14:23)=-35;
bar(27:34)=-35;

plot(bar,'linewidth',8,'color',[0.7,0.7,0.7])

fontsize(figure(3), 24, "points")

set(gca,'fontname','arial')

set(gca,'box','off','color','white')
set(gca,'TickDir','out');
set(gca,'LineWidth',1.6,'TickLength',[0.015 0.015]);

set(gcf,'color','w');
set(gca,'units','pix','pos',[100,100,500,300])



f = figure(3);
f.Position = [600 100 540 400];


title('035FOR') %change title here
x0=600;
y0=600;
width=700;
height=500
set(gcf,'position',[x0,y0,width,height])
%You can specify other units (inches, centimeters, normalized, points, or characters). For example:
set(gcf,'units','points','position',[x0,y0,width,height])

%saveas(gcf,'/Users/jhash1/Dropbox/Adam_Sunavsky/Dynamic_Adam_Study/results/figures_jh/035FOR.png') %add name and path here

%xtickangle(45)
%ytickangle(45) %90


%%


%FINAL FIGURE CODE unknown stat vs dynamic
%final plot

%change: 1 node number, 2. name of inputs (mean across subjects, and all subjects
%for sem, 3. name of title and 4. name of file to be saved. see if you need to change 
% xlim, ylim, explore color options to personalise

figure(2) %vlpag is 24 dlpag is 23
node2=137;
%node 80 for visual

sem_high=std(squeeze(train_high3mdyn(:,node2,:))')./sqrt(30);
dyhigh=sem_high';

x=[1:39]';
sem_low=std(squeeze(high3mdyn(:,node2,:))')./sqrt(30);
dylow=sem_low';

y1=train_highdynm';
ylow=y1(node2,:)';
plot(ylow,'linewidth',5,'color','yellow')
hold on
%plot(dldeg3,'color',[0.85 0.33 0.1],'linewidth',2)
fill([x;flipud(x)],[ylow-dylow;flipud(ylow+dylow)], [0.11 0.48 0.72],'linestyle','none','linewidth',0.001,'Facealpha', 0.7);

hold on
y=highdynm';
yhigh=y(node2,:)';
plot(yhigh,'linewidth',4,'color',[0.7,0.7,0.7])
hold on
%[0.83 0.31 0.18]reddish, [0.3010,0.7450,0.9330]
%bluish,[0.8500,0.3250,0.0980]rustish
fill([x;flipud(x)],[yhigh-dyhigh;flipud(yhigh+dyhigh)],[0.4660,0.6740,0.1880],'linestyle','none','linewidth',0.001,'Facealpha', 0.7);

ylim([-40.0 45.0])
set(gca,'YTick', -40:10:45)
xlabel('Time points')
xlim([0 37.0])
set(gca,'XTick', 0:10:37)
ylabel('% fMRI activation')

hold on
bar=zeros(1,39);
bar=bar*NaN;
bar(4:9)=-35;
bar(14:23)=-35;
bar(27:34)=-35;

plot(bar,'linewidth',8,'color',[0.7,0.7,0.7])

fontsize(figure(2), 24, "points")

set(gca,'fontname','arial')

set(gca,'box','off','color','white')
set(gca,'TickDir','out');
set(gca,'LineWidth',1.6,'TickLength',[0.015 0.015]);

set(gcf,'color','w');
set(gca,'units','pix','pos',[100,100,500,300])



f = figure(2);
f.Position = [600 100 540 400];


title('137NAccR') %change title here
x0=600;
y0=600;
width=700;
height=500
set(gcf,'position',[x0,y0,width,height])
%You can specify other units (inches, centimeters, normalized, points, or characters). For example:
set(gcf,'units','points','position',[x0,y0,width,height])

saveas(gcf,'/Users/jhash1/Dropbox/Adam_Sunavsky/Dynamic_Adam_Study/results/figures_jh/137NAccR.png') %add name and path here

%xtickangle(45)
%ytickangle(45) %90





%%


%FINAL FIGURE CODE unknown stat vs dynamic
%final plot

%change: 1 node number, 2. name of inputs (mean across subjects, and all subjects
%for sem, 3. name of title and 4. name of file to be saved. see if you need to change 
% xlim, ylim, explore color options to personalise

figure(2) %vlpag is 24 dlpag is 23
node2=112;
%node 80 for visual

sem_high=std(squeeze(test_lowcues_highdynamic(:,node2,:))')./sqrt(30);
dyhigh=sem_high';

x=[1:39]';
sem_low=std(squeeze(test_lowcues_highstatic(:,node2,:))')./sqrt(30);
dylow=sem_low';

y1=test_lowcues_highdynamicm';
ylow=y1(node2,:)';
plot(ylow,'linewidth',5,'color','yellow')
hold on
%plot(dldeg3,'color',[0.85 0.33 0.1],'linewidth',2)
fill([x;flipud(x)],[ylow-dylow;flipud(ylow+dylow)], [0.11 0.48 0.72],'linestyle','none','linewidth',0.001,'Facealpha', 0.7);

hold on
y=test_lowcues_highstaticm';
yhigh=y(node2,:)';
plot(yhigh,'linewidth',4,'color',[0.7,0.7,0.7])
hold on
%[0.83 0.31 0.18]reddish, [0.3010,0.7450,0.9330]
%bluish,[0.8500,0.3250,0.0980]rustish
fill([x;flipud(x)],[yhigh-dyhigh;flipud(yhigh+dyhigh)],[0.4660,0.6740,0.1880],'linestyle','none','linewidth',0.001,'Facealpha', 0.7);

ylim([-40.0 45.0])
set(gca,'YTick', -40:10:45)
xlabel('Time points')
xlim([0 37.0])
set(gca,'XTick', 0:10:37)
ylabel('% fMRI activation')

hold on
bar=zeros(1,39);
bar=bar*NaN;
bar(4:9)=-35;
bar(14:23)=-35;
bar(27:34)=-35;

plot(bar,'linewidth',8,'color',[0.7,0.7,0.7])

fontsize(figure(2), 24, "points")

set(gca,'fontname','arial')

set(gca,'box','off','color','white')
set(gca,'TickDir','out');
set(gca,'LineWidth',1.6,'TickLength',[0.015 0.015]);

set(gcf,'color','w');
set(gca,'units','pix','pos',[100,100,500,300])



f = figure(2);
f.Position = [600 100 540 400];



x0=600;
y0=600;
width=700;
height=500
set(gcf,'position',[x0,y0,width,height])
%You can specify other units (inches, centimeters, normalized, points, or characters). For example:
set(gcf,'units','points','position',[x0,y0,width,height])


title('13AmygL') %change title here
saveas(gcf,'/Users/jhash1/Dropbox/Adam_Sunavsky/Dynamic_Adam_Study/results/figures_jh/112SPLL.png') %add name and path here

%xtickangle(45)
%ytickangle(45) %90


%%
% hold on
%
% plot(epoch_timing,'k')
%hold on
%sem_low=std(squeeze(unkcues(:,node2,:))')./sqrt(38);
%dyunk=sem_unk';
% errorbar(squeeze(lowcueallmean(:,24))',sem_high,'k');
% hold on
%plot(epoch_timing)
%y1=unkcueallmean';
%yunk=y1(node2,:)';
%plot(yunk,'linewidth',2,'color','y')
%plot(dldeg3,'color',[0.85 0.33 0.1],'linewidth',2)
%fill([x;flipud(x)],[yunk-dyunk;flipud(yunk+dyunk)],'y','linestyle','none','linewidth',0.001,'Facealpha', 0.3);