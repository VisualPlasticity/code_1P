%% User Input
DataDir = {'\\live.rd.ucl.ac.uk\ritd-ag-project-rd0164-yjsun44\inscopix\session-20211118\analysis'}%'\\znas\Subjects' %'
SaveDir = 'E:\2PData\Preproccesed\'
LocalDir = 'E:\2PData\' 
MiceOpt = {'CR_Hippocanula1'};% {};%{}%{'EB006'};%,'CB007','CB008'}; %,'EB001','EB003'}%{'CB008'};%{'EB001'}%{'EB001','EB002','EB003','CB007','CB008'};%,'CB007','CB008'} %'CB007'
DataDir2Use = repmat(2,[1,length(MiceOpt)]);

maxsessnr = 2; %max nr. sessions on a day (doesn't need to be accurate)
MinDist2Include = 85; %Including only trials in which mouse reached at least xcm for some statistics (past all stimuli?)
pretrialtime = 2; %take up to x seconds prior trial
posttrialtime = 2; % take up to x seconds post trial
addpath(genpath(cd))
addpath(genpath('C:\Users\EnnyB\Documents\GitHub\Inscopix_Ana\Inscopix_Ana'))

%% Receptive field analysis
LoadSpikesRF_Inscopix

