SyncFile = readtable(fullfile(DataDir,'synced.csv'));
ChOfInterst = 'GPIO-2'

SyncIndx = find(ismember(SyncFile.ChannelName,ChOfInterst));
SyncTime = SyncFile.Time_s_(SyncIndx);
SR = 1/nanmedian(diff(SyncTime));

Values = SyncFile.Value(SyncIndx);
tresh = nanmax(Values(:))/2;
Values(Values<=tresh)=0;
Values(Values>tresh)=1;

findonset = find(Values==1,1);
SyncTime = SyncTime-SyncTime(findonset);
SyncTime(1:findonset-1)=[];
Values(1:findonset-1)=[];



%% Timeline
TimelineFile = dir(fullfile(DataDir,'*Timeline.mat'));
Timeline = load(fullfile(TimelineFile(1).folder,TimelineFile(1).name));
Timeline = Timeline.Timeline;
AllInputs = {Timeline.hw.inputs(:).name}; %Just redefine AllInputs in case
AllInputs = {'timestamps',AllInputs{:}}; % the first one is always timstamps
inputs = Timeline.hw.inputs;
flippersignal = Timeline.rawDAQData(:,strcmp({inputs(:).name},'flipper'));
tresh = nanmax(flippersignal(:))/2;
flippersignal(flippersignal<=tresh)=0;
flippersignal(flippersignal>tresh)=1;
Actualtime = Timeline.rawDAQTimestamps;


tmpdat = cat(2,Timeline.rawDAQTimestamps',Timeline.rawDAQData);
Timeline = tmpdat;

findonset = find(flippersignal==1,1);
Actualtime = Actualtime-Actualtime(findonset);
Actualtime(1:findonset-1)=[];
flippersignal(1:findonset-1)=[];
Timeline(1:findonset-1,:) =[];

figure; 
plot(SyncTime,Values)
hold on
plot(Actualtime,flippersignal+1,'r')

%% Match inscopixtime - have to match to actual time because otherwise it's so undersampled there's nothing in flipper
timesteps = unique(round(diff(inscopixtime)*1000))/1000;
tmSR = unique(round(diff(Actualtime)*1000))/1000;

resampleval = timesteps/tmSR;
tmSR = 1/tmSR;
newinsopixtime = interp(inscopixtime,resampleval);

newflipper = zeros(1,length(newinsopixtime));
for tp = 1:length(Values)-1
    newflipper(newinsopixtime>=SyncTime(tp)&newinsopixtime<SyncTime(tp+1))=Values(tp);
end

plot(newinsopixtime,newflipper+2,'g')

%% If this looks good, we can resample neural data as well
DataMatrixR = arrayfun(@(X) interp(DataMatrix(:,X),resampleval),1:size(DataMatrix,2),'UniformOutput',0);
DataMatrix = cat(2,DataMatrixR{:});

%% now all left is to exactly match newinsopixtime to actualtime
DataMatrixN = nan(length(Actualtime),size(DataMatrix,2));
newflipperN = nan(1,length(Actualtime)); 
for tp = 1:length(newinsopixtime)-1  
    idx = find(Actualtime>=newinsopixtime(tp)&Actualtime<newinsopixtime(tp+1));
    for id = idx
        DataMatrixN(id,:) = DataMatrix(tp,:);
        newflipperN(idx) = newflipper(tp);
    end
end

plot(Actualtime,newflipperN+3,'b')

%% Save in this format
DataMatrix = DataMatrixN;
clear DataMatrixN
clear DataMatrixR
