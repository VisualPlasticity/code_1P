%% User Input
fancymethod = 0
% DataDir = '\\live.rd.ucl.ac.uk\ritd-ag-project-rd0164-yjsun44\inscopix\session-20211118\10'
DataDir = strrep('..\inscopix\session-20211126\10','\',filesep);
% DataType = 'celltraces';
DataType = 'spikingAllTimes';

SaveDir = fullfile(DataDir,'Processed');

%% Automated
addpath(genpath('Dependencies'))

listfiles = dir(fullfile(DataDir))
listfiles(~cell2mat(cellfun(@isempty,cellfun(@(X) strfind(X,'.pickle'),{listfiles(:).name},'UniformOutput',0),'UniformOutput',0))) = [];
idx = find(~cell2mat(cellfun(@isempty,cellfun(@(X) strfind(X,'.p'),{listfiles(:).name},'UniformOutput',0),'UniformOutput',0)));
% read which x-file was used
if length(idx)>1
    disp('Too many files?!')
    keyboard
end
fileID = fopen(fullfile(listfiles(idx).folder,listfiles(idx).name));
A = fscanf(fileID,'%c');
fclose(fileID);

if any(strfind(A,'SparseNoise'))
    % Load protocol to check for contrast >0
    try
        Protocol = load(fullfile(fullfile(listfiles(idx).folder,'Protocol.mat')));
    catch ME
        disp(ME)
    end
    Protocol = Protocol.Protocol;
    
end

%% Loading data from inscopix
myNeuralData = dir(fullfile(DataDir,['*' DataType '.csv']));

if isempty('myNeuralData')
    disp('No data found!!')
    keyboard
end

%Saving directory
if ~isdir(fullfile(SaveDir))
    mkdir(fullfile(SaveDir))
end

sp = readtable(fullfile(myNeuralData(1).folder,myNeuralData(1).name));
head(sp)
try
    inscopixtime = sp.Var1;
catch
    inscopixtime = sp.Time_s_;
end
    %inscopixtime = cell2mat(sp.Var1(2:end));%exclude headline 'Time(s)/cell status
DataMatrix = table2array(sp);
DataMatrix = DataMatrix(:,2:end); % exclude time
nROI = size(DataMatrix,2);

%% synchronization data
AlignSyncWithTimeline_Inscopix % This will interpollate imaging data to be on the same scale as Actualtime

%% Load RF mapping session details
if exist(fullfile(DataDir,'RF_MappingStimuli.mat'))
    SS = load(fullfile(DataDir,'RF_MappingStimuli.mat'));
else
    disp([fullfile(DataDir,'RF_MappingStimuli.mat') ' does not yet exist.. create on rig using RFmapping_RegenerateStimuli'])
    keyboard
end

% extract trial information
SS = SS.SS;
ntrials = length(SS);
TrialDurations = cell2mat(cellfun(@(X) length(X.ImageSequence)./60,SS,'UniformOutput',0)); %in sec*10


%% Align to Trial Onset times
[starttrialidx,endtrialidx,trialid] = FindTrialOnAndOffSetsPhotoDiode(Actualtime,Timeline(:,ismember(AllInputs,'syncEcho')),TrialDurations);

if any((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations'>1/tmSR*10)
    disp(['flips slightly drifting... Average of ' num2str(nanmean((Actualtime(endtrialidx)-Actualtime(starttrialidx)) - TrialDurations')) 'sec'])
end


if ntrials ~= length(starttrialidx)
    warning('Can''t find enough trials')
    keyboard
end


%% Prepare Noise stimulus frames in time X rows X Cols
tmp = Timeline(:,ismember(AllInputs,'syncEcho'));
tmp = (tmp-nanmin(tmp(:)))./(nanmax(tmp(:))-nanmin(tmp(:)));
tmp = lowpass(tmp,70,tmSR);
tresh=0.3;%nanmedian(tmp);


tmp(tmp>=tresh)=1;
tmp(tmp<tresh) = 0;
photodiodetrace = tmp;

figure; plot(Actualtime,photodiodetrace);
box off
tmpimg = cat(3,SS{1}.ImageTextures{:});
uniquecolors = unique(tmpimg(:));
Greyvalue = uniquecolors(2);
stimFrames = repmat(Greyvalue,ceil(length(Actualtime)./(tmSR/60)),size(tmpimg,1),size(tmpimg,2));
countimagetotal = 1;
Stimulustimes = nan(1,ceil(length(Actualtime)./(tmSR/60)));
for tridx=1:ntrials
    countimg = 0;
    counttime = 1/60;
    diodestatus = photodiodetrace(1);
    missedflips = [];
    framesmissing = 0;
    for tp=starttrialidx(tridx):endtrialidx(tridx)
        % photodiodestate
        if photodiodetrace(tp)~= diodestatus %60 Hz flips
            %                         if counttime>1/60*1.6
            %                             xlim([Actualtime(tp-500) Actualtime(tp+500)])
            %                             if exist('h')
            %                                 delete(h)
            %                             end
            %                             h=line([Actualtime(tp) Actualtime(tp)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--');
            %                             disp('too slow')
            %                             missedflips = [missedflips tp];
            %                             drawnow
            %                             pause(0.5)
            %
            %                         end
            if counttime<1/60*0.75
                xlim([Actualtime(tp-500) Actualtime(tp+500)])
                if exist('h')
                    delete(h)
                end
                h=line([Actualtime(tp) Actualtime(tp)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--');
                disp('Too early!?')
                drawnow
                pause(0.5)
            else
                countimg = countimg+1; %go to next image
                %Load current image
                if countimg > length(SS{tridx}.ImageTextures)
                    xlim([Actualtime(tp-500) Actualtime(tp+500)])
                    if exist('h')
                        delete(h)
                    end
                    h=line([Actualtime(tp) Actualtime(tp)],get(gca,'ylim'),'color',[1 0 0],'LineStyle','--');
                    
                    
                    countimg = length(SS{tridx}.ImageTextures);
                    framesmissing = framesmissing+1;
                    drawnow
                    pause(0.5)
                    
                    if length(framesmissing)>1 %1 last flip is grey screen (should be)
                        disp('More Photodiode flips than frames?')
                        keyboard
                    end
                end
                
                tmpimg =  SS{tridx}.ImageTextures{SS{tridx}.ImageSequence(countimg)};
                stimFrames(countimagetotal,:,:) = tmpimg;
                Stimulustimes(countimagetotal)=Actualtime(tp);
                countimagetotal = countimagetotal+1;
                
                counttime = 0;
                diodestatus=photodiodetrace(tp);
            end
        end
        counttime = counttime+(1/tmSR);
    end
end
%interpolate nans stimulustime
Stimulustimes = fillmissing(Stimulustimes,'linear');

stimFrameDur = median(diff(Stimulustimes));


%RunningSpeed
RotarySignal = (Timeline(:,ismember(AllInputs,'rotaryEncoder')));
RotarySignalPerSec = nan(1,length(RotarySignal));
windowsz = tmSR; %per second
parfor tp = 1:length(RotarySignal)
    tptotake = [tp-windowsz./2:tp+windowsz./2];
    tptotake(tptotake<1|tptotake>length(RotarySignal))=[];
    RotarySignalPerSec(tp) = nanmean(diff(RotarySignal(uint16(tptotake))));
end
RotarySignalPerSec = (RotarySignalPerSec-nanmin(RotarySignalPerSec(:)))./(nanmax(RotarySignalPerSec(:))-nanmin(RotarySignalPerSec(:)));
RotarySignalPerSec = downsample(RotarySignalPerSec,15);
tmptime = downsample(Actualtime,15);
newsr = nanmedian(diff(tmptime));

% Prepare neural activity per unit according to StimulusTimes X nROI
dFFPerTP = nan(length(Stimulustimes),nROI);
parfor tp = 1:length(Stimulustimes)-1
    dFFPerTP(tp,:) = nanmean(DataMatrix(Actualtime>=Stimulustimes(tp)&Actualtime<Stimulustimes(tp+1),:),1);    
end

nX = size(stimFrames,2);
nY = size(stimFrames,3);
xPos = 1:nX;
yPos = 1:nY;
theta = linspace(0,2*pi,100);

if fancymethod % Fancy one but takes forever
    % Fancy  RF code - borrowed from https://github.com/sylviaschroeder/schroeder-et-al-2020/blob/master/scripts/main_mapReceptiveFields.m
    % for receptive field estimates
    % used for fitting 2 RFs (ON and OFF simultaneously), and fitting running
    % kernels and RFs simultaneously
    % lambdasStim = logspace(-4, 1, 6);
    % lambdasRun = logspace(0, 6, 7);
    lambdasStim = logspace(-4, -3, 2);
    lambdasRun = logspace(0, 1, 2);
    RFlimits = [0 2]; %This is the possible on/offset response time window
    % crossFolds = 10;
    crossFolds = 5;
    
    % parameters for running speed as predictor
    runKrnlLimits = [-5 5];
    RFtimesInFrames = floor(RFlimits(1) / stimFrameDur) : ...
        ceil(RFlimits(2) / stimFrameDur);
    NotNanIdx =find(sum(isnan(dFFPerTP),1)<0.8*size(dFFPerTP,1)); %exclude NAN
    
   
    % map RF
    [rFields, runKernels, runWin, ev, ev_run, ev_stim] = ...
        whiteNoise.getReceptiveField( ...
        dFFPerTP(:,5), Stimulustimes, stimFrames, Stimulustimes, ...
        RFtimesInFrames, RotarySignalPerSec, tmptime, runKrnlLimits, ...
        {lambdasStim, lambdasRun}, crossFolds);
    
    %
    v = squeeze(mean(ev,2)); % [neuron x lamStim x lamRun], average across cross-folds
    [maxEV, maxStimLam] = max(v,[],2);
    maxEV = squeeze(maxEV); % [neuron x lamRun];
    maxStimLam = squeeze(maxStimLam); % [neuron x lamRun];
    [maxEV, maxRunLam] = max(maxEV, [], 2); % [neuron x 1]
    indLam = sub2ind(size(maxStimLam), (1:size(maxStimLam,1))', maxRunLam);
    maxStimLam = maxStimLam(indLam); % [neuron x 1]
    
    vRun = squeeze(mean(ev_run,2)); % [neuron x lamStim x lamRun], average across cross-folds
    vStim = squeeze(mean(ev_stim,2)); % [neuron x lamStim x lamRun]
    inds = sub2ind(size(vRun), (1:size(vRun,1))', maxStimLam, maxRunLam);
    maxEVRun = vRun(inds); % [neuron x 1]
    maxEVStim = vStim(inds); % [neuron x 1]
    
    % test signficance of each RF - takes forever never mind
    % for now
    if 0
        [ev, ev_shift] = ...
            whiteNoise.receptiveFieldShiftTest( ...
            dFFPerTP(:,NotNanIdx), Stimulustimes, stimFrames, Stimulustimes, ...
            RFtimesInFrames, RotarySignalPerSec, tmptime, ...
            runKernels, runWin, rFields, maxStimLam, 50);
        pvals = sum(ev_shift > ev, 2) ./ size(ev_shift,2);
        pvals(isnan(ev)) = NaN;
    end
    
    % receptiveFields     [rows x cols x RFframes x RFtype x neuron]
    %                       containing linear regression solution for x in Ax=B
    %                       where A is stimulus [rows x cols x time] and B is
    %                       calcium response, for each neuron and stimulus
    %                       model; ridge regression is performed on all data
    %                       using the optimal lambda value found with
    %                       cross-validation
    %   runKernels          [t x neuron]; containing linear
    %                       regression kernel fitting calcium based on running
    %                       speed
    %   runWin              [1 x t]; time of run kernel relative to neural
    %                       response
    % Save Preprocessed data
    
    % Plot response of neurons to RunKernels
    areashere = Depth2AreaPerUnit.Area(ismember(Depth2AreaPerUnit.Cluster_ID,Good_IDx));
    areasheretmp = areashere(NotNanIdx);
    Good_IDtmp = Good_ID(NotNanIdx);
    % Fit 2D Gaussian
    [thisparams,thisgaussian,bestbeta,coeff,pvalcorr] = fit2dGaussRF(permute(squeeze(nanmax(rFields(:,:,:,1,:),[],3)),[2,1,3]));
    RFInfo = table;
    RFInfo.ROIID = NotNanIdx';
    RFInfo.params = nan(nROI,3);
    RFInfo.params(NotNanIdx,:) = thisparams';
    RFInfo.bestbeta  = nan(nROI,1);
    RFInfo.bestbeta(NotNanIdx) = bestbeta;
    RFInfo.Coeff = nan(nROI,1);
    RFInfo.Coeff(NotNanIdx) = coeff;
    RFInfo.Pval = nan(nROI,1);
    RFInfo.Pval(NotNanIdx) = pvalcorr;
    RFInfo.ExplainedVarianceRun = nan(nROI,1);
    RFInfo.ExplainedVarianceRun(NotNanIdx) = maxEVRun;
    %                 RFInfo.ExplainedVarianceShift = nan(nclus,1);
    %                 RFInfo.ExplainedVarianceShift(NotNanIdx) = ev_shift;
    RFInfo.ExplainedVarianceStim = nan(nROI,1);
    RFInfo.ExplainedVarianceStim(NotNanIdx) = maxEVStim;
    RFInfo.RunData = nan(nROI,size(runKernels,1));
    RFInfo.RunData(NotNanIdx,:) = runKernels';
    
    [~,sortid]=sort(maxEVStim,'descend');
    sortidnans = find(isnan(maxEVStim(sortid)));
    sortid = [sortid; sortid(sortidnans)];
    sortid(sortidnans)=[];
    
    figure;
    countid=1;
    for clusid=1:nROI
        subplot(5,2,countid)
        
        %Running
        h=plot(runWin,runKernels(:,sortid(clusid)),'k');
        xlabel('Time (s)')
        ylabel('Running fit')
        box off
        title([areasheretmp{sortid(clusid)} ' ,ID:' num2str(Good_IDtmp(sortid(clusid))) ', R2= ' num2str(maxEVRun(sortid(clusid))*100)])
        countid=countid+1;
        
        %Receptive field
        subplot(5,2,countid)
        h = imagesc(nanmax(rFields(:,:,:,1,sortid(clusid)),[],3));
        title(['x: ' num2str(thisparams(1,sortid(clusid))) ', y: ' num2str(thisparams(2,sortid(clusid))) ',fwhm: ' num2str(thisparams(3,sortid(clusid))) ', p='  num2str(pvalcorr(sortid(clusid))) ', R2= ' num2str(round(maxEVStim(sortid(clusid))*1000)/10)])
        colormap gray
        freezeColors
        hold on
        set(gca,'YDir','normal')
        plot(thisparams(1,sortid(clusid))+  thisparams(3,sortid(clusid))*cos(theta), thisparams(2,sortid(clusid))+ thisparams(3,sortid(clusid))*sin(theta),'r-','LineWidth',3)
        
        countid=countid+1;
        
        if countid>10
            saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'RF_FitExamples.fig'))
            saveas(gcf,fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'RF_FitExamples.bmp'))
            
            break
            figure;
            countid=1;
        end
        drawnow
    end
    
    
    %% Save Preprocessed data
    save(fullfile(SaveDir,MiceOpt{midx},thisdate,thisses,thisprobe,'SpikeData.mat'),'RFInfo','runWin','Depth2AreaPerUnit','clusinfo','rFields','-v7.3')
    
else
    
    %                 Images = SS{trid}.ImageTextures;
    
    %Take window
    windowaroundframe = [0 2];%0 to 0.5 seconds surrounding frame on
    
    rfMapOn = zeros(nX*nY,nROI);
    rfMapAll = rfMapOn;
    rfMapOff= rfMapOn;
    
    stimFrames = reshape(stimFrames,[],nX*nY);
    greyval = nanmedian(unique(stimFrames(:)));
    whiteval = nanmax(unique(stimFrames(:)));
    blackval = nanmin(unique(stimFrames(:)));
    parfor rfid = 1:nX*nY
        %All
        takeframes = find(stimFrames(:,rfid)~=greyval);
        takeframes(logical([0; diff(takeframes)==1]))=[]; %Remove multiple frames (so only take onsets)
        taketimepoints = Stimulustimes(takeframes);
        if isempty(takeframes) || sum(takeframes)==0
            warning('sparseNoiseRF: stimulus at %d was never shown', rfid);
        else
            tmp = arrayfun(@(X,Y) nanmean(dFFPerTP(Stimulustimes>=X&Stimulustimes<=Y,:),1),taketimepoints+windowaroundframe(1),taketimepoints+windowaroundframe(2),'UniformOutput',false);
            tmp = cat(1,tmp{:});
            rfMapAll(rfid,:) =  nanmean(tmp,1);
        end
        
        %ON RFs
        takeframes = find(stimFrames(:,rfid)==whiteval);
        takeframes(logical([0; diff(takeframes)==1]))=[]; %Remove multiple frames (so only take onsets)
        taketimepoints = Stimulustimes(takeframes);
        if isempty(takeframes) || sum(takeframes)==0
            warning('sparseNoiseRF: stimulus at %d was never shown', rfid);
        else
            tmp = arrayfun(@(X,Y) nanmean(dFFPerTP(Stimulustimes>=X&Stimulustimes<=Y,:),1),taketimepoints+windowaroundframe(1),taketimepoints+windowaroundframe(2),'UniformOutput',false);
            tmp = cat(1,tmp{:});
            rfMapOn(rfid,:) =  nanmean(tmp,1);
        end
        
        %off RFs
        takeframes = find(stimFrames(:,rfid)==blackval);
        takeframes(logical([0; diff(takeframes)==1]))=[]; %Remove multiple frames (so only take onsets)
        taketimepoints = Stimulustimes(takeframes);
        if isempty(takeframes) || sum(takeframes)==0
            warning('sparseNoiseRF: stimulus %d was never shown', rfid);
        else
            tmp = arrayfun(@(X,Y) nanmean(dFFPerTP(Stimulustimes>=X&Stimulustimes<=Y,:),1),taketimepoints+windowaroundframe(1),taketimepoints+windowaroundframe(2),'UniformOutput',false);
            tmp = cat(1,tmp{:});
            rfMapOff(rfid,:) =   nanmean(tmp,1);
        end
        
    end
    
    % statistical test(s)
    % - shuffle test: idea is that if you relabel each stimulus event with a different
    % position, on what proportion of relabelings do you get a peak as big as
    % the one actually observed? or can calculate this analytically from the
    % distribution of all spike counts.
    % - simplest: just the z-score of the peak relative to the whole population
    maxZ = (nanmax(rfMapAll,[],1)-nanmean(rfMapAll,1))./nanstd(rfMapAll,[],1);
    minZ = (nanmin(rfMapAll,[],1)-nanmean(rfMapAll,1))./nanstd(rfMapAll,[],1);
    peakZscore = max(abs([minZ; maxZ]),[],1);
    
    % Remove Nan clusters
    NotNanIdx=find(sum(isnan(rfMapAll),1)<size(rfMapAll,1));
    rfMapAll = rfMapAll(:,NotNanIdx);
    rfMapAll = reshape(rfMapAll,nX,nY,[]);
    parfor clusid=1:size(rfMapAll,3)
        rfMapAll(:,:,clusid) = smooth2a(rfMapAll(:,:,clusid),2);
    end
    
    % Fit 2D Gaussian
    [thisparams,thisgaussian,bestbeta,coeff,pvalcorr] = fit2dGaussRF(rfMapAll);
    theta = linspace(0,2*pi,100);
    %%
    figure;
    countid=1;
    figurecount = 1;
    for clusid=1:length(NotNanIdx)
        subplot(3,3,countid)
        tmp = rfMapAll(:,:,clusid);
        lims = [quantile(tmp(:),0.05) quantile(tmp(:),0.95)];
        try
            imagesc(tmp,lims)
            hold on
            set(gca,'YDir','normal')
            colormap gray
            freezeColors
            
            plot(thisparams(1,clusid)+thisparams(3,clusid)*cos(theta),thisparams(2,clusid)+thisparams(3,clusid)*sin(theta),'r-')
            title(['ROI ' num2str(clusid) ', RF [' num2str(thisparams(1,clusid)) ';' num2str(thisparams(2,clusid)) '], r=' num2str(round(thisparams(3,clusid)*100)/100), 'beta = ' num2str(bestbeta(clusid))])
            countid=countid+1;
        catch
            continue
        end
        
        if countid>9
            countid=1;
            saveas(gcf,fullfile(SaveDir,['RFMaps_Set' num2str(figurecount) DataType '.fig']))
            figure
            figurecount = figurecount+1;
        end
        drawnow
    end
    
    RFInfo = table;
    RFInfo.ROIID = NotNanIdx';
    RFInfo.params = nan(nROI,3);
    RFInfo.params(NotNanIdx,:) = thisparams';
    RFInfo.bestbeta  = nan(nROI,1);
    RFInfo.bestbeta(NotNanIdx) = bestbeta;
    RFInfo.Coeff = nan(nROI,1);
    RFInfo.Coeff(NotNanIdx) = coeff;
    RFInfo.Pval = nan(nROI,1);
    RFInfo.Pval(NotNanIdx) = pvalcorr;
    RFInfo.PeakZScore = peakZscore';
    %% Save Preprocessed data
    save(fullfile(SaveDir,'RFData.mat'),'RFInfo','rfMapAll','rfMapOff','rfMapOn','-v7.3')
    
end



clear sp
clear SpikeRatePerPos
clear SpikeRatePerTP
clear spikeTimesCorrected
clear allP
clear ampBins
clear cdfs
clear clusinfo
clear dataArray
clear FlipperGLX
clear Flippertimeline
clear removevec
clear rewardevents
clear spikeAmps
clear spikeCluster
clear spikeDepths
clear spikeSites
clear spikesThisTrial
clear spikeTimes
clear spikeTimestmp
clear SpikeTrialID
clear SpikeTrialTime
clear tempsUnW
clear Actualtime
clear allPowerEst
clear allPowerVar
clear F
close all
clear starttime
clear h1

%         end
%     end
% end