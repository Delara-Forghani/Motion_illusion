function ROC(condition) % 0 as control condition/1 as main condition

close all
clc

%number of trials
numOfTrials=48; %should be 48 in normal cases

%Read trial details from file
filename = './results/mofo2.txt';
delimiterIn = '\t';
%number of lines to read
headerlinesIn = numOfTrials+1;
%get the lines from file
lines = importdata(filename,delimiterIn,headerlinesIn);

%physical or illusory conditions (main or control)
picked_trials=[];   
%considering main conditions only
for i=2:numOfTrials+1
    perLine=strsplit(lines{i,1},'\t'); %Each line
    if str2double(perLine{2})==condition  %main or control condition
        picked_trials=[picked_trials,i-1];
        
    end
end


address=["./results/mofo2.edf";"./results/mofo3.edf"];

all_corr1=[];
all_corr2=[];

all_smoothen1=[];
all_smoothen2=[];

for i=1:size(address,1)
    edfFile=[char(address(i,1))];
    myFile=Edf2Mat(edfFile);
    
    events=myFile.Events.Messages;   %find events of eyelink
    samplesTime=myFile.Samples.time;  %get the time each sample was recorded
    
    %contains events.time correspond to TRIALID n
    %According to number of trials sometimes instead of 24 the array size should be 12 or 8(for each sub-trial)
    startTrial=zeros(1,24);
    %contains events.time correspond to trial OK
    %According to number of trials sometimes instead of 24 the array size should be 12 or 8(for each sub-trial)
    endTrial=zeros(1,24);
    
    %collecting start time of each trial
    %("TRIALID n" in events.info corresponds to a time in events.time)
    startTrial=events.time(1,find(contains(events.info,'TRIALID')));
    
    %"find" function used to get the index of an element in Samples.time which
    %equal to an element in startTrial
    start_indx_tmp=arrayfun(@(x) find(x==samplesTime,1),startTrial,'un',0);
    
    %press_indx=cell2mat(press_indx_tmp);
    % convert cell array which is output of ('un',0) to a vector
    startIndx = cell2mat(start_indx_tmp);
    
    %collecting end time of each trial
    %("trial OK" in events.info corresponds to a time in events.time)
    endTrial=events.time(1,find(contains(events.info,'trial OK')));
    
    %"find" function used to get the index of an element in Samples.time which
    %equal to an element in endTrial
    end_indx_tmp=arrayfun(@(x) find(x==samplesTime,1),endTrial,'un',0);
    
    %convert cell array which is output of ('un',0) to a vector
    endIndx = cell2mat(end_indx_tmp);
    
    
    %collecting start time of each trial
    %("PressTime" in events.info corresponds to a time in events.time)
    press_time=events.time(1,find(contains(events.info,'PressTime')|contains(events.info,'pressTime')));
    press_indx_tmp=arrayfun(@(x) find(x==samplesTime,1),press_time,'un',0);
    pressIndx = cell2mat(press_indx_tmp);
    
    
    %collecting start time of each trial
    %("ReleaseTime" in events.info corresponds to a time in events.time)
    release_time=events.time(1,find(contains(events.info,'ReleaseTime')));
    release_indx_tmp=arrayfun(@(x) find(x==samplesTime,1),release_time,'un',0);
    releaseIndx = cell2mat(release_indx_tmp);
    
    
    %load matlab files containing acceptable data for calculating
    %the microsaccade rates
    final_result=load('./matlab_files/mofo_final_microsaccades.mat', 'final_result');
    %collect all microsaccades in all trials
    final_result_cell = struct2cell(final_result);
    %each index contains microsaccades in a trial
    all_microsacc=final_result_cell{1,1};
    
    %loop over all of the trials
    for j=1:length(endTrial)
        % check trials conditions in both files and in the text file
        if (any(picked_trials(:)==j)&(i==1)|any(picked_trials(:)==j+24)&(i==2))
            j
            trialTime=samplesTime(startIndx(j):endIndx(j)); %time samples for the current trial
            
            tmpArr=pressIndx>startIndx(j) & pressIndx<endIndx(j);
            trialPress=samplesTime(pressIndx(tmpArr)); %time samples for press events in current trial
            
            tmpArr=releaseIndx>startIndx(j) & releaseIndx<endIndx(j);
            trialRelease=samplesTime(releaseIndx(tmpArr)); %time samples for release events in current trial
            
            events=[trialPress;trialRelease];
            startInterval=events-2000;
            endInterval=events+1000;
            eventInterval=[startInterval endInterval];
            
            trialMicrosacc=all_microsacc{j};
            microsaccStartTime=trialMicrosacc(:,1); %start time of microsaccades
            
            %loop over all events in a trial
            for z=1:size(eventInterval,1)
                windowTime=[eventInterval(z,1):eventInterval(z,2)]; %<< Padding should be symmetric
                window=zeros(length(windowTime),3);
                
                %having the start of a button press as 1
                window_press_index=arrayfun(@(x) find(x==windowTime,1),trialPress,'un',0);
                window_press_index = cell2mat(window_press_index(find(~cellfun(@isempty,window_press_index))));
                window(window_press_index,1)=1;
                
                %having the start of a button release as 1
                window_release_index=arrayfun(@(x) find(x==windowTime,1),trialRelease,'un',0);
                window_release_index =cell2mat(window_release_index(find(~cellfun(@isempty,window_release_index))));
                window(window_release_index,2)=1;
                
                %Assigning 1 to microsaccade start times
                microsacc_start_indx=arrayfun(@(x) find(x==[eventInterval(z,1)-trialTime(1):eventInterval(z,2)-trialTime(1)],1), microsaccStartTime,'un',0);
                tmp=[microsacc_start_indx{:}];
                window(tmp,3)=1;
                
                if z<=size(trialPress,1)
                    %correlation between button presses and microsaccades
                    %onsets
                    [correlation_1,lags]=xcorr(window(:,1),window(:,3),1500);
                    rate=correlation_1/length(correlation_1)*500; %sample rate (500 Hz)
                    smoothen=sgolayfilt(rate,1,151);
                    all_corr1=[all_corr1;smoothen'];
                    %size(correlation_1)
                else
                    %correlation between button releases and
                    %microsaccades onset
                    [correlation_2,lags]=xcorr(window(:,2),window(:,3),1500);
                    rate=correlation_2/length(correlation_2)*500; %sample rate (500 Hz)
                    smoothen=sgolayfilt(rate,1,151);
                    all_corr2=[all_corr2;smoothen'];
                    %size(correlation_2)
                end
                
            end
        end
    end
    
end

%average of the microsaccade rates (true positive rates and false posivtive rates)
false_positive_rate=nanmean(all_corr1);
true_positive_rate=nanmean(all_corr2);


%ROC analysis
auc=[];
w=250;
for k=-2000:2:1000
    if k-w>=-2000 & k+w<1000
        rot=true_positive_rate(k-w+2001:k+w+2001);
        stat=false_positive_rate(k-w+2001:k+w+2001);
        lable=zeros(1,length(stat));
        lable(stat>rot)=1;
        try
            [X,Y,T,AUC] = perfcurve(lable,stat,1);
            auc=[auc,AUC];
            auc=[auc,AUC];
        catch
            auc=[auc,0];
            auc=[auc,0];
        end        
    else
        auc=[auc,0];
        auc=[auc,0];
    end
end


figure;
x = -2000:1000 ;

%plot of microsaccade rate around an event
subplot(2,1,1)
plot(x,false_positive_rate,'Color','[0.6350 0.0780 0.1840]') % > press: towards stationary (red)
hold on
plot(x,true_positive_rate,'Color','[0 0.4470 0.7410]') % > release: towards rotation (blue)
xlim([-2000,1000]);
xlabel('Window Time around an event');
ylabel('Microsaccade rate');
%legend('towards stationary','towards rotation');
hold off

subplot(2,1,2)
plot(x,auc(1:end-1),'Color','[0.1 0.7 0.1410]') % > AUC
xlim([-2000,1000]);
xlabel('Window Time around an event');
ylabel('AUC');



