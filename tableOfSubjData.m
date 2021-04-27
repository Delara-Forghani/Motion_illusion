close all

%number of trials
numOfTrials=48;

%Read trial details from file
filename = './results/mofo2.txt';
delimiterIn = '\t';
%number of lines to read
headerlinesIn = numOfTrials+1;
%get the lines from file
lines = importdata(filename,delimiterIn,headerlinesIn);

%result data for this subject (1.Number of microsaccades,2.microsaccade rate,3.microsaccade magnitude
%4.microsaccade duration,5.microsaccade peak velocity,6.number of transition to rotation,7.duration of rotating periods
%8.duration of first rotating period after stimulus onset, 9.time spent in rotating periods,
%10.number of transitions to stationary, 11.duration of stationary periods,12.time spent in stationary periods)
finalResult=[];

startIndx=[];
endIndx=[];

address=["./results/mofo2.edf";"./results/mofo3.edf"]
% Read EDF and convert it to mat
for i=1:size(address,1)
    edfFile=[char(address(i,1))];
    myFile=Edf2Mat(edfFile);

    events=myFile.Events.Messages;   %find events of eyelink

    %contains events.time correspond to TRIALID n
    %According to number of trials sometimes instead of 24 the array size should be 12 or 8(for each sub-trial)
    startTrial=zeros(1,24);
    %contains events.time correspond to trial OK
    %According to number of trials sometimes instead of 24 the array size should be 12 or 8(for each sub-trial)
    endTrial=zeros(1,24);

    samples=myFile.Samples.time;  %get the time each sample was recorded

    %collecting start time of each trial
    %("TRIALID n" in events.info corresponds to a time in events.time)
    startTrial=events.time(1,find(contains(events.info,'TRIALID')));

    %"find" function used to get the index of an element in Samples.time which
    %equal to an element in startTrial
    start_indx_tmp=arrayfun(@(x) find(x==samples,1),startTrial,'un',0);


    start_indx=[];

    %press_indx=cell2mat(press_indx_tmp);
    % convert cell array which is output of ('un',0) to a vector
    start_indx = cell2mat(start_indx_tmp);
    
    startIndx=[startIndx,start_indx];
    %collecting end time of each trial
    %("trial OK" in events.info corresponds to a time in events.time)
    endTrial=events.time(1,find(contains(events.info,'trial OK')));


    %"find" function used to get the index of an element in Samples.time which
    %equal to an element in endTrial
    end_indx_tmp=arrayfun(@(x) find(x==samples,1),endTrial,'un',0);


    end_indx=[];

    %convert cell array which is output of ('un',0) to a vector
    end_indx = cell2mat(end_indx_tmp);
    endIndx=[endIndx,end_indx];
end

SAMPLING=500; %As it is recommended(sample rate)
TYPE=2;  %velocity type: TYPE=2 recommended


%total number of microsaccades
totalMicrosacc=0;
%total microsaccade rate
microsaccRate=[];
%mean microsaccade rate
meanMicrosaccRate=0;
%mean microsaccade magnitude
meanMagnitude=0;
%microsaccde magnitude for each trial
microsaccMagnitude=[];
%microsaccde duration for each trial
microsaccDuration=[];
%mean microsaccade duration
meanDuration=0;
%peak velocity for each trial
peakVelocity=[];
%mean peak velocity
meanPeakVelocity=0;

for i=1:length(endTrial)
    
    %get posx of both left and right eye during each trial
    tmpx=myFile.Samples.posX(startIndx(i):endIndx(i),:);
    
    %get posy of both left and right eye during each trial
    tmpy=myFile.Samples.posY(startIndx(i):endIndx(i),:);
    
    %tmpx,tmpy counter
    k=1;
    
    %This variable is used to calculate miccrosaccade rate (rate=Num/window size(s))
    trialMicrosacc=0;
    
    %trial microsaccade rate
    trialRate=0;
    
    while k<=size(tmpx,1)
        notNanX=[];
        notNanY=[];
        %accumulates indices with numeral value(not nan)
        while (~isnan(tmpx(k,1))&&~isnan(tmpx(k,2)))
            notNanX=[notNanX;tmpx(k,:)];
            notNanY=[notNanY;tmpy(k,:)];
            k=k+1;
            if(k>size(tmpx,1))
                break
            end
        end
        if(~isempty(notNanX) && ~isempty(notNanY))
            %Engbert algorithm
            if(size(notNanX,1)>=10&&size(notNanY,1)>=10)
                xx_left=[notNanX(:,1),notNanY(:,1)];
                v1 = vecvel(xx_left,SAMPLING,TYPE);
                
                xx_right=[notNanX(:,2),notNanY(:,2)];
                v2 = vecvel(xx_right,SAMPLING,TYPE);
                
                VFAC=6;%2;
                MINDUR=6; %2;
                [sac_left, radius_left] = microsacc(xx_left,v1,VFAC,MINDUR);
                [sac_right, radius_right]=microsacc(xx_right,v2,VFAC,MINDUR);
                %sac is binocular microsaccades
                [sac, monol, monor] = binsacc(sac_left,sac_right);
                
                result=saccpar(sac);
                
                %adding to the number of microsaccades
                trialMicrosacc=trialMicrosacc+size(sac,1);
                
                if ~isempty(result)
                    %microsaccade rate in this trial
                    trialRate=0;
                    for l=1:size(result,1)
                        trialRate=trialRate+(1/(result(l,2)-result(l,1)));
                    end
                    microsaccRate=[microsaccRate,(trialRate/size(result,1))];
                    
                    %adding microsaccade magnitude
                    trialMag=0;
                    for l=1:size(result,1)
                        trialMag=trialMag+(result(l,8)/(result(l,2)-result(l,1)));
                    end
                    microsaccMagnitude=[microsaccMagnitude,(trialMag/size(result,1))];
                    
                    %adding microsaccade duration
                    trialDur=0;
                    for l=1:size(result,1)
                        trialDur=trialDur+result(l,3);
                    end
                    microsaccDuration=[microsaccDuration,(trialDur/size(result,1))];
                    
                    %adding microsaccade peak velocity
                    trialPeak=0;
                    for l=1:size(result,1)
                        trialPeak=trialPeak+result(l,5);
                    end
                    peakVelocity=[peakVelocity,(trialPeak/size(result,1))];
                    
                end
            end
            
        end
        if(k>size(tmpx,1))
            break
        end
        while (isnan(tmpx(k,1))||isnan(tmpx(k,2)))
            if k+1 <=size(tmpx,1)
                k=k+1;
            else
                k=k+1;
                break
            end
        end
    end
    
    totalMicrosacc=totalMicrosacc+trialMicrosacc;
    meanMicrosaccRate=mean(microsaccRate);
    meanMagnitude=mean(microsaccMagnitude);
    meanDuration=mean(microsaccDuration);
    meanPeakVelocity=mean(peakVelocity);
end

%adding to finalResult
finalResult=[finalResult;totalMicrosacc;meanMicrosaccRate;meanMagnitude;meanDuration;meanPeakVelocity];

%duration of rotation periods
rotationTime=[];
%duration of stationary periods
stationaryTime=[];
%Number of transitions to rotation
numOfRotation=0;
%Number of transitions to stationary
numOfStationary=0;

durationOfFirstRotatingPeriod=0;
%find main condition trials
for i=2:numOfTrials+1
    perLine=strsplit(lines{i,1},'\t'); %Each line
    
    if size(perLine,2)>=6 & perLine{6}~='*'
        %number of transitions to rotation
        tmpRelease=strsplit(perLine{6},' ');
        tmpRelease(end) = [];
        userRelease=str2double(tmpRelease); %user release time
        %adding numbers of transitions to rotation of the current trial
        numOfRotation=numOfRotation+size(userRelease,2);
        
        %number of transitions to stationary
        tmpPress=strsplit(perLine{5},' ');
        tmpPress(end) = [];
        userPress=str2double(tmpPress); % user press time
        
        %adding numbers of transitions to stationary of the current trial
        numOfStationary=numOfStationary+size(userPress,2);
        
        %first rotation duration from onset before pressing
        firstRotationDuration=userPress(1);
        rotationTime=[rotationTime,firstRotationDuration];
        
        % first rotation periods
        durationOfFirstRotatingPeriod=durationOfFirstRotatingPeriod+firstRotationDuration;
        
        %rotation duration in the current trial
        for j=2:size(userPress,2)
            rotationTime=[rotationTime,userPress(j)-userRelease(j-1)];
        end
        %last rotation time til the end of the trial
        if userRelease(size(userRelease,2))<30
            rotationTime=[rotationTime,30-userRelease(size(userRelease,2))];
        end
        
        %stationary duration in the current trial
        for j=1:size(userPress,2)
            stationaryTime=[stationaryTime,userRelease(j)-userPress(j)];
        end
    end
end
%average over first rotation periods
durationOfFirstRotatingPeriod=durationOfFirstRotatingPeriod/24;

%Time spent in rotating periods(in percent)
rotatingPercent=(sum(rotationTime)/(sum(rotationTime)+sum(stationaryTime)))*100;

%Time spent in rotating periods(in percent)
stationaryPercent=(sum(stationaryTime)/(sum(rotationTime)+sum(stationaryTime)))*100;

%adding data to final result
finalResult=[finalResult;numOfRotation;mean(rotationTime);durationOfFirstRotatingPeriod;rotatingPercent;numOfStationary;mean(stationaryTime);stationaryPercent];
disp(finalResult);

% save(mat2str(finalResult),'finalResult');