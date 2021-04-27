close all

%number of trials
numOfTrials=48; %should be 48 in normal cases

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

%total number of microsaccades
totalMicrosacc=0;
%microsaccade rate over all trials
microsaccRate=[];
%mean microsaccade rate
meanMicrosaccRate=0;
%mean microsaccade magnitude
meanMagnitude=0;
%microsaccde magnitude for over all trials
microsaccMagnitude=[];
%microsaccde duration for over all trials
microsaccDuration=[];
%mean microsaccade duration
meanDuration=0;
%peak velocity over all trials
peakVelocity=[];
%mean peak velocity
meanPeakVelocity=0;
%amplitude of the microsaccades
acc_result=[];
%accepted microsaccades
total_result=[];

%collect number of main sequences
main_sequences=[];
for i=2:numOfTrials+1
    perLine=strsplit(lines{i,1},'\t'); %Each line
    if perLine{2}=='1'
        main_sequences=[main_sequences,i-1];
    end
end

address=["./results/mofo2.edf";"./results/mofo3.edf"];
% Read EDF and convert it to mat
experimentPart=0;
for i=1:size(address,1)
    if i==1
        experimentPart=1;
    elseif i==2
        experimentPart=2;
    end
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
    startIndx = cell2mat(start_indx_tmp);
    
    
    %collecting end time of each trial
    %("trial OK" in events.info corresponds to a time in events.time)
    endTrial=events.time(1,find(contains(events.info,'trial OK')));
    
    
    %"find" function used to get the index of an element in Samples.time which
    %equal to an element in endTrial
    end_indx_tmp=arrayfun(@(x) find(x==samples,1),endTrial,'un',0);
    
    
    end_indx=[];
    
    %convert cell array which is output of ('un',0) to a vector
    endIndx = cell2mat(end_indx_tmp);
    
    
    SAMPLING=500; %As it is recommended(sample rate)
    TYPE=2;  %velocity type: TYPE=2 recommended
    
    span=100; %It should be changed to 200 for sampling=1000
    
    for j=1:length(endTrial)
        %check if the trial is main trial
        if (any(main_sequences(:)==j)&(i==1)|any(main_sequences(:)==j+24)&(i==2))
            
            %accepted microsaccades of this trial
            total_result=[];
            %trial time record
            trial_time_rec=myFile.Samples.time(startIndx(j):endIndx(j));
            %mapping of the time of each approved data (before elimination, after elimination)
            event_time_map=[1:length(trial_time_rec)];
            %get posx of both left and right eye during each trial
            tmpx=myFile.Samples.posX(startIndx(j):endIndx(j),:);
            
            %get posy of both left and right eye during each trial
            tmpy=myFile.Samples.posY(startIndx(j):endIndx(j),:);
            
            %get the pupil size of the samples in the interval
            pupilSize = myFile.Samples.pupilSize(startIndx(j):endIndx(j),:);
            
            %pupilSize arr counter(remove samples which the increase and decrease of
            %the pupil size is more than 20 units and 200 samples before and after that)
            k=1;
            areaThreshold=20; %should be 10 for Sampling=1000 and should be 20 for Sampling=500
            while k<=(size(pupilSize,1)-1)
                if(abs(pupilSize(k+1,1)-pupilSize(k,1))>areaThreshold||abs(pupilSize(k+1,2)-pupilSize(k,2))>areaThreshold)
                    %window off size 400 is flexible according to length of the array
                    if k-span>0 && k+span<size(tmpx,1)
                        pupilSize(k-span:k+span,:) = -1;
                        %the last moment before elimination
                        before_eliminate=trial_time_rec(k-span)-trial_time_rec(1);
                        %the first moment after elimination
                        after_eliminate=trial_time_rec(k+span)-trial_time_rec(1);
                        k=k+span+1;
                    elseif k+span>size(tmpx,1)
                        pupilSize(k-span:size(pupilSize,1),:) = -1;
                        before_eliminate=trial_time_rec(k-span)-trial_time_rec(1);
                        after_eliminate=trial_time_rec(length(trial_time_rec))-trial_time_rec(1);
                        k=size(tmpx,1)+1;
                    elseif k-span<0
                        %fprintf("K and size(tmpx): %d\t%d\n",k,size(tmpx));
                        pupilSize(1:k+span,:) = -1;
                        before_eliminate=0;
                        after_eliminate=trial_time_rec(k+span)-trial_time_rec(1);
                        k=k+span+1;
                    end
                    %fprintf("Pupilsize Trial Number %d/ round %d/ Before Eliminate %d/ After Eliminate %d\n",j,k, before_eliminate,after_eliminate);
                    event_time_map(before_eliminate:after_eliminate)= -1;
                else
                    k=k+1;
                end
            end
            
            %tmpx,tmpy counter
            k=1;
            %remove samples in which the data is not recorded due to the
            %subject blink and 200 samples before and after that
            while k<=size(tmpx,1)
                % the window of size 400 is flexible according to length of the
                % array
                if(isnan(tmpx(k,1))||isnan(tmpx(k,2)))
                    if k-span>0 && k+span<size(tmpx,1)
                        tmpx(k-span:k+span,:) = -1;
                        tmpy(k-span:k+span,:) = -1;
                        before_eliminate=trial_time_rec(k-span)-trial_time_rec(1);
                        after_eliminate=trial_time_rec(k+span)-trial_time_rec(1);
                        k=k+span+1;
                    elseif k+span>size(tmpx,1) && k-span>0
                        tmpx(k-span:size(tmpx,1),:) = -1;
                        tmpy(k-span:size(tmpy,1),:) = -1;
                        before_eliminate=trial_time_rec(k-span)-trial_time_rec(1);
                        after_eliminate=trial_time_rec(length(trial_time_rec))-trial_time_rec(1);
                        k=size(tmpx,1)+1;
                    elseif k+span>size(tmpx,1) && k-span<0
                        tmpx(1:size(tmpx,1),:) = -1;
                        tmpy(1:size(tmpy,1),:) = -1;
                        before_eliminate=0;
                        after_eliminate=trial_time_rec(length(trial_time_rec))-trial_time_rec(1);
                        k=size(tmpx,1)+1;
                        break;
                    elseif k-span<0
                        tmpx(1:k+span,:) = -1;
                        tmpy(1:k+span,:) = -1;
                        before_eliminate=0;
                        after_eliminate=trial_time_rec(k+span)-trial_time_rec(1);
                        k=k+span+1;
                    end
                    
                    %fprintf("Blink Trial Number %d/ round %d/ Before Eliminate %d/ After Eliminate %d\n",j,k, before_eliminate,after_eliminate);
                    event_time_map(before_eliminate:after_eliminate)= -1;
                else
                    k=k+1;
                end
            end
            
            k=1;
            while k<size(event_time_map,2)
                if (event_time_map(k)==-1)
                    %fprintf("time: %d\n",k);
                    event_time_map(k)=[];
                else
                    k=k+1;
                end
            end
            
            k=1;
            %removing rows whith values -1 in tmpx,tmpy, and pupilSize
            while k<=size(tmpx,1)
                if (tmpx(k,:)==-1 | tmpy(k,:)==-1 | pupilSize(k,:)==-1)
                    tmpx(k,:)=[];
                    tmpy(k,:)=[];
                    pupilSize(k,:)=[];
                else
                    k=k+1;
                end
            end
            
            %This variable is used to calculate miccrosaccade rate (rate=Num/window size(s))
            trialMicrosacc=0;
            
            if(~isempty(tmpx) && ~isempty(tmpy))
                %Engbert algorithm
                if(size(tmpx,1)>=6&&size(tmpy,1)>=6)
                    xx_left=[tmpx(:,1),tmpy(:,1)];
                    v1 = vecvel(xx_left,SAMPLING,TYPE);
                    
                    xx_right=[tmpx(:,2),tmpy(:,2)];
                    v2 = vecvel(xx_right,SAMPLING,TYPE);
                    
                    VFAC=3;%it was 6 in the main paper [Engbert];
                    MINDUR=3; %it was 3 in the main paper [Engbert];
                    [sac_left, radius_left] = microsacc(xx_left,v1,VFAC,MINDUR);
                    [sac_right, radius_right]=microsacc(xx_right,v2,VFAC,MINDUR);
                    %sac is binocular microsaccades
                    [sac, monol, monor] = binsacc(sac_left,sac_right);
                    
                    result=saccpar(sac);
                    
                    %overshoot correction (remove saccades which the time interval between the offset
                    %of first one and onset of second one is less than 20ms)
                    l=2; %variable for traversing the result array
                    while l<=size(result,1)
                        if((result(l,1))-result((l-1),2)<20)
                            result(l,:)=[];
                        else
                            l=l+1;
                        end
                    end
                    
                    %removing microsaccades with magnitude greater than
                    %2degress
                    l=1; %variable for traversing the result array (one degree is 36 pixels in our config)
                    while l<=size(result,1)
                        if (abs(result(l,8))/36) >2 %maximum threshold for microssacade mag is 2degrees
                            %fprintf("the magnitude: %d\n",abs(result(l,8)));
                            result(l,:)=[];
                        else
                            l=l+1;
                        end
                        
                    end
                    
                    
                    %adding to the number of microsaccades
                    trialMicrosacc=trialMicrosacc+size(result,1);
                    total_result=[total_result;result];
                    
                    if ~isempty(result)
                        
                        % collect the total microsacc info for all trials
                        total_result=[total_result;result];
                        % get the amplitude for the trial
                        acc_result=[acc_result;result(:,8)/36];
                        
                        %                     trialRate=0;
                        %                     for l=1:size(result,1)
                        %                         trialRate=trialRate+(1/(result(l,2)-result(l,1)));
                        %                     end
                        %                     microsaccRate=[microsaccRate,size(result,1)/30]; %(trialRate/size(result,1))
                        
%                         %adding microsaccade magnitude
%                         trialMag=0;
%                         for l=1:size(result,1)
%                             trialMag=trialMag+(abs(result(l,8))/36);
%                         end
%                         microsaccMagnitude=[microsaccMagnitude,(trialMag/size(result,1))];
                        
                        %adding microsaccade duration
                        trialDur=0;
                        for l=1:size(result,1)
                            trialDur=trialDur+result(l,3);
                        end
                        microsaccDuration=[microsaccDuration,(trialDur/size(result,1))];
                        
                        %adding microsaccade peak velocity
                        trialPeak=0;
                        for l=1:size(result,1)
                            trialPeak=trialPeak+(result(l,5)/36);
                        end
                        peakVelocity=[peakVelocity,(trialPeak/size(result,1))];
                        
                    end
                end
            end
            
            totalMicrosacc=totalMicrosacc+trialMicrosacc;
            %microsaccade rate in this trial
            microsaccRate=[microsaccRate,trialMicrosacc/30];
        end
    end
end

meanMicrosaccRate=mean(microsaccRate);
microsaccMagnitude=total_result(:,8)/36;
meanDuration=mean(microsaccDuration);
meanPeakVelocity=mean(peakVelocity);


%adding to finalResult
%finalResult=[finalResult;totalMicrosacc;meanMicrosaccRate;meanMagnitude;meanDuration;meanPeakVelocity];

%duration of rotation periods
rotationTime=[];
%duration of stationary periods
stationaryTime=[];
%Number of transitions to rotation
numOfRotation=0;
%Number of transitions to stationary
numOfStationary=0;

durationOfFirstRotatingPeriod=0;

%number of accepted trials trials
count=0;

for i=2:numOfTrials+1
    perLine=strsplit(lines{i,1},'\t'); %Each line
    
    if size(perLine,2)>=6 & perLine{6}~='*'
        count=count+1;
        %number of transitions to rotation
        tmpRelease=strsplit(perLine{6},' ');
        tmpRelease(end) = [];
        
        userRelease=str2double(tmpRelease); %user release time
        
        for j=1:size(userRelease,2)
            if userRelease(j)>30
                userRelease(j) = [];
            end
        end
        
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
        numOfRotation=numOfRotation+1;
        
        % first rotation periods
        durationOfFirstRotatingPeriod=durationOfFirstRotatingPeriod+firstRotationDuration;
        
        %rotation duration in the current trial
        for j=2:size(userPress,2)
            if (j-1)<=size(userRelease,2)
                rotationTime=[rotationTime,userPress(j)-userRelease(j-1)];
            end
        end
        %last rotation time til the end of the trial
        if size(userRelease,2)==size(userPress,2)
            rotationTime=[rotationTime,30-userRelease(size(userRelease,2))];
        end
        
        %stationary duration in the current trial
        for j=1:size(userRelease,2)
            stationaryTime=[stationaryTime,userRelease(j)-userPress(j)];
        end
        if size(userPress,2)>size(userRelease,2)
            stationaryTime=[stationaryTime,30-userPress(size(userPress,2))];
        end
        %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PLOT E
        figure(4);
        
        % Each trial lasts for 30 seconds
        x=[1,30];
        % a timeline for each trial
        y=[i-1,i-1];
        
        % The whole timeline is blue
        plot(x,y,'Color',[0.3010 0.7450 0.9330],'LineWidth', 2);
        %set(gca,'Color','k')
        hold on;
        
        % Stationary moments are red
        copy_release=str2double(tmpRelease);
        
        %length(userPress) is equal to length(copy_release)
        for j=1:length(userPress)
            x=[userPress(j),copy_release(j)];
            plot(x,[i-1,i-1],'Color',[0.6350 0.0780 0.1840],'LineWidth', 2);
            %set(gca,'Color','k')
        end
        grid on;
    end
    set(gca,'Color','w')
    legend('Rotating state','Stationary state')
    xlabel('Time(s)');
    ylabel('Number of Trials');
    
end
hold off;

%average over first rotation periods
durationOfFirstRotatingPeriod=durationOfFirstRotatingPeriod/count;

%Time spent in rotating periods(in percent)
rotatingPercent=(sum(rotationTime)/(sum(rotationTime)+sum(stationaryTime)))*100;

%Time spent in rotating periods(in percent)
stationaryPercent=(sum(stationaryTime)/(sum(rotationTime)+sum(stationaryTime)))*100;

%adding data to final result
%finalResult=[finalResult;numOfRotation;mean(rotationTime);durationOfFirstRotatingPeriod;rotatingPercent;numOfStationary;mean(stationaryTime);stationaryPercent];
%disp(finalResult);



%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PLOT A
figure(1);
hist(microsaccMagnitude,150);
xlabel('Microsaccade magnitude(deg)');
ylabel('Number of microsaccades');


%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PLOT B
ample=total_result(:,8)/36; %conversion from pixel to degree visual angle
% graph B
figure(2);
vpeak=total_result(:,5)/36; %conversion from pixel to degree visual angle
disp(size(total_result));
% Next three lines is to plot the fit line
coefficients = polyfit(ample, vpeak, 1);
xFit = linspace(min(ample), max(ample), 1000);
yFit = polyval(coefficients , xFit);
hold on;
scatter(ample,vpeak,'.');
numberOfGoodPoints = sum(~isnan(ample) & ~isnan(vpeak));
plot(xFit, yFit,'b');
grid on;
xl = xlim;
yl = ylim;
xt = 0.05 * (xl(2)-xl(1)) + xl(1);
yt = 0.90 * (yl(2)-yl(1)) + yl(1);
caption=sprintf('N= %d \nAvg. Velocity= %f, +/- %f\nAvg. Magnitude= %f, +/- %f\nMicro Rate= %f/sec\nFit: y = %f * x + %f', numberOfGoodPoints, mean(vpeak),std(vpeak)/sqrt(length(vpeak)),mean(ample), std(ample)/sqrt(length(ample)),meanMicrosaccRate,coefficients(1), coefficients(2));
text(xt, yt, caption, 'FontSize', 10, 'Color', 'black', 'FontWeight', 'bold');
xlabel('Microsaccade magnitude(deg)');
ylabel('Microsaccade peak velocity (deg/s)');

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PLOT C
%intersaccadic intervals(ms)
interval=[];
for i=2:size(total_result,1)
    if(total_result(i,1)>total_result(i-1,2))
        interval=[interval;total_result(i,1)-total_result(i-1,2)];
    end
end
figure(3);
histogram(interval,'BinWidth',10,'BinLimits',[0 1000]);
xlabel('Intersaccadic intervals (ms)');
ylabel('Number of intersaccadic intervals');

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PLOT D
%%change horizontal array into vertical
%stationary=stationaryTime.';

figure(5);
hi=histogram(stationaryTime,20);
hold on;
% overlay a curve on the histogram (connect the midpoint of bins to each other)
plot(conv(hi.BinEdges, [0.5 0.5], 'valid'), hi.BinCounts,'b', 'LineWidth', 1);

% The duration when the subject perceives disks as rotating
hi2=histogram(rotationTime,10);
hold on
% overlay a curve on the histogram (connect the midpoint of bins to each other)
plot(conv(hi2.BinEdges, [0.5 0.5], 'valid'), hi2.BinCounts, 'r', 'LineWidth', 1);
legend('Stationary state bins','Stationary state curve','Rotating state bins','Rotating state curve')
xlabel('Percept duration (sec)');
ylabel('Number of percepts');
hold off;



% rng('default');
% figure(5);
% hi=histfit(rotationTime,50,'kernel');
% set(gca, 'XLim', [0 30]);
%
% hold on;
% hi(1).FaceColor = [.8 .9 1];
% hi(1).FaceAlpha = 0.5;
% hi(1).EdgeColor = [0 0 1];
% hi(2).Color = [0 0 1];
%
% %mu=mean(rotationTime);
% %sigma=std(rotationTime);
% %plot([mu, mu], ylim,'--','Color', '[0 0 1]', 'LineWidth', 2);
% %hold on;
% %plot([mu + sigma, mu + sigma], ylim, 'Color', '[0 1 0.8]', 'LineWidth', 1.5);
% %hold on;
% %plot([mu - sigma, mu - sigma], ylim, 'Color', '[0 1 0.8]', 'LineWidth', 1.5);
%
% hold on;
% hi2=histfit(stationaryTime,50,'kernel');
% set(gca, 'XLim', [0 30]);
%
% hold on;
% hi2(1).FaceColor = [1 .9 .9];
% hi2(1).FaceAlpha = 0.5;
% hi2(1).EdgeColor = [1 0 0];
% hi2(2).Color = [1 0 0];
%
% %mu2=mean(staionaryTime);
% %sigma2=std(stationaryTime);
% %plot([mu2, mu2], ylim,'--','Color', '[1 0 0]', 'LineWidth', 2);
% %hold on;
% %plot([mu2 + sigma2, mu2 + sigma2], ylim, 'Color', '[1 0 1]', 'LineWidth', 1.5);
% %hold on;
% %plot([mu2 - sigma2, mu2 - sigma2], ylim, 'Color', '[1 0 1]', 'LineWidth', 1.5);
%
% grid on;
% % xlim([0 1]);
% plots=get(gca, 'Children');
% legend([plots(2, 1) plots(4, 1)], {'Stationary', 'Rotation'});
% %title('Dashed lines are the means and side lines are -/+ SD');
% xlabel('Percept duration (sec)');
% ylabel('Number of percepts');
% hold off;
