close all
subjects={'mofo2'};

%Stationay and rotation transitions
stationary=[];
rotation=[];
%just for test
count=0;
for i=1:length(subjects)
    %number of trials
    numOfTrials=48;
    %Read trial details from file
    filename = ['./results/',char(subjects(i)),'.txt'];
    delimiterIn = '\t';
    %number of lines to read
    headerlinesIn = numOfTrials+1;
    %get the lines from file
    lines = importdata(filename,delimiterIn,headerlinesIn);
    
    %Exact transition times
    rotateRec=[];
    quitRec=[];
    
    %main loop
    for i=2:numOfTrials+1
        perLine=strsplit(lines{i,1},'\t'); %Each line        
        if perLine{6}~='*' & str2double(perLine{2})==0  %control condition
            quitRec=strsplit(perLine{7},' ');
            quitRec(end) = [];
            exactPress=str2double(quitRec);
            
            tmpPress=strsplit(perLine{5},' ');
            tmpPress(end) = [];
            userPress=str2double(tmpPress); % user press time
            fprintf("size exact press and user press %d %d\n",size(exactPress,2),size(userPress,2));
            for j=1:size(exactPress,2)
                for k=1:size(userPress,2)
                    if(abs(userPress(k)-exactPress(j))<1.5)
                        stationary=[stationary,abs(userPress(k)-exactPress(j))];
                        count=count+1;
                        break
                    end
                end
                
            end
            rotateRec=strsplit(perLine{8},' ');
            rotateRec(end) = [];
            exactRelease=str2double(rotateRec);
            
            tmpRelease=strsplit(perLine{6},' ');
            tmpRelease(end) = [];
            userRelease=str2double(tmpRelease); %user release time
            for j=1:size(exactRelease,2)
                for k=1:size(userRelease,2)
                    if(abs(userRelease(k)-exactRelease(j))<1.5)
                        rotation=[rotation,abs(userRelease(k)-exactRelease(j))];
                        break
                    end
                end
                
            end
            
        end
        
    end
end
% hisotgram with added normal curve
rng('default');
figure(1);
hi=histfit(stationary,10);
hold on;
hi(1).FaceColor = [.8 .9 1];
hi(1).FaceAlpha = 0;
hi(1).EdgeColor = 'none';
hi(2).Color = [0 0 1];
mu=mean(stationary);
sigma=std(stationary);
plot([mu, mu], ylim,'--','Color', '[0 0 1]', 'LineWidth', 2);
hold on;
plot([mu + sigma, mu + sigma], ylim, 'Color', '[0 1 0.8]', 'LineWidth', 1.5);
hold on;
plot([mu - sigma, mu - sigma], ylim, 'Color', '[0 1 0.8]', 'LineWidth', 1.5);

hold on;
hi2=histfit(rotation,10);
hold on;
hi2(1).FaceColor = [1 .9 .9];
hi2(1).FaceAlpha = 0;
hi2(1).EdgeColor = 'none';
hi2(2).Color = [1 0 0];
mu2=mean(rotation);
sigma2=std(rotation);
plot([mu2, mu2], ylim,'--','Color', '[1 0 0]', 'LineWidth', 2);
hold on;
plot([mu2 + sigma2, mu2 + sigma2], ylim, 'Color', '[1 0 1]', 'LineWidth', 1.5);
hold on;
plot([mu2 - sigma2, mu2 - sigma2], ylim, 'Color', '[1 0 1]', 'LineWidth', 1.5);
grid on;
xlim([0 2]);
plots=get(gca, 'Children');
legend([plots(4, 1) plots(9, 1)], {'Rotation', 'Stationary'});
title('Dashed lines are the means and side lines are -/+ SD');
xlabel('Latency(s)');
ylabel('Number of correctly identified transitions');
