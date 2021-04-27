% % data=[79.49,74.60,78.39,40.51,95.08];
% % average=mean(data);
% % SEM = std(data)/sqrt(length(data));
% grid on;
% % Place equation in upper left of graph.
% xl = xlim;
% yl = ylim;
% xt = 0.05 * (xl(2)-xl(1)) + xl(1)
% yt = 0.90 * (yl(2)-yl(1)) + yl(1)
% caption = sprintf('y = %f * x + %f', p(1), p(2));
% text(xt, yt, caption, 'FontSize', 16, 'Color', 'r', 'FontWeight', 'bold');
a=[2,3,1,5,6,4,7,10,9,8,12,11];
% b=a-2;
% c=a+1;
% d=[b',c'];
% disp(d);
% b=[4,6,2,10,12,8,14,20,18,16,24,22];
% R = corrcoef(a,b);
%disp(R);
% c=arrayfun(@(x) find(x==a,1),b,'un',0);
% c = cell2mat(c);
% disp(c);
%a=[1,2,3,4,5];
%b=cumsum(a)/length(a);
%disp(b);



%%
%load matlab files containing acceptable data for calculating
%the microsaccade rates
final_result=load('./matlab_files/mofo_final_microsaccades.mat', 'final_result');
microsaccs=[];
for i=1:length(final_result.final_result)
    for j=1:length(final_result.final_result{1,i})
    microsaccs=[microsaccs;final_result.final_result{1,i}(j,:)];
    end
end
magnitudes=microsaccs(:,8)/36;
histogram((microsaccs(:,8))/36,150);
save('./Data/total_microsacc.mat', 'microsaccs');
xlabel('Microsaccade magnitude(deg)');
ylabel('Number of microsaccades');
%%
%final_result=load('./matlab_files/mofo_total_microsaccades.mat');