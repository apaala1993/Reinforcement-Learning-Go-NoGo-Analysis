
files=dir('*_learning_Stim_Feedback.mat');
mouse_ids=[];
fin_perf=[];
cnt=1;
plotson=0;
for i=1:length(files)
    ind=findstr(files(i).name,'_');
    if ~isempty(ind)
        mouse_id=files(i).name(1:(ind(1)-1));
        mouse_ids{cnt}=mouse_id;
        
           
        load([mouse_id '_learning_Stim_Feedback.mat'],[mouse_id '_meanLc2c']);
        load([mouse_id '_learning_Stim_Feedback.mat'],[mouse_id '_meanLc4c']);
        eval(['meanLc2c=' mouse_id '_meanLc2c;']);
        eval(['meanLc4c=' mouse_id '_meanLc4c;']);
    
        clear D*
        
        lc2=meanLc2c{3};
        c2=meanLc2c{6};
        x2=1:length(lc2);
        lc4=meanLc4c{3};
        c4=meanLc4c{6};
        x4=1:length(lc4);
        
        fits{cnt,1}=c2;
        fits{cnt,2}=c4;
        gofs{cnt,1}=meanLc2c{7};
        gofs{cnt,2}=meanLc4c{7};
        cnt=cnt+1;

        if plotson
            figure;
            sgtitle(mouse_id)
            set(gcf,'units','normalized','position',[0.2 0.1 0.4 0.8])
            subplot(2,1,1); hold on;
            plot(lc2,'ko')
            plot(c2.a./(1+exp(-c2.b*(x2-c2.c)))+c2.d,'r')
            title('odor 2')
            ylim([0 1.])
            xlim([0 max(x2(end),x4(end))])
            
            subplot(2,1,2); hold on;
            plot(lc4,'ko')
            plot(c4.a./(1+exp(-c4.b*(x4-c4.c)))+c4.d,'r')
            title('odor 4')
            ylim([0 1.])
            xlim([0 max(x2(end),x4(end))])
        end
    
        fin_perf=[fin_perf; [c2.a+c2.d c4.a+c4.d]];
    end
end

%%
perf=fin_perf;
mperf=mean(perf);

[h,p]=ttest(perf(:,1),perf(:,2)); % Paired T Test

pw=signrank(perf(:,1),perf(:,2)); % Wilcox Rank Sum Test

[p_kw, tbl_kw, stats_kw] = kruskalwallis(perf); %Kruskal Wallice

%D1-DMS (n=3) Proposal
% p = 0.8
% pw = 0.75
% p_kw = 0.5127

%D2-DMS (n=3) Proposal
% p = 0.33
% pw = 0.25
% p_kw = 0.1266

%D1-DMS (n=9) Final
% p = 0.318
% pw = 0.4258
% p_kw = 0.2697

%D1-DLS (n=5) Final
% p = 0.849
% pw = 0.8125
% p_kw = 0.9168


%% Calculating Mean, SD, SEM

FP_o2 = fin_perf(:,1);
FP_o4 = fin_perf(:,2);

avg_FP_o2 = mean(FP_o2);
sd_FP_o2 = std (FP_o2);
sem_FP_o2 = sd_FP_o2/ sqrt(length(FP_o2));

avg_FP_o4 = mean(FP_o4);
sd_FP_o4 = std (FP_o4);
sem_FP_o4 = sd_FP_o4/ sqrt(length(FP_o4));

%% Plotting a Bar plot and Grant Altman Plot

means = [avg_FP_o2, avg_FP_o4];
std_devs = [sd_FP_o2, sd_FP_o4];
std_errors = [sem_FP_o2, sem_FP_o4];

x_val = 1:2; 

bar(x_val, means, 'FaceColor', 'flat');
h = bar(x_val, means, 'FaceColor', 'flat');

h.CData(1, :) = [0, 0, 1];  
h.CData(2, :) = [1, 0, 0];  

hold on;
errorbar(x_val, means, std_errors, 'k.', 'LineWidth', 1.5);
hold off;

xlabel('Odor: 2 and Odor:4');
ylabel('Final Performance');
title('Final Performance : Combined');
legend('Odor 2', 'Odor 4');
xticks(x_val);
xticklabels({'a', 'b'});
grid on;
ylim ([0 1.5])

% Gardner Altman Plot

x = fin_perf(:, 1);
y = fin_perf(:, 2);

figure;
hold on;

for a = 1:length(x)
    if x(a) > y(a)
        plot([0, 1], [x(a), y(a)], '-o', 'Color', 'k', 'LineWidth', 2); 
    else
        plot([0, 1], [x(a), y(a)], '-o', 'Color', 'k', 'LineWidth', 2);
    end
end

scatter(0, x, 200, 'k', 'filled', 'LineWidth', 1); % Smaller size and reduced thickness, black color
scatter(1, y, 200, 'k', 'filled', 'LineWidth', 1); % Smaller size and reduced thickness, black color

hold off;

xlabel('Odor: 2 and Odor:4');
ylabel('Final Performance');
title('Final Performance : Pre-Learning Combined');
legend('Odor 2', 'Odor 4');
xticks([0 1]);
xticklabels({'Odor 2', 'Odor 4'});
grid on;
xlim([-0.5, 1.5]);
ylim([0, 1.5]);

%%
clc
for i=1:3
    display(['file ' num2str(i) ', ' mouse_ids{i}])
    display('  ')
    display('odor2')
    display(fits{i,1})
    display('  ')
    display('odor4')
    display(fits{i,2})
    display('  ')
    display('  ')
end

%%  FP in pastel colors for SFN

% Define means, standard deviations, and standard errors
means = [avg_FP_o2, avg_FP_o4];
std_devs = [sd_FP_o2, sd_FP_o4];
std_errors = [sem_FP_o2, sem_FP_o4];

% Define x-axis values
x_val = 1:2; 

% Create bar plot with specific colors for each bar
h = bar(x_val, means, 'FaceColor', 'flat');
h.CData(1, :) = [167/255, 199/255, 231/255];  % Color #A7C7E7 for Odor 2
h.CData(2, :) = [250/255, 160/255, 160/255];  % Color #FAA0A0 for Odor 4

% Add error bars to the bar plot
hold on;
errorbar(x_val, means, std_errors, 'k.', 'LineWidth', 1.5);
hold off;

% Label and format the plot
xlabel('Odor: 2 and Odor: 4');
ylabel('Final Performance');
title('Final Performance : Combined');
legend('Odor 2', 'Odor 4');
xticks(x_val);
xticklabels({'a', 'b'});
grid on;
ylim([0, 1.5]);

% Gardner Altman Plot
x = fin_perf(:, 1);
y = fin_perf(:, 2);

figure;
hold on;

% Plot lines connecting each pair of points
for a = 1:length(x)
    plot([0, 1], [x(a), y(a)], '-o', 'Color', 'k', 'LineWidth', 2);
end

% Scatter plots for Odor 2 and Odor 4
scatter(0, x, 200, 'k', 'filled', 'LineWidth', 1);
scatter(1, y, 200, 'k', 'filled', 'LineWidth', 1);

hold off;

% Label and format the plot
xlabel('Odor: 2 and Odor: 4');
ylabel('Final Performance');
title('Final Performance : Pre-Learning Combined');
legend('Odor 2', 'Odor 4');
xticks([0 1]);
xticklabels({'Odor 2', 'Odor 4'});
grid on;
xlim([-0.5, 1.5]);
ylim([0, 1.5]);
