clc
clearvars

%% Loading the folder

[path]=uigetdir;
cd(path)
files=dir('*.mat');
%%

bintime=30;  % in seconds
cnt=1;
for nfile=1:length(files)
    filename=files(nfile).name;
    if filename(1)~='.'
        load(filename);
        stimside=mean(target_position(target_position(:,4)==1,2));
        nonstimside=mean(target_position(target_position(:,4)==0,2));
        refside=mean(target_position(target_position(:,4)==2,2));
        left_selected=round(mean(mapImage(:,100)));
        
        if ~isnan(refside)
            stim='ref';
            sideinfo=zeros(size(target_position,1),1);
            for i=1:size(target_position,1)
                if mapImage(round(target_position(i,1)),round(target_position(i,2)))==left_selected
                    sideinfo(i)=1;
                end
            end
        else
            if stimside>nonstimside
                stim='right';
            else
                stim='left';
            end
        end
        stimside=stim;
                
        clear stim nonstimside refside         
        t=cumsum(target_position(:,3))/1000;
        nbins=ceil(max(t)/bintime);
        
        bintimes=zeros(nbins,2);
        fleft=zeros(nbins,1);
        
        for i=1:nbins
            inds=find(t>(i-1)*bintime & t<=bintime*i);
            indmin=min(inds);
            indmax=max(inds);
            bintimes(i,:)=[t(indmin) t(indmax)];
        
            fstim=mean(target_position(indmin:indmax,4));
            if strcmp(stimside,'left')
                fleft(i)=fstim;
            end
            if strcmp(stimside,'right')
                fleft(i)=1-fstim;
            end
            if strcmp(stimside,'ref')
                fleft(i)=mean(sideinfo(indmin:indmax));
            end
        end
        
        data{cnt,1}=filename;
        data{cnt,2}=bintimes;        %y and x axis coordinates for each dataset 
        data{cnt,3}=bintime;
        data{cnt,4}=fleft;     %fleft is the fraction of time spent on the left side
        data{cnt,5}=1-fleft;
        data{cnt,6}=stimside;
    
        cnt=cnt+1;
    
        %clear indmin indmax i fstim sideinfo 
    end
end

eval(['data' num2str(bintime) '=data;']);

%clearvars -except data files
% Data column 4 always shows the time spent on left side
% Data column 5 always shows the time spent on right side
%% Assigning stimulation times
% Assigning data Ref and Left Side (Stim and No-Stim)

times_ref_l = data {1,4};
times_R11_l = data {2,4};
times_R21_l = data {3,4};
times_R12_l = data {4,4};
times_R22_l = data {5,4};

%% Calculating mean, standard deviation nad standard error for left side
% Confidence Interval 95% (* plus-minus 1.96 SEM)

mean_ref_l = mean (times_ref_l);
std_dev_ref_l = std (times_ref_l);
std_err_ref_l = (std_dev_ref_l)/(sqrt(length(times_ref_l)));
x1=1:length(times_ref_l);
std_ref_low = mean_ref_l - std_dev_ref_l;
std_ref_high = mean_ref_l + std_dev_ref_l;
std_err_ref_low = ((mean_ref_l) - (1.96 * std_err_ref_l));
std_err_ref_high = ((mean_ref_l) + (1.96 * std_err_ref_l));

mean_R11_l = mean (times_R11_l);
std_dev_R11_l = std (times_R11_l);
std_err_R11_l = (std_dev_R11_l)/(sqrt(length(times_R11_l)));
x2=1:length(times_R11_l);
std_R11_low = mean_R11_l - std_dev_R11_l;
std_R11_high = mean_R11_l + std_dev_R11_l;
std_err_R11_low = ((mean_R11_l) - (1.96 * std_err_R11_l));
std_err_R11_high = ((mean_R11_l) + (1.96 * std_err_R11_l));

mean_R21_l = mean (times_R21_l);
std_dev_R21_l = std (times_R21_l);
std_err_R21_l = (std_dev_R21_l)/(sqrt(length(times_R21_l)));
x3=1:length(times_R21_l);
std_R21_low = mean_R21_l - std_dev_R21_l;
std_R21_high = mean_R21_l + std_dev_R21_l;
std_err_R21_low = ((mean_R21_l) - (1.96 * std_err_R21_l));
std_err_R21_high = ((mean_R21_l) + (1.96 * std_err_R21_l));

mean_R12_l = mean (times_R12_l);
std_dev_R12_l = std (times_R12_l);
std_err_R12_l = (std_dev_R12_l)/(sqrt(length(times_R12_l)));
x4=1:length(times_R12_l);
std_R12_low = mean_R12_l - std_dev_R12_l;
std_R12_high = mean_R12_l + std_dev_R12_l;
std_err_R12_low = ((mean_R12_l) - (1.96 * std_err_R12_l));
std_err_R12_high = ((mean_R12_l) + (1.96 * std_err_R12_l));

mean_R22_l = mean (times_R22_l);
std_dev_R22_l = std (times_R22_l);
std_err_R22_l = (std_dev_R22_l)/(sqrt(length(times_R22_l)));
x5=1:length(times_R22_l);
std_R22_low = mean_R22_l - std_dev_R22_l;
std_R22_high = mean_R22_l + std_dev_R22_l;
std_err_R22_low = ((mean_R22_l) - (1.96 * std_err_R22_l));
std_err_R22_high = ((mean_R22_l) + (1.96 * std_err_R22_l));

%% Plotting Ref and Left Stim Only

figure;
subplot(1,5,1)
p1 = bar(x1,times_ref_l);
set(p1,'FaceColor','cyan');
grid on;
title( 'Reference')
xlabel ('Bins')
ylabel ('Time Spent in Reference (in seconds)')
ylim ([0 1.2])
yline(0.5,'-','Threshold');
yline (mean_ref_l,'-','Mean')
% yline (std_dev_ref_l,'-','St. Dev')
% yline (std_err_ref_l,'-','St. Err')

hold on;

subplot(1,5,2)
p2 = bar(times_R11_l);
set(p2,'FaceColor','green');
grid on;
title( 'Region:1 - Left Stim')
xlabel ('Bins')
ylabel ('Time Spent in Region:1 (in seconds)')
ylim ([0 1.2])
yline(0.5,'-','Threshold');
yline (mean_R11_l,'-','Mean')
% yline (std_dev_R11_l,'-','St. Dev')
% yline (std_err_R11_l,'-','St. Err')

hold on;

subplot(1,5,3)
p3 = bar(times_R21_l);
set(p3,'FaceColor','magenta');
grid on;
title( 'Region:2 - Left No Stim')
xlabel ('Bins')
ylabel ('Time Spent in Region:2 (in seconds)')
ylim ([0 1.2])
yline(0.5,'-','Threshold');
yline (mean_R21_l,'-','Mean')
% yline (std_dev_R21_l,'-','St. Dev')
% yline (std_err_R21_l,'-','St. Err')

hold on;

subplot(1,5,4)
p4 = bar(times_R12_l);
set(p4,'FaceColor','green');
grid on;
title( 'Region:1 - Left Stim')
xlabel ('Bins')
ylabel ('Time Spent in Region:1 (in seconds)')
ylim ([0 1.2])
yline(0.5,'-','Threshold');
yline (mean_R12_l,'-','Mean')
% yline (std_dev_R12_l,'-','St. Dev')
% yline (std_err_R12_l,'-','St. Err')

hold on;

subplot(1,5,5)
p5 = bar(times_R22_l);
set(p5,'FaceColor','magenta');
grid on;
title( 'Region:2 - Left No Stim')
xlabel ('Bins')
ylabel ('Time Spent in Region:2 (in seconds)')
ylim ([0 1.2])
yline(0.5,'-','Threshold');
yline (mean_R22_l,'-','Mean')
% yline (std_dev_R22_l,'-','St. Dev')
% yline (std_err_R22_l,'-','St. Err')

sgtitle('OPEN FIELD TEST LEFT SIDE PREFERENCE') 

% Pairwise Comparison of Ref and LS1, Ref and LNS1, Ref and LS2, Ref and LNS2
% Statistical Analysis if the time spent in the stimulated region is signifiant or not
% T-Test (Confidence Interval = 95%)
%[h,p] = ttest(x,y,'Alpha',0.05)

[h_11,p_11] = ttest(times_ref_l,times_R11_l,'Alpha',0.05); % For Reference and Left Stim R11
[h_21,p_21] = ttest(times_R11_l,times_R21_l,'Alpha',0.05); % For Left Stim R11 and Left No Stim R21
[h_12,p_12] = ttest(times_R21_l,times_R12_l,'Alpha',0.05); % For Left No Stim R21 and Left Stim R12
[h_22,p_22] = ttest(times_R12_l,times_R22_l,'Alpha',0.05); % For Left Stim R12 and Left No Stim R22

% Bar Plot for Significance 

figure;
x_val = [mean_ref_l,mean_R11_l,mean_R21_l,mean_R12_l,mean_R22_l];
left_len = 1:length(x_val);
errhigh = [std_dev_ref_l,std_dev_R11_l,std_dev_R21_l,std_dev_R21_l,std_dev_R22_l];

m = bar (left_len,x_val);
m.FaceColor = 'flat';
m.CData(2,:) = [.9 0 .0];
m.CData(3,:) = [.4 0 .5];
m.CData(4,:) = [.9 0 .0];
m.CData(5,:) = [.4 0 .5];
grid on;
title( 'LEFT SIDE STIMULATION VS. NO-STIMULATION')
xlabel ('Reference,LS, LNS, LS, LNS')
ylabel ('Mean and Standard Deviation')
ylim ([0 1.5])

hold on;

er = errorbar(left_len,x_val,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

% Pooling the stimulated and Non Stimulated Data together

ref_pool = data {1,4};
stim_pool = [times_R11_l ; times_R12_l];
nostim_pool = [times_R21_l ; times_R22_l]; 

mean_ref_pool = mean (ref_pool);
mean_stim_pool = mean (stim_pool);
mean_nostim_pool = mean (nostim_pool);

std_ref_pool = std (ref_pool);
std_stim_pool = std (stim_pool);
std_nostim_pool = std (nostim_pool);

se_ref_pool = std (ref_pool)/ sqrt (length(ref_pool));
se_stim_pool = std (stim_pool)/ sqrt (length(stim_pool));
se_nostim_pool = std (nostim_pool)/ sqrt (length(nostim_pool));

% Plotting

figure;
y_val = [mean_ref_pool,mean_stim_pool,mean_nostim_pool];
left_len_y = 1:length(y_val);
errhigh_pool = [std_ref_pool,std_stim_pool,std_nostim_pool];

n = bar (left_len_y,y_val);
n.FaceColor = 'flat';
n.CData(2,:) = [.9 0 .0];
n.CData(3,:) = [.4 0 .5];

grid on;

title( 'LEFT SIDE STIMULATION VS. NO-STIMULATION')
xlabel ('Reference,LS, LNS')
ylabel ('Mean and Standard Deviation')
ylim ([0 1.2])

hold on;

er_y = errorbar(left_len_y,y_val,errhigh_pool);    
er_y.Color = [0 0 0];                            
er_y.LineStyle = 'none';  

hold on;
legend ('Location', 'northeast')

hold off

% Statistics (TTest) : if the time spent in the stimulated region is signifiant or not
% T-Test (Confidence Interval = 95%)
%[h,p] = ttest(x,y,'Alpha',0.05)

[h_y1,p_y1] = ttest2(ref_pool,stim_pool,'Alpha',0.05); % For Reference and Left Stim Pool
[h_y2,p_y2] = ttest2(stim_pool,nostim_pool,'Alpha',0.05); % For Stim Pool and NoStim Pool
[h_y3,p_y3] = ttest2(ref_pool,nostim_pool,'Alpha',0.05); % For Stim Pool and NoStim Pool

%% REMOVING FIRST ONE MINUTE FROM RFERENCE, STIM and NOSTIM

data_new = data;
%bin_exclude = 2;

data_new {1,7} = [times_ref_l(3:end)];
data_new {2,7} = [times_R11_l(3:end)];
data_new {3,7} = [times_R21_l(3:end)];
data_new {4,7} = [times_R12_l(3:end)];
data_new {5,7} = [times_R22_l(3:end)];

% Calculatimng mean, sd, sem

new_ref_mean = mean (data_new {1,7});
new_R11_mean = mean (data_new {2,7});
new_R21_mean = mean (data_new {3,7});
new_R12_mean = mean (data_new {4,7});
new_R22_mean = mean (data_new {5,7});

new_ref_std = std (data_new {1,7});
new_R11_std = std (data_new {2,7});
new_R21_std = std (data_new {3,7});
new_R12_std = std (data_new {4,7});
new_R22_std = std (data_new {5,7});

new_ref_sem = (std (data_new {1,7}))/sqrt(length(data_new {1,7}));
new_R11_sem = (std (data_new {2,7}))/sqrt(length(data_new {2,7}));
new_R21_sem = (std (data_new {3,7}))/sqrt(length(data_new {3,7}));
new_R12_sem = (std (data_new {4,7}))/sqrt(length(data_new {4,7}));
new_R22_sem = (std (data_new {5,7}))/sqrt(length(data_new {5,7}));

x6=1:length(data_new {1,7});
x7=1:length(data_new {2,7});
x8=1:length(data_new {3,7});
x9=1:length(data_new {4,7});
x10=1:length(data_new {5,7});

% Keeping the first two value (1 minute) of each condition in a new column
% of data_new (col:8)

data_new {1,8} = [times_ref_l(1:2)];
data_new {2,8} = [times_R11_l(1:2)];
data_new {3,8} = [times_R21_l(1:2)];
data_new {4,8} = [times_R12_l(1:2)];
data_new {5,8} = [times_R22_l(1:2)];

% Calculating mean, sd, sem

new_ref_mean_nex = mean (data_new {1,8});
new_R11_mean_nex = mean (data_new {2,8});
new_R21_mean_nex = mean (data_new {3,8});
new_R12_mean_nex = mean (data_new {4,8});
new_R22_mean_nex = mean (data_new {5,8});

new_ref_std_nex = std (data_new {1,8});
new_R11_std_nex = std (data_new {2,8});
new_R21_std_nex = std (data_new {3,8});
new_R12_std_nex = std (data_new {4,8});
new_R22_std_nex = std (data_new {5,8});

new_ref_sem_nex = (std (data_new {1,8}))/sqrt(length(data_new {1,8}));
new_R11_sem_nex = (std (data_new {2,8}))/sqrt(length(data_new {2,8}));
new_R21_sem_nex = (std (data_new {3,8}))/sqrt(length(data_new {3,8}));
new_R12_sem_nex = (std (data_new {4,8}))/sqrt(length(data_new {4,8}));
new_R22_sem_nex = (std (data_new {5,8}))/sqrt(length(data_new {5,8}));

%% PLOTTING AFTER REMOVING FIRST ONE MINUTE

figure;
subplot(1,5,1)
p6 = bar(data_new {1,7});
set(p6,'FaceColor','cyan');
grid on;
title( 'Reference')
xlabel ('Bins')
ylabel ('Time Spent in Reference (in seconds)')
ylim ([0 1.2])
yline(0.5,'-','Threshold');
yline (mean_ref_l,'-','Mean')

hold on;

subplot(1,5,2)
p7 = bar(data_new {2,7});
set(p7,'FaceColor','green');
grid on;
title( 'Region:1 - Left Stim')
xlabel ('Bins')
ylabel ('Time Spent in Region:1 (in seconds)')
ylim ([0 1.2])
yline(0.5,'-','Threshold');
yline (mean_R11_l,'-','Mean')

hold on;

subplot(1,5,3)
p8 = bar(data_new {3,7});
set(p8,'FaceColor','magenta');
grid on;
title( 'Region:2 - Left No Stim')
xlabel ('Bins')
ylabel ('Time Spent in Region:2 (in seconds)')
ylim ([0 1.2])
yline(0.5,'-','Threshold');
yline (mean_R21_l,'-','Mean')

hold on;

subplot(1,5,4)
p9 = bar(data_new {4,7});
set(p9,'FaceColor','green');
grid on;
title( 'Region:1 - Left Stim')
xlabel ('Bins')
ylabel ('Time Spent in Region:1 (in seconds)')
ylim ([0 1.2])
yline(0.5,'-','Threshold');
yline (mean_R12_l,'-','Mean')

hold on;

subplot(1,5,5)
p10 = bar(data_new {5,7});
set(p10,'FaceColor','magenta');
grid on;
title( 'Region:2 - Left No Stim')
xlabel ('Bins')
ylabel ('Time Spent in Region:2 (in seconds)')
ylim ([0 1.2])
yline(0.5,'-','Threshold');
yline (mean_R22_l,'-','Mean')

sgtitle('OPEN FIELD TEST LEFT SIDE PREFERENCE') 

% Bar Plot for Significance after Removing 1 minute

figure;
x_val_new = [new_ref_mean,new_R11_mean,new_R21_mean,new_R12_mean,new_R22_mean];
left_len_new = 1:length(x_val_new);
errhigh_new = [new_ref_std,new_R11_std,new_R21_std,new_R12_std,new_R22_std];

a = bar (left_len_new,x_val_new);
a.FaceColor = 'flat';
a.CData(2,:) = [.9 0 .0];
a.CData(3,:) = [.4 0 .5];
a.CData(4,:) = [.9 0 .0];
a.CData(5,:) = [.4 0 .5];
grid on;
title( 'LEFT SIDE STIMULATION VS. NO-STIMULATION')
xlabel ('Reference,LS, LNS, LS, LNS')
ylabel ('Mean and Standard Deviation')
ylim ([0 1.5])

hold on;

er = errorbar(left_len_new,x_val_new,errhigh_new);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off

% Pooled Bar Plot for Significance after Removing 1 minute

l1=data_new {2,7};
l2=data_new {3,7};
l3=data_new {4,7};
l4=data_new {5,7};

ref_pool_new = data_new {1,7};
stim_pool_new = [l1;l3];
nostim_pool_new =[l2;l4];

mean_ref_pool_new = mean (ref_pool_new);
mean_stim_pool_new = mean (stim_pool_new);
mean_nostim_pool_new = mean (nostim_pool_new);

std_ref_pool_new = std (ref_pool_new);
std_stim_pool_new = std (stim_pool_new);
std_nostim_pool_new = std (nostim_pool_new);

se_ref_pool_new = std (ref_pool_new)/ sqrt (length(ref_pool_new));
se_stim_pool_new = std (stim_pool_new)/ sqrt (length(stim_pool_new));
se_nostim_pool_new = std (nostim_pool_new)/ sqrt (length(nostim_pool_new));

% Plotting

figure;
y_val_new = [mean_ref_pool_new,mean_stim_pool_new,mean_nostim_pool_new];
left_len_y_new = 1:length(y_val_new);
errhigh_pool_new = [std_ref_pool_new,std_stim_pool_new,std_nostim_pool_new];

b = bar (left_len_y_new,y_val_new);
b.FaceColor = 'flat';
b.CData(2,:) = [.9 0 .0];
b.CData(3,:) = [.4 0 .5];

grid on;

title( 'LEFT SIDE STIMULATION VS. NO-STIMULATION AFTER REMOVING FIRST ONE MINUTE')
xlabel ('Reference,LS, LNS')
ylabel ('Mean and Standard Deviation')
ylim ([0 1.2])

hold on;

er_y_new = errorbar(left_len_y_new,y_val_new,errhigh_pool_new);    
er_y_new.Color = [0 0 0];                            
er_y_new.LineStyle = 'none';  

hold off

% Statistics (TTest) : if the time spent in the stimulated region is signifiant or not
% T-Test (Confidence Interval = 95%)
%[h,p] = ttest(x,y,'Alpha',0.05)

[h_y4,p_y4] = ttest2(ref_pool_new,stim_pool_new,'Alpha',0.05); % For Reference and Left Stim Pool
[h_y5,p_y5] = ttest2(stim_pool_new,nostim_pool_new,'Alpha',0.05); % For Stim Pool and NoStim Pool
[h_y6,p_y6] = ttest2(ref_pool_new,nostim_pool_new,'Alpha',0.05); % For Ref Pool and NoStim Pool

%% 
clc
clearvars

%% Excel File Import as a numeric matrix: Saved at Desktop
% D1-DMS (Col:5)
% D1-DLS (Col:10)
% D1-NAc (Col:15)

excel_file = 'Feb_OFT_DATA_SUMMARY_Pvalues.xlxs';
[num_data, txt_data, raw_data] = xlsread(excel_file);

%% Plotting ANOVA

data = FebOFTDATASUMMARYPvalues (:,3);
data_valid = data(~isnan(data));

% One Way ANOVA

[p, tbl, stats] = anova1(data_valid);

% Display ANOVA results

disp(tbl);
disp(['p-value: ', num2str(p)]);

%% Kruskal-Wallis Test

[p, tbl, stats] = kruskalwallis(data_valid);

% Kruskal-Wallis test results

disp(tbl);
disp(['p-value: ', num2str(p)]);

%% %% General Trace Plot (DO NOT RUN)

clc
clearvars

%%

axis_y=target_position(:,1);
axis_x=target_position(:,2);
figure;
plot(axis_x,axis_y)
title('TRACE GENERAL')

