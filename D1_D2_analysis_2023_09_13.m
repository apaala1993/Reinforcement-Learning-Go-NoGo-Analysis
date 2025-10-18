%% Step 0: Get all the raw data for a mouse

clearvars;
clear all;

[folders]=uigetdir2;  %Directory selection dialog box which remembers the last directory selected, open the uigetdir2 function

% Generate a list of files by session dates 

nfolders=length(folders);     %nfolders goes through all the sub-folders in the chosen folder
cnt=1;                        %cnt variable keeps a count of the number of sub-folders in the chosen folder
odorset_dates=cell(1);
for nfolder=1:nfolders        %this loop goes through each subfolder in the chosen folder
    folder=folders{nfolder};  
    cd(folder)                %change directory
    odorfolders=dir;          %lists the files in a folder
    odorfolders=odorfolders(3:end);

    nodorfolders=length(odorfolders); %number of odor folders present in the sub-folder (odorwise division) 
    for nodorfolder=1:nodorfolders    %odorfolder take each of the folders seperated by odors 
        odorfolder=odorfolders(nodorfolder).name; 
        if odorfolder (1) ~='.' && odorfolders(nodorfolder).isdir % first character . , for mac skip it
            cd(odorfolder); %if not ., it goes to that folder
            fnames=dir;  
            fnamecnt=1;   % number of files, dummy variable to skip the files
            fname=fnames(fnamecnt).name;
            while fname(1)=='.'                 % if it is an odor folder, then we get the file names only if the first char is not a .
                fname=fnames(fnamecnt).name;    
                fnamecnt=fnamecnt+1;            %fnamecnt is now the number of files without .
            end
            fdate=split(fname,'_');             % '_' is being used as a seperator to get the file dates
            mouseID=fdate{1};
            fdate=str2num(fdate{2})';           % this converts the charcters to a string
            odorset_dates{cnt,1}=odorfolder;    % cnt is the number of odor folders
            
            odorset_dates{cnt,2}=fdate;
            odorset_dates{cnt,3}=[folder '/' odorfolder];  % record the position of teh odor sets
            cnt=cnt+1;                                     % cnt has the total number of folders
            cd('..');
        end
    end
    cd('..');
end

alldata=sortrows(odorset_dates,2); % sorts the files according to the dates according to teh dates they were done
clearvars -except alldata mouseID
info_alldata{1,1}='odor set';
info_alldata{2,1}='date';
info_alldata{3,1}='file location';

% Get session data

nfolders=size(alldata,1);
for i=1:nfolders
    cd(alldata{i,3});
    fnames=dir('*.txt');
    data=[];
    for j=1:size(fnames,1)
        fdate=split(fnames(j).name,'_');             
        fdates{j,1}=str2num(fdate{2})'; 
        temp{j,1}=load(fnames(j).name);
        temp{j,1}(:,end)=j;
        data=[data; temp{j,1}];
    end
    alldata{i,2}=fdates;
    alldata{i,4}=temp;
    data(:,1)=1:size(data,1);
    alldata{i,5}=data;
    clear fdates temp
end
info_alldata{4,1}='Session data';
info_alldata{5,1}='Combined session data';
clearvars -except alldata info_alldata mouseID
cd('..');

nodorsets=size(alldata,1);
cnt=1;
for i=1:nodorsets
    dates=alldata{i,2};
    data=alldata{i,4};
    for j=1:length(dates)
        allsessions{cnt,1}=alldata{i,1};
        allsessions{cnt,2}=dates{j};
        allsessions{cnt,3}=data{j};
        cnt=cnt+1;
    end
end

info{1,1}='odor set';
info{2,1}='date';
info{3,1}='session data';

 nsessions=size(allsessions,1);

 for i=1:nsessions
    data=allsessions{i,3};
    data=data(data(:,3)~=0 | data(:,4)~=0,:);
    allsessions{i,4}=data;
    for j=1:4
        odordata{j,1}=data(data(:,2)==j,:);
    end
    allsessions{i,5}=odordata;
    clear odordata data
 end

 info{4,1}='incomplete trials removed';
 info{5,1}='odors separated after removing the incomplete trials';

clearvars -except alldata info_alldata info mouseID allsessions
save([mouseID '_alldata.mat'], 'alldata', 'allsessions', 'info_alldata', 'info', 'mouseID');

%% Step1 (Run to load the _alldata file): 
% Calculate the learning curves with sliding window,find indices of
% learning condition and get reaction times
% Data will be saved as _processed_data
 
winslide=11;
learn_cond=[3 3];  % three correct in three trials

[fname,pname]=uigetfile( '*.mat');
cd(pname);
load(fname);

nsessions=size(allsessions,1);
for i=1:nsessions
    for j=1:4
        data=allsessions{i,5}{j};
        learning_data{j,1}=[data(:,1) movmean(data(:,3),winslide)];
        rttemp=(data(:,9)-data(:,8)-700)/1000;
        rttemp(rttemp<0)=NaN;
        rt{j,1}=[data(:,1) rttemp];
        %rt{j,1}=[data(:,1) movmean(rttemp,winslide,'omitnan')];
        ncond{j,1}=find(movsum(data(:,3),learn_cond(2))==learn_cond(1));
    end
    allsessions{i,6}=learning_data;
    allsessions{i,7}=ncond;
    allsessions{i,8}=rt;
    clear learning_data ncond rt rttemp
end

info{6,1}='learning curves';
info{7,1}='n (learning trial number)';
info{8,1}='reaction times';

clearvars -except alldata info_alldata info mouseID allsessions
save([mouseID '_processed_data.mat'], 'alldata', 'allsessions', 'info_alldata', 'info', 'mouseID');


%% Check inclusion criteria (applies to NoGo odors (odors 2 and 4) only) and align the learning curves
% 1. Mouse did not know the odors: less than 20% correct in the first 11 trials
% 2. Mouse learned the odor set (3 consecutive correct trials)
% 3. After learning the performance did not drop below 20% for at least 11 trials

pcor=0.2;
winslide=11;
nfirst=(winslide-1)/2;

nsessions=size(allsessions,1);
for i=1:nsessions
    % check criteria 
    allsessions{i,9}=0;
    if ~isempty(allsessions{i,7}{2}) && ~isempty(allsessions{i,7}{4})
        lc2=allsessions{i,6}{2}(:,2);
        lc4=allsessions{i,6}{4}(:,2);
        if isempty(find(lc2(1:nfirst)>pcor)) && isempty(find(lc4(1:nfirst)>pcor))
            n2=min(allsessions{i,7}{2});
            n4=min(allsessions{i,7}{4});
            n2end=min(length(lc2),n2+winslide);
            n4end=min(length(lc4),n4+winslide);
            if isempty(find(lc2(1:nfirst)>pcor)) && isempty(find(lc4(1:nfirst)>pcor))
                if isempty(find(lc2(n2end:end)<pcor)) && isempty(find(lc4(n4end:end)<pcor))
                    allsessions{i,9}=1;
                end
            end
        end
    end
end

clearvars -except perfdata info mouseID allsessions
        
% Align learning curves and reaction times

nsessions=size(allsessions,1);
lcmat2=nan(nsessions,1000);
lcmat4=nan(nsessions,1000);
rtmat2=nan(nsessions,1000);
rtmat4=nan(nsessions,1000);
for i=1:nsessions
    if allsessions{i,9}==1
        n2=min(allsessions{i,7}{2});
        n4=min(allsessions{i,7}{4});
        lc2=allsessions{i,6}{2}(:,2);
        lc4=allsessions{i,6}{4}(:,2);
        rt2=allsessions{i,8}{2}(:,2);
        rt4=allsessions{i,8}{4}(:,2);
        startind2=500-n2+1;
        stopind2=500-n2+length(lc2);
        startind4=500-n4+1;
        stopind4=500-n4+length(lc4);
        lcmat2(i,startind2:stopind2)=lc2;
        lcmat4(i,startind4:stopind4)=lc4;
        rtmat2(i,startind2:stopind2)=rt2;
        rtmat4(i,startind4:stopind4)=rt4;
    end
end

clearvars -except allsessions info mouseID lcmat2 lcmat4 rtmat2 rtmat4

save([mouseID '_learning.mat'],'allsessions', 'info', 'mouseID', 'lcmat2', 'lcmat4', 'rtmat2', 'rtmat4');

%%

figure;
subplot(2,1,1); hold on;
plot(nanmean(lcmat2),'b');
plot(nanmean(lcmat2)+nanstd(lcmat2),'b--');
plot(nanmean(lcmat2)-nanstd(lcmat2),'b--');
plot(nanmean(lcmat4),'r');
plot(nanmean(lcmat4)+nanstd(lcmat4),'r--');
plot(nanmean(lcmat4)-nanstd(lcmat4),'r--');

subplot(2,1,2); hold on;
plot(nanmean(rtmat2),'b');
% plot(nanmean(rtmat2)+nanstd(rtmat2),'b--');
% plot(nanmean(rtmat2)-nanstd(rtmat2),'b--');
plot(nanmean(rtmat4),'r');
% plot(nanmean(rtmat4)+nanstd(rtmat4),'r--');
% plot(nanmean(rtmat4)-nanstd(rtmat4),'r--');
ylim([0.3 0.8])