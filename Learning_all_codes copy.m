clc
clearvars       

%This code is run on the folder named PHD_ALLDATA_OZDEN ->
%LEARNING_FINAL_DATA -> D1-D2 Project -> D1-DMS/D2-DMS

%% PART:1 Run this on the folder containing all animals (This sorts the data into individual animals matrix)
%This program collects all odorset information from all mice and saves them as _mouseID_allodors.mat

stimConditionFolderName='Stim_Feedback';
[mainfolder]=uigetdir; 
cd(mainfolder)               
animalfolders=dir; 
animalfolders=animalfolders(3:end);
inds=[];
for i=1:length(animalfolders)
    if ~animalfolders(i).isdir
        inds=[inds i];
    end
end
animalfolders(inds)=[];
nanimalfolders=length(animalfolders); 
for nanimalfolder=1:nanimalfolders
    cnt=1;
    animalfolder=animalfolders(nanimalfolder).name;
    cd(animalfolder);
    cd(stimConditionFolderName);
    odorfolders=dir;          %lists the files in a folder
    odorfolders=odorfolders(3:end);
    odorset_dates=cell(1);
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
            odorset_dates{cnt,3}=[animalfolder '/' stimConditionFolderName '/' odorfolder];  % record the position of teh odor sets
            cnt=cnt+1;                                     % cnt has the total number of folders
            cd('..');
        end
    end
    allodors=sortrows(odorset_dates,2); % sorts the files according to the dates according to teh dates they were done
    info_allodors{1,1}='odor set';
    info_allodors{2,1}='date';
    info_allodors{3,1}='file location';
    clear odorset_dates
    cd('..');
    cd('..');
    nfolders=size(allodors,1);
    for nfolder=1:nfolders
        cd(allodors{nfolder,3});
        fnames=dir('*.txt');
        data=[];
        for j=1:size(fnames,1)
            fdate=split(fnames(j).name,'_');             
            fdates{j,1}=str2num(fdate{2})'; 
            temp{j,1}=load(fnames(j).name);
            temp{j,1}(:,end)=j;
            if temp{j,1}(2,29)==10
                datatemp=temp{j,1};
                inds1=find(datatemp(:,2)==1);
                inds2=find(datatemp(:,2)==2);
                inds3=find(datatemp(:,2)==3);
                inds4=find(datatemp(:,2)==4);
                temp{j,1}(inds1,2)=3;
                temp{j,1}(inds2,2)=4;
                temp{j,1}(inds3,2)=1;
                temp{j,1}(inds4,2)=2;
            end
            data=[data; temp{j,1}];
        end
        clear datatemp
        allodors{nfolder,2}=fdates;
        allodors{nfolder,4}=temp;
        %data(:,1)=1:size(data,1);
        allodors{nfolder,5}=data;
        odorsep=cell(4,1);
        odorsep{1,1}=data(data(:,2)==1,:);
        odorsep{2,1}=data(data(:,2)==2,:);
        odorsep{3,1}=data(data(:,2)==3,:);
        odorsep{4,1}=data(data(:,2)==4,:);
        allodors{nfolder,7}=odorsep;

        data=data(data(:,3)~=0 | data(:,4)~=0,:);
        allodors{nfolder,6}=data;
        odorsep=cell(4,1);
        odorsep{1,1}=data(data(:,2)==1,:);
        odorsep{2,1}=data(data(:,2)==2,:);
        odorsep{3,1}=data(data(:,2)==3,:);
        odorsep{4,1}=data(data(:,2)==4,:);
        allodors{nfolder,8}=odorsep;
        clear fdates temp
        cd('..');
        cd('..');
        cd('..')
    end
    info_allodors{4,1}='session data';
    info_allodors{5,1}='combined session data';
    info_allodors{6,1}='combined session data with inclomplete trials removed';
    info_allodors{7,1}='combined session data separated according to odors';
    info_allodors{8,1}='combined session data with inclomplete trials removed separated according to odors';
    nodorsets=size(allodors,1);
    cnt=1;
    for i=1:nodorsets
        dates=allodors{i,2};
        data=allodors{i,4};
        for j=1:length(dates)
            allsessions{cnt,1}=allodors{i,1};
            allsessions{cnt,2}=dates{j};
            allsessions{cnt,3}=allodors{i,3};
            allsessions{cnt,4}=data{j};
            cnt=cnt+1;
        end
    end
    
    info_allsessions{1,1}='odor set';
    info_allsessions{2,1}='date';
    info_allsessions{3,1}='file_location';
    info_allsessions{4,1}='session data';
    
     nsessions=size(allsessions,1);
    
     for i=1:nsessions
        data=allsessions{i,4};
        data=data(data(:,3)~=0 | data(:,4)~=0,:);
        allsessions{i,5}=data;
        for j=1:4
            odordata{j,1}=data(data(:,2)==j,:);
        end
        allsessions{i,6}=odordata;
        clear odordata data
     end
    
     info_allsessions{5,1}='incomplete trials removed';
     info_allsessions{6,1}='odors separated after removing the incomplete trials';
    animalID=findstr(animalfolder,'_');
    animalID=animalfolder(1:(animalID-1));
    eval([animalID '_allodors=allodors;']);
    eval([animalID '_allsessions=allsessions;']);
    save([animalID '_allodors_' stimConditionFolderName '.mat'], [animalID '_allodors'], [animalID '_allsessions'],'info_allodors', 'info_allsessions', 'animalID');
    clear allodors allsessions odordata dates odorsep
end

clearvars -except *allodors *allsessions info_allsessions info_allodors animalID
clear allodors allsessions

%% PART:2 Step:1 Run this on the same Folder
%This program collects all odorset information from all mice and saves them as _mouseID_allodors.mat

% Run to load the _allodors files and apply inclusion criteria: 
% Calculate the learning curves with sliding window,find indices of
% learning condition and get reaction times
% Data will be saved as _processed_data

% Check inclusion criteria (applies to NoGo odors (odors 2 and 4) only) and align the learning curves (Runs on the _processed_data.mat)
% 1. Mouse did not know the odors: less than 20% correct in the first 11 trials
% 2. Mouse learned the odor set (3 consecutive correct trials)
% 3. After learning the performance did not drop below 20% for at least 11 trials

winslide=11;
learn_cond=[3 5];  % three correct in three trials

nfirst=(winslide-1)/2;
stimConditionFolderName='Stim_Feedback';

[mainfolder]=uigetdir; 
cd(mainfolder)               
allodorsfiles=dir(['*' 'allodors_' stimConditionFolderName '.mat']);
nfiles=length(allodorsfiles);

for nfile=1:nfiles
    fname=allodorsfiles(nfile).name;
    load(fname);
    eval(['allsessions=' animalID '_allsessions;']);
    eval(['allodors=' animalID '_allodors;']);
    
    nsessions=size(allsessions,1);
    for i=1:nsessions
        for j=1:4
            data=allsessions{i,6}{j};
            learning_data{j,1}=[data(:,1) movmean(data(:,3),winslide)];
            rttemp=(data(:,9)-data(:,8)-700)/1000;
            rttemp(rttemp<0)=NaN;
            rt{j,1}=[data(:,1) rttemp];
            %rt{j,1}=[data(:,1) movmean(rttemp,winslide,'omitnan')];
            trcond{j,1}=movsum(data(:,3),learn_cond(2));
            ncond{j,1}=find(trcond{j,1}==learn_cond(1));
        end
        allsessions{i,7}=learning_data;
        allsessions{i,8}=ncond;
        allsessions{i,9}=rt;
        allsessions{i,10}=trcond;
        clear learning_data ncond rt rttemp trcond
    end 
    eval([animalID '_allsessions=allsessions;']);
    
    nsessions=size(allodors,1);
    for i=1:nsessions
        for j=1:4
            data=allodors{i,8}{j};
            learning_data{j,1}=[data(:,1) movmean(data(:,3),winslide)];
            rttemp=(data(:,9)-data(:,8)-700)/1000;
            rttemp(rttemp<0)=NaN;
            rt{j,1}=[data(:,1) rttemp];
            %rt{j,1}=[data(:,1) movmean(rttemp,winslide,'omitnan')];
            trcond{j,1}=movsum(data(:,3),learn_cond(2));
            ncond{j,1}=find(trcond{j,1}==learn_cond(1));
        end
        allodors{i,9}=learning_data;
        allodors{i,10}=ncond;
        allodors{i,11}=rt;
        allodors{i,12}=trcond;
        clear learning_data ncond rt rttemp trcond
    end
    eval([animalID '_allodors=allodors;']);

    info_allodors{9,1}='learning curves with sliding window';
    info_allodors{10,1}='prelearning';
    info_allodors{11,1}='reaction times';
    info_allodors{12,1}='learning condition trace';
    
    info_allsessions{7,1}='learning curves with sliding window';
    info_allsessions{8,1}='prelearning';
    info_allsessions{9,1}='reaction times';
    info_allsessions{10,1}='learning condition trace';

    save([animalID '_processed_data_' stimConditionFolderName '.mat'], [animalID '_allodors'], [animalID '_allsessions'], 'info_allodors', 'info_allsessions', 'animalID', 'winslide','learn_cond');
    clear allodors allsessions info_allodors info
end

clearvars -except *allodors* *allsessions* info animalID

%% PART: 2 Step2 (Apply inclusion criteria and calculate learning curves that passed the criteria): 
% Check inclusion criteria (applies to NoGo odors (odors 2 and 4) only) and align the learning curves (Runs on the _processed_data.mat)
% 1. Mouse did not know the odors: less than 20% correct in the first 11 trials
% 2. Mouse learned the odor set (3 consecutive correct trials)
% 3. After learning the performance did not drop below 20% for at least 11 trials

pcor=0.4;
winslide=11;
nmin=10;
stimConditionFolderName='Stim_Feedback';

[mainfolder]=uigetdir; 
cd(mainfolder)               
alldatafiles=dir(['*_processed_data*' stimConditionFolderName '.mat']);
nfiles=length(alldatafiles);

for nfile=1:nfiles
    fname=alldatafiles(nfile).name;
    load(fname);
    eval(['allsessions=' animalID '_allsessions;']);
    eval(['allodors=' animalID '_allodors;']);

    nsessions=size(allsessions,1);
    lc_nc_al2=nan(nsessions,1000);
    lc_nc_nal2=nan(nsessions,1000);
    lc_nc_al4=nan(nsessions,1000);
    lc_nc_nal4=nan(nsessions,1000);
    rt_nc_al2=nan(nsessions,1000);
    rt_nc_nal2=nan(nsessions,1000);
    rt_nc_al4=nan(nsessions,1000);
    rt_nc_nal4=nan(nsessions,1000);
   
    for i=1:nsessions
        % check criteria 
        allsessions{i,11}=0;
        if ~isempty(allsessions{i,8}{2}) && ~isempty(allsessions{i,8}{4})
            lc2=allsessions{i,7}{2}(:,2);
            lc4=allsessions{i,7}{4}(:,2);
            rt2=allsessions{i,9}{2}(:,2);
            rt4=allsessions{i,9}{4}(:,2);
            nlc2=length(lc2);
            nlc4=length(lc4);
            lc_nc_nal2(i,1:nlc2)=lc2;
            lc_nc_nal4(i,1:nlc4)=lc4;

            rt_nc_nal2(i,1:nlc2)=rt2;
            rt_nc_nal4(i,1:nlc4)=rt4;
            n2=allsessions{i,8}{2};
            n4=allsessions{i,8}{4};
            % if isempty(find(lc2(1:nfirst)>pcor)) && isempty(find(lc4(1:nfirst)>pcor))
            if mean(lc2(1:nmin))<pcor && mean(lc4(1:nmin))<pcor
                nval2=[];
                if ~isempty(n2)
                    for nval=1:length(n2)
                        %if isempty(find(lc2(n2(nval):end)<pcor))
                        if mean(lc2(n2(nval):end))>pcor
                            nval2=n2(nval);
                            break
                        end
                    end
                end
                nval4=[];
                if ~isempty(n4)
                    for nval=1:length(n4)
                        %if isempty(find(lc4(n4(nval):end)<pcor))
                        if mean(lc4(n4(nval):end))>pcor
                            nval4=n4(nval);
                            break
                        end
                    end
                end
                allsessions{i,8}{2}=nval2;
                allsessions{i,8}{4}=nval4;

                if ~isempty(nval2) && ~isempty(nval4)
                    allsessions{i,11}=1;
                    startind2=500-nval2+1;
                    stopind2=500-nval2+length(lc2);
                    startind4=500-nval4+1;
                    stopind4=500-nval4+length(lc4);
                    lc_nc_al2(i,startind2:stopind2)=lc2;
                    lc_nc_al4(i,startind4:stopind4)=lc4;
                    rt_nc_al2(i,startind2:stopind2)=rt2;
                    rt_nc_al4(i,startind4:stopind4)=rt4;
                end
            end
        end
    end

    nsessions=size(allodors,1);
    nNan=5000;
    lc_c_al2=nan(nsessions,nNan);
    lc_c_nal2=nan(nsessions,nNan);
    lc_c_al4=nan(nsessions,nNan);
    lc_c_nal4=nan(nsessions,nNan);
    rt_c_al2=nan(nsessions,nNan);
    rt_c_nal2=nan(nsessions,nNan);
    rt_c_al4=nan(nsessions,nNan);
    rt_c_nal4=nan(nsessions,nNan);
    for i=1:nsessions
        % check criteria 
        allodors{i,13}=0;
        if ~isempty(allodors{i,8}{2}) && ~isempty(allodors{i,8}{4})
            lc2=allodors{i,9}{2}(:,2);
            lc4=allodors{i,9}{4}(:,2);
            nlc2=length(lc2);
            nlc4=length(lc4);
            rt2=allodors{i,11}{2}(:,2);
            rt4=allodors{i,11}{4}(:,2);
            lc_c_nal2(i,1:nlc2)=lc2;
            lc_c_nal4(i,1:nlc4)=lc4;
            rt_c_nal2(i,1:nlc2)=rt2;
            rt_c_nal4(i,1:nlc4)=rt4;
            n2=allodors{i,10}{2};
            n4=allodors{i,10}{4};
            % if isempty(find(lc2(1:nfirst)>pcor)) && isempty(find(lc4(1:nfirst)>pcor))
            if mean(lc2(1:nmin))<pcor && mean(lc4(1:nmin))<pcor
                nval2=[];
                if ~isempty(n2)
                    for nval=1:length(n2)
                        %if isempty(find(lc2(n2(nval):end)<pcor))
                        if mean(lc2(n2(nval):end))>pcor
                            nval2=n2(nval);
                            break
                        end
                    end
                end
                nval4=[];
                if ~isempty(n4)
                    for nval=1:length(n4)
                        %if isempty(find(lc4(n4(nval):end)<pcor))
                        if mean(lc4(n4(nval):end))>pcor
                            nval4=n4(nval);
                            break
                        end
                    end
                end
                allodors{i,10}{2}=nval2;
                allodors{i,10}{4}=nval4;
                if ~isempty(nval2) && ~isempty(nval4)
                    allodors{i,13}=1;
                    startind2=nNan/2-nval2+1;
                    stopind2=nNan/2-nval2+length(lc2);
                    startind4=nNan/2-nval4+1;
                    stopind4=nNan/2-nval4+length(lc4);
                    lc_c_al2(i,startind2:stopind2)=lc2;
                    lc_c_al4(i,startind4:stopind4)=lc4;
                        y2=lc_c_al2(i,:);
                        y4=lc_c_al4(i,:);
                        Ny2=length(y2);
                        Ny4=length(y4);
                        inds=~isnan(y2);
                        ind21=find(inds==1,1);
                        ind22=find(inds==1, 1, 'last' );
                        inds=~isnan(y4);
                        ind41=find(inds==1,1);
                        ind42=find(inds==1, 1, 'last' );
                        
                        a=0.75;
                        b=1;
                        d=0.05;
                        c2=Ny2/2-ind21; 
                        x2=1:(ind22-ind21+1);
                        y2=y2(ind21:ind22);

                        c4=Ny4/2-ind41; 
                        x4=1:(ind42-ind41+1);
                        y4=y4(ind41:ind42);

                        [coeffs2,~]=io_sigmoid_fit(x2,y2,[a b c2 d]);
                        [coeffs4,~]=io_sigmoid_fit(x4,y4,[a b c4 d]);
                        coeffss{1,1}=coeffs2;
                        coeffss{2,1}=coeffs4;
                    
                    allodors{i,14}=coeffss;
                    rt_c_al2(i,startind2:stopind2)=rt2;
                    rt_c_al4(i,startind4:stopind4)=rt4;
                end
            end
        end
    end

    info_allodors{13,1}='1: passed the learning criteria';
    info_allodors{14,1}='Sigmoid fit parameters a, b and c (a/(1+exp(-b*(x-c)))';
    info_allsessions{11,1}='1: passed the learning criteria';
    totalPassedCombined=sum(nanmean(lc_c_al2')>0);
    totalPassedNotCombined=sum(nanmean(lc_nc_al2')>0);


    pinds=sum(~isnan(lc_c_al2'))>0;
    nsessions=sum(pinds);
    cond2=sum(~isnan(lc_c_al2));
    ninit2=find(cond2>2,1);
    nfin2=find(cond2>2, 1,'last');
    cond4=sum(~isnan(lc_c_al4));
    ninit4=find(cond4>2,1);
    nfin4=find(cond4>2, 1,'last');
    Nl=size(lc_c_al2,2);
   
    cond2=cond2(ninit2:nfin2);
    cond4=cond4(ninit4:nfin4);
    meanLc2c{1,1}=lc_c_al2(pinds,ninit2:nfin2);
    meanLc2c{2,1}=Nl/2-ninit2;
    meanLc2c{3,1}=nanmean(meanLc2c{1,1});
    meanLc2c{4,1}=nanstd(meanLc2c{1,1});
    meanLc2c{5,1}=meanLc2c{4,1}./sqrt(cond2);
    [cs,gof]=io_sigmoid_fit(1:(nfin2-ninit2+1),meanLc2c{3,1},[0.75 1 meanLc2c{2,1} 0.05]);
    meanLc2c{6,1}=cs;
    meanLc2c{7,1}=gof;
    
    meanLc4c{1,1}=lc_c_al4(pinds,ninit4:nfin4);
    meanLc4c{2,1}=Nl/2-ninit4;
    meanLc4c{3,1}=nanmean(meanLc4c{1,1});
    meanLc4c{4,1}=nanstd(meanLc4c{1,1});
    meanLc4c{5,1}=meanLc4c{4,1}./sqrt(cond4);
    [cs,gof]=io_sigmoid_fit(1:(nfin4-ninit4+1),meanLc4c{3,1},[0.75 1 meanLc4c{2,1} 0.05]);
    meanLc4c{6,1}=cs;
    meanLc4c{7,1}=gof;
    
    eval([animalID '_allsessions=allsessions;']);
    eval([animalID '_allodors=allodors;']);
    eval([animalID '_lc_nc_al2=lc_nc_al2;']);
    eval([animalID '_lc_nc_al4=lc_nc_al4;']);
    eval([animalID '_lc_nc_nal2=lc_nc_nal2;']);
    eval([animalID '_lc_nc_nal4=lc_nc_nal4;']);
    eval([animalID '_lc_c_al2=lc_c_al2;']);
    eval([animalID '_lc_c_al4=lc_c_al4;']);
    eval([animalID '_lc_c_nal2=lc_c_nal2;']);
    eval([animalID '_lc_c_nal4=lc_c_nal4;']);

    eval([animalID '_rt_nc_al2=rt_nc_al2;']);
    eval([animalID '_rt_nc_al4=rt_nc_al4;']);
    eval([animalID '_rt_nc_nal2=rt_nc_nal2;']);
    eval([animalID '_rt_nc_nal4=rt_nc_nal4;']);
    eval([animalID '_rt_c_al2=rt_c_al2;']);
    eval([animalID '_rt_c_al4=rt_c_al4;']);
    eval([animalID '_rt_c_nal2=rt_c_nal2;']);
    eval([animalID '_rt_c_nal4=rt_c_nal4;']);

    eval([animalID '_meanLc2c=meanLc2c;']);
    eval([animalID '_meanLc4c=meanLc4c;']);


    eval([animalID '_totalPassedCombined=totalPassedCombined;']);
    eval([animalID '_totalPassedNotCombined=totalPassedNotCombined;']);
  
    save([animalID '_learning_' stimConditionFolderName '.mat'], [animalID '_allsessions'], [animalID '_allodors'], 'info_allsessions', 'info_allodors', 'animalID', ...
        [animalID '_lc_nc_al2'],[animalID '_lc_nc_al4'], [animalID '_rt_nc_al2'],[animalID '_rt_nc_al4'],...
        [animalID '_lc_nc_nal2'],[animalID '_lc_nc_nal4'], [animalID '_rt_nc_nal2'],[animalID '_rt_nc_nal4'],...
        [animalID '_lc_c_al2'],[animalID '_lc_c_al4'], [animalID '_rt_c_al2'],[animalID '_rt_c_al4'],...
        [animalID '_lc_c_nal2'],[animalID '_lc_c_nal4'], [animalID '_rt_c_nal2'],[animalID '_rt_c_nal4'],...
        [animalID '_totalPassedCombined'],[animalID '_totalPassedNotCombined'], [animalID '_meanLc2c'],[animalID '_meanLc4c']);
    clear allsessions info animalID lc_nc_al2 rt_nc_al2 lc_nc_al4 rt_nc_al4 lc_nc_nal2 rt_nc_nal2 lc_nc_nal4 rt_nc_nal4
    clear lc_c_al2 rt_c_al2 lc_c_al4 rt_c_al4 lc_c_nal2 rt_c_nal2 lc_c_nal4 rt_c_nal4 y2 y4 x2 x4 Ny2 NY4 ind21 ind22 ind41 ind42 inds
end

clearvars -except D* info_allodors info_allsessions

%% PART:2 Step3 (Apply inclusion criteria and calculate learning curves): 
% Check inclusion criteria (applies to NoGo odors (odors 2 and 4) only) and align the learning curves (Runs on the _processed_data.mat)
% 1. Mouse did not know the odors: less than 20% correct in the first 11 trials
% 2. Mouse learned the odor set (3 consecutive correct trials)
% 3. After learning the performance did not drop below 20% for at least 11 trials

clear all

[mainfolder]=uigetdir; 
cd(mainfolder)               
alldatafiles=dir(['*_learning_Stim_Feedback*' '.mat']);
nfiles=length(alldatafiles);

prelearningCombined=cell(0);
prelearningNotCombined=cell(0);
PLc=[];
PLnc=[];

for nfile=1:nfiles
    fname=alldatafiles(nfile).name;
    load(fname);
    eval(['allodors = ' animalID '_allodors;']);
    eval(['allsessions = ' animalID '_allsessions;']);
    
    nC=[];
    nNc=[];
    prelearningCombined{nfile,1}=fname;
    prelearningNotCombined{nfile,1}=fname;

    nodors=size(allodors,1);
    for i=1:nodors
        if allodors{i,13}==1
            nC=[nC; [allodors{i,10}{2} allodors{i,10}{4}]];
        end
    end
    if ~isempty(nC)
        prelearningCombined{nfile,2}=nC;
        if size(nC,1)>1
            prelearningCombined{nfile,3}=mean(nC);
            PLc=[PLc; mean(nC)];
            prelearningCombined{nfile,4}=std(nC);
        else
            prelearningCombined{nfile,3}=nC;
            PLnc=[PLc; nC];
            prelearningCombined{nfile,4}=NaN;
        end
        [h,p]=ttest(nC(:,1),nC(:,2),'Alpha',0.05);
        prelearningCombined{nfile,5}=p;
    end


    nsessions=size(allsessions,1);
    for i=1:nsessions
        if allsessions{i,11}==1
            nNc=[nNc; [allsessions{i,8}{2} allsessions{i,8}{4}]];
        end
    end
    if ~isempty(nNc)
        prelearningNotCombined{nfile,2}=nNc;
        if size(nNc,1)>1
            prelearningNotCombined{nfile,3}=mean(nNc);
            PLnc=[PLnc; mean(nNc)];
            prelearningNotCombined{nfile,4}=std(nNc);
        else
            prelearningNotCombined{nfile,3}=nNc;
            PLnc=[PLnc; nNc];
            prelearningNotCombined{nfile,4}=NaN;
        end
        [h,p]=ttest(nNc(:,1),nNc(:,2),'Alpha',0.05);
        prelearningNotCombined{nfile,5}=p;
    end

end

clearvars -except prelearningCombined prelearningNotCombined PLc PLnc

% CALCULATING THE PRELEARNING SECTION OF THE LEARNING CURVES by Trial Number

mean_PLc_nostim = mean(PLc(:, 1));
mean_PLc_stim = mean(PLc(:, 2));

std_PLc_nostim = std(PLc(:, 1));
std_PLc_stim = std(PLc(:, 2));

sem_PLc_nostim = std_PLc_nostim / sqrt(length(PLc(:, 1)));
sem_PLc_stim = std_PLc_stim / sqrt(length(PLc(:, 1)));

means = [mean_PLc_nostim, mean_PLc_stim];
std_devs = [std_PLc_nostim, std_PLc_stim];
std_errors = [sem_PLc_nostim, sem_PLc_stim];


x_val = 1:2; 

h = bar(x_val, means, 'FaceColor', 'flat');
h.CData(1, :) = [91/255, 141/255, 184/255];  % Color #A7C7E7 for no stim
h.CData(2, :) = [229/255, 115/255, 115/255];  % Color #FAA0A0 for stim

hold on;
errorbar(x_val, means, std_errors, 'k.', 'LineWidth', 1.5);
hold off;

xlabel('No Stim and Stim');
ylabel('Performance');
title('Pre Learning, C_AL : No Stim vs. Stim');
legend('No Stim', 'Stim');
xticks(x_val);
xticklabels({'No Stim', 'Stim'});
%grid on;
ylim([0, 200]);

[h_PLc,p_PLc]=ttest(PLc(:,1),PLc(:,2)); % Paired T Test

% Gardner-Altman plot for Pre-learning (PLc)

x = PLc(:, 1);  % No Stim
y = PLc(:, 2);  % Stim

noStimColor = [91/255, 141/255, 184/255];   
stimColor   = [229/255, 115/255, 115/255];

figure; hold on;

hLines = gobjects(length(x), 1);
for a = 1:length(x)
    hLines(a) = plot([0, 1], [x(a), y(a)], 'Color', 'k', 'LineWidth', 2);
end

hNoStim = scatter(zeros(size(x)), x, 400, noStimColor, 'filled', 'LineWidth', 1);
hStim   = scatter(ones(size(y)), y, 400, stimColor, 'filled', 'LineWidth', 1);

hold off;

xlabel('No Stim and Stim');
ylabel('Performance');
title('Pre Learning, C_AL : Gardner-Altman');
legend([hLines(1), hNoStim, hStim], {'Connection', 'No Stim', 'Stim'}, 'Location', 'northeast');
xticks([0 1]);
xticklabels({'No Stim', 'Stim'});
%grid on;
xlim([-0.5, 1.5]);

ylim([0, 200]);
[h_PLCGA, p_PLCGA] = ttest(x, y);
disp(['Paired t-test p-value: ', num2str(p_PLCGA)]);


%% === Precentage Correct: Pre-learning (Odor 2 vs Odor 4) ===

[mainfolder] = uigetdir([], 'Select folder with *_learning_Stim_Feedback*.mat');
cd(mainfolder)
alldatafiles = dir('*_learning_Stim_Feedback*.mat');

PLpct_session_all = [];
animal_ids = {};

for f = 1:numel(alldatafiles)
    fname = alldatafiles(f).name;
    load(fname);
    eval(['allsessions = ' animalID '_allsessions;']);

    nsess = size(allsessions,1);
    for i = 1:nsess
        if isempty(allsessions{i,11}) || allsessions{i,11} ~= 1
            continue
        end

        nval2 = allsessions{i,8}{2};   % Odor 2 (NoStim)
        nval4 = allsessions{i,8}{4};   % Odor 4 (Stim)
        if isempty(nval2) || isempty(nval4) || nval2 <= 1 || nval4 <= 1
            continue
        end

        data2 = allsessions{i,6}{2};
        data4 = allsessions{i,6}{4};

        n2 = min(nval2-1, size(data2,1));
        n4 = min(nval4-1, size(data4,1));
        if n2 < 1 || n4 < 1, continue; end

        c2 = data2(1:n2, 3);
        c4 = data4(1:n4, 3);
        pct2 = 100 * mean(c2 == 1, 'omitnan');
        pct4 = 100 * mean(c4 == 1, 'omitnan');

        PLpct_session_all = [PLpct_session_all; pct2, pct4];
        animal_ids{end+1,1} = animalID;
    end
end

if isempty(PLpct_session_all)
    warning('No qualifying sessions found for pre-learning %% correct.');
else
    g1 = PLpct_session_all(:,1);
    g2 = PLpct_session_all(:,2);

    m1 = mean(g1,'omitnan');  m2 = mean(g2,'omitnan');
    s1 = std(g1,'omitnan');   s2 = std(g2,'omitnan');
    n1 = sum(~isnan(g1));     n2 = sum(~isnan(g2));
    e1 = s1 / sqrt(max(n1,1)); e2 = s2 / sqrt(max(n2,1));

    fprintf('\n=== Pre-learning %% Correct (sessions pooled) ===\n');
    fprintf('Odor 2 (NoStim): mean = %.1f%%, n = %d\n', m1, n1);
    fprintf('Odor 4 (Stim)  : mean = %.1f%%, n = %d\n', m2, n2);
    [~, p_unpaired] = ttest2(g1, g2, 'Vartype','unequal');
    fprintf('Unpaired t-test: p = %.4f\n', p_unpaired);

    % Bar plot + SEM with legend

    figure;
    b = bar(1:2, [m1 m2], 'FaceColor','flat'); hold on;
    b.CData(1,:) = [113 125 146]/255;  % NoStim
    b.CData(2,:) = [166 156 125]/255; % Stim
    
    er = errorbar(1:2, [m1 m2], [e1 e2], 'k.', 'LineWidth', 1.5);
    
    xticks([1 2]);
    xticklabels({'No Stim', 'Stim'});
    ylabel('% correct (pre-learning window)');
    title('Pre-learning Performance');
    legend(b, {'No Stim', 'Stim'}, 'Location', 'northeast');
    
    grid off;
    box on; % keeps border on all four sides
    set(gca, 'LineWidth', 1.5, 'XColor', 'k', 'YColor', 'k'); % black border
    ylim([0 25]);


    % Gardner Altman Plot

    if ~isempty(animal_ids)
        T = table(animal_ids, PLpct_session_all(:,1), PLpct_session_all(:,2), ...
                  'VariableNames', {'animal','pct2','pct4'});
        u = unique(T.animal);
        A = nan(numel(u), 2);
        for k = 1:numel(u)
            rows = strcmp(T.animal, u{k});
            A(k,1) = mean(T.pct2(rows), 'omitnan');
            A(k,2) = mean(T.pct4(rows), 'omitnan');
        end
        A = A(all(~isnan(A),2),:);
        if size(A,1) >= 2
            [~, p_paired] = ttest(A(:,1), A(:,2));
            fprintf('Per-animal paired means: n = %d, p = %.4f\n', size(A,1), p_paired);

            figure; hold on;
            hLines = gobjects(size(A,1),1);
            for r = 1:size(A,1)
                hLines(r) = plot([0 1], [A(r,1) A(r,2)], 'k-', 'LineWidth', 1.5);
            end
            hNoStim = scatter(zeros(size(A,1),1), A(:,1), 400, [113 125 146]/255, 'filled', 'MarkerEdgeColor','k');
            hStim   = scatter(ones(size(A,1),1),  A(:,2), 400, [166 156 125]/255, 'filled', 'MarkerEdgeColor','k');

            xticks([0 1]);
            xticklabels({'No Stim', 'Stim'});
            xlim([-0.5 1.5]); ylim([0 25]);
            grid off; box off;
            ylabel('% correct (pre-learning window)');
            title('Pre-learning % Correct — Per-animal Paired Means');
            legend([hLines(1), hNoStim, hStim], {'Connection', 'No Stim', 'Stim'}, 'Location', 'northeast');
        end
    end
end

% [91 141 184]/255   old colors
% [229 115 115]/255

%% PART:3 PLOTTING LEARNING CURVES : Condition Combined and Aligned (c_al) for all animals in the folder
%Step:1

condToPlot = 'c_al';
ind = findstr(condToPlot, '_');
dataFiles = dir('*_learning_Stim_Feedback.mat');

stimColor = [166 156 125]/255;   
noStimColor  = [113 125 146]/255;

% noStimColor = [91/255, 141/255, 184/255];   
% stimColor   = [229/255, 115/255, 115/255];

for fileIdx = 1:numel(dataFiles)
    fname = dataFiles(fileIdx).name;
    load(fname);
    [~, filename, ~] = fileparts(dataFiles(fileIdx).name);
    animalID = strsplit(filename, '_');
    animalID = animalID{1};

    eval(['lcmat2 = ' animalID '_lc_' condToPlot '2;']);
    eval(['lcmat4 = ' animalID '_lc_' condToPlot '4;']);
    eval(['rtmat2 = ' animalID '_rt_' condToPlot '2;']);
    eval(['rtmat4 = ' animalID '_rt_' condToPlot '4;']);

    nsessions = sum(sum(~isnan(lcmat2')) > 0);
    cond = sum(~isnan(lcmat2));
    ninit = find(cond > 1, 1);
    nfin = find(cond > 1, 1, 'last');
    cond = sum(~isnan(lcmat4));
    ninit = max(ninit, find(cond > 0, 1));
    nfin = min(nfin, find(cond > 1, 1, 'last'));

    if nsessions > 1
        figure('WindowState', 'maximized');

        subplot(2,1,1); hold on;
        shadedErrorBar(1:size(lcmat2, 2), nanmean(lcmat2), nanstd(lcmat2), ...
            'lineProps', {'Color', noStimColor, 'LineWidth', 5});
        shadedErrorBar(1:size(lcmat4, 2), nanmean(lcmat4), nanstd(lcmat4), ...
            'lineProps', {'Color', stimColor, 'LineWidth', 5});
        ylim([0 1.0]);
        xlim([ninit nfin]);
        legend('Non Stimulated', 'Stimulated');
        title(['Learning Curves for ' animalID ', ' ...
            condToPlot(1:ind-1) '\_' condToPlot(ind+1:end) ', ' ...
            num2str(nsessions) ' sessions']);

        subplot(2,1,2); hold on;
        shadedErrorBar(1:size(rtmat2, 2), nanmean(rtmat2), nanstd(rtmat2), ...
            'lineProps', {'Color', noStimColor, 'LineWidth', 5});
        shadedErrorBar(1:size(rtmat4, 2), nanmean(rtmat4), nanstd(rtmat4), ...
            'lineProps', {'Color', stimColor, 'LineWidth', 5});
        ylim([0 2]);
        xlim([ninit nfin]);
        legend('Non Stimulated', 'Stimulated');
        title(['Reaction Times for ' animalID ', ' ...
            condToPlot(1:ind-1) '\_' condToPlot(ind+1:end) ', ' ...
            num2str(nsessions) ' sessions']);
    end
end


%% Plotting Mean and SEM - Shaded Error Bars
%Step:2

condToPlot = 'c_al';

dataFiles = dir('*_learning_Stim_Feedback.mat');

for fileIdx = 1:numel(dataFiles)
    load(dataFiles(fileIdx).name);
    [~, filename, ~] = fileparts(dataFiles(fileIdx).name);
    animalID = strsplit(filename, '_');
    animalID = animalID{1};

    eval(['lcmat2 = ' animalID '_lc_' condToPlot '2;']);
    eval(['lcmat4 = ' animalID '_lc_' condToPlot '4;']);
    eval(['rtmat2 = ' animalID '_rt_' condToPlot '2;']);
    eval(['rtmat4 = ' animalID '_rt_' condToPlot '4;']);

    nsessions=sum(sum(~isnan(lcmat2'))>0);
    cond2=sum(~isnan(lcmat2));
    ninit=find(cond2>1,1);
    nfin=find(cond2>1, 1,'last');
    cond4=sum(~isnan(lcmat4));
    ninit=max(ninit,find(cond4>0, 1));
    nfin=min(nfin,find(cond4>1,1,'last'));
    
    figure;
    subplot(2,1,1); hold on;
    shadedErrorBar(1:size(lcmat2, 2), nanmean(lcmat2), nanstd(lcmat2)./sqrt(cond2), 'lineprops', {'b', 'linewidth', 1});
    shadedErrorBar(1:size(lcmat4, 2), nanmean(lcmat4), nanstd(lcmat4)./sqrt(cond4), 'lineprops', {'r', 'linewidth', 1});
    xlim([ninit nfin]);
    title(['Learning Curves for ' animalID ' ' condToPlot ', ' num2str(nsessions) ' sessions']);

    subplot(2,1,2); hold on;
    shadedErrorBar(1:size(rtmat2, 2), nanmean(rtmat2), nanstd(rtmat2)./sqrt(cond2), 'lineprops', {'b', 'linewidth', 1});
    shadedErrorBar(1:size(rtmat4, 2), nanmean(rtmat4), nanstd(rtmat4)./sqrt(cond4), 'lineprops', {'r', 'linewidth', 1});
    ylim([0.01 2]);
    xlim([ninit nfin]);
    title(['Reaction Times for ' animalID ' ' condToPlot  ', ' num2str(nsessions) ' sessions']);
end

%% T-Test on the Stimulated vs. Non Stimulated LC and RT
%Step:3
% Unpaired T-Test
% Confidence Interval = 95%

condToPlot = 'c_al';

dataFiles = dir('*_learning_Stim_Feedback.mat');
numFiles = numel(dataFiles);

%for fileIdx = 1:numFiles
    fileIdx=2;
    load(dataFiles(fileIdx).name);
    [~, filename, ~] = fileparts(dataFiles(fileIdx).name);
    animalID = strsplit(filename, '_');
    animalID = animalID{1};

    eval(['lcmat2 = ' animalID '_lc_' condToPlot '2;']);
    eval(['lcmat4 = ' animalID '_lc_' condToPlot '4;']);

    nsessions=sum(sum(~isnan(lcmat2'))>0);
    cond2=sum(~isnan(lcmat2));
    ninit=find(cond2>1,1);
    nfin=find(cond2>1, 1,'last');
    cond4=sum(~isnan(lcmat4));
    ninit=max(ninit,find(cond4>0, 1));
    nfin=min(nfin,find(cond4>1,1,'last'));
    
    figure;
    subplot(2,1,1); hold on;
    shadedErrorBar(1:size(lcmat2, 2), nanmean(lcmat2), nanstd(lcmat2)./sqrt(cond2), 'lineprops', {'b', 'linewidth', 1});
    shadedErrorBar(1:size(lcmat4, 2), nanmean(lcmat4), nanstd(lcmat4)./sqrt(cond4), 'lineprops', {'r', 'linewidth', 1});
    xlim([ninit nfin]);
    title(['Learning Curves for ' animalID ' ' condToPlot ', ' num2str(nsessions) ' sessions']);

    d=nanmean(lcmat2)-nanmean(lcmat4);

    [h, p] = ttest(d(ninit:nfin))
%end

%% Pooled P-Values
%Step:4

num_animals = numel(dataFiles);
ttest_results_all = cell(num_animals, 3);

for fileIdx = 1:num_animals
    data = load(dataFiles(fileIdx).name);
    [~, filename, ~] = fileparts(dataFiles(fileIdx).name);
    animalID = strsplit(filename, '_');
    animalID = animalID{1};

    eval(['lcmat2 = ' animalID '_lc_' condToPlot '2;']);
    eval(['lcmat4 = ' animalID '_lc_' condToPlot '4;']);

    
    [h, p] = ttest2(lcmat2(:), lcmat4(:));
    
    ttest_results_all{fileIdx, 1} = animalID;
    ttest_results_all{fileIdx, 2} = h;
    ttest_results_all{fileIdx, 3} = p;
end

%% T-Test on the Stimulated vs. Non Stimulated LC and RT
%Step:5
% Unpaired T-Test
% Confidence Interval = 95%

condToPlot = 'c_al';

noStimColor = [91/255, 141/255, 184/255];   
stimColor   = [229/255, 115/255, 115/255];

dataFiles = dir('*_learning_Stim_Feedback.mat');
numFiles = numel(dataFiles);

    fileIdx = 2;
    load(dataFiles(fileIdx).name);
    [~, filename, ~] = fileparts(dataFiles(fileIdx).name);
    animalID = strsplit(filename, '_');
    animalID = animalID{1};

    eval(['lcmat2 = ' animalID '_lc_' condToPlot '2;']);
    eval(['lcmat4 = ' animalID '_lc_' condToPlot '4;']);

    nsessions = sum(sum(~isnan(lcmat2')) > 0);
    cond2 = sum(~isnan(lcmat2));
    ninit = find(cond2 > 1, 1);
    nfin = find(cond2 > 1, 1, 'last');
    cond4 = sum(~isnan(lcmat4));
    ninit = max(ninit, find(cond4 > 0, 1));
    nfin = min(nfin, find(cond4 > 1, 1, 'last'));
    
    figure;
    subplot(2,1,1); hold on;
    shadedErrorBar(1:size(lcmat2, 2), nanmean(lcmat2), nanstd(lcmat2)./sqrt(cond2), ...
        'lineProps', {'Color', noStimColor, 'LineWidth', 1});
    shadedErrorBar(1:size(lcmat4, 2), nanmean(lcmat4), nanstd(lcmat4)./sqrt(cond4), ...
        'lineProps', {'Color', stimColor, 'LineWidth', 1});
    xlim([ninit nfin]);
    title(['Learning Curves for ' animalID ' ' condToPlot ', ' num2str(nsessions) ' sessions']);

    d = nanmean(lcmat2) - nanmean(lcmat4);
    [h, p] = ttest(d(ninit:nfin))

%% Pooled P-Values
%Step:6

num_animals = numel(dataFiles);
ttest_results_all = cell(num_animals, 3);

for fileIdx = 1:num_animals
    data = load(dataFiles(fileIdx).name);
    [~, filename, ~] = fileparts(dataFiles(fileIdx).name);
    animalID = strsplit(filename, '_');
    animalID = animalID{1};

    eval(['lcmat2 = ' animalID '_lc_' condToPlot '2;']);
    eval(['lcmat4 = ' animalID '_lc_' condToPlot '4;']);

    [h, p] = ttest2(lcmat2(:), lcmat4(:));
    
    ttest_results_all{fileIdx, 1} = animalID;
    ttest_results_all{fileIdx, 2} = h;
    ttest_results_all{fileIdx, 3} = p;
end

%% Plotting Only Learning Curve.
%Step:7

condToPlot = 'c_al';
ind = findstr(condToPlot, '_');
dataFiles = dir('*_learning_Stim_Feedback.mat');
 
stimColor = [166 156 125]/255;   
noStimColor  = [113 125 146]/255;

for fileIdx = 1:numel(dataFiles)
    fname = dataFiles(fileIdx).name;
    load(fname);
    [~, filename, ~] = fileparts(dataFiles(fileIdx).name);
    animalID = strsplit(filename, '_');
    animalID = animalID{1};

    eval(['lcmat2 = ' animalID '_lc_' condToPlot '2;']);
    eval(['lcmat4 = ' animalID '_lc_' condToPlot '4;']);
    eval(['rtmat2 = ' animalID '_rt_' condToPlot '2;']);
    eval(['rtmat4 = ' animalID '_rt_' condToPlot '4;']);

    nsessions = sum(sum(~isnan(lcmat2')) > 0);
    cond = sum(~isnan(lcmat2));
    ninit = find(cond > 1, 1);
    nfin = find(cond > 1, 1, 'last');
    cond = sum(~isnan(lcmat4));
    ninit = max(ninit, find(cond > 0, 1));
    nfin = min(nfin, find(cond > 1, 1, 'last'));

    if nsessions > 1
        figure('WindowState', 'maximized');
        hold on;
        shadedErrorBar(1:size(lcmat2, 2), nanmean(lcmat2), nanstd(lcmat2), ...
            'lineProps', {'Color', noStimColor, 'LineWidth', 5});
        shadedErrorBar(1:size(lcmat4, 2), nanmean(lcmat4), nanstd(lcmat4), ...
            'lineProps', {'Color', stimColor, 'LineWidth', 5});
        ylim([0 1.2]); 
        xlim([ninit nfin]);
        legend('Non Stimulated', 'Stimulated');
        title(['Learning Curves for ' animalID ', ' condToPlot(1:ind-1) '\_' ...
            condToPlot(ind+1:end) ', ' num2str(nsessions) ' sessions']);
    end
end

%% PART:4 Extract N numbers for each animal/session

[mainfolder] = uigetdir; 
cd(mainfolder);

alldatafiles = dir(['*_learning_Stim_Feedback*.mat']);
nfiles = length(alldatafiles);

N_numbers_combined = cell(nfiles, 1);  % for allodors
N_numbers_notCombined = cell(nfiles, 1);  % for allsessions
N_numbers_combined = cell(nfiles, 1);  % for allodors
N_numbers_notCombined = cell(nfiles, 1);  % for allsessions

for nfile = 1:nfiles
    fname = alldatafiles(nfile).name;
    load(fname);
    eval(['allodors = ' animalID '_allodors;']);
    eval(['allsessions = ' animalID '_allsessions;']);

    N_combined = []; 
    N_notCombined = []; 

   
    nodorsets = size(allodors, 1);
    for i = 1:nodorsets
        if allodors{i, 13} == 1  
            nval2 = allodors{i, 10}{2};
            nval4 = allodors{i, 10}{4};
            N_combined = [N_combined; [nval2, nval4]];
        end
    end

   
    nsessions = size(allsessions, 1);
    for i = 1:nsessions
        if allsessions{i, 11} == 1  
            nval2 = allsessions{i, 8}{2};
            nval4 = allsessions{i, 8}{4};
            N_notCombined = [N_notCombined; [nval2, nval4]];
        end
    end

    N_numbers_combined{nfile, 1} = {fname, N_combined};
    N_numbers_notCombined{nfile, 1} = {fname, N_notCombined};
end

% Plotting
% Pooled Average N Number (Learning Trial) with Mean, SD, SEM 

pooled_combined = [];
pooled_notCombined = [];

for i = 1:nfiles
    N_combined = N_numbers_combined{i}{2};
    if ~isempty(N_combined)
        pooled_combined = [pooled_combined; N_combined];
    end
    N_notCombined = N_numbers_notCombined{i}{2};
    if ~isempty(N_notCombined)
        pooled_notCombined = [pooled_notCombined; N_notCombined];
    end
end

mean_combined = mean(pooled_combined, 1, 'omitnan');
std_combined  = std(pooled_combined, [], 1, 'omitnan');
sem_combined  = std_combined ./ sqrt(sum(~isnan(pooled_combined), 1));

mean_notCombined = mean(pooled_notCombined, 1, 'omitnan');
std_notCombined  = std(pooled_notCombined, [], 1, 'omitnan');
sem_notCombined  = std_notCombined ./ sqrt(sum(~isnan(pooled_notCombined), 1));

[h_combined, p_combined] = ttest(pooled_combined(:,1), pooled_combined(:,2));
[h_notCombined, p_notCombined] = ttest(pooled_notCombined(:,1), pooled_notCombined(:,2));

disp(['Combined: h = ' num2str(h_combined) ', p = ' num2str(p_combined)])
disp(['Not Combined: h = ' num2str(h_notCombined) ', p = ' num2str(p_notCombined)])

% nostim_color = [91/255, 141/255, 184/255];
% stim_color   = [229/255, 115/255, 115/255];

nostim_color = [113 125 146]/255;
stim_color   = [166 156 125]/255;

x = [1 2];
figure;
b = bar(x, mean_combined, 'FaceColor', 'flat');
b.CData(1,:) = nostim_color;
b.CData(2,:) = stim_color;
hold on;
errorbar(x, mean_combined, sem_combined, 'k.', 'LineWidth', 2);
hold off;
xticks(x);
xticklabels({'Odor 2 (NoStim)', 'Odor 4 (Stim)'});
ylabel('Trial Number of Learning (N)');
title('Pooled N Number - Combined');
h1 = patch(NaN, NaN, nostim_color);
h2 = patch(NaN, NaN, stim_color);
legend([h1 h2], {'NoStim', 'Stim'}, 'Location', 'northeast');
ylim ([0 200])
%grid on;

figure;
b = bar(x, mean_notCombined, 'FaceColor', 'flat');
b.CData(1,:) = nostim_color;
b.CData(2,:) = stim_color;
hold on;
errorbar(x, mean_notCombined, sem_notCombined, 'k.', 'LineWidth', 2);
hold off;
xticks(x);
xticklabels({'Odor 2 (NoStim)', 'Odor 4 (Stim)'});
ylabel('Trial Number of Learning (N)');
title('Pooled N Number - Not Combined');
h1 = patch(NaN, NaN, nostim_color);
h2 = patch(NaN, NaN, stim_color);
legend([h1 h2], {'NoStim', 'Stim'}, 'Location', 'northeast');
ylim ([0 60])
%grid on;

%% N-Number GA Plot


[mainfolder] = uigetdir;
cd(mainfolder);

alldatafiles = dir(['*_learning_Stim_Feedback*.mat']);
nfiles = length(alldatafiles);

N_numbers_combined = cell(nfiles, 1);
N_numbers_notCombined = cell(nfiles, 1);

for nfile = 1:nfiles
    fname = alldatafiles(nfile).name;
    load(fname);
    eval(['allodors = ' animalID '_allodors;']);
    eval(['allsessions = ' animalID '_allsessions;']);

    N_combined = [];
    N_notCombined = [];

    % === Use FINAL entry that passed in allodors ===
    valid_indices = find(cell2mat(allodors(:,13)) == 1);
    if ~isempty(valid_indices)
        last_idx = valid_indices(end);  % Only final one
        nval2 = allodors{last_idx, 10}{2};
        nval4 = allodors{last_idx, 10}{4};
        N_combined = [N_combined; [nval2, nval4]];
    end

    % === Use FINAL entry that passed in allsessions ===
    valid_indices_sessions = find(cell2mat(allsessions(:,11)) == 1);
    if ~isempty(valid_indices_sessions)
        last_idx = valid_indices_sessions(end);  % Only final one
        nval2 = allsessions{last_idx, 8}{2};
        nval4 = allsessions{last_idx, 8}{4};
        N_notCombined = [N_notCombined; [nval2, nval4]];
    end

    N_numbers_combined{nfile, 1} = {fname, N_combined};
    N_numbers_notCombined{nfile, 1} = {fname, N_notCombined};
end

%% === Collect pooled data ===

pooled_combined = [];
pooled_notCombined = [];

for i = 1:nfiles
    N_combined = N_numbers_combined{i}{2};
    if ~isempty(N_combined)
        pooled_combined = [pooled_combined; N_combined];
    end
    N_notCombined = N_numbers_notCombined{i}{2};
    if ~isempty(N_notCombined)
        pooled_notCombined = [pooled_notCombined; N_notCombined];
    end
end

%% === Paired Line Plot for N-Num Combined Data ===

PLc = pooled_combined;
x = PLc(:, 1);  % No Stim
y = PLc(:, 2);  % Stim

% noStimColor = [91/255, 141/255, 184/255];
% stimColor   = [229/255, 115/255, 115/255];

noStimColor = [113 125 146]/255;
stimColor   = [166 156 125]/255;

figure; hold on;
hLines = gobjects(length(x), 1);
for a = 1:length(x)
    hLines(a) = plot([0, 1], [x(a), y(a)], 'Color', 'k', 'LineWidth', 2);
end

hNoStim = scatter(zeros(size(x)), x, 400, noStimColor, 'filled', 'LineWidth', 1);
hStim   = scatter(ones(size(y)), y, 400, stimColor, 'filled', 'LineWidth', 1);

hold off;
xlabel('No Stim and Stim');
ylabel('N');
title('Learning Trial Numbers (N): Gardner-Altman Paired Plot');
legend([hLines(1), hNoStim, hStim], {'Connection', 'No Stim', 'Stim'}, 'Location', 'northeast');
xticks([0 1]);
xticklabels({'No Stim', 'Stim'});
%grid on;
xlim([-0.5, 1.5]);
ylim([0, max([x; y]) + 20]);

% % Statistics
% [h_N_GA, p_N_GA] = ttest(x, y);
% disp(['Paired t-test (Combined N): p = ', num2str(p_N_GA)]);

%% PART:5 FINAL PERFORMANCE for COMBINED ALIGNED 
% Step:1 Sorting the Final Performance 

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

%% PART:4 FINAL PERFORMANCE CALCULATION AND PLOTTING
% Step:1 Sorting the Final Performance

perf=fin_perf;
mperf=mean(perf);

[h,p]=ttest(perf(:,1),perf(:,2)); % Paired T Test
pw=signrank(perf(:,1),perf(:,2)); % Wilcox Rank Sum Test
[p_kw, tbl_kw, stats_kw] = kruskalwallis(perf); %Kruskal Wallice

% Calculating Mean, SD, SEM

FP_o2 = fin_perf(:,1);
FP_o4 = fin_perf(:,2);

avg_FP_o2 = mean(FP_o2);
sd_FP_o2 = std (FP_o2);
sem_FP_o2 = sd_FP_o2/ sqrt(length(FP_o2));

avg_FP_o4 = mean(FP_o4);
sd_FP_o4 = std (FP_o4);
sem_FP_o4 = sd_FP_o4/ sqrt(length(FP_o4));

%%  Final Performance in Pastel Colors 

means = [avg_FP_o2, avg_FP_o4];
std_devs = [sd_FP_o2, sd_FP_o4];
std_errors = [sem_FP_o2, sem_FP_o4];

% nostim_color = [91/255, 141/255, 184/255];  % Darker pastel blue
% stim_color   = [229/255, 115/255, 115/255]; % Darker pastel red

 stim_color  = [166 156 125]/255;  % Darker pastel brown
 nostim_color  = [113 125 146]/255; % Darker pastel purple

x_val = 1:2;

h = bar(x_val, means, 'FaceColor', 'flat');
h.CData(1, :) = nostim_color;
h.CData(2, :) = stim_color;

hold on;
errorbar(x_val, means, std_errors, 'k.', 'LineWidth', 1.5);
hold off;

xlabel('Odor: 2 and Odor: 4');
ylabel('Final Performance');
title('Final Performance : Combined');

h1 = patch(NaN, NaN, nostim_color);
h2 = patch(NaN, NaN, stim_color);
legend([h1 h2], {'Odor 2', 'Odor 4'}, 'Location', 'northeast');

xticks(x_val);
xticklabels({'a', 'b'});
%grid on;
ylim([0, 1.5]);

% Gardner Altman Plot showing Pairwise Connections

x = fin_perf(:, 1);
y = fin_perf(:, 2);

figure;
hold on;

hLines = gobjects(length(x), 1);
for a = 1:length(x)
    hLines(a) = plot([0, 1], [x(a), y(a)], 'Color', 'k', 'LineWidth', 2);
end

hOdor2 = scatter(zeros(size(x)), x, 400, nostim_color, 'filled', 'LineWidth', 1);
hOdor4 = scatter(ones(size(y)), y, 400, stim_color, 'filled', 'LineWidth', 1);

hold off;

xlabel('Odor: 2 and Odor: 4');
ylabel('Final Performance');
title('Final Performance : Combined');

legend([hLines(1), hOdor2, hOdor4], ...
    {'Connection', 'Odor 2', 'Odor 4'}, ...
    'Location', 'best');

xticks([0 1]);
xticklabels({'Odor 2', 'Odor 4'});
%grid on;
xlim([-0.5, 1.5]);
ylim([0, 1.5]);

% Paired t-test for Final Performance

[h_FP, p_FP] = ttest(fin_perf(:, 1), fin_perf(:, 2));

disp(['Final Performance t-test: h_FP = ' num2str(h_FP) ', p_FP = ' num2str(p_FP)]);


%% Generalized Linear Model for Animals (Odor:2 and Odor:4)

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

%% SIGMOID CURVE FITTING FOR THE LEARNING CURVE
% RUN THIS FOR D1 

clc

load('D21_learning_Stim_Feedback.mat');

D1_mean_2 = D21_meanLc2c{3, 1};
coeffs_2  = D21_meanLc2c{6, 1};
gof_2     = D21_meanLc2c{7, 1};

D1_mean_4 = D21_meanLc4c{3, 1};
coeffs_4  = D21_meanLc4c{6, 1};
gof_4     = D21_meanLc4c{7, 1};

% noStimColor = [91/255, 141/255, 184/255];   % Deep pastel blue
% stimColor   = [229/255, 115/255, 115/255];  % Deep pastel red

 noStimColor = [166 156 125]/255;  % Darker pastel brown
 stimColor   = [113 125 146]/255; % Darker pastel purple

figure;

scatter(1:length(D1_mean_2), D1_mean_2, 30, 'o', 'MarkerEdgeColor', noStimColor);
hold on;
scatter(1:length(D1_mean_4), D1_mean_4, 30, 'o', 'MarkerEdgeColor', stimColor);

x2 = 1:length(D1_mean_2);
x4 = 1:length(D1_mean_4);

plot(x2, coeffs_2.a ./ (1 + exp(-coeffs_2.b * (x2 - coeffs_2.c))) + coeffs_2.d, ...
    '-', 'Color', noStimColor, 'LineWidth', 3);

plot(x4, coeffs_4.a ./ (1 + exp(-coeffs_4.b * (x4 - coeffs_4.c))) + coeffs_4.d, ...
    '-', 'Color', stimColor, 'LineWidth', 3);

xlabel('Trial Number');
ylabel('Learning Performance');

legend('Odor:2 No-Stim', 'Odor:4 Stim', ...
       'Sigmoid Fit for Odor:2', 'Sigmoid Fit for Odor:4', ...
       'Location', 'southeast');

title('Sigmoid Curve Fitting for Learning Curve');
%grid on; 
box on;


%% D2 Sigmoid: ONLY RUN THIS WHILE DOING D2

clc

load('D74_learning_Stim_Feedback.mat');

D2_mean_2 = D74_meanLc2c{3, 1};
coeffs_2  = D74_meanLc2c{6, 1};
gof_2     = D74_meanLc2c{7, 1};

D2_mean_4 = D74_meanLc4c{3, 1};
coeffs_4  = D74_meanLc4c{6, 1};
gof_4     = D74_meanLc4c{7, 1};

noStimColor = [91/255, 141/255, 184/255];
stimColor   = [229/255, 115/255, 115/255];

figure;

hScatter1 = scatter(1:length(D2_mean_2), D2_mean_2, 30, 'o', 'MarkerEdgeColor', noStimColor);
hold on;
hScatter2 = scatter(1:length(D2_mean_4), D2_mean_4, 30, 'o', 'MarkerEdgeColor', stimColor);

x2 = 1:length(D2_mean_2);
x4 = 1:length(D2_mean_4);

hFit2 = plot(x2, coeffs_2.a ./ (1 + exp(-coeffs_2.b * (x2 - coeffs_2.c))) + coeffs_2.d, ...
    '-', 'Color', noStimColor, 'LineWidth', 3);
hFit4 = plot(x4, coeffs_4.a ./ (1 + exp(-coeffs_4.b * (x4 - coeffs_4.c))) + coeffs_4.d, ...
    '-', 'Color', stimColor, 'LineWidth', 3);

xlabel('Trial Number');
ylabel('Learning Performance');

legend([hScatter1, hScatter2, hFit2, hFit4], ...
    {'Odor:2 No-Stim', 'Odor:4 Stim', ...
     'Sigmoid Fit for Odor:2', 'Sigmoid Fit for Odor:4'}, ...
    'Location', 'southeast');

title('D2 Animal - Sigmoid Curve Fitting for Learning Curve');
grid on; box on;



%% D1_DMS for second project

clc

load('D116_learning_Stim_Feedback.mat');

% Extract mean learning curves and sigmoid coefficients
D1_mean_2 = D116_meanLc2c{3, 1};
coeffs_2  = D116_meanLc2c{6, 1};
gof_2     = D116_meanLc2c{7, 1};

D1_mean_4 = D116_meanLc4c{3, 1};
coeffs_4  = D116_meanLc4c{6, 1};
gof_4     = D116_meanLc4c{7, 1};

% Define deep pastel colors
noStimColor = [91/255, 141/255, 184/255];   % Deep pastel blue
stimColor   = [229/255, 115/255, 115/255];  % Deep pastel red

% Plot raw scatter points
figure;

scatter(1:length(D1_mean_2), D1_mean_2, 30, 'o', 'MarkerEdgeColor', noStimColor);
hold on;
scatter(1:length(D1_mean_4), D1_mean_4, 30, 'o', 'MarkerEdgeColor', stimColor);

% X values for fits
x2 = 1:length(D1_mean_2);
x4 = 1:length(D1_mean_4);

% Plot sigmoid fits
plot(x2, coeffs_2.a ./ (1 + exp(-coeffs_2.b * (x2 - coeffs_2.c))) + coeffs_2.d, ...
    '-', 'Color', noStimColor, 'LineWidth', 3);

plot(x4, coeffs_4.a ./ (1 + exp(-coeffs_4.b * (x4 - coeffs_4.c))) + coeffs_4.d, ...
    '-', 'Color', stimColor, 'LineWidth', 3);

% Labels, legend, title
xlabel('Trial Number');
ylabel('Learning Performance');

legend('Odor:2 No-Stim', 'Odor:4 Stim', ...
       'Sigmoid Fit for Odor:2', 'Sigmoid Fit for Odor:4', ...
       'Location', 'southeast');

title('Sigmoid Curve Fitting for Learning Curve');
grid on; box on;

%% D1_DLS

clc

load('D91_learning_Stim_Feedback.mat');

% Extract mean learning curves and sigmoid coefficients
D1_mean_2 = D91_meanLc2c{3, 1};
coeffs_2  = D91_meanLc2c{6, 1};
gof_2     = D91_meanLc2c{7, 1};

D1_mean_4 = D91_meanLc4c{3, 1};
coeffs_4  = D91_meanLc4c{6, 1};
gof_4     = D91_meanLc4c{7, 1};

% % Define deep pastel colors
% noStimColor = [91/255, 141/255, 184/255];   % Deep pastel blue
% stimColor   = [229/255, 115/255, 115/255];  % Deep pastel red

% Define deep pastel colors
noStimColor =  [113 125 146]/255;   % Deep pastel blue
stimColor   = [166 156 125]/255;  % Deep pastel red

% Plot raw scatter points
figure;

scatter(1:length(D1_mean_2), D1_mean_2, 30, 'o', 'MarkerEdgeColor', noStimColor);
hold on;
scatter(1:length(D1_mean_4), D1_mean_4, 30, 'o', 'MarkerEdgeColor', stimColor);

% X values for fits
x2 = 1:length(D1_mean_2);
x4 = 1:length(D1_mean_4);

% Plot sigmoid fits
plot(x2, coeffs_2.a ./ (1 + exp(-coeffs_2.b * (x2 - coeffs_2.c))) + coeffs_2.d, ...
    '-', 'Color', noStimColor, 'LineWidth', 3);

plot(x4, coeffs_4.a ./ (1 + exp(-coeffs_4.b * (x4 - coeffs_4.c))) + coeffs_4.d, ...
    '-', 'Color', stimColor, 'LineWidth', 3);

% Labels, legend, title
xlabel('Trial Number');
ylabel('Learning Performance');

legend('Odor:2 No-Stim', 'Odor:4 Stim', ...
       'Sigmoid Fit for Odor:2', 'Sigmoid Fit for Odor:4', ...
       'Location', 'southeast');

title('Sigmoid Curve Fitting for Learning Curve');
grid on; box on;
