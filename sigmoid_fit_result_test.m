
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
perf=fin_perf
mperf=mean(perf)
[h,p]=ttest(perf(:,1),perf(:,2))
pw=signrank(perf(:,1),perf(:,2))

%%
clc
for i=1:11
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