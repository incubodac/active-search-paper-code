  %% Create Summary files to use with Mensner TFCE
%This script creates .m summary files from 
%fun_makeStructToERP(ERPstruct,condVec,iSuj,folder,subjName) output,
%to be use with Mensen's TFCE. 
%Damian Care 2/02/21: 

folderAux =  session_path.aux_data;
Nsuj = numel(ERPstruct);
data1 = [];
data2 = [];
Summary=[];

for iSuj = [1:17 19:Nsuj]
             
   % data = ERPstruct(iSuj).tar_O;
   data = ERPstruct(iSuj).(condVec{3}(2:end-4));

   %data = ERPstruct(iSuj).dis_O;
   data1= cat(3,data1,data);
    
end
data1 = permute(data1,[3 1 2]);
Summary = data1;
fileOut1 = [folderAux  'summary250' condVec{3}(2:end-4)];

save(fileOut1,'Summary')
Summary=[];
for iSuj = [1:17 19:Nsuj]
    
   % data = ERPstruct(iSuj).tar_F;
   data = ERPstruct(iSuj).(condVec{4}(2:end-4));

   %data = ERPstruct(iSuj).dis_F;
    data2= cat(3,data2,data);
    
end
data2 = permute(data2,[3 1 2]);
Summary = data2;
fileOut2 = [folderAux  'summary250' condVec{4}(2:end-4)];

save(fileOut2,'Summary')
%% Create location file for TFE (loading some analysis EEG first is required)
%e_loc file containing electrode locations is used by TFCE function
%Damian Care 2/02/21: 

folderAux =  session_path.aux_data;
fileOut = [folderAux  'locationFile'];

e_loc = EEG.chanlocs;

save(fileOut,'e_loc')
%% Simple t-test (load first summary1 and summary2)

ROIchans = {'O2','PO4','PO8'};
ventana  = [190 211];
indchans = find(ismember({chanlocs.labels},ROIchans));
cond1 = Summary1(:,indchans,:);
cond1 = mean(mean(cond1, 2),3);

cond2 = Summary2(:,indchans,:);
cond2 = mean(mean(cond2, 2),3);

[h,p] = ttest(cond1,cond2,'Alpha',0.05) 



%% Main TFCE function
Results = fp.ept_TFCE()

%How to compile mex file on MAC

%Mensen commented on 14 Jun 2018 • 
%edited 
%Hey!

%So matlab finally fixed the issue, you can read about it here

%You can easily fix quickly by using the -compatibleArrayDims when mexing the 
%file (use the older API), but I'm going to fix the mex files themselves when I get a chance... then it should all work out with a new update as well.

%mex ept_mex_TFCE.c -compatibleArrayDims
%% TOPOS +  TFCE significant results for channel 

y=sum(Results.P_Values(:,times>175 & times<225)<0.05,2);

figure; topoplot(y,EEG.chanlocs)


figure; topoplot(y,EEG.chanlocs); colormap bone; colbar

%figure; topoplot(y,EEG.chanlocs,'maplimits'); colormap bone; colorbar

% Topo with amount of significative samples 
figure; topoplot(y,EEG.chanlocs,'maplimits',[0 20]); colormap bone; colorbar
% sighnificative inverted colors
%poster

c = summer

a=figure;topoplot(y,EEG.chanlocs,'maplimits',[0 20]); colormap(flipud(c)); h = colorbar; set(gcf,'Color','w') 
% ylabel(h, 'Amount of Significant Values')

set(h,'TickLabels',[])

figure; topoplot(1-(y>0),EEG.chanlocs,'maplimits',[0 1]); colormap bone; colorbar


%% TFCE  plots poster

%figure(); imagesc(times,1:64,Results.P_Values<0.05); colormap bone; colorbar
%figure(); imagesc(times,1:64,Results.P_Values(Results.P_Values<0.05)); colormap bone; colorbar
%figure(); imagesc(times,1:64,Results.P_Values); colormap bone; colorbar
%figure(); imagesc(times,1:64,Results.P_Values); colormap summer; colorbar
%figure(); imagesc(times,1:64,Results.P_Values([1:32 64:-1:33],:)); colormap summer; colorbar

Y = Results.P_Values; Y(Y>0.05)=1;
%figure(); imagesc(times(-100<times<durMinFix),1:64,Y([1:32 64:-1:33],find(times<250))); colormap summer; colorbar
%figure(); imagesc(times,1:64,Y([1:32 64:-1:33],:),[0 0.05]); colormap summer; colorbar


% Poster plot VSS
figure(); imagesc(1:64,times(-100<times & times<1000*durMinFix),Y([1:32 64:-1:33],find(-100<times & times<1000*durMinFix))',[0 0.05]); colormap hot; colorbar% NTO
%SAN unfold
figure(); imagesc(1:64,times(-.1<times & times<.25),Y([1:32 64:-1:33],find(-.1<times & times<.25))',[0 0.05]); colormap hot; colorbar% NTO

%figure(); imagesc(1:64,times,Y(:,:)',[0 0.05]); colormap hot; colorbar

title('TFCE result','fontsize',22)
ylabel('Time [ms]','fontsize',22)
xlabel('Channels','fontsize',22)
xticks([7 20 27 33 40 55])
xticklabels(nuChans([7 20 27 33 40 55]))%int8(linspace(1,64,6))
%y=sum(Y(:,times>175 & times<225)<0.05,2);
%figure; topoplot(y,EEG.chanlocs,'maplimits',[0 0.05]); colormap bone; colorbar

%y=mean(Y(:,times>175 & times<225),2);
%figure; topoplot(y,EEG.chanlocs,'maplimits',[0 0.05]); colormap bone; colorbar

%y=sum(Y(:,times>175 & times<225)<0.05,2);
%figure; topoplot(y,EEG.chanlocs,'maplimits',[0 20]); colormap bone; colorbar
set(gcf,'Color','w')
set(gca,'fontsize',19)
%set(gca,'XTickLabel',nuChans([7 25 27 33 35 55]),'fontsize',8)
ax = gca
%ax.XAxis.FontSize = 10;


%%
nuChans = {EEG.chanlocs([1:32 64:-1:33]).labels};
%% Poster version
ventanas={[85 115],[160 190],[190 220]}

   
    figure(100); clf
        for i=1:3
            y=[];
           
            ventana = find(times>ventanas{i}(1) & times<ventanas{i}(2));
            ventana
            y=sum(Y(:, ventana)<0.05,2);

            fh(i)=subplot(1,4,i)
            hold on
            str1 = sprintf('%d ms' ,[ventanas{i}(1)])
            str2 = sprintf('%d ms' ,[ventanas{i}(2)])
            title([str1 '-' str2],'fontsize',21)%[ventanas{i}(1) ' ms' ;ventanas{i}(2) ' ms'] )
            topoplot(y,EEG.chanlocs,'maplimits',[0 15]); 
            %colormap bone;
            colormap(flipud(bone));

            if i==3;
                cb=colorbar;
                cb.Position = cb.Position + [.1, -.05, 0.01, .16];
                ylabel(cb, '# of significative samples','FontSize',12)

            end;
            

            hold off
        end
                set(gcf,'Color','w')
%% quantile CI 95%
x=rand(1,100)

quantile(x,[0.025 0.975])

x=rand(13,100)
quantile(x,[0.025 0.975])
%% TFCE for SAN one sample t test 


folderAux =  session_path.aux_data;
subj = fields(unfoldSubjStruct)
data1 = [];
Summary=[];

for iSuj = {subj{:}};
             
   % data = ERPstruct(iSuj).tar_O;
  data = unfoldSubjStruct.(iSuj{1}).beta_dc(:,:,3);
   
  data1= cat(3,data1,data);
    
end
data1 = permute(data1,[3 1 2]);
Summary = data1;
fileOut1 = [folderAux  'SANsecondIntentEasy' ];

save(fileOut1,'Summary')


times=unfoldSubjStruct.E01.times

%% Create location file for TFE (loading some analysis EEG first is required)
folderAux =  session_path.aux_data;
fileOut = [folderAux  'locationFile'];

e_loc = EEG.chanlocs;

save(fileOut,'e_loc')
%%

for i={pp{:}}
    i
end
