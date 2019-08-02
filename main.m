clear all
% dir_path='F:\ORC\2018_11_08_orc\WAB1_cropped\';
% dir_path='E:\Mustafa\ORC_mitosis\WAB1_cropped\';
dir_path='E:\Mustafa\ORC_mitosis\WT_cropped\'; % directory where data is saved

redo_L=0; % set to recalculate label matrices
level=0.42; % adaptive thresholding level
minA=2000; % minimum area(#pixels, filter for binary mask)6

max_linking_distance=50; %max linking distancs for tracking
max_gap_closing=0; %number of gap frames allowed for tracking

redo_track=0;
dir_list=dir(dir_path); dirFlags=[dir_list.isdir]; % find the directories
dir_list=dir_list(dirFlags);  dir_list(1:2)=[]; % extract the directory list and remove the . and .. listings

redo_plots=0; 

dir_names={dir_list.name};d=[]; d=find(strcmp(dir_names,'data')); if ~isempty(d); dir_list(d)=[]; end %remove data folder from list

%%
for cfolder=1:length(dir_list) % loop through directories

    cfolder
    cpath=strcat(dir_path,dir_list(cfolder).name,'\'); % folder path
    sav_path=[dir_path,'data\',dir_list(cfolder).name];
    dat_path=[sav_path,'_data\'];
    lpath=[sav_path,'_Ltif\']; %folder for saving initial label images
    tracked_path=[sav_path,'_L_tracktif\']; %folder for saving tracked label mimages; 
    
 
    cfilesA=dir(strcat(cpath,'*CamA*.tif'));
    cfilesB=dir(strcat(cpath,'*CamB*.tif'));
    mipA_path=[dir_path,dir_list(cfolder).name,'_mip_CamA.tif']; cmipA=loadtiff(mipA_path);
    mipB_path=[dir_path,dir_list(cfolder).name,'_mip_CamB.tif']; cmipB=loadtiff(mipB_path);
    
    t_s=str2double(cfilesA(1).name(end-20:end-11)); %starting time
    for ct=1:size(cfilesA,1) % go through file names
        t_vec(ct,1)=str2double(cfilesA(ct).name(end-20:end-11))-t_s; % make time vector
    end
    t_vec=t_vec/1000; %convert to seconds 
    
    
    if redo_L || ~isdir(lpath) % if redo_L or segmentation isn't done
        segment_images(cfilesA,level, minA, cpath,sav_path)
        
    end
    'segmentation done'
    if redo_track || ~isdir(tracked_path) %% if redo_track or trackign isn't done
        track_nucs(cpath,sav_path,max_linking_distance,max_gap_closing)
    end
    'tracking done'

    


%% calculate intensities

cL=dir(strcat(tracked_path,'*CamA*.tif'));
clear points nucs_list
parfor ct=1:size(cL,1)  % go through each time
%     ct=10
    cstackL=loadtiff([tracked_path,cL(ct).name]);  % load in label matrix
    cstackA=loadtiff([cpath,cfilesA(ct).name]);     
    cstackB=loadtiff([cpath,cfilesB(ct).name]);     
    propsA=regionprops(cstackL,cstackA,'MeanIntensity'); % calculate region props
    propsB=regionprops(cstackL,cstackB,'MeanIntensity');
    propsL=regionprops(cstackL,'Centroid','Area','PixelIdxList','BoundingBox');
    point_temp=[];
    % calculate background from mipsA
    mipAtmp=cmipA(:,:,ct); bgA(ct,1)=prctile(mipAtmp(:),2);
    mipBtmp=cmipB(:,:,ct); bgB(ct,1)=prctile(mipBtmp(:),2);
    ctr=1;
    for cp=1:size(propsA,1)  
        if propsL(cp).Area>0
        point_temp(ctr,1)=t_vec(ct); % time
        point_temp(ctr,2)=propsA(cp).MeanIntensity; % channelA intensity
        point_temp(ctr,3)=propsB(cp).MeanIntensity; % channelB intensity
        point_temp(ctr,4)=mean(cstackL(propsL(cp).PixelIdxList)); % label
        point_temp(ctr,5:7)=propsL(cp).Centroid; % centroid position     
        point_temp(ctr,8:13)=propsL(cp).BoundingBox;
        point_temp(ctr,14)=ct; % frame
        point_temp(ctr,15)=bgA(ct,1); % backgroundA
        point_temp(ctr,16)=bgB(ct,1); % backgroundB
        ctr=ctr+1;
        end
    end 
    nucs_list{ct}=point_temp; % make a nucs list cell array 
end
nucs_list= vertcat(nucs_list{:});
%% make nucs data list
clear data
nucs=unique(nucs_list(:,4));
min_en=1.06; % filter out short trajectories that end before anaphase based on dynamic range of change
max_l=15;
ctr=1;
for k=1:length(nucs)
    cnuc=nucs_list(nucs_list(:,4)==nucs(k),:);
    if max(cnuc(:,3))/min(cnuc(:,3)) > min_en && size(cnuc,1)>max_l
        data(ctr).tvec=cnuc(:,1); 
        data(ctr).Frame=cnuc(:,14); 
        data(ctr).ORC=cnuc(:,3);
        data(ctr).H2B=cnuc(:,2);
        data(ctr).Cent=cnuc(:,5:7);
        data(ctr).BBox=cnuc(:,8:13);
        data(ctr).cpath=cpath;
        data(ctr).mipA_path=mipA_path;
        data(ctr).mipB_path=mipB_path;
        data(ctr).dat_path=dat_path;
        data(ctr).label=nucs(k); % label
        data(ctr).ORCbg=cnuc(:,16); % BackgroundB
        data(ctr).H2Bbg=cnuc(:,15); % BackgroundA
        ctr=ctr+1;
    end        
end
save([sav_path,'_data.mat'],'data');
%%
if redo_plots || ~isdir(dat_path)
    mkdir(dat_path)
    'saving plots'
for k=1:size(data,2)
    
    figure ('Visible','off','Position',[25,512,2000,600])

%     k=12;
    [m,I]=max(data(k).ORC); max_fr=data(k).Frame(I);    
    tb=20; tm=5; tnorm=12; 
    if (I+tm)> size(data(k).Frame,1);  tm=size(data(k).Frame,1)-I;end;  if (I-tb)<1;  tb=I-1;end
    if (I-tnorm)<1;  normI=1; else;  normI=I-tnorm; end
    
    ctr=1;
 
    for cI=I-tb:1:I+tm
       cfr=data(k).Frame(cI);
       subplot(3,size(I-tb:1:I+tm,2),ctr),...
        imagesc(cmipA((round(data(k).BBox(cI,2))-1):(round(data(k).BBox(cI,2))+data(k).BBox(cI,5)+1),...
        (round(data(k).BBox(cI,1))-1):(round(data(k).BBox(cI,1))+data(k).BBox(cI,4)+1),...
        data(k).Frame(cI))), axis image
        colormap gray
        box off
        axis off       
        if cI==I;title('max'); elseif cI==normI; title('norm'); else; title(data(k).tvec(cI)); end
        
        subplot(3,size(I-tb:1:I+tm,2),size(I-tb:1:I+tm,2)+ctr),...
        imagesc(cmipB((round(data(k).BBox(cI,2))-1):(round(data(k).BBox(cI,2))+data(k).BBox(cI,5)+1),...
        (round(data(k).BBox(cI,1))-1):(round(data(k).BBox(cI,1))+data(k).BBox(cI,4)+1),...
        data(k).Frame(cI))), axis image
        colormap gray
        box off
        axis off
        ctr=ctr+1;
    end
    
    subplot(3,size(I-tb:1:I+tm,2),(ctr*2-1):(ctr*2-1)+size(I-tb:1:I+tm,2)-1),plot(data(k).tvec,data(k).ORC./data(k).ORC(normI),'b'); %/max(data(k).ORC));
    hold on
    plot(data(k).tvec,data(k).H2B./data(k).H2B(normI),'r');
    hold off
    xlim([data(k).tvec(I-tb),data(k).tvec(I+tm)])
    print([dat_path,'nuc_',num2str(data(k).label),'data.bmp'],'-dbmp');
    close all
end

end
end
