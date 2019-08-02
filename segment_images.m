function segment_images(cfilesA,level, minA, cpath,sav_path)
    % generates label matrices from the H2B channel and saves them as 16
    % bit tif stacks
%     sav_path=cpath(1:end-1);
     parfor ct=1:size(cfilesA,1) 
        %%
%         ct=17;
        ct
        cstackA=loadtiff([cpath,cfilesA(ct).name]);      
        cstackAmed=medfilt3(cstackA); %apply median filter 
        cstackAgauss = imgaussfilt3(cstackAmed); %apply gaussian filter
        
        T = adaptthresh(cstackAgauss,level); %adaptive thresholding
        bw1= imbinarize(cstackAgauss,T); %bianrize stack

        bw2=imclearborder(bw1); % clear border
        bw3 = bwareaopen(bw2, minA); %remove small areas
        CC = bwconncomp(bw3); %creates connected component matrix 
        L = labelmatrix(CC); %creates label matrix fromCC 

%         figure,imagesc(cstackA1(:,:,41)), axis image       
%         figure,imagesc(L(:,:,41)), axis image
        if ~ isdir([sav_path,'_Ltif\']); mkdir([sav_path,'_Ltif\']); end
        write3Dtiff(uint16(L),[sav_path,'_Ltif\',cfilesA(ct).name])
        
        
       
    end
end