%%%%%%%%
%%% Processing the calcium imaging stacks
%%% Have two major functions: 
%%% 1) segment bone and 2) segment Rhode5N
%%% requirements: 1) Bone segmentation doesn't include osteocytes
%%% 2) Rhode5N segmentation is smooth and narrow
%%%%%%%


clc; close all hidden; fclose('all'); clearvars;


path = '/Users/AllisonYeh/Desktop/Calcium manuscript/R4';      % the folder that stores all the images 
cd(path)
stacks=dir('*.tif*');    
len=length(stacks);


Diameter_um_min=5.0; 
XY_PxLength=400/1000; 
Z_PxLength=1.0;

Red_BG = 15;        %%% the background of red channels
Green_BG = 23;      %%% the background of green channels
Blue_BG = 27;      %%% the background of green channels


Osx_Thresh=[9;1]; %%% used for Otsu thresholding of osx
Col_Thresh=[9;2]; %%% used for Otsu thresholding of col1
SHG_Thresh = [5;1]; %%% used for SHG thresholding

NOOD = ones(3,3);
SE=strel('disk',3,4);
SE2=strel('disk',2,4);
SE3=strel('disk',10,4);
boneblur = strel('disk',10,8);


%% Calculate px size for kernel size
% Kernel size should be ~1/2 of feature diameter, w/ variance

Diameter_px_min=Diameter_um_min/XY_PxLength;
Kernel_px_min=floor(Diameter_px_min/2-(Diameter_px_min/2/4));

if mod(Kernel_px_min,2)==0  % Must be odd
    Kernel_px_min=Kernel_px_min+1;
end

tophatMask = strel('disk',20);      %%% masks that used to subtract the unevenness of the image.

blur = 1; %%% guassian filter for bone smoothness
Osteo_size = [100,800];



%% test to read in the entire stack 


NumofStacks = len;

size_x  = 500;
size_y  = 1000;

Ori_Red = zeros (size_x,size_y,NumofStacks);
Ori_Green = zeros (size_x,size_y,NumofStacks);

Osx_img = zeros (size_x,size_y,NumofStacks);
Shg_img = zeros (size_x,size_y,NumofStacks);

for ind = 1:NumofStacks
    
    fname = stacks(ind).name;
    image_Input = imread(fname);
    
    %%% read in and preprocess Red channel for Rhod5N
    Ori_Red(:,:,ind) = uint8(image_Input(:,:,1));
    Ori_Green(:,:,ind) = uint8(image_Input(:,:,2));
    
    Osx_img_tempt = single(image_Input(:,:,1));
    Osx_img_tempt = imgaussfilt(Osx_img_tempt,1);       %%% smooth to remove the noise use for SNARF and don't use for Rhod5N
    Osx_img_tempt = Osx_img_tempt - Red_BG;
    Osx_img_tempt(Osx_img_tempt<0)=0;
    %Osx_img_tempt = imgaussfilt(Osx_img_tempt,1);
    %Osx_img_tempt = imtophat(Osx_img_tempt,tophatMask);
    %Osx_img_tempt = Osx_img_tempt+1;
    %log_Osx = log(Osx_img_tempt);
    
    %log_Osx = log_Osx-min(log_Osx(:));
    %log_Osx(log_Osx<0)=0;
    %log_Osx = 255.*log_Osx./(max(log_Osx(:)));
    
      
    OsxImg= Osx_img_tempt;
    %%% to remove the big bright blobs which are most likly to be osteoids and OB
    %OsxImg = Osx_img_tempt - OsxImg;
    %OsxImg(OsxImg<=0) = 0; 
    
    OsxImg = OsxImg./max(OsxImg(:));
    %OsxImg(OsxImg>0.3) = 0.3;             %%% compress the contrast
    %OsxImg = OsxImg./0.3;
    %OsxImg = double(OsxImg);
    
    
    
    Osx_img(:,:,ind) = OsxImg;
    
    %%% read in and preprocess SHG channel to get a SHG mask
    
    SHG = uint8(image_Input(:,:,3));
    SHG = SHG-Blue_BG;
    SHG(SHG<0)=0;
    SHG = double(SHG);
    SHG = imgaussfilt(SHG,blur);    %%% blur the SHG to reduce noise
    
    
    
    
    if max(SHG(:))==0
        SHG = zeros(size_x,size_y);
    else
        SHG = SHG./max(SHG(:));
    end
    
    
%     SHG(SHG>0.7) = 0.7;             %%% compress the contrast
%     SHG = SHG./0.7;
%     SHG = double(SHG);
    %SHG = imdilate(SHG,boneblur);
    %SHG = imerode(SHG,boneblur);
    
    Shg_img(:,:,ind) = SHG;
    
    
end

Osx_img=Osx_img-0.0;
Osx_img(Osx_img<0) = 0;
Osx_img(Osx_img>0.8) = 0.8;     %% compress the contrast           
Osx_img = Osx_img./0.8;
Osx_img = double(Osx_img);

%Rhod5N = Osx_img+0.001;    
%Rhod_log = log(Rhod5N);        %% equalize the histogram with log transform
%Rhod_log = Rhod_log - min( Rhod_log(:));
%Rhod_eq = Rhod_log./max(Rhod_log(:));

Rhod_eq = histeq(Osx_img);       %% backup histogram equalizer
Rhod_eq = Rhod_eq - min(Rhod_eq(:));
Rhod_eq = Rhod_eq./max(Rhod_eq(:));

Shg_img=Shg_img-0;
Shg_img(Shg_img<0) = 0;



%% adjust the bone decay
clc
AdjustedSHG = Mask3DBone (Shg_img);
AdjustedSHG = AdjustedSHG - 0.01;
AdjustedSHG(AdjustedSHG<0) = 0;
AdjustedSHG(AdjustedSHG>0.7) = 0.7;             %%% compress the contrast
AdjustedSHG = AdjustedSHG./0.7;
AdjustedSHG = double(AdjustedSHG);
%AdjustedSHG = Shg_img;

%% Bone  

% for z=1:NumofStacks %%% can consider do z= BMStart:zslice
%             
%             j=figure(1); clf;
%             set(gcf,'Position',[225 420 1015 415]);        
%             subplot(1,2,1);
%             titlelabel=sprintf(' plane %d / %d', z,NumofStacks);
%             imagesc(AdjustedSHG(:,:,z));
%             title(titlelabel);
%             subplot(1,2,2);
%             imshow(Shg_img(:,:,z),[]);
%             title(titlelabel);
%             
%             pause(0.01);
% end
% 
% promptMessage = sprintf('Do you want to continue with this,\n or try with different parameters?');
% button = questdlg(promptMessage, 'Continue?', 'Continue', 'Restart', 'Continue');

%Shg_img = imgaussfilt3(Shg_img,1);

bone_mask = zeros(size_x,size_y,NumofStacks);

bonestart = 3;   %%% the frame when bone is in focus and takes more than 10% of the frame
bone_size1 = 5000;
bone_size2 = 20000;

for ind = bonestart:NumofStacks
       
    bone = AdjustedSHG(:,:,ind);
    %%% generate masks for SHG
    
    %bone_ThreshLevel=multithresh(AdjustedSHG(:), SHG_Thresh(1)); % ***** 2017/09/01
    %bone_ThreshLevel=bone_ThreshLevel(1,SHG_Thresh(2)); % *****  2017/09/01
    %SHGmask=imbinarize(bone,bone_ThreshLevel);
    
    SHGmask=bone>0.01;
    SHGmask = imdilate(SHGmask,SE);
    SHGmask = imerode(SHGmask,SE); 
    
    if ind<=bonestart
        SHGmask = bwareaopen(SHGmask, bone_size1);
    else      
        SHGmask = bwareaopen(SHGmask, bone_size2);
    end
      
   
    %%%fill out the holes that might be osteocytes
    Ost = ~SHGmask;
    Ost = bwareaopen(Ost, 2000);
    %figure, imshow(Ost,[])
    
    SHGmask = ~Ost;
    bone_mask(:,:,ind) = SHGmask;
    %figure, imshow(SHGmask,[]);
    
    
   
end
%figure, imshow(bone_mask(:,:,21),[])



%% Rhode5N segmentation
close all
Rhod5N_mask = zeros(size_x,size_y,NumofStacks);
Rhod5N_mask2 = zeros(size_x,size_y,NumofStacks);
Vess_mask = zeros(size_x,size_y,NumofStacks);
%figure, imshow(Osx_img(:,:,36),[0 30])
offset = 0;
%INITPSF = ones(8,8);
vessErode = strel('disk',5);
bErode = strel('disk',30);

%sig_start = 36;

for ind = 1:NumofStacks
       
    Rhod5N = Rhod_eq(:,:,ind);
    
    %%% get rid of the signals inside the bones if we analyze interstitial    
    %Rhod5N_clean = Rhod5N.*(~temp);
    Rhod5N_clean = Rhod5N;  
%     Rhod5N_c = Rhod5N_clean+0.001;
%     
%     Rhod5N_log = log(Rhod5N_c);
%     Rhod5N_log = Rhod5N_log - min( Rhod5N_log(:));
%     Rhod5N_log = Rhod5N_log./max(Rhod5N_log(:));
    

    
    %%% Interstitial space masking
    %%% 1) BG subtraction
    %figure, imshow(Rhod5N_clean,[]);
    
    %Rhod5N_BGf = imtophat(Rhod5N_clean,tophatMask);
    BG = imopen(Rhod5N_clean,tophatMask);
    BG = imgaussfilt(BG,50);
    Rhod5N_BGf = Rhod5N_clean - 0.8*BG;
    Rhod5N_BGf(Rhod5N_BGf<0) = 0;
    %figure, imshow(Rhod5N_BGf,[]);
    %%% 2) contrast enhacement
    
    Rhod5N_en = Rhod5N_BGf./max(Rhod5N_BGf(:));
    Rhod5N_en = Rhod5N_en - 0.0;
    Rhod5N_en(Rhod5N_en<0) = 0;
    Rhod5N_en(Rhod5N_en>0.7) = 0.7;
    Rhod5N_en = Rhod5N_en./0.7;
    
    %figure, imshow(Rhod5N_en,[]);
    
    T = adaptthresh(Rhod5N_en, 0.4,'NeighborhoodSize',[5 5],'Statistic','gaussian');
    %figure, imshow(T,[]);
    mask1 = imbinarize(Rhod5N_en,T);
    
    %%% 3) filter out the masks
    
    mask = imgaussfilt(single(mask1),1);
    %figure, imshow(mask,[]);
    mask_den = im2bw(mask,0.4);
    %figure, imshow(mask_den,[]);
    %mask_den = bwareaopen(mask_den, 10, 4);
    
    %Rhod5Nmask = bwareaopen(Rhod5Nmask, 10, 4);
    Rhod5Nmask = imdilate(mask1,SE2);
    Rhod5Nmask = imerode(Rhod5Nmask,SE2);
    %Rhod5Nmask = imerode(Rhod5Nmask,NOOD);
    %Rhod5Nmask = bwareaopen(Rhod5Nmask, 200,8);
    Rhod5Nmask = bwareaopen(Rhod5Nmask, 50, 4);
    
    Temp = ~Rhod5Nmask;
    Temp = bwareaopen(Temp, 30,4);    
    Rhod5Nmask_fin = ~Temp;
    %Rhod5Nmask_fin = imerode(Rhod5Nmask_fin,[1,1;1,1]);

    
    %mask2 = imdilate(mask2,NOOD);
    %mask2 = imerode(mask2,NOOD);
    %mask2 = bwareaopen(mask2, 5,8);
    
    %bw = activecontour(J, mask2, 100,'Chan-Vese','SmoothFactor', 0.0);
    %bw = activecontour(Rhod5N_t, Rhod5Nmask_fin, 100,'Chan-Vese','SmoothFactor', 0.0);
    %bw = bwareaopen(bw, 5,8);
    
    %figure, imshow(Rhod5Nmask_fin,[]);
    temp1 = bone_mask(:,:,ind);
    temp = imerode(temp1,SE3);
    Rhod5Nmask_fin = Rhod5Nmask_fin.*(~temp);
    Rhod5N_mask(:,:,ind) = Rhod5Nmask_fin; 
    Rhod5N_mask2(:,:,ind) = mask_den.*(~temp);
    
    
   %%%%%%%% find the vessels

    Rhod_ThreshLevel=multithresh(Rhod5N_en, 9); % ***** 2017/09/01
    Rhod_ThreshLevel=Rhod_ThreshLevel(1,5); % *****  2017/09/01
    Vesselmask=imbinarize(Rhod5N_en,Rhod_ThreshLevel);
    
    Vesselmask = imdilate(Vesselmask,SE2);
    Vesselmask = imerode(Vesselmask,SE2);
    
    
    %figure, imshow(Vesselmask);
    
    Vess_ini = imerode(Vesselmask, vessErode);
    %figure, imshow(Vess_ini)
    Vess_ini = bwareaopen (Vess_ini, 800,4);
    %figure, imshow(Vess_ini)
    
    Vess_dil = imdilate(Vess_ini, vessErode);
    %figure, imshow(Vess_mask)
    
    
    boneM = imdilate(temp,bErode);
    Vess_mask(:,:,ind)= Vess_dil.*(~boneM);
    %figure, imshow(Vess_mask(:,:,ind))
end


%% save masks
% SaveZSARFilePath=strcat(path,'RhodAdj','/');  %%% save rhod5N masks   
% 
%     if exist(SaveZSARFilePath, 'dir')
%         cd(SaveZSARFilePath);
%     else
%         mkdir(SaveZSARFilePath);
%         cd(SaveZSARFilePath);
%     end
%     
% for k = 1:NumofStacks
%     
%         imwrite(Rhod5N_mask(:,:,k),strcat('Rhod5N_',num2str(k,'%03.0f'),'.tiff'),'Compression','none');
%                
% end    
% 
% 
% SaveZSARFilePath=strcat(path,'RhodAdj2','/');  %%% save rhod5N masks   
% 
%     if exist(SaveZSARFilePath, 'dir')
%         cd(SaveZSARFilePath);
%     else
%         mkdir(SaveZSARFilePath);
%         cd(SaveZSARFilePath);
%     end
%     
% for k = 1:NumofStacks
%     
%         imwrite(Rhod5N_mask2(:,:,k),strcat('Rhod5N_',num2str(k,'%03.0f'),'.tiff'),'Compression','none');
%                
% end 
%     
    
%%     
SaveZSARFilePath=strcat(path,'bone','/');   %%% save bone masks

    if exist(SaveZSARFilePath, 'dir')
        cd(SaveZSARFilePath);
    else
        mkdir(SaveZSARFilePath);
        cd(SaveZSARFilePath);
    end
    
for k = 1:NumofStacks
    
        imwrite(bone_mask(:,:,k),strcat('Bone_',num2str(k,'%03.0f'),'.tiff'),'Compression','none');
 
               
end     

%% vessel
SaveZSARFilePath=strcat(path,'Vessel','/');   %%% save bone masks

    if exist(SaveZSARFilePath, 'dir')
        cd(SaveZSARFilePath);
    else
        mkdir(SaveZSARFilePath);
        cd(SaveZSARFilePath);
    end
    
for k = 1:NumofStacks
    
        imwrite(Vess_mask(:,:,k),strcat('Vess_',num2str(k,'%03.0f'),'.tiff'),'Compression','none');
 
               
end   


   

