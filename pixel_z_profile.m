% This code generates the thickness and the depth map
% These maps serve as input for depth_correction.m
% The expected running time in a normal computer is within 10 seconds. 

clc;clear all; close all;
%% Settings
FileName='Bone_mask.tif'; %segmented bone stack
zstep=3; %step size (microns)
%% Load files
tiff_info = imfinfo(FileName); % return tiff structure, one element per image
tiff_stack = imread(FileName, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(FileName, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end
nTifImages=size(tiff_info,1);% the number of z slices (z-1)
nRows = getfield(tiff_info, 'Height'); nCols = getfield(tiff_info, 'Width'); %xy dimension


%% This gets the thickness from bone surface to endosteum
tStart=tic;
map_thickness2 = zeros(nRows,nCols,'uint8');% bone thickness map
for i = 1:nRows;% loop over rows
     for j = 1:nCols;% loop over columns
         idx2=find(tiff_stack(i,j,:)==255);%index of segmented bone
         tf=isempty(idx2);% for regions only with bone
         if tf==1;
             idx2=[1:1:nTifImages];
         end
         
         indexF = idx2(1);%first slice
         for k = 2:length(idx2);
         if idx2(k)- idx2(k-1)>1; 
             idx2(k:end)=0; 
             
         end
         end
         indexL=idx2(length(find(idx2>0)));% stop before entering bone marrow
         map_thickness2(i,j)= (indexL-indexF)*zstep;% 3 um step % no need to minus 1
     end
end
%figure; imshow(map_thickness2);title('bone thickness2');
imwrite(map_thickness2,'bone thickness map.tif')

%% This get the slice number for the endosteum
% The first slice needs to be all black
map_endo = zeros(nRows,nCols,'uint8');
for i = 1:nRows;% loop over rows
     for j = 1:nCols;% loop over columns
         idx3=find(tiff_stack(i,j,:)==0);% bone marrow or black space
         
         tf=isempty(idx3);% for regions only with bone
         if tf==1;
             idx3=[1:1:nTifImages];
         end
     
         for k = 2:length(idx3);
             if idx3(k)- idx3(k-1)>1;% entering endosteum
             idx3(k+1:end)=0; % k+1 is the slice number when the bone marrow just appear
             
             end
             
         end
         
      
         index_endo=idx3(length(find(idx3>0)));% stop before entering bone marrow
         map_endo(i,j)= index_endo;% the slice entering BM
     end
end
% figure; imshow(map_endo);title('endosteum');
% imwrite(map_endo,'endo.tif')

%% Get the depth below endosteum
tiff_stack2 = zeros(nRows,nCols,'uint8');
for m=1:nRows;% e.g. 175
   for n = 1:nCols;%e.g. 525
        for z=1:nTifImages
           
            depth= (z-map_endo(m,n))+1;% 1*3 is the offset and needs to be subtracted
            if depth <= 0;
               depth =1;% only accept positive integers
            end
          
            tiff_stack2(m,n,z)= depth*zstep;
        end
   end
end


%% Mask the depth map to show only from the bone marrow region
BM_mask=imcomplement(tiff_stack);
depth_BM=zeros(nRows,nCols,nTifImages);
tiff_stack2=double(tiff_stack2);
for k=1:nTifImages;
depth_BM(:,:,k)=tiff_stack2(:,:,k) & BM_mask (:,:,k);
depth_BM(:,:,k)=tiff_stack2(:,:,k).*depth_BM(:,:,k);
end

depth_BM=depth_BM-3; % subtract offset

% Save depth stack
I = zeros(nRows,nCols,'uint8');
        for d = 1:nTifImages; %loop over image
           
            
            I(:,:,d)= depth_BM(:,:,d);
           
          
            StackName='BM depth map.tif';
           
          
            StackImage = I(:,:,d);
            imwrite(StackImage,StackName,'WriteMode', 'append','Compression','none');
            
        end
        
        
display(sprintf('The files were saved. Elapsed time : %.3f s.', toc(tStart)))