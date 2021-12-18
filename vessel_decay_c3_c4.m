% The output files: the 1st row of R intensity decay.csv and G intensity decay.csv are used to obtain C3 and C4, the attenuation coefficients in the bone marrow.
% The output files also include the thickness correcte vessel stacks for visualization:
% G_thickness corrected vessels.tif and R_thickness corrected vessels.tif

clear all; close all; clc;
%% settings
FileName_depth='BM depth map.tif';% depth map
%ss=21; %starting slice that has signals;
StackName_G='G_BM_thickness corrected.tif';% thickness corrected stack
StackName_R='R_BM_thickness corrected.tif';% thickness corrected stack
FileName_mask='Vessel_mask.tif';% segmented vessels
FileName_GG='G intensity decay';% saved file name
FileName_RR='R intensity decay';
MaxDepth=63;% max. depth in the BM

%% Load depth map
tiff_info4 = imfinfo(FileName_depth); % return tiff structure, one element per image
tiff_stack4 = imread(FileName_depth, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info4, 1)
    temp_tiff = imread(FileName_depth, ii);
    tiff_stack4 = cat(3 , tiff_stack4, temp_tiff);
end
stack_depth=double(tiff_stack4);
nTifImages=size(tiff_info4,1);
nRows = getfield(tiff_info4, 'Height'); nCols = getfield(tiff_info4, 'Width'); %xy dimension

%% Load thickness corrected stacks
tiff_info6 = imfinfo(StackName_G); % return tiff structure, one element per image
tiff_stack6 = imread(StackName_G, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info6, 1)
    temp_tiff = imread(StackName_G, ii);
    tiff_stack6 = cat(3 , tiff_stack6, temp_tiff);
end
stack_G=double(tiff_stack6);

tiff_info8 = imfinfo(StackName_R); % return tiff structure, one element per image
tiff_stack8 = imread(StackName_R, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info8, 1)
    temp_tiff = imread(StackName_R, ii);
    tiff_stack8 = cat(3 , tiff_stack8, temp_tiff);
end
stack_R=double(tiff_stack8);

%% load vessel mask and generate thickness corrected vessel stacks
mask_tiff_info = imfinfo(FileName_mask);
mask_stack = imread(FileName_mask, 1) ;
for m = 2 : size(mask_tiff_info, 1)
    temp_tiff = imread(FileName_mask, m);
    mask_stack = cat(3 , mask_stack, temp_tiff);% load mask stack
end

% Red channel: only show the vessel using the mask
stack_R_vessels=zeros(nRows,nCols,nTifImages);
for k=1:nTifImages;
stack_R_vessels(:,:,k)=stack_R(:,:,k) & mask_stack (:,:,k);
stack_R_vessels(:,:,k)=stack_R_vessels(:,:,k).*stack_R(:,:,k);
end
stack_R=stack_R_vessels;

% Greened channel: only show the vessel using the mask
stack_G_vessels=zeros(nRows,nCols,nTifImages);
for k=1:nTifImages;
stack_G_vessels(:,:,k)=stack_G(:,:,k) & mask_stack (:,:,k);
stack_G_vessels(:,:,k)=stack_G_vessels(:,:,k).*stack_G(:,:,k);
end
stack_G=stack_G_vessels;
% save the thickness corrected vessel stacks
IG_vessels = zeros(nRows,nCols,'uint8');% save Red channel (BM only) corrected for bone thickness
        for ff = 1:nTifImages; %loop over image
            IG_vessels(:,:,ff)=stack_G(:,:,ff);
            StackName='G_thickness corrected vessels.tif';
            StackImage = IG_vessels(:,:,ff);
            imwrite(StackImage,StackName,'WriteMode', 'append','Compression','none');      
        end
        IR_vessels = zeros(nRows,nCols,'uint8');% save Red channel (BM only) corrected for bone thickness
        for ff = 1:nTifImages; %loop over image
            IR_vessels(:,:,ff)=stack_R(:,:,ff);
            StackName='R_thickness corrected vessels.tif';
            StackImage = IR_vessels(:,:,ff);
            imwrite(StackImage,StackName,'WriteMode', 'append','Compression','none');      
        end


%% Generate the vessel intensity decay with depth
GG=zeros(2,MaxDepth/3);

for dd=1:MaxDepth/3
Temp=stack_depth==dd*3;
depth_plot_G=stack_G(Temp);
GG(1,dd)=mean(nonzeros(depth_plot_G)); 
GG(2,dd)=std(nonzeros(depth_plot_G)); 
end
xlswrite(FileName_GG, GG); 

RR=zeros(2,MaxDepth/3);

for dd=1:MaxDepth/3
Temp=stack_depth==dd*3;
depth_plot_R=stack_R(Temp);
RR(1,dd)=mean(nonzeros(depth_plot_R)); 
RR(2,dd)=std(nonzeros(depth_plot_R)); 
end
xlswrite(FileName_RR, RR); 

display(sprintf('The files were saved.'))


