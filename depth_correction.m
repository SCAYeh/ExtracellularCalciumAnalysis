
% the thickness corrected files will be used as input for obtaining the attenuation coefficents (C3, C4) in the bone marrow.
% C3, C4, derived as described in detail in Methods, will be used to run section II (Line 170)
% the final depth-corrected files shown below will be used as input for clearance correction:
% G_BM_depth corrected for clearance analyses_T1.tif;
% R_BM_depth corrected for clearance analyses_T1.tif;
% G_BM_depth corrected for clearance analyses_T2.tif;
% R_BM_depth corrected for clearance analyses_T2.tif;

% the final depth-corrected files need to be converted to 32-bit then multiplied by the TotalScaleFactor (saved in the excel output: Scale factor T1 or Scale factor T2)        

clc;clear all; close all;
%% Settings
FileName='R4_T1.Tif';%raw data
FileName_T='T1';% time point of a given ROI 
FileName_BGR='R_BG.Tif';%for subtracting background from the red channel 
FileName_BGG='G_BG.Tif';%for subtracting background from the green channel
FileName_thickness='bone thickness map.tif';% thickness map
FileName_depth='BM depth map.tif';% depth map
FileName_mask='RedMask.tif';
C1=0.02454695; % bone thickness correction constant for red channel 
C2=0.02628571;% bone thickness correction constant for green channel
C3=0.035420099;% depth correction for red channel
C4=0.038666667;% depth correction for green channel
zstep=3; %step size (microns)
leakage_corr=[0.01 0.12];% G leak to R and R leak to G, respectively

%% Load original files
tiff_info = imfinfo(FileName); % return tiff structure, one element per image
tiff_stack = imread(FileName, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(FileName, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end
nTifImages=size(tiff_info,1)/3;
nRows = getfield(tiff_info, 'Height'); nCols = getfield(tiff_info, 'Width'); %xy dimension

R=tiff_stack(:,:,1:3:end);% red channel
G=tiff_stack(:,:,2:3:end);% green channel
R=double(R);
G=double(G);
%% Load depth map
tiff_info4 = imfinfo(FileName_depth); % return tiff structure, one element per image
tiff_stack4 = imread(FileName_depth, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info4, 1)
    temp_tiff = imread(FileName_depth, ii);
    tiff_stack4 = cat(3 , tiff_stack4, temp_tiff);
end
stack_depth=double(tiff_stack4);

%% Load background map
tiff_info6 = imfinfo(FileName_BGG); % return tiff structure, one element per image
tiff_stack6 = imread(FileName_BGG, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info6, 1)
    temp_tiff = imread(FileName_BGG, ii);
    tiff_stack6 = cat(3 , tiff_stack6, temp_tiff);
end
stack_BGG=double(tiff_stack6);

tiff_info7 = imfinfo(FileName_BGR); % return tiff structure, one element per image
tiff_stack7 = imread(FileName_BGR, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info7, 1)
    temp_tiff = imread(FileName_BGR, ii);
    tiff_stack7 = cat(3 , tiff_stack7, temp_tiff);
end
stack_BGR=double(tiff_stack7);


%% subtract DC from the whole image
R_noDC= zeros(nRows,nCols, nTifImages);

for dd=1:nTifImages
R_noDC(:,:,dd) = R(:,:,dd) - stack_BGR(:,:,dd);
end

G_noDC= zeros(nRows,nCols, nTifImages);

for dd=1:nTifImages
G_noDC(:,:,dd) = G(:,:,dd) - stack_BGG(:,:,dd);
end

%% leakage correction to save stacks before depth correction
R_noDC1=R_noDC-G_noDC*leakage_corr(1);% correct for signal crosstalk, 1% G leaked to R
G_noDC1=G_noDC-R_noDC*leakage_corr(2);% correct for signal crosstalk, 12% R leaked to G

%% Save R and G stacks (DC subtracted, leakage corrected, before depth correction)            
        IR = zeros(nRows,nCols,'uint8');
        for d = 1:nTifImages; %loop over image
            IR(:,:,d)= R_noDC1(:,:,d);          
            StackName=['R_before correction_' FileName];
            StackImage = IR(:,:,d);
            imwrite(StackImage,StackName,'WriteMode', 'append','Compression','none');            
        end
        
        IG = zeros(nRows,nCols,'uint8');
        for d = 1:nTifImages; %loop over image             
            IG(:,:,d)= G_noDC1(:,:,d);          
            StackName=['G_before corretion_' FileName];        
            StackImage = IG(:,:,d);
            imwrite(StackImage,StackName,'WriteMode', 'append','Compression','none');             
        end
        
%% Section I: Correction for bone thickness (correction based on images before clearance and leakage correction)
Bone_Thickness= imread(FileName_thickness);
Bone_Thickness = double(Bone_Thickness);
A1=Bone_Thickness.*(-C1);
corr_R=exp(A1);
%IR_corr=zeros(nRows,nCols,nTifImages, 'uint8');
IR_corr=zeros(nRows,nCols,nTifImages);
for s=1:nTifImages;
    IR_corr(:,:,s)=R_noDC(:,:,s)./corr_R;
end

% load mask stack
mask_tiff_info = imfinfo(FileName_mask);
mask_stack = imread(FileName_mask, 1) ;
for m = 2 : size(mask_tiff_info, 1)
    temp_tiff = imread(FileName_mask, m);
    mask_stack = cat(3 , mask_stack, temp_tiff);% load mask stack
end

% only show the BM using the mask
IR_corr_BM=zeros(nRows,nCols,nTifImages);
for k=1:nTifImages;
IR_corr_BM(:,:,k)=IR_corr(:,:,k) & mask_stack (:,:,k);
IR_corr_BM(:,:,k)=IR_corr_BM(:,:,k).*IR_corr(:,:,k);
end
        
% for green channel
A2=Bone_Thickness.*(-C2);
corr_G=exp(A2);
IG_corr=zeros(nRows,nCols,nTifImages);
for s=1:nTifImages;
    IG_corr(:,:,s)=G_noDC(:,:,s)./corr_G;
end

% only show the BM using the mask
IG_corr_BM=zeros(nRows,nCols,nTifImages);
for k=1:nTifImages;
IG_corr_BM(:,:,k)=IG_corr(:,:,k) & mask_stack (:,:,k);
IG_corr_BM(:,:,k)=IG_corr_BM(:,:,k).*IG_corr(:,:,k);
end

max1=max(max(max(IR_corr_BM)));%rescale max value to 255
max2=max(max(max(IG_corr_BM)));%rescale max value to 255
m1=[max1 max2];m1=max(m1);
factor=m1/255;
IR_corr_BM=IR_corr_BM./factor;
IG_corr_BM=IG_corr_BM./factor;

% save the corrected stack

IR_new = zeros(nRows,nCols,'uint8');% save Red channel (BM only) corrected for bone thickness
        for d = 1:nTifImages; %loop over image
            IR_new(:,:,d)= IR_corr_BM(:,:,d);
            StackName='R_BM_thickness corrected.tif';
            StackImage = IR_new(:,:,d);
            imwrite(StackImage,StackName,'WriteMode', 'append','Compression','none');      
        end

IG_new = zeros(nRows,nCols,'uint8');% save Red channel (BM only) corrected for bone thickness
        for d = 1:nTifImages; %loop over image
            IG_new(:,:,d)= IG_corr_BM(:,:,d);
            StackName='G_BM_thickness corrected.tif';
            StackImage = IG_new(:,:,d);
            imwrite(StackImage,StackName,'WriteMode', 'append','Compression','none');      
        end
        
        
        
        
%% Section II: Correction for depth to endosteum (run this section after obtaining C3, C4)

% correct red channel
IR_corr_BM = double (IR_corr_BM);
IR_corr_depth=zeros(nRows,nCols,nTifImages);
depth_corr_R=zeros(nRows,nCols,nTifImages);
A3=zeros(nRows,nCols,nTifImages);
for u=1:nTifImages;
    A3(:,:,u)=stack_depth(:,:,u).*(-C3);
    depth_corr_R(:,:,u)=exp(A3(:,:,u));
    IR_corr_depth(:,:,u)=IR_corr_BM(:,:,u)./depth_corr_R(:,:,u);%IR_corr_BM is thickness corrected already
  
end
        
% correct green channel
IG_corr_BM = double (IG_corr_BM);
IG_corr_depth=zeros(nRows,nCols,nTifImages);
depth_corr_G=zeros(nRows,nCols,nTifImages);
A4=zeros(nRows,nCols,nTifImages);
for v=1:nTifImages;
    A4(:,:,v)=stack_depth(:,:,v).*(-C4);
    depth_corr_G(:,:,v)=exp(A4(:,:,v));
    IG_corr_depth(:,:,v)=IG_corr_BM(:,:,v)./depth_corr_G(:,:,v);
end
   
%% leakage correction for depth corrected stacks
IR_corr_depth1=IR_corr_depth-IG_corr_depth*leakage_corr(1);% correct for signal crosstalk, 1% G leaked to R
IG_corr_depth1=IG_corr_depth-IR_corr_depth1*leakage_corr(2);% correct for signal crosstalk, 12% R leaked to G

%% save leakage and depth corrected files for clearance processing later

max1=max(max(max(IR_corr_depth1)));%rescale max value to 255
max2=max(max(max(IG_corr_depth1)));%rescale max value to 255
m1=[max1 max2];m1=max(m1);
factor1=m1/255;
IR_corr_depth1=IR_corr_depth1./factor1;
IG_corr_depth1=IG_corr_depth1./factor1;


IR_new_depth1 = zeros(nRows,nCols,'uint8');% save Red channel (BM only) corrected for bone thickness
        for dd = 1:nTifImages; %loop over image
            IR_new_depth1(:,:,dd)= IR_corr_depth1(:,:,dd);
            StackName=['R_BM_depth corrected for clearance analyses_' FileName_T '.tif'];
            StackImage = IR_new_depth1(:,:,dd);
            imwrite(StackImage,StackName,'WriteMode', 'append','Compression','none');      
        end

IG_new_depth1 = zeros(nRows,nCols,'uint8');% save Red channel (BM only) corrected for bone thickness
        for ff = 1:nTifImages; %loop over image
            IG_new_depth1(:,:,ff)= IG_corr_depth1(:,:,ff);
            StackName=['G_BM_depth corrected for clearance analyses_' FileName_T '.tif'];
            StackImage = IG_new_depth1(:,:,ff);
            imwrite(StackImage,StackName,'WriteMode', 'append','Compression','none');      
        end

TotalScaleFactor=factor*factor1 % the scale factor used to scale the intensity to 255 
ScaleFactor=[factor factor1 TotalScaleFactor]; % the saved 8-bit file needs to x TotalScaleFactor to get the original intensity
FileName_factor=['Scale factor' FileName_T ];
xlswrite(FileName_factor, ScaleFactor); 

        
display(sprintf('The files were saved.'))