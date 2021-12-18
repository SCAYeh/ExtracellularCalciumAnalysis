%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% do a 25x25 moving block clearance decay correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; close all hidden; fclose('all'); clearvars;

size_x  = 500;
size_y  = 1000;

Red_BG = 18;        %%% the background of red channels
Green_BG = 23;      %%% the background of green channels
Blue_BG = 20;      %%% the background of green channels





%% read in images at time point 1
path_t1 = '/Users/allisonyeh/Desktop/Calcium manuscript/Re__Your_manuscript_NCOMMS-21-05892_-_please_provide_further_information';      % the folder that stores all the images at time point 1 
cd(path_t1)
Gt1='G_BM_depth corrected for clearance analyses_T1.tif';
Gt2='G_BM_depth corrected for clearance analyses_T2.tif';

Rt1='R_BM_depth corrected for clearance analyses_T1.tif';
Rt2='R_BM_depth corrected for clearance analyses_T2.tif';



InfoMov=imfinfo(Gt1);
NumofStacks = length(InfoMov);

Redt1 = zeros (size_x,size_y,NumofStacks);
Greent1 = zeros (size_x,size_y,NumofStacks);
Redt2 = zeros (size_x,size_y,NumofStacks);
Greent2 = zeros (size_x,size_y,NumofStacks);


for ind = 1:NumofStacks
    
       
    Redt1(:,:,ind) = single(imread(Rt1, ind));
    Redt2(:,:,ind) = single(imread(Rt2, ind));
    
    %%% normalize to t1
    %Redt1(:,:,ind) = Redt1(:,:,ind)./Redt1(:,:,ind);
    %Redt1(:,:,ind) = Redt1(:,:,ind)./Redt1(:,:,ind);
    
    
    Greent1(:,:,ind) = single(imread(Gt1, ind));
    Greent2(:,:,ind) = single(imread(Gt2, ind));
    
    
end






%% Align the two time points using the green channel. the green channel is more stable
%%%% this section uses the function provided by the Matlab

% [optimizer,metric] = imregconfig('multimodal');
% 
% 
% optimizer.InitialRadius = 1.5e-6;
% optimizer.Epsilon = 1.5e-6;
% optimizer.GrowthFactor = 1.01;
% optimizer.MaximumIterations = 300;
% 
% moving = Greent2(:,:,28);
% fixed = Greent1(:,:,25);
% 
% tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
% 
% movingRegistered = imwarp(moving,tform);
% 
% 
% figure
% imshowpair(fixed, movingRegistered,'Scaling','joint')


%% read the images and calculate the decay coefficient

rollRadius = 20;
xblock = size_x/rollRadius;
yblock = size_y/rollRadius;


t1 = 1;     %% the time the first stack is taken after injection
t2 = 8;     %% the time the second stack was taken after injection
tinj = 0.16;   %% the time we fit it to

t = [t1, t2];

CorrRt1 = zeros (size_x,size_y,NumofStacks);
CorrGt1 = zeros (size_x,size_y,NumofStacks);
CorrRt2 = zeros (size_x,size_y,NumofStacks);
CorrGt2 = zeros (size_x,size_y,NumofStacks);
coeffR = zeros (size_x,size_y);
coeffG = zeros (size_x,size_y);

for ind = 16:NumofStacks
    
       
    currentRt1 = Redt1(:,:,ind);
    currentRt2 = Redt2(:,:,ind);
    
    currentGt1 = Greent1(:,:,ind);
    currentGt2 = Greent2(:,:,ind);
    
    %%% delete the pixels that are close to the image edges
       
    currentRt1 (1:rollRadius,:) = 0;
    currentRt1 (:,1:rollRadius) = 0;
    currentRt1 (size_x-rollRadius:end,:) = 0;
    currentRt1 (:,size_y-rollRadius:end) = 0;
    
    currentGt1 (1:rollRadius,:) = 0;
    currentGt1 (:,1:rollRadius) = 0;
    currentGt1 (size_x-rollRadius:end,:) = 0;
    currentGt1 (:,size_y-rollRadius:end) = 0;
    
    currentRt2 (1:rollRadius,:) = 0;
    currentRt2 (:,1:rollRadius) = 0;
    currentRt2 (size_x-rollRadius:end,:) = 0;
    currentRt2 (:,size_y-rollRadius:end) = 0;
    
    currentGt2 (1:rollRadius,:) = 0;
    currentGt2 (:,1:rollRadius) = 0;
    currentGt2 (size_x-rollRadius:end,:) = 0;
    currentGt2 (:,size_y-rollRadius:end) = 0;
    
    %%%% only sample the pixels that is not 0 in both channels
    R1nz = currentRt1~=0;
    R2nz = currentRt2~=0;
    mask = R1nz&R2nz;
    mask = bwareaopen(mask, 400);
    
    currentRt1 = currentRt1.*single(mask);
    currentRt2 = currentRt2.*single(mask);
    currentGt1 = currentGt1.*single(mask);
    currentGt2 = currentGt2.*single(mask);
    
    %%% down-sample the image to calculate the coefficienct
    
    newRt1 = zeros(xblock,yblock);
    newRt2 = zeros(xblock,yblock);
    newGt1 = zeros(xblock,yblock);
    newGt2 = zeros(xblock,yblock);
    
    for i = 1:xblock            
        for j = 1:yblock
            
            startx = (i-1)*rollRadius+1;
            starty = (j-1)*rollRadius+1;
            
            
            R1ij = currentRt1(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1));
            R2ij = currentRt2(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1));
            G1ij = currentGt1(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1));
            G2ij = currentGt2(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1));
            
            if nnz(R1ij)==0
                newRt1(i,j) = 0;
                newRt2(i,j) = 0;
                newGt1(i,j) = 0;
                newGt2(i,j) = 0;
                
            else
                
                newRt1(i,j) = mean(R1ij(R1ij~=0));
                newRt2(i,j) = mean(R2ij(R1ij~=0));
                newGt1(i,j) = mean(G1ij(R1ij~=0));
                newGt2(i,j) = mean(G2ij(R1ij~=0));
                
            end
            
            
            
        end              
    end
    
    
    
    [row, col] = find(newRt1~=0);   %%% find the pixels that is non zero
    NoPix = size(row,1);
    
    
    for i = 1:NoPix
        
        R1 = newRt1(row(i),col(i));
        R2 = newRt2(row(i),col(i));
        G1 = newGt1(row(i),col(i));
        G2 = newGt2(row(i),col(i));  
        
        yR = [R1, R2];
        yG = [G1, G2];
        
        
        Rf = fit(t.',yR.','exp1');
        Gf = fit(t.',yG.','exp1');
        
        
        
        coRt1 = (Rf.a*exp(Rf.b*tinj))./(Rf.a*exp(Rf.b*t1));
        coRt2 = (Rf.a*exp(Rf.b*tinj))./(Rf.a*exp(Rf.b*t2));
        
        coGt1 = (Gf.a*exp(Gf.b*tinj))./(Gf.a*exp(Gf.b*t1));
        coGt2 = (Gf.a*exp(Gf.b*tinj))./(Gf.a*exp(Gf.b*t2));
        
        
        
        startx = (row(i)-1)*rollRadius+1;
        starty = (col(i)-1)*rollRadius+1;
        
        CorrRt1(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1),ind) = coRt1.*currentRt1(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1));
        CorrRt2(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1),ind) = coRt2.*currentRt2(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1));

        CorrGt1(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1),ind) = coGt1.*currentGt1(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1));
        CorrGt2(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1),ind) = coGt2.*currentGt2(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1));
       
        
        coeffR(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1)) = Rf.b;
        coeffG(startx:(startx+rollRadius-1),starty:(starty+rollRadius-1)) = Gf.b;
        
        
    end
    
       
    
    
end


%% save the corrected images

SaveZSARFilePath=strcat(path_t1,'Rt1_Corr','/');   %%% save R @ t1

    if exist(SaveZSARFilePath, 'dir')
        cd(SaveZSARFilePath);
    else
        mkdir(SaveZSARFilePath);
        cd(SaveZSARFilePath);
    end
    
for k = 1:NumofStacks
    
        imwrite(uint16(CorrRt1(:,:,k)),strcat('CorrRT1_',num2str(k,'%02.0f'),'.tiff'),'Compression','none');
 
               
end 


SaveZSARFilePath=strcat(path_t1,'Rt2_Corr','/');   %%% save R @ t2

    if exist(SaveZSARFilePath, 'dir')
        cd(SaveZSARFilePath);
    else
        mkdir(SaveZSARFilePath);
        cd(SaveZSARFilePath);
    end
    
for k = 1:NumofStacks
    
        imwrite(uint16(CorrRt2(:,:,k)),strcat('CorrRT2_',num2str(k,'%02.0f'),'.tiff'),'Compression','none');
 
               
end 

SaveZSARFilePath=strcat(path_t1,'Gt1_Corr','/');   %%% save G @ t1

    if exist(SaveZSARFilePath, 'dir')
        cd(SaveZSARFilePath);
    else
        mkdir(SaveZSARFilePath);
        cd(SaveZSARFilePath);
    end
    
for k = 1:NumofStacks
    
        imwrite(uint16(CorrGt1(:,:,k)),strcat('CorrGT1_',num2str(k,'%02.0f'),'.tiff'),'Compression','none');
 
               
end 


SaveZSARFilePath=strcat(path_t1,'Gt2_Corr','/');   %%% save G @ t2

    if exist(SaveZSARFilePath, 'dir')
        cd(SaveZSARFilePath);
    else
        mkdir(SaveZSARFilePath);
        cd(SaveZSARFilePath);
    end
    
for k = 1:NumofStacks
    
        imwrite(uint16(CorrGt2(:,:,k)),strcat('CorrGT2_',num2str(k,'%02.0f'),'.tiff'),'Compression','none');
 
               
end




