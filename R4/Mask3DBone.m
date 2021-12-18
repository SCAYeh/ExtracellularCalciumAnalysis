
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% function modified from Juwell's algorithm for bone segmentation
%%%%% this function is used for determine osteocytes from 3D stacks
%%%%%%%%%%%%%%%%%%%%%%%%%

function AdjustedStack = Mask3DBone (ImgStack)

% Interactive signal recovery by Z plane for bone stack
% (Empirical; alternative to increase by exponential fit)
% Users are asked when to start and end signal recovery
% (Hint: Estimate these values using Enhanced Contrast in ImageJ)
% 20180416: *** WORKING VERSION ***

%% User Input

% === INPUT: ImgStack to be Filtered
% InputImg Tif vs MAT Query:
% TIF: img from microscope (ex. 1st pass)
% MAT: 1 channel single precision matrix (ex. 2nd pass; have done XY, now XZ)

Img_Height=500;
Img_Width = 1000;
NumImgSlices = size(ImgStack,3);
Z_PxLength = 1;

% Flat vs Non-Flat Query
% Non-flat slightly smoother outcome
FlatStrel_Query='n';
strel_Radius=3;     % Radius of flat element (disk)
offsetstrel_Radius=3;       % Radius of non-flat element (ball), >3 for effect 
offsetstrel_MaxOffset=0.05;  % Max offset of non-flat element (ball), larger offset->image brighter

%% 20180417: Sharpen along Z axis (5-pixel 1D Gaussian)

%fprintf('Sharpening along Z axis by Gaussian deconvolution...\n');

ImgStack_Input=ImgStack;

% YZ Reslice
ImgStack_Perm=double(permute(ImgStack_Input,[3 2 1]));  

% Sharpen along Z axis using Gaussian 1D Filter
Gauss1D_sigma=0.5;
Gauss1D_size=5;    % Length of Gauss 1D Filter vector
Gauss1D_colnum=linspace(-Gauss1D_size/2,Gauss1D_size/2,Gauss1D_size);
Gauss1D_filter=exp(-Gauss1D_colnum .^ 2/(2 * Gauss1D_sigma ^ 2));
Gauss1D_filter=Gauss1D_filter./sum(Gauss1D_filter); % Normalize

ZSharpenFilter=(padarray(Gauss1D_filter,(Gauss1D_size-1)/2,'both'))';

ImgStack_PermZSharp=zeros(size(ImgStack_Perm));
for k=1:size(ImgStack_Perm,3)   % parfor faster than for 
    [ImgStack_PermZSharp(:,:,k),~]=deconvblind(ImgStack_Perm(:,:,k),ZSharpenFilter,10);
end

ImgStack_ZSharp=double(permute(ImgStack_PermZSharp,[3 2 1]));   % Do not rescale intensity
%ImgStack_ZSharp = ImgStack_ZSharp./max(ImgStack_ZSharp(:));
ImgStack_ZSharp(ImgStack_ZSharp>=0.99)=0;
ImgStack_ZSharp = ImgStack_ZSharp./max(ImgStack_ZSharp(:));


%% intensity profile by z

ImgStack_Input=ImgStack_ZSharp;

ImgStack_IP=ImgStack_Input;
ImgStack_IP(ImgStack_IP>0.999)=0;   % Do not consider saturation (255 in 8-bit scale)
ImgStack_IP(ImgStack_IP<0.05)=0;    % Do not consider v. low signal

% Match Intensity Profile 
IntensityProfile=zeros(size(ImgStack_IP,3),1);

for k=1:size(ImgStack_IP,3)
    
    Img=reshape(ImgStack_IP(:,:,k),Img_Height*Img_Width,1);  % column vector
    Img(Img<eps)=[]; % Remove zeros
    [IPSort,IPSortIdx]=sort(Img,'descend');
    IPSortCtCutoff=round(Img_Height*Img_Width*0.00125);  % if #pixels available less than this number, skip
    
    % Exclude brightest and darkest intensities for mean calculation
    if ~isempty(IPSort)        
        IPSortCalc=horzcat(IPSort,IPSortIdx);
        % IPSortCalc(IPSortCalc(:,1)>0.9*IPSort(1,1),:)=NaN;
        % IPSortCalc(IPSortCalc(:,1)<1.1*IPSort(end,1),:)=NaN;
        IPSortCalc(isnan(IPSortCalc(:,1)),:)=[];
        
        if length(IPSortCalc(:,1))>IPSortCtCutoff 
            IntensityProfile(k,1)=mean(IPSortCalc(:,1),1);
        else
            IntensityProfile(k,1)=NaN; % Insufficient data pts to give meaningful value
        end
    else
        IntensityProfile(k,1)=NaN;
    end
end

% Patch occasional holes in intensity profile
% Create moving averages of 3 and 5 pixel span;
% ie. if NaN segment is < 5 elements long, fill by approximate


IPMovMean3=movmean(IntensityProfile,3,'omitnan');
IntensityProfile_MM3=IntensityProfile;
IntensityProfile_MM3(find(isnan(IntensityProfile)),1)=IPMovMean3(find(isnan(IntensityProfile)),1);

IPMovMean5=movmean(IntensityProfile,5,'omitnan');
IntensityProfile_MM5=IntensityProfile_MM3;
IntensityProfile_MM5(find(isnan(IntensityProfile_MM3)),1)=IPMovMean5(find(isnan(IntensityProfile_MM3)),1);



%% 20180412: Brightness-Contrast Adjustment by Z Intensity Profile

%fprintf('Brightness-Contrast Adjustment by Z Intensity Profile...\n');

ImgStack_Input=ImgStack;
%ImgStack_Input=ImgStack_ZSharp;

IntensityProfile_Input=IntensityProfile_MM5;

% == Finalize IntensityScaleFactor
IntensityScaleFactor=1./IntensityProfile_Input;
IntensityScaleFactor=IntensityScaleFactor/min(IntensityScaleFactor(:));

%  = Intensity Scale Factor 1: Fill NaN with constants, no extrapolation
% * Before 1st non-NAN  number: fill with 1st non-NaN number
% * After last non-NAN  number: fill with last non-NaN number
% * The rest fill with largest non-NAN number
IntensityScaleFactor1=IntensityScaleFactor;
IntensityScaleFactor1_FirstNum=IntensityScaleFactor1(find(~isnan(IntensityScaleFactor1),1,'first'),1);
IntensityScaleFactor1_LastNum=IntensityScaleFactor1(find(~isnan(IntensityScaleFactor1),1,'last'),1);
IntensityScaleFactor1(1:(find(~isnan(IntensityScaleFactor1),1,'first')-1),1)=...
    repmat(IntensityScaleFactor1_FirstNum,find(~isnan(IntensityScaleFactor1),1,'first')-1,1);
IntensityScaleFactor1((find(~isnan(IntensityScaleFactor1),1,'last')+1:end),1)=...
    repmat(IntensityScaleFactor1_LastNum,(size(IntensityScaleFactor1,1)-find(~isnan(IntensityScaleFactor1),1,'last')),1);
IntensityScaleFactor1(isnan(IntensityScaleFactor1),1)=...
    repmat(max(IntensityScaleFactor1,[],1),size(IntensityScaleFactor1(isnan(IntensityScaleFactor1),1),1),1);

% = Request user input for start and end Z planes
% For Intensity Scale Factor 1: planes out of user specified range set to 0
% For Intensity Scale Factor 2: planes out of user specified range, or factor>5, set to 0
% IntensityProfile_ZStart=input('\nPlease enter Start Z-plane based on intensity profile:');
% IntensityProfile_ZEnd=input('Please enter End Z-plane based on intensity profile:');
% fprintf('\n');

% IntensityScaleFactor1(1:IntensityProfile_ZStart,1)=0;
% IntensityScaleFactor1(IntensityProfile_ZEnd:end,1)=0;

ImgStack_IPBC=single(ImgStack_Input);
IntensityScaleFactorTemp=IntensityScaleFactor1;   

ImgStack_IPBC=ImgStack_IPBC.*...
        repmat(reshape(IntensityScaleFactorTemp,1,1,numel(IntensityScaleFactorTemp)),Img_Height,Img_Width,1);
ImgStack_IPBC=ImgStack_IPBC./max(ImgStack_IPBC(:));


%AdjustedStack=ImgStack_IPBC.*255;




%% 17Jun16: Normalize ImgStack w/ Maximum Projection along Z-axis 

% [ImgStack_ZMaxNorm]=ImgStackZMaxNorm_20170328(ImgStack,Z_PxLength,ImgMask_Query,ImgMaskStack);

%fprintf('Normalizing by maximum Z-projection...\n');

ImgStack_Input=ImgStack_IPBC;

ImgStack_ZMaxSort=ImgStack_Input;


ImgStack_ZMaxSort(ImgStack_ZMaxSort>=1)=0;   % Do not consider SatVal
ImgStack_ZMaxSort=sort(ImgStack_ZMaxSort,3,'descend');    % Matrix should not have NaN

 % Normalize by largest values along Z for each [Row,Col]
 % If [Row,Col] have < 5-10 um Z (~ 1-2 cell diameter) with non-zero value, 
 % remove from threshold consideration; otherwise take average of non-zero values 
 % and set as max
ImgStack_ZMaxSortZlim=round(5/Z_PxLength); 
ImgStack_ZMaxSortZlimZeroIdx=find(ImgStack_ZMaxSort(:,:,ImgStack_ZMaxSortZlim)==0); % List of [Row,Col] in linear index
for i=1:ImgStack_ZMaxSortZlim    % Small i count, do not use parfor
    ImgStack_ZMaxSort(ImgStack_ZMaxSortZlimZeroIdx+Img_Height*Img_Width*(i-1))=0;
end
ImgStack_ZMax=single(mean(ImgStack_ZMaxSort(:,:,1:ImgStack_ZMaxSortZlim),3));
ImgStack_ZMax(ImgStack_ZMax < 0.1)=1; % Ensure 0/0 does not happen; smallest val=eps; 
ImgStack_ZMax=repmat(ImgStack_ZMax,1,1,NumImgSlices);

% Divide to normalize
ImgStack_ZMaxNorm=single(ImgStack_Input./ImgStack_ZMax);

ImgStack_ZMaxNorm(ImgStack_ZMaxNorm>1)=1;

%AdjustedStack=ImgStack_ZMaxNorm.*255;


%% 17Jun16: Average Filter in XY plane

%fprintf('Running XY average filter...\n');

%ImgStack_Input=ImgStack_ZMaxNorm;
ImgStack_Input=ImgStack_IPBC;       %%% skip Zmax Norm procedure




% Average filter, n-pixel radius circle
if FlatStrel_Query=='y'
    h_av=fspecial('disk',strel_Radius);
else
    h_av=fspecial('disk',offsetstrel_Radius);
end

% Average filter on each Z slice.
% imfilter() is okay with NaN; if neighbourhood touches NaN pixels, 
% Center pixel = NaN 
ImgStack_XYAv=single(zeros(size(ImgStack_Input)));
for k=1:NumImgSlices   % for faster than parfor
    ImgStack_XYAv(:,:,k)=single(imfilter(ImgStack_Input(:,:,k),h_av));
end

%figure, imagesc(ImgStack_XYAv(:,:,22));


%% 17Jun16: Image closing -OR- Maximum Filter in XY plane
% Maximum filter includes more bone area ~2pixel width); Image closing
% looks more natural OPTION 1: Imge Closing

%fprintf('Performing XY image closing...\n');

ImgStack_Input=ImgStack_XYAv;

if FlatStrel_Query=='y'
    SE=strel('disk',strel_Radius,8);
else
    SE=offsetstrel('ball',offsetstrel_Radius,offsetstrel_MaxOffset,8);
end

% Image closing on each Z slice. Image closing is dilation, followed by
% erotion; fill in small holes MATLAB's imdilate is the same as max filter
% for grayscale images MATLAB's imerode is the same as min filter for
% grayscale images imdilate does not take NaN; set NaN to 0 (largely does
% not affect max)
ImgStack_XYClose=single(ImgStack_Input);
ImgStack_XYClose(isnan(ImgStack_XYClose))=0;

for k=1:NumImgSlices   % for faster than parfor
    ImgStack_XYClose(:,:,k)=imdilate(ImgStack_XYClose(:,:,k),SE);
    ImgStack_XYClose(:,:,k)=imerode(ImgStack_XYClose(:,:,k),SE);
    
end

%figure, imagesc(ImgStack_XYClose(:,:,22));

AdjustedStack=ImgStack_XYClose;

%% 17Jun16: Rescale ImgStack_XYNormSlideAvMax to fill Greyscale [0,1]

% % [ImgStack_HistRescale]=ImgStackHistRescale_20170328(ImgStack_XYClose,ImgMask_Query,ImgMaskStack,FlatStrel_Query,offsetstrel_MaxOffset);
% 
% %fprintf('Rescaling ImgStack Greyscale Intensity to [0,1]...\n');
% 
% % ImgStack_Input=ImgStack_ZMaxNorm;
% ImgStack_Input=ImgStack_XYClose;
% % ImgStack_Input=ImgStack_XYAv;
% 
% ImgStack_HistRescale_Prep=ImgStack_Input;
% 
% % Stack histogram to find scaling factor: 
% % New saturation intensity determined as follows:
% % 1) calculate moving av of ImgStack_XYNormSlideAvMax_Ct (smooth curve)
% % 2) Determine largest intensity bin where bin intensity > HistCtperBinEst,
% % which is the estimated bin count if counts are equally distributed among
% % histogram bins
% NumHistBins=(Img_Height*Img_Width*NumImgSlices)/1E3;   
% HistCtEdge=min(ImgStack_Input(:)):(1/NumHistBins):max(ImgStack_Input(:)); 
% HistCt_BinLoc=0.5*(HistCtEdge(1:end-1)+HistCtEdge(2:end));
% [ImgStack_HistCt,~] = histcounts(ImgStack_Input,HistCtEdge);
% 
% 
% HistCtperBinEst=(Img_Height*Img_Width*NumImgSlices)/NumHistBins;
% MovAv_kernel=ones(101,1)/4;  % 101-pixel moving average (50 pts on each side); pixel count should not include central pixel
% MovAv_kernel(51)=0;   % central pixel set to 0
% ImgStack_XYNormSlideAvMax_Ct_MovAv=conv(ImgStack_HistCt,MovAv_kernel,'same');
% 
% ImgStack_SatIndex=find(ImgStack_XYNormSlideAvMax_Ct_MovAv'>1*HistCtperBinEst,1,'last');
% ImgStack_Sat=HistCt_BinLoc(ImgStack_SatIndex);
% 
% % Plot Histogram counts of intensity range to be rescaled to full range
% % If flat structural element, smallest (bkgd) value after max filter = 0
% % If non-flat structural element, smallest (bkgd) value after max filter ~
% % offset_MaxOffset
% %figure;
% if FlatStrel_Query=='y'
%     %plot(ImgStack_HistCt(2:ImgStack_SatIndex));
%     ImgStack_Bkgd=0;
% else
%     % Locate bin of offset_MaxOffset (osMO);
%     % Then check for 1st index of positive HistCt slope, after moving average
%     [osMO_Ct,~]=histcounts(offsetstrel_MaxOffset,HistCtEdge);
%     ImgStack_osMO_BkgdIndex=find(osMO_Ct);
%     ImgStack_BkgdIndex_HistCtDiff=single(diff(ImgStack_HistCt(ImgStack_osMO_BkgdIndex:ImgStack_osMO_BkgdIndex+round(NumHistBins/10))))';
%     ImgStack_BkgdIndex_HistCtDiff_MovAv=single(conv(ImgStack_BkgdIndex_HistCtDiff,MovAv_kernel,'same'));
%     ImgStack_BkgdIndex_HistCtDiff_MovAv(ImgStack_BkgdIndex_HistCtDiff_MovAv<0)=0;
%     ImgStack_BkgdIndex_HistCtDiff_MovAv(ImgStack_BkgdIndex_HistCtDiff_MovAv>0)=1;
%     ImgStack_BkgdIndex=find(ImgStack_BkgdIndex_HistCtDiff_MovAv,1,'first')+ImgStack_osMO_BkgdIndex;    
%     ImgStack_Bkgd=HistCt_BinLoc(ImgStack_BkgdIndex);
%     %plot(ImgStack_HistCt(ImgStack_BkgdIndex:ImgStack_SatIndex));
% end
% 
% % Rescale
% ImgStack_HistRescale=single((ImgStack_HistRescale_Prep-ImgStack_Bkgd).*1/(ImgStack_Sat-ImgStack_Bkgd)); 
% ImgStack_HistRescale(ImgStack_HistRescale>1)=1;
% ImgStack_HistRescale(ImgStack_HistRescale<0)=0;
% 
% figure, imagesc(ImgStack_HistRescale(:,:,22));
% 
% AdjustedStack=ImgStack_HistRescale.*255;



end




