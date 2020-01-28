%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%this script builds on demo for Cytoplasmic Actin Intensity extraction
%changes are as following: 
%intensities are per cell count cell count from max projections of dapi 
%nuclear mask uses same structuring element as cytosolicMask
%volumes are counted as sum of pixels in the area per cell count
%includes region prop analysis of blobs per slice
%previous cortical varibales were updated to cytosolic
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%##################################################################
% use following conversion if have different resolution images
% %pixel conversion 1024x1024: .124um/px
% pxsize= 8.06; %px/um 
% bestradius = 4.96; %40px÷9px/um = um 
% radiconvertpx = bestradius*pxsize; 
%################################################################## 



clear
%Image Inputs
addpath('/Users/Mojde/Documents/Cilia/confocal/tomography scripts/matlab pics/std/Taz/E_F1R_4/');
file = 'E_F1R_4';
pref.test = file;
pref.ext = '.tif';


%stack specifics
endslice = 112;
startslice = 1;
numslice = floor((endslice-startslice+1)/2);

%variable preallocation
nucAIS= zeros(1,numslice);  cytosolicAIS=zeros(1,numslice);
cellBodyAIS=zeros(1,numslice); dapiBstack=[];
nucvol=0; cellvol=0;


for i=startslice:numslice
    %8bit grayscale background20pxsubtract autocontrast 2pxmedian filter
    slided = num2str(10000+i);%DAPI
    slidea = num2str(10000+i+ numslice); %AF488
    actin = imread([pref.test,slidea(2:end),pref.ext]);
    dapi= imread([pref.test,slided(2:end),pref.ext]);
    OGactin= imread([pref.test,'OG',slidea(2:end),pref.ext]); %Raw Intensity
   
    figure(1); subplot(2,2,1), imshow(dapi); title('fixed dapi');
    subplot(2,2,2), imshow(OGactin); title('OGactin');
    
    %============(1) Nuclear Mask ============
    dapiRaw = dapi;
    dapi = imbinarize(dapi); %create logical mask
    dapi = imfill(dapi,'holes'); %fill
    dapi = bwmorph(dapi,'shrink'); %erode
    radius = 40; %pxs rolling ball
    structuringElement = strel('ball', radius,0);
	dapi = imclose(uint8(dapi), structuringElement);
    dapi = uint8(bwareafilt(imbinarize(dapi), [50,1500],8)); %area filter
    nucactin = OGactin.*uint8(dapi); %convolute mask & actin
%     figure(5); imshow(nucactin); title('dapi appd');
    subplot(2,2,3), imshow(nucactin); title('dapi appd');
    dapiBstack= cat(3,dapi,dapiBstack);
    
    
	
	%============(2) Cell Body Mask ============
	actinbinary=imbinarize(actin); 
    smallestArea=240; %nuclear area
	cytosolicMask=uint8(bwareafilt(actinbinary,[smallestArea,2500])); %area filter
	% Smooth the border by morphological closing 
    radius = 40; %pxs rolling ball
    structuringElement = strel('ball', radius,0);
	cytosolicMask = imclose(cytosolicMask, structuringElement);
	% Fill holes since likely part of cell body
	cytosolicMask = uint8(imfill(cytosolicMask, 'holes'));
    cytosolicMask=uint8(bwareafilt(imbinarize(cytosolicMask),[smallestArea,2500])); %area filter

    
    %convolusions w actin channel
    C = OGactin.*uint8(cytosolicMask);
%     figure(6); imshow(C); title('cytosolic Mask applied');
    Co = C.*uint8(~dapi);
%     figure(7);imshow(Co); title('cytosolicMask applied rmvd nuclei');
    subplot(2,2,4),imshow(Co); title('cytosolicMask applied rmvd nuclei');
    %pause(1);
    %°°°°°° exert saving patch here from demo.m if want to save imgs °°°°°    
    
    

    %============(3) Actin Intensities Sum ============
    cytosolicAIS(i) = sum(Co(:)); %add slice AI to stack
    cellBodyAIS(i) = sum(C(:)); %total intensity of cell body
    nucAIS(i)= sum(nucactin(:)); %add slice nuclear AI to stack
    nucvol = nucvol +nnz(nucactin); %nuclear pixel num
    cellvol= cellvol +nnz(C); %cellbody pixel num
    
    
    %============(4) Mean Cell AI blob ============
    slicestats1=regionprops('table',imbinarize(dapi),OGactin,'Area','MeanIntensity');
    nucMAI = sum(slicestats1.MeanIntensity);
    
    slicestats2=regionprops('table',cytosolicMask,OGactin,'Area','MeanIntensity');
    cellBodyMAI = sum(slicestats2.MeanIntensity);
    
    slicestats3=regionprops('table',~imbinarize(dapi),C,'Area','MeanIntensity');
    cytosolicMAI = sum(slicestats3.MeanIntensity);

end

%Max Intensity Projections
mipd = max(dapiBstack,[],3);
mipd= bwpropfilt(imbinarize(mipd),'Area',[300,2500]);
dd = bwconncomp(mipd); %find blobs/objs
n=dd.NumObjects; %cell count




cortex = sum(cytosolicAIS)/n;
cellbody = sum(cellBodyAIS)/n;
nuclear = sum(nucAIS)/n;
%strange sightings: nuc+cort > cell
nucvolume= nucvol/n; %per cell
cellvolume= cellvol/n;

mcortex = sum(cytosolicMAI)/n;
mcellbody = sum(cellBodyMAI)/n;
mnuclear = sum(nucMAI)/n;

results = table(cellbody,cortex,nuclear,nucvolume,cellvolume,n)
figure(); imshow(mipd);