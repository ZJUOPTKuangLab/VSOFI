clear;clc;close all;
addpath('.\utils');
addpath('.\func');
addpath('.\deBKGD');

%% setting
filePath = '.\img\2025_2_15_0001_Confocal_Intensity_1';
% filePath = '.\img\2025_2_26\0007_Confocal_Intensity';
% filePath = '.\img\2025_2_15_0001_Confocal_Intensity_4';
% filePath = '.\img\2025_2_15_0003_Confocal_Intensity';
% cut frames
startNum = 1;
endNum = 20;
% deBKGD
interpNum = 2;
binning = 2;
dark = 1;
% PSF
flagGpu = 1;
lamdaEX = 640;
lamdaEM = 690;
pinsize = 80*1000/(100*300/70*60/50)/2;
% SOFI order
ACorder = 2;
% virtual order
virtual_order = [10,15,20];
% savePath
if interpNum == 1
    savePath = fullfile(filePath,[num2str(endNum),'result']);
else
    savePath = fullfile(filePath,[num2str(endNum),'result-interpNum',num2str(interpNum)]);
end
if ~exist(savePath, 'dir')
    mkdir(savePath);
end
%% load data
disp(['Load data from DriftCorrection mat...']);
matName = 'DriftCorrectionmat\images';
matFile = fullfile(filePath, [matName, '.mat']);
data = load(matFile);
cutflag  = 1;

images = data.images;
txtname = fullfile(filePath,'info.txt');
delimiterIn = ' ';
temp = importdata(txtname,delimiterIn);
info = temp.data;
infoName = temp.textdata;
exposeTime = info(3);
pixelsize = info(4);
outPower485 = info(5);
outPower640 = info(6);
clear data;

detNum = size(images,3);
imgRes = size(images,1);
imgNum = size(images,4);
%% cut image
if cutflag == 1
    if ~(startNum == 1 && endNum == imgNum)
        disp([num2str(endNum-startNum+1),'frames...']);
        images = images(:,:,:,startNum+1:endNum+1);
        imgNum = size(images,4);
    end
end
%% Fourier interpolation
disp(['Fourier interpolation...'])
for j = 1:size(images,3)
    temp = squeeze(images(:,:,j,:));
    temp = fourierInterpolation(temp,[interpNum,interpNum,1],'lateral');
    temp(temp < 0) = 0;
    images_interp(:,:,j,:) = temp;
    disp(['    -- Fourier interpolation ',num2str(j),' / ',num2str(19),'...']);
end
images = images_interp;
clear temp images_interp;
pixelsize = pixelsize/interpNum;
%% Photon Reassign
detNum = size(images,3);
imgRes = size(images,1);
imgNum = size(images,4);
images = single(images);
imgAiry = sum(images,4);
Confocal = sum(imgAiry,3);
[yoffset,xoffset] = calcOffset(imgAiry);
imgPR = photonReassign(imgAiry,round(yoffset),round(xoffset));
%% deBKGD-->mask
disp(['deBKGD...']);
imgISM = sum(imgPR,3);
[imgIF,mask] = supBKGD_mask(imgPR,binning,dark,savePath);
figure;imshow(mask,[]);title('mask');
mask19 = repmat(mask, [1, 1, 19]);
clear imgPR;
imgOUT = imgISM-imgIF;
figure;
subplot(2,2,1);imshow(Confocal,[]);colormap(hot),title('Confocal');
subplot(2,2,2);imshow(imgISM,[]);colormap(hot),title('ISM');
subplot(2,2,3);imshow(imgIF,[]);colormap(hot),title('IF');
subplot(2,2,4);imshow(imgOUT,[]);colormap(hot),title('OUT');
imwrite(mat2gray(imgIF),[savePath,'\imgIF.tif']);
imwrite(mat2gray(imgOUT),[savePath,'\imgOUT.tif']);
imwrite(mat2gray(Confocal),[savePath,'\Confocal.tif']);
imwrite(mat2gray(imgISM),[savePath,'\imgISM.tif']);
imwrite(mat2gray(mask),[savePath,'\mask.tif']);
%% Apply in-focus mask
images_nodeBKGD = images;
clear images;
disp(['Apply in-focus mask...']);
for i = 1:imgNum
    temp19 = photonReassign(squeeze(images_nodeBKGD(:,:,:,i)),round(yoffset),round(xoffset));
    imgPRs(:,:,:,i) = temp19.*mask19;
    imgAirys(:,:,:,i) = photonReassign(imgPRs(:,:,:,i),(-1)*round(yoffset),(-1)*round(xoffset));
    clear temp19;
    disp(['  Apply mask ' num2str(i) ' | ' num2str(imgNum)]);
end
% save([savePath,'/imgAirys.mat'],'imgAirys','-v7.3');
%% PSF
disp(['PSF calculation...']);
na = 1.45;
psize = pixelsize;
PSFex = calcPSF(lamdaEX,psize,na,flagGpu);
padSize = (imgRes-length(PSFex))/2;
PSFex = padarray(PSFex,[padSize,padSize],0,'both');
PSFex = PSFex/sum(PSFex(:));
OTFex = fftshift(fft2(PSFex));
PSFem = calcPSF(lamdaEM,psize,na,flagGpu);
PSFem = padarray(PSFem,[padSize,padSize],0,'both');
[xx,yy] = meshgrid(-(imgRes-1)/2:(imgRes-1)/2,-(imgRes-1)/2:(imgRes-1)/2);
freqUnit = sqrt(xx.^2+yy.^2)/imgRes/psize;
Fpinhole = besselj(1,freqUnit*2*pi*pinsize)./freqUnit;
OTFem = fftshift(fft2(PSFem)).*Fpinhole;
OTFem = OTFem/max(OTFem(:));
PSFISM = PSFex.*abs(ifft2(ifftshift(OTFem)));
PSFISM = PSFISM/sum(PSFISM(:));
OTFISM = fftshift(fft2(PSFISM));
%% Virtual modulation-demodulation
disp(['Virtual modulation-demodulation...']);
k0 = 0.4;
freqRound = 3;
kc = na*4*pi/lamdaEM;
ydet = zeros(imgRes,imgRes,detNum);
xdet = zeros(imgRes,imgRes,detNum);
for ii = 2:detNum
    ydet(:,:,ii) = 2*yoffset(ii-1);
    xdet(:,:,ii) = 2*xoffset(ii-1);
end
k = k0*kc*psize;
for i = 1:imgNum
    disp(['  NO. ' num2str(i) ' | ' num2str(imgNum)]);
    temp = squeeze(imgAirys(:,:,:,i));
    tempISM = sum(squeeze(imgPRs(:,:,:,i)),3);
    [vsd,patterneff] = virtual_modulation(temp,abs(OTFem),k,ydet,xdet,freqRound);
    for j = 1:numel(virtual_order)
        disp(['  --iter',num2str(virtual_order(j))]);
        Iobj(:,:,i,j) = virtual_recovery(vsd,abs(OTFex),abs(OTFISM),patterneff,tempISM,virtual_order(j));
    end
end
savePathtemp = [savePath,'\deBKGD-virtual'];
if ~exist(savePathtemp, 'dir')
    mkdir(savePathtemp);
end
for i = 1:imgNum
    imwrite(mat2gray(Iobj(:,:,i,3)),[savePathtemp,'\',num2str(i),'.tif']);
end
% save([savePath,'/Iobj.mat'],'Iobj','-v7.3');
%% VSOFI
disp(['VSOFI...']);
for j = 1:numel(virtual_order)
    disp(['  --iter',num2str(virtual_order(j))]);
    rawSOFI = squeeze(Iobj(:,:,:,j));
    nan_images = all(isnan(rawSOFI), [1, 2]);
    cleaned = rawSOFI(:,:,~nan_images);
    meanSOFI = mean(cleaned,3);
    subSOFI = cleaned - meanSOFI;
    subSOFI = abs(subSOFI);
    VSOFI(:,:,j) = cumulant(subSOFI,ACorder);
end
clear Confocal imgIF imgOUT temp;
for i = 1:numel(virtual_order)
    imwrite(mat2gray(VSOFI(:,:,i)),[savePath,'\VSOFI-iter',num2str(virtual_order(i)),'.tif']);
end




