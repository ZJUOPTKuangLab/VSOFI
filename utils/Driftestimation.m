function imagesout = Driftestimation(images,filePath,DriftestimationSaveflag)
imgRes = size(images,1);
imgNum = size(images,4);
detNum = size(images,3);
for i = 1:detNum
    tic;
    stacka = squeeze(images(:,:,i,:));
    stackb = stacka;
    t = size(stacka,3);
    for j = 2:t
        [stackb(:,:,j),Shift]=register(stacka(:,:,j),stacka(:,:,1));
    end
    disp(['  --detNum',num2str(i)]);
    images(:,:,i,:) = stackb;
    toc;
end
images(imgRes-5:end,imgRes-5:end,:,:) = 0;
images(1:5,1:5,:,:) = 0;
%% Drift estimation save
DriftimagesFolderPath = fullfile(filePath, 'DriftCorrectionimages');
if ~exist(DriftimagesFolderPath, 'dir')
    mkdir(DriftimagesFolderPath);
end
if DriftestimationSaveflag == 1
    DriftmatFolderPath = fullfile(filePath, 'DriftCorrectionmat');
    if ~exist(DriftmatFolderPath, 'dir')
        mkdir(DriftmatFolderPath);
    end
end

for i = 1:imgNum
    temp = images(:,:,:,i);
    WFtemp(:,:,i) = mat2gray(sum(temp, 3));
    tiffFilePath = fullfile(DriftimagesFolderPath, ['images-', num2str(i), '.tif']);
    imwrite(WFtemp(:,:,i), tiffFilePath);
end
imagesout = images;
if DriftestimationSaveflag == 1
    DriftmatFilePath = fullfile(DriftmatFolderPath, 'images.mat');
    disp('  --save drift correction mat...');
    tic;
    save(DriftmatFilePath,'images','-v7.3');
    toc;
end





