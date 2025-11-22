function [imgIF,mask] = supBKGD_mask(imgPR,binning,dark,savePath)

imgISM = sum(imgPR,3);
imgRes = size(imgPR,1);
rdets = zeros(imgRes,imgRes,19);
rdets(:,:,2:7) = 0.25;%(2*r0/rairy)^2;
rdets(:,:,8:2:18) = 0.75;%(sqrt(3)*2*r0/rairy)^2;
rdets(:,:,9:2:19) = 1;%;(2*2*r0/rairy)^2;

if binning > 0
    kernel = fspecial('disk',binning);
    for ii = 1:19
        imgPR1(:,:,ii) = imfilter(imgPR(:,:,ii),kernel);
    end
else
    imgPR1 = imgPR;
end
imgDark = mean(imgPR1,3)>dark;
figure;imshow(imgDark);title('imgDark');

gB = sum(imgPR1.*cos(2*pi*rdets),3)./sum(imgPR1,3);
sB = sum(imgPR1.*sin(2*pi*rdets),3)./sum(imgPR1,3);

imgDark1 = max(imgPR1,[],3)>dark;
Coeffs = fit(sB(imgDark1),gB(imgDark1),'poly1');
a = 1/Coeffs.p1;
b = -Coeffs.p2/Coeffs.p1;
plotPhasor1(sB,gB,a,b,imgDark);

frame = getframe(gcf);
img = frame2im(frame);
imwrite(img, [savePath,'\phasor.tif']);

empmatr = ones(imgRes,imgRes);
x1 = 0;
y1 = a*x1+b;
xmatr = (gB+a*sB-a*b*empmatr)/(1+a^2);
ymatr = (a*gB+a^2*sB+b*empmatr)/(1+a^2);
dpie = sqrt((x1*empmatr-xmatr).^2+(y1*empmatr-ymatr).^2);
dpie = dpie-min(dpie(imgDark));
dpie(~imgDark) = max(dpie(imgDark));

dpie = mat2gray(dpie);

[lowOut,highOut] = phasor2hist1_save2(dpie(imgDark),savePath);

mask = 1-imadjust(dpie,[lowOut,highOut],[0,1]);
if binning > 0
    mask = filterGauss(mask,imgRes/3/binning);
end

imgIF = imgISM.*mask;

end