function Iout = filterGauss(Iin,sigma)
% Frequency domain Gaussian filter

imgRes = size(Iin,1);
[xx,yy]=meshgrid(-imgRes/2:imgRes/2-1,-imgRes/2:imgRes/2-1);
rr = sqrt(xx.^2+yy.^2);
fGauss = exp(-rr.^2/2/sigma^2);
Iout = ifft2(ifftshift(fftshift(fft2(Iin)).*fGauss));
Iout = abs(Iout);
end