function [yoffset,xoffset] = calcOffset(imgAiry)
detNum = size(imgAiry,3);
xoffset = zeros(1,detNum-1);
yoffset = zeros(1,detNum-1);

FI1 = fft2(imgAiry(:,:,1));
FRI1I1 = FI1.*conj(FI1);
RI1 = ifft2(FRI1I1);
RI1 = fftshift(RI1);
numI1 = find(RI1==max(RI1(:)));
[iI1,jI1] = ind2sub(size(RI1), numI1);

for k=2:detNum
    FIn = fft2(imgAiry(:,:,k));
    FRI1In = FI1.*conj(FIn);
    R = ifft2(FRI1In);
    R = fftshift(R);
    num = find(R==max(R(:)));
    [i,j] = ind2sub(size(R), num);
    if length(i) == 1
        yoffset(k-1) = iI1-i;
        xoffset(k-1) = jI1-j;
    end
end
end