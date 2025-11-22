 function Iobj = virtual_recovery(raw_data,OTFex,OTFISM,patterneff,widefield,max_iter)

raw_data = gpuArray(single(raw_data));
OTFex = gpuArray(single(OTFex));
patterneff = gpuArray(single(patterneff));
widefield = gpuArray(single(widefield));

num=size(raw_data,3);
widefield = widefield/max(widefield(:));
Iobj = widefield;
imgRes = size(raw_data,1);

n=1;

while n <= max_iter
%     tic;
    for i = 1 : num
        Iop = Iobj .* patterneff(:,:,i);
        fop = fftshift(fft2(Iop));
        fop1 = fop+conj(OTFex).* (fftshift(fft2(raw_data(:,:,i)))-OTFex.*fop);
        Iop1 = ifft2(ifftshift(fop1));
        Iobj1 = Iobj+patterneff(:,:,i).*(Iop1-Iop)./(max(max(patterneff(:,:,i)))^2);
        Iobj1(Iobj1<0) = 0;      
        Iobj = abs(Iobj1);
    end
%     ttime = toc ;
%     disp(['  iter ' num2str(n) ' | ' num2str(max_iter) ', took ' num2str(ttime) ' secs']);
    n=n+1;
end

Iobj = gather(double(Iobj));

end