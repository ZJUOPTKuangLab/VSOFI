function [lowOut,highOut] = phasor2hist1_save2(img_phasor,savePath)

phasorLine = img_phasor(:);
phasorLine(phasorLine == 0) = [];
phasorLine(phasorLine == 1) = [];
Nbin = 200;
[phasorH,bins] = hist(phasorLine,Nbin);
options = fitoptions('gauss2', 'Lower', zeros(1,6),'upper', [+Inf,1,1,+Inf,1,1]);
phasorFit = fit(bins',phasorH','gauss2',options);

figure;
hold on;
histogram(phasorLine, Nbin,'FaceColor', 'yellow');
xlabel('bins');
ylabel('counts');
legend('phasorLine Histogram','box','off');
box on;
hold off;

frame = getframe(gcf);
img = frame2im(frame);
imwrite(img, [savePath,'\hist.tif']);

if phasorFit.b1 < phasorFit.b2
    ifp = phasorFit.b1;
    ifs = phasorFit.c1/sqrt(2);
    ofp = phasorFit.b2;
    ofs = phasorFit.c2/sqrt(2);
else
    ifp = phasorFit.b2;
    ifs = phasorFit.c2/sqrt(2);
    ofp = phasorFit.b1;
    ofs = phasorFit.c1/sqrt(2);
end

lowOut = bins(100);
highOut = bins(180);

end