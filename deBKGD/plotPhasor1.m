function plotPhasor1(sB,gB,a,b,imgDark)
yshift = 0.2;
sB = sB+yshift;
b = b+yshift;

mapResolution = 500;
PhasorMap = hist3([sB(:),gB(:)],'Edges', {1/mapResolution:1/mapResolution:1,1/mapResolution:1/mapResolution:1});
cmp = [zeros(3,33) jet(256)']';
cmp(1:33,1) = 1:-1/32:0;
cmp(1:33,2) = 1:-1/32:0;
cmp(1:33,3) = 1:-0.5/32:0.5;
figure;imagesc(PhasorMap(1:0.6*mapResolution,1:0.6*mapResolution)); colormap(cmp);
lx = (0.6-b)/a:0.001:-b/a;
hold on;plot(lx*mapResolution,(a*lx+b)*mapResolution,'--k','LineWidth',2);

hold on;quiver(-0.03*mapResolution,0.2*mapResolution,0.700*mapResolution,0,'color','k','LineWidth',2);
hold on;quiver(0,-0.03*mapResolution,0,0.7*mapResolution,'color','k','LineWidth',2);
text(-0.02*mapResolution,(yshift-0.04)*mapResolution,'0', 'FontSize', 16,'HorizontalAlignment','right');

xylabel = (-0.01:0.001:0.01)*mapResolution;
for ii= 0.2:0.2:0.4
    hold on;plot(ones(1,length(xylabel))*ii*mapResolution,xylabel+yshift*mapResolution,'-k','LineWidth',2);
    text(ii*mapResolution,(yshift-0.04)*mapResolution,num2str(ii), 'FontSize', 16,'HorizontalAlignment','center');
end
ii = 0.6;
text((ii-0.005)*mapResolution,(yshift-0.035)*mapResolution,'g', 'FontSize', 25,'HorizontalAlignment','right');

ii = 0;
hold on;plot(xylabel,ones(1,length(xylabel))*ii*mapResolution,'-k','LineWidth',2);
text((-0.02)*mapResolution,(ii+0.005)*mapResolution,num2str(ii-yshift), 'FontSize', 16,'HorizontalAlignment','right');
for ii = 0.4:0.2:0.4
    hold on;plot(xylabel,ones(1,length(xylabel))*ii*mapResolution,'-k','LineWidth',2);
    text((-0.02)*mapResolution,(ii+0.005)*mapResolution,num2str(ii-yshift), 'FontSize', 16,'HorizontalAlignment','right');
end
ii = 0.6;
text((-0.025)*mapResolution,(ii-0.015)*mapResolution,'s', 'FontSize', 25,'HorizontalAlignment','right');
axis xy; axis image;axis off;

sB = sB-yshift;
b = b-yshift;

x1 = 0;
y1 = a*x1+b;
xmatr = (gB+a*sB-a*b)/(1+a^2);
ymatr = (a*gB+a^2*sB+b)/(1+a^2);
dpie = sqrt((x1-xmatr).^2+(y1-ymatr).^2);

angle=atan(a);

gm = xmatr(dpie == min(dpie(imgDark)))-sin(-angle)*0.0;
sm = ymatr(dpie == min(dpie(imgDark)))-cos(-angle)*0.0;

phasorLine = dpie(imgDark);
phasorLine = phasorLine(:)-min(phasorLine(:));
phasorLine(phasorLine == 0) = [];
phasorLine(phasorLine == 1) = [];
Nbin = 120;
[phasorH,bins] = hist(phasorLine,Nbin-1);
options = fitoptions('gauss2', 'Lower', zeros(1,6));
phasorFit = fit(bins',phasorH','gauss2',options);
bins = [0,bins]+gm;
phasorH = [0,-phasorH/max(phasorH(:))*0.225]+sm;

if phasorFit.b1 < phasorFit.b2
    ifp = phasorFit.b1+gm;
    ifs = phasorFit.c1/sqrt(2);
    ofp = phasorFit.b2+gm;
    ofs = phasorFit.c2/sqrt(2);
else
    ifp = phasorFit.b2+gm;
    ifs = phasorFit.c2/sqrt(2);
    ofp = phasorFit.b1+gm;
    ofs = phasorFit.c1/sqrt(2);
end

    lowOut = ifp;

if ifp+2*ifs < ofp
    highOut = ofp;
else
    highOut = ofp+ofs;
end

c1 = [0.6,0.05,0.1];
c2 = [0.05,0.4,0.6];
color = hsv(Nbin);
lowPos = find(bins<lowOut,1,'last');
highPos = find(bins>highOut,1);
startPos = find(phasorH<sm-4e-3,1);
stopPos = find(phasorH<sm-4e-3,1,'last');

Rx0=gm;
Ry0=sm;
for i=1:length(phasorH)
    x0=(bins(i)-Rx0)*cos(angle)-(phasorH(i)-Ry0)*sin(angle)+Rx0;
    y0=(bins(i)-Rx0)*sin(angle)+(phasorH(i)-Ry0)*cos(angle)+Ry0;
    bins(i)=x0;
    phasorH(i)=y0+yshift;
end

hold on;
plot(bins(startPos:lowPos)*mapResolution, phasorH(startPos:lowPos)*mapResolution,'LineWidth',3.5,'Color',c1);
hold on;
plot(bins(highPos:stopPos)*mapResolution, phasorH(highPos:stopPos)*mapResolution,'LineWidth',3.5,'Color',c2);

a =repmat((1:highPos-lowPos)'/(highPos-lowPos+1),1,3);
c = repmat(c1,highPos-lowPos,1).*(1-a)+repmat(c2,highPos-lowPos,1).*a;

for ii = lowPos:highPos-1
  hold on
  plot(bins(ii:ii+1)*mapResolution, phasorH(ii:ii+1)*mapResolution,'LineWidth',3.5,'Color',c(ii-lowPos+1,:));
end

end