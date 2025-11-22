function imgPR = photonReassign(imgAiry,yoffset,xoffset)
imageRes = size(imgAiry,1);
det_num = size(imgAiry,3);

ystart = -yoffset.*(yoffset < 0) + 1;
yend = imageRes - yoffset.*(yoffset >= 0);
xstart = -xoffset.*(xoffset < 0) + 1;
xend = imageRes - xoffset.*(xoffset >= 0);
imgPR = zeros(imageRes,imageRes,det_num);
imgPR(:,:,1) = imgAiry(:,:,1);
for k=1:det_num-1
    y_index = ystart(k):yend(k);
    x_index = xstart(k):xend(k);
    imgPR(y_index,x_index,k+1) = imgAiry(y_index+yoffset(k),x_index+xoffset(k),k+1);
end
end