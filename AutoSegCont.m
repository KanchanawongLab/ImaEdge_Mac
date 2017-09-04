% Copyright (c) 2017.
% All rights reserved. Please read the 'license.txt' for license terms.
% 
% Developers: Zhen Zhang, Dr Fumio Motegi, Dr Pakorn Kanchanawong
% Contact: 
% Dr Pakorn Kanchanawong (biekp@nus.edu.sg)
% Dr Fumio Motegi (fmotegi@tll.org.sg)

function CurrBW = AutoSegCont(eComposite,MedianFilter,MaskRegion,ThreshFactor)
[Gmag,Gdir] = imgradient(eComposite,'sobel');
ThreshFactor = ThreshFactor/100;
AllAverGrad = [];
for zzz = 1:length(ThreshFactor);
    Aver = eComposite;
    Aver = uint8(Aver*255);
    AverWholeMask = zeros((size(eComposite,1)),(size(eComposite,2)));
    for i = 1:size(MaskRegion,3)
        [xr yr] = find(MaskRegion(:,:,i)==1);
        temp = Aver(sub2ind(size(Aver),xr,yr));
        temp1 = zeros(256,1);
        c = hist(double(temp),max(temp));
        temp1(1:length(c)) = c;
        level = ThreshFactor(zzz)*SubregionOtsu(temp1);
        k = find(temp>=level);
        AverWholeMask(sub2ind(size(Aver),xr(k),yr(k))) = 1;
    end
    temp = bwfill(AverWholeMask,'hole');
    [temp, number] = bwlabel(temp,8);
    unqA=unique(temp(:));
    pix_count=histc(temp(:),unqA);
    unqA = unqA(2:end);pix_count=pix_count(2:end);
    pix_max = find(pix_count==max(pix_count));
    temp(temp==find(pix_count==max(pix_count)))=pix_max(1);
    pix_count(pix_count==max(pix_count))=0;
    temp(temp~=pix_max)=0;
    temp=im2bw(temp,0.5);
    temp = bwmorph(temp,'spur');
    if MedianFilter~=0
        temp = medfilt2(temp,[MedianFilter MedianFilter]);
    end
    temp = bwfill(AverWholeMask,'hole');
    [temp, number] = bwlabel(temp,8);
    unqA=unique(temp(:));
    pix_count=histc(temp(:),unqA);
    unqA = unqA(2:end);pix_count=pix_count(2:end);
    pix_max = find(pix_count==max(pix_count));
    temp(temp==find(pix_count==max(pix_count)))=pix_max(1);
    pix_count(pix_count==max(pix_count))=0;
    temp(temp~=pix_max)=0;
    temp=im2bw(temp,0.5);
    temp = bwmorph(temp,'spur');
    SmoothImg = temp;
    B=bwboundaries(SmoothImg);
    B = B{1};
    
    AllAverGrad = [AllAverGrad  sum(Gmag(sub2ind(size(Gmag),B(:,1),B(:,2))))/size(B,1)];
end


PkIndx = find(AllAverGrad==max(AllAverGrad));
PkIndx = PkIndx(1);
OptThreshFactor = ThreshFactor(PkIndx);

Aver = eComposite;
Aver = uint8(Aver*255);
AverWholeMask = zeros((size(eComposite,1)),(size(eComposite,2)));
for i = 1:size(MaskRegion,3)
    [xr yr] = find(MaskRegion(:,:,i)==1);
    temp = Aver(sub2ind(size(Aver),xr,yr));
    temp1 = zeros(256,1);
    c = hist(double(temp),max(temp));
    temp1(1:length(c)) = c;
    level = OptThreshFactor*SubregionOtsu(temp1);
    if level>=max(temp(:))
        msgbox('The threshold is too high that segmentation errors may occur !');
    end
    k = find(temp>=level);
    AverWholeMask(sub2ind(size(Aver),xr(k),yr(k))) = 1;
end
temp = bwfill(AverWholeMask,'hole');
[temp, number] = bwlabel(temp,8);
unqA=unique(temp(:));
pix_count=histc(temp(:),unqA);
unqA = unqA(2:end);pix_count=pix_count(2:end);
pix_max = find(pix_count==max(pix_count));
temp(temp==find(pix_count==max(pix_count)))=pix_max(1);
pix_count(pix_count==max(pix_count))=0;
temp(temp~=pix_max)=0;
temp=im2bw(temp,0.5);
temp = bwmorph(temp,'spur');
if MedianFilter~=0
    temp = medfilt2(temp,[MedianFilter MedianFilter]);
end
temp = bwfill(temp,'hole');
[temp, number] = bwlabel(temp,8);
unqA=unique(temp(:));
pix_count=histc(temp(:),unqA);
unqA = unqA(2:end);pix_count=pix_count(2:end);
pix_max = find(pix_count==max(pix_count));
temp(temp==find(pix_count==max(pix_count)))=pix_max(1);
pix_count(pix_count==max(pix_count))=0;
temp(temp~=pix_max)=0;
temp=im2bw(temp,0.5);
temp = bwmorph(temp,'spur');
CurrBW = temp;


