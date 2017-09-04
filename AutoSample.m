% Copyright (c) 2017.
% All rights reserved. Please read the 'license.txt' for license terms.
% 
% Developers: Zhen Zhang, Dr Fumio Motegi, Dr Pakorn Kanchanawong
% Contact: 
% Dr Pakorn Kanchanawong (biekp@nus.edu.sg)
% Dr Fumio Motegi (fmotegi@tll.org.sg)

function  inner = AutoSample(eCompositeImg,Mask,BoundMin,BoundMax)

[Gmag,Gdir] = imgradient(eCompositeImg,'sobel');

B=bwboundaries(Mask);
B = B{1};

temp=Mask;
temp=im2bw(mat2gray(temp),0.5);
mask = im2single(temp);
mask(mask==0)=-1;
disF=bwdist(bwmorph(temp,'remove'));
disF = disF.*mask;
inner = disF;
inner(disF>=BoundMin) = 1;
inner(disF<BoundMin) = 0;

A=bwboundaries(inner);
A = A{1};
AllDepth = [];
AllGrad = [];
for i = BoundMin:BoundMax
    AllGrad = [AllGrad  sum(Gmag(sub2ind(size(Gmag),A(:,1),A(:,2))))/length(A)];
    AllDepth = [AllDepth  i];
    inner = disF;
    inner(disF>=i) = 1;
    inner(disF<i) = 0;
    A=bwboundaries(inner);
    A = A{1};
end

PkIndx = find(AllGrad==max(AllGrad));
PkIndx = PkIndx(1);
OptDepth = AllDepth(PkIndx);

inner = disF;
inner(disF>=OptDepth) = 1;
inner(disF<OptDepth) = 0;
% figure;plot([BoundMin:BoundMax],AllGrad);


