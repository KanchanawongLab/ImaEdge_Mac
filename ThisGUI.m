% Copyright (c) 2017.
% All rights reserved. Please read the 'license.txt' for license terms.
% 
% Developers: Zhen Zhang, Dr Fumio Motegi, Dr Pakorn Kanchanawong
% Contact: 
% Dr Pakorn Kanchanawong (biekp@nus.edu.sg)
% Dr Fumio Motegi (fmotegi@tll.org.sg)

function varargout = ThisGUI(varargin)
% THISGUI MATLAB code for ThisGUI.fig
%      THISGUI, by itself, creates a new THISGUI or raises the existing
%      singleton*.
%
%      H = THISGUI returns the handle to a new THISGUI or the handle to
%      the existing singleton*.
%
%      THISGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THISGUI.M with the given input arguments.
%
%      THISGUI('Property','Value',...) creates a new THISGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ThisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ThisGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ThisGUI

% Last Modified by GUIDE v2.5 23-Nov-2016 20:26:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ThisGUI_OpeningFcn, ...
    'gui_OutputFcn',  @ThisGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ThisGUI is made visible.
function ThisGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ThisGUI (see VARARGIN)

% Choose default command line output for ThisGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ThisGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ThisGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Nchannels.
function Nchannels_Callback(hObject, eventdata, handles)
% hObject    handle to Nchannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Nchannels contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Nchannels


% --- Executes during object creation, after setting all properties.
function Nchannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nchannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Channel4Ana.
function Channel4Ana_Callback(hObject, eventdata, handles)
% hObject    handle to Channel4Ana (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Channel4Ana contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Channel4Ana


% --- Executes during object creation, after setting all properties.
function Channel4Ana_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Channel4Ana (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function FactorOtsu_Callback(hObject, eventdata, handles)
% hObject    handle to FactorOtsu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FactorOtsu as text
%        str2double(get(hObject,'String')) returns contents of FactorOtsu as a double


% --- Executes during object creation, after setting all properties.
function FactorOtsu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FactorOtsu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ThickOut_Callback(hObject, eventdata, handles)
% hObject    handle to ThickOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThickOut as text
%        str2double(get(hObject,'String')) returns contents of ThickOut as a double


% --- Executes during object creation, after setting all properties.
function ThickOut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThickOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ThickIn_Callback(hObject, eventdata, handles)
% hObject    handle to ThickIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ThickIn as text
%        str2double(get(hObject,'String')) returns contents of ThickIn as a double


% --- Executes during object creation, after setting all properties.
function ThickIn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ThickIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nwindows_Callback(hObject, eventdata, handles)
% hObject    handle to Nwindows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nwindows as text
%        str2double(get(hObject,'String')) returns contents of Nwindows as a double


% --- Executes during object creation, after setting all properties.
function Nwindows_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nwindows (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LoadImg.
function LoadImg_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
close(figure(1));clc;
Nchan = get(handles.Nchannels,'Value');

for j = 1:Nchan
    ChanPath=imgetfile;
    InfoImage=imfinfo(ChanPath);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
    h = waitbar(0,['Reading Channel ',num2str(j),'/',num2str(Nchan),'.']);
    for i=1:NumberImages
        waitbar(i/NumberImages,h);
        OriAllChans(:,:,i,j) = imread(ChanPath,i);
    end
    close(h);
end
OriAllChans = mat2gray(OriAllChans);
mkdir data;
mkdir verification;
save data/OriAllChans.mat OriAllChans;


% --- Executes on button press in PrevRadial.
function PrevRadial_Callback(hObject, eventdata, handles)
% hObject    handle to PrevRadial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(figure(1));clc;
load data/ROI;
load data/AllMerged;
Nregions = get(handles.NsubRegions4Seg,'Value');

MaskRegion = zeros(size(AllMerged,1), size(AllMerged,2), Nregions);

AngleRange = [zeros(1,Nregions) + 360/Nregions*((1:Nregions)-ones(1,Nregions))  360];
[x y] = find(ROI==1);
mx = mean(x);
my = mean(y);
Angles = atan2((mx*ones(length(x),1)-x), (y-my*ones(length(x),1)))*180/pi;
Angles(find(Angles<0)) = Angles(find(Angles<0)) + 360;
h = waitbar(0,'generate radial masks ...');
for i = 1:Nregions
    waitbar(i/Nregions,h);
    k1 = find(Angles>=AngleRange(i));
    k2 = find(Angles(k1)<AngleRange(i+1));
    k = k1(k2);
    xt = x(k); yt = y(k);
    temp = zeros(size(AllMerged));
    temp(sub2ind(size(temp),xt,yt)) = 1;
    MaskRegion(:,:,i) = temp;
end
close(h);
save data/MaskRegion.mat MaskRegion;
figure(1);imshow(AllMerged);axis off;hold on;
for i = 1:Nregions
    B=bwboundaries(MaskRegion(:,:,i));
    B = B{1};
    figure(1);plot(B(:,2),B(:,1),'color','b');
end



% --- Executes on button press in showROI.
function showROI_Callback(hObject, eventdata, handles)
% hObject    handle to showROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(figure(1));clc;
Nchan = get(handles.Nchannels,'Value');
load data/OriAllChans;

for i = 1:Nchan
    Merged(:,:,((1:size(OriAllChans,3)) + (i-1)*size(OriAllChans,3))) = OriAllChans(:,:,:,i);
end
all = sum(Merged,3)/size(Merged,3);
all = imadjust(all);
temp = im2bw(all,graythresh(all));
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

margin = str2num(get(handles.Margin,'String'));

mask = im2single(temp);
mask(mask==0)=-1;
disF=bwdist(bwmorph(temp,'remove'));
disF = disF.*mask;

ROI = disF;
ROI(disF>=(-margin)) = 1;
ROI(disF<(-margin)) = 0;

B=bwboundaries(ROI);
B = B{1};
figure(1);imshow(all);hold on;axis off;plot(B(:,2),B(:,1));
save data/ROI.mat ROI;
AllMerged = all;
save data/AllMerged.mat AllMerged;






% --- Executes on button press in Seg.
function Seg_Callback(hObject, eventdata, handles)
% hObject    handle to Seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
close(figure(1));clc;
load data/OriAllChans;
load data/MaskRegion;
CompositeChanList = get(handles.Chans4Seg,'Value');
ThreshFactor = str2num(get(handles.FactorOtsu,'String'));
MedianFilter = str2num(get(handles.MedFilter,'String'));
CompositeImg = [];
eCompositeImg = [];
switch CompositeChanList
    case 1
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I1 = OriAllChans(:,:,i,1);
            I2 = zeros(size(I1));
            I3 = zeros(size(I1));
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
    case 2
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I2 = OriAllChans(:,:,i,2);
            I1 = zeros(size(I2));
            I3 = zeros(size(I2));
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
    case 3
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I3 = OriAllChans(:,:,i,3);
            I2 = zeros(size(I3));
            I1 = zeros(size(I3));
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
    case 4
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I1 = OriAllChans(:,:,i,1);
            I2 = OriAllChans(:,:,i,2);
            I3 = zeros(size(I1));
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
    case 5
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I1 = OriAllChans(:,:,i,1);
            I2 = zeros(size(I1));
            I3 = OriAllChans(:,:,i,3);
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
    case 6
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I2 = OriAllChans(:,:,i,2);
            I1 = zeros(size(I2));
            I3 = OriAllChans(:,:,i,3);
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
    case 7
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I1 = OriAllChans(:,:,i,1);
            I2 = OriAllChans(:,:,i,2);
            I3 = OriAllChans(:,:,i,3);
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
end

save data/CompositeImg.mat CompositeImg;
save data/eCompositeImg.mat eCompositeImg;

% start radial segmentation
AllMask = [];
mkdir verification/ContourVerify;
h = waitbar(0,'generating contours ...');
FinalContour = CompositeImg;
for z = 1:size(CompositeImg,3);
    waitbar(z/size(CompositeImg,3),h);
    Aver = eCompositeImg(:,:,z);
    Aver = uint8(Aver*255);
    AverWholeMask = zeros(size(I2));
    for i = 1:size(MaskRegion,3)
        [xr yr] = find(MaskRegion(:,:,i)==1);
        temp = Aver(sub2ind(size(Aver),xr,yr));
        temp1 = zeros(256,1);
        c = hist(double(temp),max(temp));
        temp1(1:length(c)) = c;
        level = ThreshFactor*SubregionOtsu(temp1);
        if level>=max(temp1(:))
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
        SmoothImg = medfilt2(temp,[MedianFilter MedianFilter]);
    end
    AllMask(:,:,z) = SmoothImg;
    B=bwboundaries(SmoothImg);
    B = B{1};
    temp2 = FinalContour(:,:,z);
    temp2(sub2ind(size(temp2),B(:,1),B(:,2))) = 255 ;
    FinalContour(:,:,z) = temp2;
    imwrite(temp2,['verification/ContourVerify/Contour',num2str(z),'.tif']);
end
close(h);
save data/AllMask.mat AllMask;
msgbox('Please verify the contours segmented in the folder verification -> ContourVerify');





% --- Executes on button press in Sample.
function Sample_Callback(hObject, eventdata, handles)
% hObject    handle to Sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
close(figure(1));clc;
load data/AllMask.mat;
load data/OriAllChans.mat;
Nsample = str2num(get(handles.Nwindows,'String'));
proj_depth = str2num(get(handles.ThickIn,'String'));
proj_out = str2num(get(handles.ThickOut,'String'));

h = waitbar(0,'Analyzing Each Frame...');
for i=1:size(AllMask,3)
    waitbar(i/size(AllMask,3),h,'Analyzing Each Frame...');
    display(['Currently Processing Frame ',num2str(i),'/',num2str(size(AllMask,3)),'.']);
    temp=AllMask(:,:,i);
    temp=im2bw(mat2gray(temp),0.5);
    mask = im2single(temp);
    mask(mask==0)=-1;
    disF=bwdist(bwmorph(temp,'remove'));
    disF = disF.*mask;
    inner = disF;
    inner(disF>=proj_depth) = 1;
    inner(disF<proj_depth) = 0;
    
    outer = disF;
    outer(disF>=(-proj_out)) = 1;
    outer(disF<(-proj_out)) = 0;
    
    B=bwboundaries(outer);
    B = B{1};
    B(1,3) = 0;
    TotalLength = 0;
    for j = 2:length(B)
        TotalLength = TotalLength + pdist([B(j-1,1:2);B(j,1:2)]);
        B(j,3) = TotalLength;
    end
    SampleInterval = TotalLength/Nsample;
    AllCalcuLengths = ((1:Nsample)-1)*SampleInterval;
    C = abs(bsxfun(@minus,B(:,3),AllCalcuLengths));
    [dummy,idx] = min(C(:,1:size(C,2)));
    SelectPts = B(idx,1:2);
    SelectPts(1:Nsample-1,3:4) = SelectPts(2:Nsample,1:2);
    SelectPts(Nsample,3:4) = SelectPts(1,1:2);
    
    a=bwboundaries(inner);
    A=a{1};
    inner_idx = dsearchn(A,SelectPts(:,1:2));
    inner_SelectPts = A(inner_idx,1:2);
    inner_SelectPts(1:Nsample-1,3:4) = inner_SelectPts(2:Nsample,1:2);
    inner_SelectPts(Nsample,3:4) = inner_SelectPts(1,1:2);
    AllCentroids = [(inner_SelectPts(:,1) + inner_SelectPts(:,3) + SelectPts(:,1) + SelectPts(:,3))/4, ...
        (inner_SelectPts(:,2) + inner_SelectPts(:,4) + SelectPts(:,2) + SelectPts(:,4))/4];
    
    if i~=1
        align_idx = dsearchn(AllCentroids,CentroidRecord(1,:,i-1));
        if align_idx~=1
            tempPts = inner_SelectPts(1:(align_idx-1),:);
            inner_SelectPts(1:(size(inner_SelectPts,1)-align_idx+1),:) = inner_SelectPts(align_idx:end,:);
            inner_SelectPts((size(inner_SelectPts,1)-align_idx+2):end,:) = tempPts;
            
            tempPts = SelectPts(1:align_idx-1,:);
            SelectPts(1:(size(SelectPts,1)-align_idx+1),:) = SelectPts(align_idx:end,:);
            SelectPts((size(SelectPts,1)-align_idx+2):end,:) = tempPts;
            
            tempPts = AllCentroids(1:align_idx-1,:);
            AllCentroids(1:(size(AllCentroids,1)-align_idx+1),:) = AllCentroids(align_idx:end,:);
            AllCentroids((size(AllCentroids,1)-align_idx+2):end,:) = tempPts;
        end
        
    end
    
    % show and save sampled image below
    mkdir verification/SampleVerify;
    for zz = 1:size(OriAllChans,4)
        close(figure(1));close(figure(2));close(figure(3));
        figure(zz);imshow(OriAllChans(:,:,i,zz));hold on;axis off;
        plot(B(:,2),B(:,1),'color','c');
        plot(A(:,2),A(:,1),'color','c');
        for zzz = 1:size(SelectPts,1)
            plot([SelectPts(zzz,2)  inner_SelectPts(zzz,2)],[SelectPts(zzz,1)  inner_SelectPts(zzz,1)],'color','c');
        end
        saveas(figure(zz),['verification/SampleVerify/Channel',num2str(zz),'-frame',num2str(i),'.tif']);
        figure(figure(zz));
    end
    close(figure(1));close(figure(2));close(figure(3));
    % show and save sampled image above
    CentroidRecord(:,:,i) = AllCentroids;
    for j = 1:Nsample
        dummy = zeros(size(temp));
        c = [SelectPts(j,2)  SelectPts(j,4)  inner_SelectPts(j,2)  inner_SelectPts(j,4)];
        r = [SelectPts(j,1)  SelectPts(j,3)  inner_SelectPts(j,1)  inner_SelectPts(j,3)];
        BW_window = roipoly(dummy,c,r);
        [x y] = find(BW_window==1);
        if ~isempty(x)
            for z = 1:size(OriAllChans,4)
                temp = OriAllChans(:,:,i,z);
                CortexIntensity(j,z,i) = mean(temp(sub2ind(size(temp),x,y)));
            end
        else
            for z = 1:size(OriAllChans,4)
                temp = OriAllChans(:,:,i,z);
                CortexIntensity(j,z,i) = mean(temp(sub2ind(size(temp),r,c)));
            end
        end
    end
end
close(h);
save data/CentroidRecord.mat CentroidRecord;
save data/CortexIntensity.mat CortexIntensity;



% --- Executes on button press in Analysis.
function Analysis_Callback(hObject, eventdata, handles)
% hObject    handle to Analysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
close(figure(1));clc;

load data/CortexIntensity.mat;
sigmaX = str2num(get(handles.FilterX,'String'));
sigmaY = str2num(get(handles.FilterY,'String'));
ChannelIdx = get(handles.Channel4Ana,'Value');
ChanofInterest = [];
for i = 1:size(CortexIntensity,3)
    ChanofInterest = [ChanofInterest  CortexIntensity(:,ChannelIdx,i)];
end
if sigmaX~=0 && sigmaY~=0
    H = fspecial('gaussian',[sigmaY  sigmaX],1);
    ChanofInterest = imfilter(ChanofInterest,H);
end
ChanofInterest = mat2gray(ChanofInterest);
figure(1);pcolor(ChanofInterest); shading flat;colorbar;colormap(jet);
xlabel('Frame Number');ylabel('Window Index');title(['Channel ',num2str(ChannelIdx),' Normalized Intenisty in Cell Cortex']);
mkdir result;
saveas(figure(1),['result/Channel',num2str(ChannelIdx),'.tif']);



function FilterX_Callback(hObject, eventdata, handles)
% hObject    handle to FilterX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FilterX as text
%        str2double(get(hObject,'String')) returns contents of FilterX as a double


% --- Executes during object creation, after setting all properties.
function FilterX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FilterX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FilterY_Callback(hObject, eventdata, handles)
% hObject    handle to FilterY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FilterY as text
%        str2double(get(hObject,'String')) returns contents of FilterY as a double


% --- Executes during object creation, after setting all properties.
function FilterY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FilterY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CurrWin.
function CurrWin_Callback(hObject, eventdata, handles)
% hObject    handle to CurrWin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
close(figure(1));clc;
load data/CentroidRecord.mat;
load data/CompositeImg.mat;

figure(1);imshow(CompositeImg(:,:,1));hold on;
for j = 1:5:size(CentroidRecord,1)
    text(CentroidRecord(j,2,1),CentroidRecord(j,1,1),['(',int2str(j),')'],'color','r');
end


function NewPos_Callback(hObject, eventdata, handles)
% hObject    handle to NewPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NewPos as text
%        str2double(get(hObject,'String')) returns contents of NewPos as a double


% --- Executes during object creation, after setting all properties.
function NewPos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NewPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Align.
function Align_Callback(hObject, eventdata, handles)
% hObject    handle to Align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
close(figure(1));clc;

load data/CentroidRecord.mat;
load data/CompositeImg.mat;
load data/CortexIntensity.mat;

newPosition = str2num(get(handles.NewPos,'String'));
Nsample = str2num(get(handles.Nwindows,'String'));

if newPosition~=1
    tempM = CortexIntensity(newPosition:end,:,:);
    CortexIntensity((Nsample-newPosition+2):end,:,:) = CortexIntensity(1:(newPosition-1),:,:);
    CortexIntensity(1:(Nsample-newPosition+1),:,:) = tempM;
    
    tempM = CentroidRecord(newPosition:end,:,:);
    CentroidRecord((Nsample-newPosition+2):end,:,:) = CentroidRecord(1:(newPosition-1),:,:);
    CentroidRecord(1:(Nsample-newPosition+1),:,:) = tempM;
end

save data/CortexIntensity.mat CortexIntensity;
save data/CentroidRecord.mat CentroidRecord;

mkdir verification/ImgWinsVerify;
for i = 1:size(CompositeImg,3)
    display(['Writing Image with Window Index: ',int2str(i),'/',int2str(size(CompositeImg,3))]);
    figure(1);
    imshow(CompositeImg(:,:,i));hold on;axis off;
    for j = 1:5:Nsample
        text(CentroidRecord(j,2,i),CentroidRecord(j,1,i),['(',int2str(j),')'],'color','r');
    end
    saveas(figure(1),['verification/ImgWinsVerify/ImgWinIdx',int2str(i),'.tif']);
    close(figure(1));
end
msgbox({'Window Alignment Done !',' ','You Can Find the Images with Window Index in Folder verification -> ImgWinsVerify.'});




% --- Executes on selection change in Chans4Seg.
function Chans4Seg_Callback(hObject, eventdata, handles)
% hObject    handle to Chans4Seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Chans4Seg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Chans4Seg


% --- Executes during object creation, after setting all properties.
function Chans4Seg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Chans4Seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Margin_Callback(hObject, eventdata, handles)
% hObject    handle to Margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Margin as text
%        str2double(get(hObject,'String')) returns contents of Margin as a double


% --- Executes during object creation, after setting all properties.
function Margin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in NsubRegions4Seg.
function NsubRegions4Seg_Callback(hObject, eventdata, handles)
% hObject    handle to NsubRegions4Seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns NsubRegions4Seg contents as cell array
%        contents{get(hObject,'Value')} returns selected item from NsubRegions4Seg


% --- Executes during object creation, after setting all properties.
function NsubRegions4Seg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NsubRegions4Seg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PreviewSample.
function PreviewSample_Callback(hObject, eventdata, handles)
% hObject    handle to PreviewSample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
close(figure(1));clc;
load data/AllMask.mat;
load data/OriAllChans.mat;

% parameters
Nsample = str2num(get(handles.Nwindows,'String'));
proj_depth = str2num(get(handles.ThickIn,'String'));
proj_out = str2num(get(handles.ThickOut,'String'));
% parameters
ImgShow = OriAllChans(:,:,1,1);


temp=AllMask(:,:,1);
temp=im2bw(mat2gray(temp),0.5);
mask = im2single(temp);
mask(mask==0)=-1;
disF=bwdist(bwmorph(temp,'remove'));
disF = disF.*mask;
inner = disF;
inner(disF>=proj_depth) = 1;
inner(disF<proj_depth) = 0;

outer = disF;
outer(disF>=(-proj_out)) = 1;
outer(disF<(-proj_out)) = 0;


B=bwboundaries(outer);
B = B{1};
B(1,3) = 0;
TotalLength = 0;
for j = 2:length(B)
    TotalLength = TotalLength + pdist([B(j-1,1:2);B(j,1:2)]);
    B(j,3) = TotalLength;
end
SampleInterval = TotalLength/Nsample;
AllCalcuLengths = ((1:Nsample)-1)*SampleInterval;
C = abs(bsxfun(@minus,B(:,3),AllCalcuLengths));
[dummy,idx] = min(C(:,1:size(C,2)));
SelectPts = B(idx,1:2);
SelectPts(1:Nsample-1,3:4) = SelectPts(2:Nsample,1:2);
SelectPts(Nsample,3:4) = SelectPts(1,1:2);

a=bwboundaries(inner);
A=a{1};
inner_idx = dsearchn(A,SelectPts(:,1:2));
inner_SelectPts = A(inner_idx,1:2);
inner_SelectPts(1:Nsample-1,3:4) = inner_SelectPts(2:Nsample,1:2);
inner_SelectPts(Nsample,3:4) = inner_SelectPts(1,1:2);
AllCentroids = [(inner_SelectPts(:,1) + inner_SelectPts(:,3) + SelectPts(:,1) + SelectPts(:,3))/4, ...
    (inner_SelectPts(:,2) + inner_SelectPts(:,4) + SelectPts(:,2) + SelectPts(:,4))/4];

figure(1);imshow(imadjust(ImgShow));hold on;axis off;
plot(B(:,2),B(:,1),'color','c');
plot(A(:,2),A(:,1),'color','c');
for i = 1:size(SelectPts,1)
    plot([SelectPts(i,2)  inner_SelectPts(i,2)],[SelectPts(i,1)  inner_SelectPts(i,1)],'color','c');
end


% --- Executes on selection change in MultiCoreList.
function MultiCoreList_Callback(hObject, eventdata, handles)
% hObject    handle to MultiCoreList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MultiCoreList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MultiCoreList


% --- Executes during object creation, after setting all properties.
function MultiCoreList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MultiCoreList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AutoSeg.
function AutoSeg_Callback(hObject, eventdata, handles)
% hObject    handle to AutoSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

warning off;
close(figure(1));clc;
load data/OriAllChans;
load data/MaskRegion;
CompositeChanList = get(handles.Chans4Seg,'Value');
MedianFilter = str2num(get(handles.MedFilter,'String'));
MaskRegion = MaskRegion;
p1 = str2num(get(handles.AutoSegStart,'String'));
p2 = str2num(get(handles.AutoSegStep,'String'));
p3 = str2num(get(handles.AutoSegEnd,'String'));
ThreshFactor = p1:p2:p3;
CompositeImg = [];
eCompositeImg = [];
switch CompositeChanList
    case 1
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I1 = OriAllChans(:,:,i,1);
            I2 = zeros(size(I1));
            I3 = zeros(size(I1));
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
    case 2
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I2 = OriAllChans(:,:,i,2);
            I1 = zeros(size(I2));
            I3 = zeros(size(I2));
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
    case 3
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I3 = OriAllChans(:,:,i,3);
            I2 = zeros(size(I3));
            I1 = zeros(size(I3));
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
    case 4
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I1 = OriAllChans(:,:,i,1);
            I2 = OriAllChans(:,:,i,2);
            I3 = zeros(size(I1));
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
    case 5
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I1 = OriAllChans(:,:,i,1);
            I2 = zeros(size(I1));
            I3 = OriAllChans(:,:,i,3);
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
    case 6
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I2 = OriAllChans(:,:,i,2);
            I1 = zeros(size(I2));
            I3 = OriAllChans(:,:,i,3);
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
    case 7
        h = waitbar(0,'Generating CompositeImg image ...');
        for i = 1:size(OriAllChans,3)
            waitbar(i/size(OriAllChans,3),h);
            I1 = OriAllChans(:,:,i,1);
            I2 = OriAllChans(:,:,i,2);
            I3 = OriAllChans(:,:,i,3);
            temp = cat(3,I1,I2,I3);
            CompositeImg(:,:,i) = mat2gray(rgb2gray(temp));
            eCompositeImg(:,:,i) = imadjust(mat2gray(rgb2gray(temp)));
        end
        close(h);
end

save data/CompositeImg.mat CompositeImg;
save data/eCompositeImg.mat eCompositeImg;
mkdir verification/ContourVerify;

MultiCore = get(handles.MultiCoreList,'Value');
if MultiCore==1
    h = waitbar(0,'Automated Segmentation in Progress ...');
    for i = 1:size(eCompositeImg,3)
        tic;
        waitbar(i/size(eCompositeImg,3),h);
        AllMask(:,:,i) = AutoSegCont(eCompositeImg(:,:,i),MedianFilter,MaskRegion,ThreshFactor);
        display(['Currently Processing: Frame ',num2str(i),'/',num2str(size(eCompositeImg,3)),'. Time:',num2str(toc),'s']);
    end
    close(h);
else
    tic;
    delete(gcp);
    if MultiCore==2
        parpool ('local',round(feature('numCores')/2));
    else
        parpool ('local',feature('numCores'));
    end
    % parellel computing is used and it may take a few minutes for large data set
    parfor i = 1:size(eCompositeImg,3)
        AllMask(:,:,i) = AutoSegCont(eCompositeImg(:,:,i),MedianFilter,MaskRegion,ThreshFactor);
    end
    delete(gcp);
    toc;
end


FinalContour = CompositeImg;
h = waitbar(0,'Writing Contours ...');
for i = 1:size(AllMask,3)
    waitbar(i/size(AllMask,3),h);
    B=bwboundaries(AllMask(:,:,i));
    B = B{1};
    temp2 = FinalContour(:,:,i);
    temp2(sub2ind(size(temp2),round(B(:,1)),round(B(:,2)))) = 255;
    
    imwrite(temp2,['verification/ContourVerify/Contour',num2str(i),'.tif']);
end
close(h);
save data/AllMask.mat AllMask;
msgbox('Please verify the contours segmented in the folder verification -> ContourVerify');



function MedFilter_Callback(hObject, eventdata, handles)
% hObject    handle to MedFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MedFilter as text
%        str2double(get(hObject,'String')) returns contents of MedFilter as a double


% --- Executes during object creation, after setting all properties.
function MedFilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MedFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AutoSample.
function AutoSample_Callback(hObject, eventdata, handles)
% hObject    handle to AutoSample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
close(figure(1));clc;
load data/AllMask.mat;
load data/OriAllChans.mat;
load data/eCompositeImg;
Nsample = str2num(get(handles.Nwindows,'String'));
BoundMin = str2num(get(handles.AutoSamMinBound,'String'));
BoundMax = str2num(get(handles.AutoSamMaxBound,'String'));

MultiCore = get(handles.MultiCoreAutoSample,'Value');
if MultiCore==1
    h = waitbar(0,'Automated Sampling in Progress ...');
    for i = 1:size(AllMask,3)
        tic;
        waitbar(i/size(AllMask,3),h);
        inners(:,:,i) = AutoSample(eCompositeImg(:,:,i),AllMask(:,:,i),BoundMin,BoundMax);
        display(['Currently Processing: Frame ',num2str(i),'/',num2str(size(eCompositeImg,3)),'. Time:',num2str(toc),'s']);
    end
    close(h);
else
    tic;
    delete(gcp);
    if MultiCore==2
        parpool ('local',round(feature('numCores')/2));
    else
        parpool ('local',feature('numCores'));
    end
    % parellel computing is used and it may take a few minutes for large data set
    parfor i = 1:size(AllMask,3)
        inners(:,:,i) = AutoSample(eCompositeImg(:,:,i),AllMask(:,:,i),BoundMin,BoundMax);
    end
    delete(gcp);
    toc;
end

outers = AllMask;
h = waitbar(0,'Analyzing Each Frame...');
for i=1:size(AllMask,3)
    waitbar(i/size(AllMask,3),h,'Analyzing Each Frame...');
    display(['Currently Processing Frame ',num2str(i),'/',num2str(size(AllMask,3)),'.']);
    
    outer = outers(:,:,i);
    inner = inners(:,:,i);
    
    B=bwboundaries(outer);
    B = B{1};
    B(1,3) = 0;
    TotalLength = 0;
    for j = 2:length(B)
        TotalLength = TotalLength + pdist([B(j-1,1:2);B(j,1:2)]);
        B(j,3) = TotalLength;
    end
    SampleInterval = TotalLength/Nsample;
    AllCalcuLengths = ((1:Nsample)-1)*SampleInterval;
    C = abs(bsxfun(@minus,B(:,3),AllCalcuLengths));
    [dummy,idx] = min(C(:,1:size(C,2)));
    SelectPts = B(idx,1:2);
    SelectPts(1:Nsample-1,3:4) = SelectPts(2:Nsample,1:2);
    SelectPts(Nsample,3:4) = SelectPts(1,1:2);
    
    a=bwboundaries(inner);
    A=a{1};
    inner_idx = dsearchn(A,SelectPts(:,1:2));
    inner_SelectPts = A(inner_idx,1:2);
    inner_SelectPts(1:Nsample-1,3:4) = inner_SelectPts(2:Nsample,1:2);
    inner_SelectPts(Nsample,3:4) = inner_SelectPts(1,1:2);
    AllCentroids = [(inner_SelectPts(:,1) + inner_SelectPts(:,3) + SelectPts(:,1) + SelectPts(:,3))/4, ...
        (inner_SelectPts(:,2) + inner_SelectPts(:,4) + SelectPts(:,2) + SelectPts(:,4))/4];
    
    if i~=1
        align_idx = dsearchn(AllCentroids,CentroidRecord(1,:,i-1));
        if align_idx~=1
            tempPts = inner_SelectPts(1:(align_idx-1),:);
            inner_SelectPts(1:(size(inner_SelectPts,1)-align_idx+1),:) = inner_SelectPts(align_idx:end,:);
            inner_SelectPts((size(inner_SelectPts,1)-align_idx+2):end,:) = tempPts;
            
            tempPts = SelectPts(1:align_idx-1,:);
            SelectPts(1:(size(SelectPts,1)-align_idx+1),:) = SelectPts(align_idx:end,:);
            SelectPts((size(SelectPts,1)-align_idx+2):end,:) = tempPts;
            
            tempPts = AllCentroids(1:align_idx-1,:);
            AllCentroids(1:(size(AllCentroids,1)-align_idx+1),:) = AllCentroids(align_idx:end,:);
            AllCentroids((size(AllCentroids,1)-align_idx+2):end,:) = tempPts;
        end
        
    end
    
    % show and save sampled image below
    mkdir verification/SampleVerify;
    for zz = 1:size(OriAllChans,4)
        close(figure(1));close(figure(2));close(figure(3));
        figure(zz);imshow(OriAllChans(:,:,i,zz));hold on;axis off;
        plot(B(:,2),B(:,1),'color','c');
        plot(A(:,2),A(:,1),'color','c');
        for zzz = 1:size(SelectPts,1)
            plot([SelectPts(zzz,2)  inner_SelectPts(zzz,2)],[SelectPts(zzz,1)  inner_SelectPts(zzz,1)],'color','c');
        end
        saveas(figure(zz),['verification/SampleVerify/Channel',num2str(zz),'-frame',num2str(i),'.tif']);
        close(figure(zz));
    end
    close(figure(1));close(figure(2));close(figure(3));
    % show and save sampled image above
    CentroidRecord(:,:,i) = AllCentroids;
    for j = 1:Nsample
        dummy = zeros(size(AllMask,1), size(AllMask,2));
        c = [SelectPts(j,2)  SelectPts(j,4)  inner_SelectPts(j,2)  inner_SelectPts(j,4)];
        r = [SelectPts(j,1)  SelectPts(j,3)  inner_SelectPts(j,1)  inner_SelectPts(j,3)];
        BW_window = roipoly(dummy,c,r);
        [x y] = find(BW_window==1);
        if ~isempty(x)
            for z = 1:size(OriAllChans,4)
                temp = OriAllChans(:,:,i,z);
                CortexIntensity(j,z,i) = mean(temp(sub2ind(size(temp),x,y)));
            end
        else
            for z = 1:size(OriAllChans,4)
                temp = OriAllChans(:,:,i,z);
                CortexIntensity(j,z,i) = mean(temp(sub2ind(size(temp),r,c)));
            end
        end
    end
end
close(h);
save data/CentroidRecord.mat CentroidRecord;
save data/CortexIntensity.mat CortexIntensity;


% --- Executes on selection change in MultiCoreAutoSample.
function MultiCoreAutoSample_Callback(hObject, eventdata, handles)
% hObject    handle to MultiCoreAutoSample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MultiCoreAutoSample contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MultiCoreAutoSample


% --- Executes during object creation, after setting all properties.
function MultiCoreAutoSample_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MultiCoreAutoSample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function AutoSegStart_Callback(hObject, eventdata, handles)
% hObject    handle to AutoSegStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AutoSegStart as text
%        str2double(get(hObject,'String')) returns contents of AutoSegStart as a double


% --- Executes during object creation, after setting all properties.
function AutoSegStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AutoSegStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AutoSegStep_Callback(hObject, eventdata, handles)
% hObject    handle to AutoSegStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AutoSegStep as text
%        str2double(get(hObject,'String')) returns contents of AutoSegStep as a double


% --- Executes during object creation, after setting all properties.
function AutoSegStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AutoSegStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AutoSegEnd_Callback(hObject, eventdata, handles)
% hObject    handle to AutoSegEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AutoSegEnd as text
%        str2double(get(hObject,'String')) returns contents of AutoSegEnd as a double


% --- Executes during object creation, after setting all properties.
function AutoSegEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AutoSegEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AutoSamRef.
function AutoSamRef_Callback(hObject, eventdata, handles)
% hObject    handle to AutoSamRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off;
load data/eCompositeImg;
load data/AllMask.mat;
close(figure(1));
h = waitbar(0,'Checking Gradients for All Frames...');
colorlist = jet(size(eCompositeImg,3));
for i = 1:size(eCompositeImg,3)
    waitbar(i/size(eCompositeImg,3),h,'Checking Gradients for All Frames...');
    [Gmag,Gdir] = imgradient(eCompositeImg(:,:,i),'sobel');
    Mask = AllMask(:,:,i);
    B=bwboundaries(Mask);
    B = B{1};
    iniPos = str2num(get(handles.AutoSamMinBound,'String')); % unit: pixels
    endPos = str2num(get(handles.AutoSamMaxBound,'String')); % unit: pixels
    CurrStep = 1;
    temp=Mask;
    temp=im2bw(mat2gray(temp),0.5);
    mask = im2single(temp);
    mask(mask==0)=-1;
    
    disF=bwdist(bwmorph(temp,'remove'));
    disF = disF.*mask;
    inner = disF;
    inner(disF>=iniPos) = 1;
    inner(disF<iniPos) = 0;
    
    A=bwboundaries(inner);
    A = A{1};
    AllDepth = iniPos;
    AllGrad = [];
    AllDepth = [];
    for ibound = iniPos:CurrStep:endPos
        AllGrad = [AllGrad  sum(Gmag(sub2ind(size(Gmag),A(:,1),A(:,2))))/length(A)];
        AllDepth = [AllDepth  ibound];
        inner = disF;
        inner(disF>=ibound) = 1;
        inner(disF<ibound) = 0;
        A=bwboundaries(inner);
        A = A{1};
    end
    iniPos(end) = [];
    
    figure(1);plot(AllDepth, AllGrad, 'color',colorlist(i,:));hold on;
    xlabel('Distance to Cell Edge (pixel)');ylabel('Image Gradient');title('Trend in Image Gradient');
end
close(h);


function AutoSamMaxBound_Callback(hObject, eventdata, handles)
% hObject    handle to AutoSamMaxBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AutoSamMaxBound as text
%        str2double(get(hObject,'String')) returns contents of AutoSamMaxBound as a double


% --- Executes during object creation, after setting all properties.
function AutoSamMaxBound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AutoSamMaxBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AutoSamMinBound_Callback(hObject, eventdata, handles)
% hObject    handle to AutoSamMinBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AutoSamMinBound as text
%        str2double(get(hObject,'String')) returns contents of AutoSamMinBound as a double


% --- Executes during object creation, after setting all properties.
function AutoSamMinBound_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AutoSamMinBound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
