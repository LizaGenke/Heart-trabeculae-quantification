function varargout = project(varargin)
%-------------------------------------------------------------------------
%                                                              MAIA Master
%                                                          Medical sensors  
%                                      Quantification of trabeculae inside 
%                                the heart from MRI usingf ractal analysis
%                                                 Professor: Alain Lalande
% Authors:
% Daria Zotova
% Elizaveta Genke
% 
% January 2018
%-------------------------------------------------------------------------

% Begin initialization code 
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @project_OpeningFcn, ...
                   'gui_OutputFcn',  @project_OutputFcn, ...
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


% --- Executes just before project is made visible.
function project_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = project_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


% --- Executes on button press in openImage.
function openImage_Callback(hObject, eventdata, handles)
set(handles.drawEpicardialBorder, 'enable', 'on');
set(handles.clearEpicardiumBorder, 'enable', 'on');
set(handles.selectLeftVentricular, 'enable', 'off');
set(handles.setContractionBias, 'enable', 'off');
set(handles.selectPapillaryMuscles, 'enable', 'off');
set(handles.staticThresholdPM, 'enable', 'off');
set(handles.setThresholdForPM, 'enable', 'off');
set(handles.clearPapillaryMuscles, 'enable', 'off');
set(handles.calculateAreas, 'enable', 'off');
set(handles.calculateRatio, 'enable', 'off');
clearEpicardiumBorder_Callback(hObject, eventdata, handles);
handles = guidata(hObject);

axesHandlesToChildObjects = findobj(handles.axes1, 'Type', 'image');
if ~isempty(axesHandlesToChildObjects)
	delete(axesHandlesToChildObjects);
end
handles.currentImage = [];
guidata(hObject, handles);
axes(handles.axes1);
imshow(handles.currentImage,[]);

[f,p] = uigetfile('*.dcm','Select the image in DICOM format');
path = fullfile(p, f);
imageInfo = dicominfo(path);
I = dicomread(imageInfo);
I = im2double(I);
handles.currentImage = I;
handles.imageInfo = imageInfo;
guidata(hObject, handles);

 renderData(handles);


% --- Executes on button press in drawEpicardialBorder.
function drawEpicardialBorder_Callback(hObject, eventdata, handles)
clearEpicardiumBorder_Callback(hObject, eventdata, handles);

guidata(hObject,handles);
hFH = imfreehand();
BW = hFH.createMask();

I = handles.currentImage;
epicardialArea = I; 
epicardialArea(BW==0) = 0;

handles.epicardialArea = epicardialArea;
handles.epicardialAreaRoi = hFH;
guidata(hObject,handles);

renderData(handles);

set(handles.selectLeftVentricular, 'enable', 'on');
set(handles.staticThresholdLV, 'enable', 'on');
set(handles.setThresholdForLV, 'enable', 'on');
set(handles.setContractionBias, 'enable', 'on');



% --- Executes on button press in selectLeftVentricular.
function selectLeftVentricular_Callback(hObject, eventdata, handles)

epicardialArea = handles.epicardialArea;
ncAreaThreshold = handles.ncAreaThreshold;
valueContractionBias = handles.valueContractionBias;

trabeculaeArea = regiongrowing(epicardialArea,ncAreaThreshold);
trabeculaeArea = imfill(trabeculaeArea,'holes');
iterations = 200;
endocardialArea = activecontour(epicardialArea,trabeculaeArea,iterations,'Chan-Vese','SmoothFactor',3,'ContractionBias',handles.valueContractionBias);
level = graythresh(epicardialArea);
epicardialAreabw = im2bw(epicardialArea,level);
handles.epicardialAreabw = epicardialAreabw;
handles.trabeculaeArea = trabeculaeArea;
handles.endocardialArea = endocardialArea;

set(handles.selectPapillaryMuscles, 'enable', 'on');
set(handles.staticThresholdPM, 'enable', 'on');
set(handles.setThresholdForPM, 'enable', 'on');
set(handles.clearPapillaryMuscles, 'enable', 'on');
set(handles.calculateAreas, 'enable', 'on');
set(handles.calculateRatio, 'enable', 'on');

renderData(handles);

guidata(hObject,handles);


% --- Executes on button press in clearEpicardiumBorder.
function clearEpicardiumBorder_Callback(hObject, eventdata, handles)
if isfield(handles, 'epicardialArea')
    hold off;
    handles = rmfield(handles, 'epicardialArea');
    guidata(hObject, handles);
    axes(handles.axes1);
    imshow(handles.currentImage,[]);
end
if isfield(handles, 'endocardialArea')
    hold off;
    handles = rmfield(handles, 'endocardialArea');
    guidata(hObject, handles);
    axes(handles.axes1);
    imshow(handles.currentImage,[]);
end
if isfield(handles, 'trabeculaeArea')
    hold off;
    handles = rmfield(handles, 'trabeculaeArea');
    guidata(hObject, handles);
    axes(handles.axes1);
    imshow(handles.currentImage,[]);
end
if isfield(handles, 'papillaryMuscles')
    hold off;
    handles = rmfield(handles, 'papillaryMuscles');
    guidata(hObject, handles);
    axes(handles.axes1);
    imshow(handles.currentImage,[]);
end


% --- Executes on button press in resetImage.
function resetImage_Callback(hObject, eventdata, handles)

set(handles.textResults, 'visible','off');
set(handles.textNC, 'visible','off');
set(handles.textC, 'visible','off');
set(handles.textRatio, 'visible','off');
set(handles.showNC, 'visible','off');
set(handles.showC, 'visible','off');
set(handles.showRatio, 'visible','off');
set(handles.textCBSA, 'visible','off');
set(handles.textNCBSA, 'visible','off');
set(handles.showCBSA, 'visible','off');
set(handles.showNCBSA, 'visible','off');

set(handles.drawEpicardialBorder, 'enable', 'off');
set(handles.clearEpicardiumBorder, 'enable', 'off');
set(handles.selectLeftVentricular, 'enable', 'off');
set(handles.staticThresholdLV, 'enable', 'off');
set(handles.setThresholdForLV, 'enable', 'off');
set(handles.setContractionBias, 'enable', 'off');
set(handles.selectPapillaryMuscles, 'enable', 'off');
set(handles.staticThresholdPM, 'enable', 'off');
set(handles.setThresholdForPM, 'enable', 'off');
set(handles.clearPapillaryMuscles, 'enable', 'off');
set(handles.calculateAreas, 'enable', 'off');
set(handles.calculateRatio, 'enable', 'off');

set(handles.cropImage, 'enable', 'off');
set(handles.drawBorderF, 'enable', 'off');
set(handles.clearBorderF, 'enable', 'off');
set(handles.calculateFD, 'enable', 'off');
set(handles.calculateAllFD, 'enable', 'off');

clearEpicardiumBorder_Callback(hObject, eventdata, handles);
handles = guidata(hObject);

axesHandlesToChildObjects = findobj(handles.axes1, 'Type', 'image');
if ~isempty(axesHandlesToChildObjects)
	delete(axesHandlesToChildObjects);
end
handles.currentImage = [];
guidata(hObject, handles);
axes(handles.axes1);
imshow(handles.currentImage,[]);



function setThresholdForLV_Callback(hObject, eventdata, handles)

input = str2double(get(hObject,'String'));
if (isnan(input) || isempty(input) || input<0 || input>1)
    setThresholdForLV_CreateFcn(hObject, eventdata, handles)
else
    handles.ncAreaThreshold = input;
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function setThresholdForLV_CreateFcn(hObject, eventdata, handles)

handles.ncAreaThreshold = 0.001;
set(hObject,'String',handles.ncAreaThreshold)
set(hObject,'Value',handles.ncAreaThreshold)
guidata(hObject, handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function setContractionBias_Callback(hObject, eventdata, handles)

input = str2double(get(hObject,'String'));
if (isnan(input) || isempty(input) || input<-1 || input>1)
    setContractionBias_CreateFcn(hObject, eventdata, handles)
else
    handles.valueContractionBias = input;
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function setContractionBias_CreateFcn(hObject, eventdata, handles)

handles.valueContractionBias = -0.3;
set(hObject,'String',handles.valueContractionBias)
set(hObject,'Value',handles.valueContractionBias)
guidata(hObject, handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selectPapillaryMuscles.
function selectPapillaryMuscles_Callback(hObject, eventdata, handles)

I = handles.currentImage;
epicardialArea = handles.epicardialArea;
trabeculaeArea = handles.trabeculaeArea;
endocardialArea = handles.endocardialArea;
papillaryMusclesThreshold = handles.papillaryMusclesThreshold;

papillaryMuscles = regiongrowing(epicardialArea,handles.papillaryMusclesThreshold);
papillaryMuscles = imfill(papillaryMuscles,'holes');

if ~isfield(handles, 'papillaryMuscles')
    handles.papillaryMuscles = {};
end
handles.papillaryMuscles{end + 1} = papillaryMuscles;
guidata(hObject, handles);

renderData(handles);



function setThresholdForPM_Callback(hObject, eventdata, handles)

input = str2double(get(hObject,'String'));
if (isnan(input) || isempty(input) || input<0 || input>1)
    setThresholdForPM_CreateFcn(hObject, eventdata, handles)
else
    handles.papillaryMusclesThreshold = input;
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function setThresholdForPM_CreateFcn(hObject, eventdata, handles)

handles.papillaryMusclesThreshold = 0.001;
set(hObject,'String',handles.papillaryMusclesThreshold)
guidata(hObject, handles);
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearPapillaryMuscles.
function clearPapillaryMuscles_Callback(hObject, eventdata, handles)

I = handles.currentImage;
epicardialArea = handles.epicardialArea;
trabeculaeArea = handles.trabeculaeArea;
endocardialArea = handles.endocardialArea;

handles = rmfield(handles, 'papillaryMuscles');
guidata(hObject, handles);

renderData(handles);

function renderData(handles)
axes(handles.axes1);
I = handles.currentImage;
imshow(I,[]);
if ~I
    return
end
hold on;
if isfield(handles, 'epicardialArea')
    visboundaries(handles.epicardialArea,'Color','b');
end
if isfield(handles, 'trabeculaeArea')
    visboundaries(handles.trabeculaeArea,'Color','g');
end
if isfield(handles, 'endocardialArea')
    visboundaries(handles.endocardialArea,'Color','r');
end
if isfield(handles, 'papillaryMuscles')
    for k=1:length(handles.papillaryMuscles)
        data = handles.papillaryMuscles{1,k};
        visboundaries(data,'Color','y');
    end
end


% --- Executes on button press in calculateAreas.
function calculateAreas_Callback(hObject, eventdata, handles)

imageInfo = handles.imageInfo;
I = handles.currentImage;
trabeculaeArea = handles.trabeculaeArea;
endocardialArea = handles.endocardialArea;
epicardialAreabw = handles.epicardialAreabw;

if ~isfield(handles, 'papillaryMuscles')
    sumPapillaryMuscles = 0;
    papillaryMusclesAreaInMm = 0;
else
    papillaryMuscles = handles.papillaryMuscles;
    num_white_pixels_pmArea = nnz(sum([papillaryMuscles{:}]));
    spacing = imageInfo.PixelSpacing;
    per_pixel_area = spacing(1)*spacing(2);
    papillaryMusclesAreaInMm = num_white_pixels_pmArea*per_pixel_area; %area of PM1+PM2+PM#

    n = size(papillaryMuscles,2);
    matSize = size(papillaryMuscles{1},1);
    B = reshape(cell2mat(papillaryMuscles),matSize,[],n);
    sumPapillaryMuscles = sum(B,3);
end

spacing = imageInfo.PixelSpacing;
per_pixel_area = spacing(1)*spacing(2); %size of pixel in mm
num_white_pixels_epicardialAreabw = nnz(epicardialAreabw);
epicardialAreaInMm = num_white_pixels_epicardialAreabw*per_pixel_area; %area of Blue
disp(epicardialAreaInMm);

num_white_pixels_endocardialArea = nnz(endocardialArea);
endocardialAreaInMm = num_white_pixels_endocardialArea*per_pixel_area; %area of Red
disp(endocardialAreaInMm);

num_white_pixels_trabeculaeArea = nnz(trabeculaeArea);
trabeculaeAreaInMm = num_white_pixels_trabeculaeArea*per_pixel_area; %area of Green

nonC_white_pixels=nnz(endocardialArea&not(trabeculaeArea)&not(sumPapillaryMuscles));
nonCompArea = nonC_white_pixels*per_pixel_area;
compArea = epicardialAreaInMm - endocardialAreaInMm + papillaryMusclesAreaInMm;

if ~isfield(handles, 'arrayNonCompArea')
    handles.arrayNonCompArea = {};
end
handles.arrayNonCompArea{end + 1} = nonCompArea;
guidata(hObject, handles);

if ~isfield(handles, 'arrayCompArea')
    handles.arrayCompArea = {};
end
handles.arrayCompArea{end + 1} = compArea;
guidata(hObject, handles);

if ~isfield(handles, 'numberOfImages')
    handles.numberOfImages = 0;
end
handles.numberOfImages = handles.numberOfImages+1;
guidata(hObject, handles);
set(handles.setNumberOfImages, 'String', handles.numberOfImages);



% --- Executes on button press in calculateAreas.
function calculateRatio_Callback(hObject, eventdata, handles)

imageInfo = handles.imageInfo;
arrayCompArea = handles.arrayCompArea;
arrayNonCompArea = handles.arrayNonCompArea;
density = 0.00105;
if isfield(handles.imageInfo, 'SpacingBetweenSlices')
    thickness = imageInfo.SpacingBetweenSlices;
else thickness = 10;
end
ratio = sum([arrayNonCompArea{:}])*100/sum([arrayCompArea{:}]);
handles.ratio = ratio;
guidata(hObject, handles);

ratio = sprintf('%.2f', ratio);
nonCompArea = sprintf('%.2f', sum([arrayNonCompArea{:}])*density*thickness);
compArea = sprintf('%.2f', sum([arrayCompArea{:}])*density*thickness);

calculateBSA = 0.007184*((imageInfo.PatientSize*100)^0.725)*(imageInfo.PatientWeight^0.425);
BSA = sprintf('%.2f', calculateBSA);
calcnonCompAreaBSA = sum([arrayNonCompArea{:}])*density*thickness/calculateBSA;
nonCompAreaBSA = sprintf('%.2f', calcnonCompAreaBSA);
calccompAreaBSA = sum([arrayCompArea{:}])*density*thickness/calculateBSA;
compAreaBSA = sprintf('%.2f', calccompAreaBSA);

set(handles.textResults, 'visible','on');
set(handles.textNC, 'visible','on');
set(handles.textC, 'visible','on');
set(handles.textRatio, 'visible','on');
set(handles.textCBSA, 'visible','on');
set(handles.textNCBSA, 'visible','on');
set(handles.showNC, 'visible','on');
set(handles.showNC, 'String',nonCompArea);
set(handles.showC, 'visible','on');
set(handles.showNCBSA, 'visible','on');
set(handles.showCBSA, 'visible','on');

set(handles.showC, 'String',compArea);
set(handles.showRatio, 'visible','on');
set(handles.showRatio, 'String',ratio);
set(handles.showNCBSA, 'String',nonCompAreaBSA);
set(handles.showCBSA, 'String',compAreaBSA);




% --- Executes on button press in resetResults.
function resetResults_Callback(hObject, eventdata, handles)

resetImage_Callback(hObject, eventdata, handles);
set(handles.setNumberOfImages, 'String', '0');
handles.arrayCompArea = {};
handles.arrayNonCompArea = {};
handles.fractalDimension = [];
handles.numberOfImages = 0;
guidata(hObject, handles);


% --- Executes on button press in chooseSemiAuto.
function chooseSemiAuto_Callback(hObject, eventdata, handles)

resetResults_Callback(hObject, eventdata, handles);
set(handles.openImage, 'enable', 'on');
set(handles.resetImage, 'enable', 'on');
set(handles.openImageF, 'enable', 'off');
set(handles.resetImageF, 'enable', 'off');
set(handles.drawBorderF, 'enable', 'off');
set(handles.clearBorderF, 'enable', 'off');
set(handles.calculateFD, 'enable', 'off');
set(handles.calculateAllFD, 'enable', 'off');


% --- Executes on button press in chooseFractal.
function chooseFractal_Callback(hObject, eventdata, handles)

resetResults_Callback(hObject, eventdata, handles);
set(handles.textResults, 'visible','off');
set(handles.textNC, 'visible','off');
set(handles.textC, 'visible','off');
set(handles.textRatio, 'visible','off');
set(handles.showNC, 'visible','off');
set(handles.showC, 'visible','off');
set(handles.showRatio, 'visible','off');

set(handles.openImageF, 'enable', 'on');
set(handles.resetImageF, 'enable', 'on');
set(handles.drawBorderF, 'enable', 'off');
set(handles.clearBorderF, 'enable', 'off');
set(handles.calculateFD, 'enable', 'off');
set(handles.calculateAllFD, 'enable', 'off');

set(handles.openImage, 'enable', 'off');
set(handles.resetImage, 'enable', 'off');
set(handles.drawEpicardialBorder, 'enable', 'off');
set(handles.clearEpicardiumBorder, 'enable', 'off');
set(handles.selectLeftVentricular, 'enable', 'off');
set(handles.staticThresholdLV, 'enable', 'off');
set(handles.setThresholdForLV, 'enable', 'off');
set(handles.setContractionBias, 'enable', 'off');
set(handles.selectPapillaryMuscles, 'enable', 'off');
set(handles.staticThresholdPM, 'enable', 'off');
set(handles.setThresholdForPM, 'enable', 'off');
set(handles.clearPapillaryMuscles, 'enable', 'off');
set(handles.calculateAreas, 'enable', 'off');
set(handles.calculateRatio, 'enable', 'off');

% --- Executes on button press in openImageF.
function openImageF_Callback(hObject, eventdata, handles)

resetImage_Callback(hObject, eventdata, handles);
set(handles.cropImage, 'enable', 'on');
set(handles.drawBorderF, 'enable', 'off');
set(handles.clearBorderF, 'enable', 'off');

[f,p] = uigetfile('*.dcm','Select the image in DICOM format');
path = fullfile(p, f);
imageInfo = dicominfo(path);
I = dicomread(imageInfo);

handles.currentImage = I;
handles.imageInfo = imageInfo;
guidata(hObject, handles);
axes(handles.axes1);
imshow(handles.currentImage,[]);


% --- Executes on button press in cropImage.
function cropImage_Callback(hObject, eventdata, handles)

set(handles.drawBorderF, 'enable', 'on');
set(handles.clearBorderF, 'enable', 'on');
I = handles.currentImage;
croppedImage = imcrop(I,[]);
imshow(croppedImage,[]);
handles.croppedImage = croppedImage;
guidata(hObject, handles);



% --- Executes on button press in drawBorderF.
function drawBorderF_Callback(hObject, eventdata, handles)

hFH = imfreehand();
BW = hFH.createMask();

I = handles.croppedImage;
LV = I; 
LV(BW==0) = 0;

handles.LV = LV;
guidata(hObject,handles);
axes(handles.axes1);
imshow(handles.LV,[]);

set(handles.calculateFD, 'enable', 'on');
set(handles.calculateAllFD, 'enable', 'on');



% --- Executes on button press in clearBorderF.
function clearBorderF_Callback(hObject, eventdata, handles)

axes(handles.axes1);
imshow(handles.croppedImage,[]);


% --- Executes on button press in calculateFD.
function calculateFD_Callback(hObject, eventdata, handles)

LV = handles.LV;
LV8 = uint8(LV);
level = graythresh(LV8);
LVbw = im2bw(LV8,level);
LVb = bwmorph(LVbw,'remove');
axes(handles.axes1);
imshow(LVb,[]);
fractalDimension = boxcount(LVb);

if ~isfield(handles, 'fractalDimension')
    handles.fractalDimension = [];
end
handles.fractalDimension(end+1) = fractalDimension;

if ~isfield(handles, 'numberOfImages')
    handles.numberOfImages = 0;
end
handles.numberOfImages = handles.numberOfImages+1;
guidata(hObject, handles);
set(handles.setNumberOfImages, 'String', handles.numberOfImages);

guidata(hObject,handles);


% --- Executes on button press in calculateAllFD.
function calculateAllFD_Callback(hObject, eventdata, handles)

fractalDimension = handles.fractalDimension;
x = 1:length(fractalDimension);
figure
plot(x, fractalDimension);
xlabel('Base--------------->Mid--------------->Apex');
title('Fractal dimension');
