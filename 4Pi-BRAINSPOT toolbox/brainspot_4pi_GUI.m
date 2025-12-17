function varargout = brainspot_4pi_GUI(varargin)
% brainspot_4pi_GUI MATLAB code for brainspot_4pi_GUI.fig
%
% (C) Copyright                     The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu
%

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @brainspot_4pi_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @brainspot_4pi_GUI_OutputFcn, ...
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





% --- Executes just before brainspot_4pi_GUI is made visible.
function brainspot_4pi_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to brainspot_4pi_GUI (see VARARGIN)

% Choose default command line output for brainspot_4pi_GUI

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


global INSPR4pi; 

INSPR4pi = [];

% setup parameters
INSPR4pi.setup.workspace = '';
INSPR4pi.setup.Pixelsize = 0.129; %um
INSPR4pi.setup.RefractiveIndex = 1.512;
INSPR4pi.setup.nMed = 1.352;
INSPR4pi.setup.Lambda = 0.69;   %um
INSPR4pi.setup.NA = 1.32;
INSPR4pi.setup.offset = 100.0;
INSPR4pi.setup.gain = 2.0;

INSPR4pi.setup.is_sCMOS = 0;    %0 for EMCCD camera; 1 for sCMOS camera
INSPR4pi.setup.sCMOS_input = []; 
INSPR4pi.setup.is_imgsz = 1;    % 0 for image size smaller than 100 x 100 pixels
INSPR4pi.setup.is_bg = 0;       %default is no background subtraction


% pupil parameters
INSPR4pi.pupil.init_z_bot = zeros(1,21);
INSPR4pi.pupil.init_z_bot(1) = 1.2;
INSPR4pi.pupil.init_z_top = zeros(1,21);
INSPR4pi.pupil.init_z_top(1) = -1.2;
INSPR4pi.pupil.ZernikeorderN = 7;  
INSPR4pi.pupil.Zshift_mode = 1; 
INSPR4pi.pupil.Zrange_low = -0.8;
INSPR4pi.pupil.Zrange_mid = 0.05;
INSPR4pi.pupil.Zrange_high = 0.8;
INSPR4pi.pupil.Z_pos = [-0.8:0.05:0.8];
INSPR4pi.pupil.bin_lowerBound = 30; %the lower bound of each group images
INSPR4pi.pupil.iter = 4;
INSPR4pi.pupil.min_similarity = 0.6; %min similarity in NCC calculation

INSPR4pi.pupil.iter_mode = 0;  

INSPR4pi.pupil.Phasediff = -1.58;   % phase difference between s- and p-polarizations
INSPR4pi.pupil.Iratio = 0.92;       % transmission ratio between bottom and top emission path
INSPR4pi.pupil.Phi0 = 0;            % cavity phase
INSPR4pi.pupil.ModulationDepth = 0.8; 
INSPR4pi.pupil.Zoffset = 0;
INSPR4pi.pupil.min_photon = 1000;
INSPR4pi.pupil.blur_sigma = 2;

% cavity phase and objective misalignment parameters
INSPR4pi.cavity.isNewdata = 0;   %cavity parameters
INSPR4pi.cavity.isNewpupil = 0;
INSPR4pi.cavity.isSeg = 1;
INSPR4pi.cavity.isObj = 0;

INSPR4pi.cavity.seg_thresh_low = 25; %segmentation threshold
INSPR4pi.cavity.seg_thresh_high = 40;

% reconstruction parameters
INSPR4pi.recon.isNewdata = 0;   %Reconstruction parameters
INSPR4pi.recon.isNewpupil = 0;
INSPR4pi.recon.isNewphi0 = 0;

INSPR4pi.recon.isSeg = 1;
INSPR4pi.recon.isObj = 0;   %default is no objective misalignment
INSPR4pi.recon.isRej = 1;
INSPR4pi.recon.isDC = 1;
INSPR4pi.recon.isGPU = 1;   %default is GPU version
INSPR4pi.recon.is_bg = 0;   %default is no background subtraction

INSPR4pi.recon.seg_thresh_low = 25; %segmentation threshold
INSPR4pi.recon.seg_thresh_high = 40;

INSPR4pi.recon.rej.min_photon = 1500;   %rejection parameters
INSPR4pi.recon.rej.llthreshold = 2400;
INSPR4pi.recon.rej.loc_uncerxy_max = 15;  % nm
INSPR4pi.recon.rej.loc_uncerz_max = 10; 
INSPR4pi.recon.rej.zmask_low = -0.5;    % um
INSPR4pi.recon.rej.zmask_high = 0.5;

INSPR4pi.recon.dc.frmpfile = 1999; 
INSPR4pi.recon.dc.z_offset = 2000;
INSPR4pi.recon.dc.step_ini = 400;

% display parameters
INSPR4pi.display.Recon_color_highb = 1;
INSPR4pi.display.imagesz = 144;
INSPR4pi.display.zm = 10;

%default parameters, load default parameters


%% display usage information in the mouseover event  

% setup module
handles.setup_workspace_pushbutton.TooltipString = 'set default input and output paths';
handles.setup_import_pushbutton.TooltipString = 'import general setting parameters';
handles.setup_set_pushbutton.TooltipString = 'modify setting parameters by users';
handles.setup_export_pushbutton.TooltipString = 'export setting parameters';


% data module
handles.data_bg_checkbox.TooltipString = '4Pi-brainspot supports the background subtraction option in cases with high background. During background subtraction, the statistical properties of the raw detected camera counts will be no longer maintained, it may decrease localization precisions';
handles.data_import_pushbutton.TooltipString = 'import the 4Pi single-molecule dataset';
handles.data_show_img_pushbutton.TooltipString = 'display the projection images of the 4Pi single-molecule dataset';
handles.data_show_bg_pushbutton.TooltipString = 'display the background of the 4Pi single-molecule dataset';
handles.data_show_merge_pushbutton.TooltipString = 'show the overlaid projection images from 4 channels compared to first channel';

% channel registration module
handles.registration_import_pushbutton.TooltipString = 'import the alignment calibration file to align the images from the rest channels to first channel';
handles.registration_cal_pushbutton.TooltipString = 'calculate the position relationship between 4 detection channels using affine transformation';
handles.registration_process_pushbutton.TooltipString = 'align the images from 4 channels based on the alignment calibration file';
handles.registration_export_pushbutton.TooltipString = 'export the alignment calibration file';
handles.registration_merge_pushbutton.TooltipString = 'show the overlaid projection images from 4 channels after alignment calibration';


% segmentation module
handles.text2.TooltipString = 'sub-region size of the cropped sub-regions from the single molecule dataset (unit: pixels). Recommended value is from 28 to 40';
handles.text3.TooltipString = 'distance threshold to make sure that each selected sub-region contains only one molecule (unit: pixels). Recommend value is around Box size divide by sqrt(2)';
handles.text4.TooltipString = 'initial intensity threshold to obtain the candidate sub-regions';
handles.text5.TooltipString = 'segmentation threshold to select sub-regions with higher photon counts compared to initial candidate sub-regions';
handles.seg_process_pushbutton.TooltipString = 'crop sub-regions from the dataset.The details can be seen in the software instruction';
handles.seg_export_pushbutton.TooltipString = 'export the cropped sub-regions';
handles.text27.TooltipString = 'frame index in the single-molecule dataset';
handles.seg_show_img_pushbutton.TooltipString = 'show the cropped sub-regions with the given frame index';

% display module
handles.display_import_pushbutton.TooltipString = 'import the 3D super-resolution reconstruction result by users';
handles.display_set_pushbutton.TooltipString = 'set display parameters';
handles.display_show_pushbutton.TooltipString = 'show the x-y view of the reconstructed super-resolution image with each molecule color-coded by its axial position';
handles.display_export_pushbutton.TooltipString = 'export 3D the super-resolution image';

% in situ model generation module
handles.pupil_data_checkbox.TooltipString = 'If this box is selected, the previous sub-regions can be imported. If not selected, 4Pi-brainspot uses the sub-regions from 4Pi-PSF segmentation module.';
handles.pupil_import_pushbutton.TooltipString = 'import the sub-regions by users';
handles.pupil_zernike_checkbox.TooltipString = 'If this box is selected, the users can modify the initial coefficients of 21 Zernike modes (Wyant order, from vertical astigmatism to tertiary spherical aberration)';
handles.pupil_change_pushbutton.TooltipString = 'modify the initial coefficients of 21 Zernike modes';
handles.text6.TooltipString = 'three input texts from left to right are minimum axial position, axial step size, and maximum axial position';
handles.text24.TooltipString = 'lateral alignment mode';
handles.text8.TooltipString = 'lateral and axial positions optimization mode for averaged sub-regions';
handles.text10.TooltipString = 'number of output Zernike modes (Wyant order)';
handles.text7.TooltipString = 'number threshold to reject an axial position group which contains fewer sub-regions than this threshold. Recommended value is from 25 to 50.';
handles.text11.TooltipString = 'iteration number of in situ model generation. This process usually converged in 4 iterations.';
handles.text12.TooltipString = 'similarity threshold to reject a sub-region with similarity lower than this threshold. Recommended value is from 0.5 to 0.6';
handles.pupil_process_pushbutton.TooltipString = 'carry out in situ 3D PSF model generation. The details can be seen in the software instruction';
handles.pupil_stop_pushbutton.TooltipString = 'stop in situ 3D PSF model generation';
handles.pupil_export_pushbutton.TooltipString = 'export the in situ 3D PSF model';
handles.pupil_showPSFs_pushbutton.TooltipString = 'show the retrieved PSFs along the axial direction';
handles.pupil_show_pupil_pushbutton.TooltipString = 'show the retrieved pupil (including its magnitude and phase)';
handles.pupil_show_zernike_pushbutton.TooltipString = 'show the decomposed Zernike coefficients from the retrieved phase';
handles.pupil_showPSFs_fixed_pushbutton.TooltipString = 'show the retrieved PSFs at the fixed axial direction';
handles.pupil_zshift_mode_popupmenu.TooltipString = 'Z shift mode = Shift, meaning the lateral and axial positions are optimized. Z shift mode = No shift, meaning the lateral and axial positions are not optimized';

% 3D localization module
handles.recon_data_checkbox.TooltipString = 'If this box is selected, the single molecule dataset can be imported. If not selected, 4Pi-brainspot uses the dataset from data import module';
handles.recon_import_data_pushbutton.TooltipString = 'import the single molecule dataset by users';
handles.recon_tform_checkbox.TooltipString = 'If this box is selected, the alignment calibration file can be imported. If not selected, 4Pi-brainspot toolbox uses the calibration file from channel alignment module';
handles.recon_import_tform_pushbutton.TooltipString = 'import the alignment calibration file by users';
handles.recon_pupil_checkbox.TooltipString = 'If this box is selected, in situ 3D PSF model can be imported. If not selected, 4Pi-brainspot toolbox uses the in situ 3D PSF model from in situ model generation module';
handles.recon_import_pupil_pushbutton.TooltipString = 'import the in situ 3D PSF model by users';
handles.recon_seg_checkbox.TooltipString = 'If this box is selected, the threshold can be set to crop sub-regions';
handles.text14.TooltipString = 'initial and segmentation thresholds. Same as the setting of 4Pi-PSF segmentation module';
handles.recon_dc_checkbox.TooltipString = 'If this box is selected, 4Pi-brainspot will carry out 3D drift correction and volume alignment';
handles.recon_rej_checkbox.TooltipString = 'If this box is selected, 4Pi-brainspot will carry out the rejection process';
handles.recon_gpu_checkbox.TooltipString = 'If this box is selected, 4Pi-brainspot will run the GPU version for pupil-based 3D localization';
handles.recon_bg_checkbox.TooltipString = 'If this box is selected, 4Pi-brainspot will estimate background using temporal median filtering, then carry out background subtraction';
handles.recon_process_pushbutton.TooltipString = 'carry out 3D super-resolution reconstruction. The details can be seen in the software instruction';
handles.recon_stop_pushbutton.TooltipString = 'stop 3D super-resolution reconstruction';
handles.recon_export_pushbutton.TooltipString = 'export super-resolution reconstruction results';
handles.recon_export_csv_pushbutton.TooltipString = 'export (x, y, z) positions of single molecules by using csv format';


% UIWAIT makes brainspot_4pi_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = brainspot_4pi_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function seg_boxsz_edit_Callback(hObject, eventdata, handles)
% hObject    handle to seg_boxsz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_boxsz_edit as text
%        str2double(get(hObject,'String')) returns contents of seg_boxsz_edit as a double


% --- Executes during object creation, after setting all properties.
function seg_boxsz_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_boxsz_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in registration_import_pushbutton.
function registration_import_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to registration_import_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi
[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select the transformation model');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.registration_export_pushbutton,'Enable','on');
   set(handles.registration_process_pushbutton,'Enable','on');

   set(handles.registration_import_edit,'String', filename);

   % load file
   tmp = load([pathname filename]);
   
   if isfield(tmp,'tform_all')
       INSPR4pi.tform_all = tmp.tform_all;
       
       msgbox('Finish importing the transformation model!');
   else
       msgbox('Please import the correct data!');
   end
   
end



function registration_import_edit_Callback(hObject, eventdata, handles)
% hObject    handle to registration_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of registration_import_edit as text
%        str2double(get(hObject,'String')) returns contents of registration_import_edit as a double


% --- Executes during object creation, after setting all properties.
function registration_import_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to registration_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in registration_cal_pushbutton.
function registration_cal_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to registration_cal_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath('.\Channel_alignment\')

global INSPR4pi

disp('Calculate the transformation model...')

if isfield(INSPR4pi,'qd1') && isfield(INSPR4pi,'qd2') && isfield(INSPR4pi,'qd3') && isfield(INSPR4pi,'qd4')
    tform_all = registration_4pi(INSPR4pi.qd1,INSPR4pi.qd2,INSPR4pi.qd3,INSPR4pi.qd4);
    INSPR4pi.tform_all = tform_all;
    
    set(handles.registration_export_pushbutton,'Enable','on');
    set(handles.registration_process_pushbutton,'Enable','on');
    msgbox('Finish calculating the transformation model!');
    
else
    msgbox('Please import the data!');
end




% --- Executes on button press in registration_export_pushbutton.
function registration_export_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to registration_export_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

[filename, pathname] = uiputfile(fullfile(INSPR4pi.setup.workspace,'*.mat'), 'save the transformation model');
if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   disp(['User selected ',fullfile(pathname,filename)])

   tform_all = INSPR4pi.tform_all;
   save(fullfile(pathname,filename),'tform_all');
end



% --- Executes on button press in data_import_pushbutton.
function data_import_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to data_import_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi
[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select the 4Pi interference data');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.data_show_img_pushbutton,'Enable','on');
   set(handles.data_show_merge_pushbutton,'Enable','on');
   set(handles.data_import_edit,'String', filename);

   % load file
   tmp = load([pathname filename]);
   
   if isfield(tmp,'qd1') && isfield(tmp,'qd2') && isfield(tmp,'qd3') && isfield(tmp,'qd4')    
       if (size(tmp.qd1,1) == size(tmp.qd2,1) && size(tmp.qd1,2) == size(tmp.qd2,2) && size(tmp.qd1,3) == size(tmp.qd2,3))
       
           if isfield(INSPR4pi,'qd1')
               set(handles.registration_merge_pushbutton,'Enable','off');
               set(handles.seg_process_pushbutton,'Enable','off');
               %            set(handles.seg_subregion_pushbutton,'Enable','off');
               set(handles.seg_export_pushbutton,'Enable','off');
               set(handles.seg_show_img_pushbutton,'Enable','off');
           end
           
           
           if INSPR4pi.setup.is_sCMOS   %sCMOS case, not finish
               % sCMOS parameters
               
               offsetim_chs = repmat(INSPR4pi.setup.sCMOS_input.qd_offset,[1 1 size(tmp.qd1,3) 1]);  
               gainim_chs = repmat(INSPR4pi.setup.sCMOS_input.qd_gain,[1 1 size(tmp.qd1,3) 1]);
               
               qd1_in = (tmp.qd1 - offsetim_chs(:,:,:,1)) ./ gainim_chs(:,:,:,1);
               qd2_in = (tmp.qd2 - offsetim_chs(:,:,:,2)) ./ gainim_chs(:,:,:,2);
               qd3_in = (tmp.qd3 - offsetim_chs(:,:,:,3)) ./ gainim_chs(:,:,:,3);
               qd4_in = (tmp.qd4 - offsetim_chs(:,:,:,4)) ./ gainim_chs(:,:,:,4);
           else    %EMCCD case
               qd1_in = (tmp.qd1 - INSPR4pi.setup.offset) /INSPR4pi.setup.gain;
               qd2_in = (tmp.qd2 - INSPR4pi.setup.offset) /INSPR4pi.setup.gain;
               qd3_in = (tmp.qd3 - INSPR4pi.setup.offset) /INSPR4pi.setup.gain;
               qd4_in = (tmp.qd4 - INSPR4pi.setup.offset) /INSPR4pi.setup.gain;
           end
           qd1_in(qd1_in<=0) = 1e-6;
           qd2_in(qd2_in<=0) = 1e-6;
           qd3_in(qd3_in<=0) = 1e-6;
           qd4_in(qd4_in<=0) = 1e-6;

           
           INSPR4pi.display.imagesz = size(tmp.qd1,1);
           
           if (size(tmp.qd1,1) > 100)
               INSPR4pi.setup.is_imgsz = 1;
           else
               INSPR4pi.setup.is_imgsz = 0;
           end
                       
           
           if INSPR4pi.setup.is_bg == 1

               bg_img_1 = median(qd1_in,3);
               bg_img_2 = median(qd2_in,3);
               bg_img_3 = median(qd3_in,3);
               bg_img_4 = median(qd4_in,3);
                
               INSPR4pi.bg_img_1 = bg_img_1;
               INSPR4pi.bg_img_2 = bg_img_2;
               INSPR4pi.bg_img_3 = bg_img_3;
               INSPR4pi.bg_img_4 = bg_img_4;
               
               subtract_img_1 = qd1_in - bg_img_1;
               subtract_img_2 = qd2_in - bg_img_2;
               subtract_img_3 = qd3_in - bg_img_3;
               subtract_img_4 = qd4_in - bg_img_4;
               subtract_img_1(subtract_img_1<=0) = 1e-6;
               subtract_img_2(subtract_img_2<=0) = 1e-6;
               subtract_img_3(subtract_img_3<=0) = 1e-6;
               subtract_img_4(subtract_img_4<=0) = 1e-6;
               
               INSPR4pi.qd1 = subtract_img_1;
               INSPR4pi.qd2 = subtract_img_2;
               INSPR4pi.qd3 = subtract_img_3;
               INSPR4pi.qd4 = subtract_img_4;
               
               msgbox('Finish subtracting background and importing 4Pi data!');
           else
               INSPR4pi.qd1 = qd1_in;
               INSPR4pi.qd2 = qd2_in;
               INSPR4pi.qd3 = qd3_in;
               INSPR4pi.qd4 = qd4_in;
               
               msgbox('Finish importing 4Pi data!');
           end
       else
           msgbox('Channels should be the same size!');
       end
       
   else
       msgbox('Please import the correct data!');
   end
end


function data_import_edit_Callback(hObject, eventdata, handles)
% hObject    handle to data_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_import_edit as text
%        str2double(get(hObject,'String')) returns contents of data_import_edit as a double


% --- Executes during object creation, after setting all properties.
function data_import_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in data_show_img_pushbutton.
function data_show_img_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to data_show_img_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi


disp('Show maximum projection images in 4 channels...')

imouts = [];
for ii = 1 : 4
    im1 = eval(['INSPR4pi.qd', num2str(ii)]);
    tmp_img = max(im1,[],3);
    imouts = cat(2,imouts,tmp_img);
end

figure; imshow(imouts,[]);
axis tight
title('Maximum projection in 4 channels (P1, S2, P2, S1)');




% --- Executes on button press in data_show_merge_pushbutton.
function data_show_merge_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to data_show_merge_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

disp('Show merge data...')

imref = [];
imouts = [];
im1 = max(INSPR4pi.qd1,[],3);
for ii = 1 : 4
    im2 = eval(['INSPR4pi.qd', num2str(ii)]);
    tmp_img = max(im2,[],3);
    
    imref = cat(2,imref,im1);
    imouts = cat(2,imouts,tmp_img);
end

figure; imshowpair(imref, imouts,'ColorChannels',[1 2 0], 'Scaling','joint');
axis tight
title('Merged images in 4 channels compared to P1 channel');




% --- Executes on button press in setup_import_pushbutton.
function setup_import_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setup_import_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select the setup configuration');


if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.setup_export_pushbutton,'Enable','on');
   set(handles.setup_import_edit,'String',filename);


   % load file
   tmp = load([pathname filename]);
   if isfield(tmp,'setup')
       
       if isfield(tmp.setup,'Pixelsize')
           INSPR4pi.setup.Pixelsize = tmp.setup.Pixelsize;
       else
           msgbox('No pixel size. Use the default!');
       end
       
       if isfield(tmp.setup,'RefractiveIndex')
           INSPR4pi.setup.RefractiveIndex = tmp.setup.RefractiveIndex;
       else
           msgbox('No Refractive Index. Use the default!');
       end
       
       if isfield(tmp.setup,'nMed')
           INSPR4pi.setup.nMed = tmp.setup.nMed;
       else
           msgbox('No Medium Index. Use the default!');
       end
       
       if isfield(tmp.setup,'Lambda')
           INSPR4pi.setup.Lambda = tmp.setup.Lambda;
       else
           msgbox('No Lambda. Use the default!');
       end
       
       if isfield(tmp.setup,'NA')
           INSPR4pi.setup.NA = tmp.setup.NA;
       else
           msgbox('No NA. Use the default!');
       end
       
       if isfield(tmp.setup,'offset')
           INSPR4pi.setup.offset = tmp.setup.offset;
       else
           msgbox('No camera offset. Use the default!');
       end
       
       if isfield(tmp.setup,'gain')
           INSPR4pi.setup.gain = tmp.setup.gain;
       else
           msgbox('No camera gain. Use the default!');
       end
       
       msgbox('Finish importing setup configuration!');
       
   else
       
       msgbox('Please import the correct data!');
   end

end



% --- Executes on button press in setup_set_pushbutton.
function setup_set_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setup_set_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi 
%%
display('Setup parameters...');
% Draw set figure
f_set = figure('Position',[600 100 400 600],'Name','Setup parameters');

%static text
uicontrol(f_set, 'Style','text','String','Pixel size (um)','FontSize',10,'Position',[50 470 100 50],...
    'HorizontalAlignment','left','TooltipString','effective pixel size on the camera');
uicontrol(f_set, 'Style','text','String','Refractive index of immersion medium','FontSize',10,'Position',[50 430 115 50],...
    'HorizontalAlignment','left','TooltipString','refractive index of the immersion medium of the objective lens');
uicontrol(f_set, 'Style','text','String','Refractive index of sample medium','FontSize',10,'Position',[50 380 110 50],...
    'HorizontalAlignment','left','TooltipString','refractive index of the imaging medium');
uicontrol(f_set, 'Style','text','String','Lambda (um)','FontSize',10,'Position',[50 320 100 50],...
    'HorizontalAlignment','left','TooltipString','emission wavelength');
uicontrol(f_set, 'Style','text','String','NA','FontSize',10,'Position',[50 270 100 50],...
    'HorizontalAlignment','left','TooltipString','numerical aperture of the objective lens');
uicontrol(f_set, 'Style','text','String','Camera offset','FontSize',10,'Position',[50 220 100 50],...
    'HorizontalAlignment','left','TooltipString','offset on the camera');
uicontrol(f_set, 'Style','text','String','Camera gain','FontSize',10,'Position',[50 170 100 50],...
    'HorizontalAlignment','left','TooltipString','gain on the camera');

%edit 
h_Pixelsize = uicontrol(f_set, 'Style','edit','String',num2str(INSPR4pi.setup.Pixelsize),...
    'FontSize',10,'Position',[200 500 150 25]);
h_RefractiveIndex = uicontrol(f_set, 'Style','edit','String',num2str(INSPR4pi.setup.RefractiveIndex),...
    'FontSize',10,'Position',[200 450 150 25]);
h_nMed = uicontrol(f_set, 'Style','edit','String',num2str(INSPR4pi.setup.nMed),...
    'FontSize',10,'Position',[200 400 150 25]);
h_Lambda  = uicontrol(f_set, 'Style','edit','String',num2str(INSPR4pi.setup.Lambda),...
    'FontSize',10,'Position',[200 350 150 25]);
h_NA = uicontrol(f_set, 'Style','edit','String',num2str(INSPR4pi.setup.NA),...
    'FontSize',10,'Position',[200 300 150 25]);
h_offset = uicontrol(f_set, 'Style','edit','String',num2str(INSPR4pi.setup.offset),...
    'FontSize',10,'Position',[200 250 150 25]);
h_gain = uicontrol(f_set, 'Style','edit','String',num2str(INSPR4pi.setup.gain),...
    'FontSize',10,'Position',[200 200 150 25]);


%save figure handle
h_set.h_Pixelsize = h_Pixelsize;
h_set.h_RefractiveIndex = h_RefractiveIndex;
h_set.h_nMed = h_nMed;
h_set.h_Lambda = h_Lambda;
h_set.h_NA = h_NA;
h_set.h_offset = h_offset;
h_set.h_gain = h_gain;

% sCMOS part 
h_sCMOS_edit = uicontrol(f_set, 'Style','edit','String','',...
    'FontSize',10,'Position',[200 110 150 25]);

h_sCMOS_import = uicontrol(f_set,'Style','pushbutton','String','Import calibration file','Position',[50 110 130 25],'FontSize',10, 'Callback',...
    {@setup_import_sCMOS_pushbutton_Callback,h_sCMOS_edit},'TooltipString','import sCMOS calibration file');

h_set.h_sCMOS_edit = h_sCMOS_edit;
h_set.h_sCMOS_import = h_sCMOS_import;
h_sCMOS_chk = uicontrol(f_set, 'Style','checkbox','String','sCMOS camera mode','Value',INSPR4pi.setup.is_sCMOS,...
    'Position',[50 140 200 50],'FontSize',10,...
    'Callback',{@setup_sCMOS_checkBox_Callback,h_set},'TooltipString','If this check box is selected, INSPR will carry out sCMOS calibration. If not selected, INSPR will use EMCCD camera mode, which uses the same offset and gain for each pixel on the camera');

h_set.h_sCMOS_chk = h_sCMOS_chk;

if INSPR4pi.setup.is_sCMOS 
    set(h_set.h_offset,'Enable','off'); %EMCCD off
    set(h_set.h_gain,'Enable','off');
else
    set(h_set.h_sCMOS_edit,'Enable','off'); %sCMOS off
    set(h_set.h_sCMOS_import,'Enable','off');
end

% reset and save button
uicontrol(f_set,'Style','pushbutton','String','Reset','Position',[50 40 100 30],'FontSize',10, 'Callback',...
    {@setup_reset_pushbutton_Callback,h_set});

uicontrol(f_set,'Style','pushbutton','String','Save','Position',[200 40 100 30],'FontSize',10, 'Callback',...
    {@setup_save_pushbutton_Callback,h_set});
%%

set(handles.setup_export_pushbutton,'Enable','on');

%Reset parameters
function setup_sCMOS_checkBox_Callback(src,event,t)
global INSPR4pi

INSPR4pi.setup.is_sCMOS = get(src,'Value');

if INSPR4pi.setup.is_sCMOS
    set(t.h_offset,'Enable','off'); %EMCCD off
    set(t.h_gain,'Enable','off');
    
    set(t.h_sCMOS_edit,'Enable','on'); %sCMOS on
    set(t.h_sCMOS_import,'Enable','on');
else
    set(t.h_offset,'Enable','on');  %EMCCD on
    set(t.h_gain,'Enable','on');
    
    set(t.h_sCMOS_edit,'Enable','off'); %sCMOS off
    set(t.h_sCMOS_import,'Enable','off');
end


function setup_import_sCMOS_pushbutton_Callback(src,event,t)

global INSPR4pi

[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select the sCMOS calibration file');


if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   % load file
   sCMOS_input = load([pathname filename]);   
   set(t,'String', filename);
   
   if isfield(sCMOS_input,'qd_offset') && isfield(sCMOS_input,'qd_var') && ...
           isfield(sCMOS_input,'qd_gain')
       INSPR4pi.setup.sCMOS_input = sCMOS_input;
       
       msgbox('Finish importing sCMOS calibration file!'); 
   else
       msgbox('sCMOS calibration file import error!');
   end
   
end


%Reset parameters
function setup_reset_pushbutton_Callback(src,event,t)
global INSPR4pi

default_cfg = load('default_cfg.mat');

%update INSPR4pi
INSPR4pi.setup.Pixelsize = default_cfg.setup.Pixelsize;
INSPR4pi.setup.RefractiveIndex = default_cfg.setup.RefractiveIndex;
INSPR4pi.setup.nMed = default_cfg.setup.nMed;
INSPR4pi.setup.Lambda = default_cfg.setup.Lambda;
INSPR4pi.setup.NA = default_cfg.setup.NA;
INSPR4pi.setup.offset = default_cfg.setup.offset;
INSPR4pi.setup.gain = default_cfg.setup.gain;
INSPR4pi.setup.is_sCMOS = default_cfg.setup.is_sCMOS;
INSPR4pi.setup.sCMOS_input = default_cfg.setup.sCMOS_input;

%update UI
set(t.h_Pixelsize,'String', num2str(INSPR4pi.setup.Pixelsize));
set(t.h_RefractiveIndex,'String', num2str(INSPR4pi.setup.RefractiveIndex));
set(t.h_nMed,'String', num2str(INSPR4pi.setup.nMed));
set(t.h_Lambda,'String', num2str(INSPR4pi.setup.Lambda));
set(t.h_NA,'String', num2str(INSPR4pi.setup.NA));
set(t.h_offset,'String', num2str(INSPR4pi.setup.offset));
set(t.h_gain,'String', num2str(INSPR4pi.setup.gain));
set(t.h_offset,'Enable','on'); 
set(t.h_gain,'Enable','on'); 

set(t.h_sCMOS_chk,'Value',INSPR4pi.setup.is_sCMOS); %sCMOS
set(t.h_sCMOS_import,'Enable','off'); 
set(t.h_sCMOS_edit,'Enable','off'); 
set(t.h_sCMOS_edit,'String', ''); 


%save set parameters
function setup_save_pushbutton_Callback(src,event,t)
global INSPR4pi 

Pixelsize = str2num( get(t.h_Pixelsize,'String') ); 
RefractiveIndex = str2num( get(t.h_RefractiveIndex,'String') ); 
nMed = str2num( get(t.h_nMed,'String') ); 
Lambda = str2num( get(t.h_Lambda,'String') ); 
NA = str2num( get(t.h_NA,'String') ); 
offset = str2num( get(t.h_offset,'String') ); 
gain = str2num( get(t.h_gain,'String') ); 


INSPR4pi.setup.Pixelsize = Pixelsize;
INSPR4pi.setup.RefractiveIndex = RefractiveIndex;
INSPR4pi.setup.nMed = nMed;
INSPR4pi.setup.Lambda = Lambda;
INSPR4pi.setup.NA = NA;
INSPR4pi.setup.offset = offset;
INSPR4pi.setup.gain = gain;




function setup_import_edit_Callback(hObject, eventdata, handles)
% hObject    handle to setup_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setup_import_edit as text
%        str2double(get(hObject,'String')) returns contents of setup_import_edit as a double


% --- Executes during object creation, after setting all properties.
function setup_import_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setup_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in setup_export_pushbutton.
function setup_export_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setup_export_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

[filename, pathname] = uiputfile(fullfile(INSPR4pi.setup.workspace,'*.mat'), 'save the setup configuration');
if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   disp(['User selected ',fullfile(pathname,filename)])

   setup = INSPR4pi.setup;
   save(fullfile(pathname,filename),'setup');
end


function seg_thresh_dist_edit_Callback(hObject, eventdata, handles)
% hObject    handle to seg_thresh_dist_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_thresh_dist_edit as text
%        str2double(get(hObject,'String')) returns contents of seg_thresh_dist_edit as a double


% --- Executes during object creation, after setting all properties.
function seg_thresh_dist_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_thresh_dist_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function seg_thresh_low_edit_Callback(hObject, eventdata, handles)
% hObject    handle to seg_thresh_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_thresh_low_edit as text
%        str2double(get(hObject,'String')) returns contents of seg_thresh_low_edit as a double


% --- Executes during object creation, after setting all properties.
function seg_thresh_low_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_thresh_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function seg_thresh_high_edit_Callback(hObject, eventdata, handles)
% hObject    handle to seg_thresh_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_thresh_high_edit as text
%        str2double(get(hObject,'String')) returns contents of seg_thresh_high_edit as a double


% --- Executes during object creation, after setting all properties.
function seg_thresh_high_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_thresh_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in seg_subregion_pushbutton.
% function seg_subregion_pushbutton_Callback(hObject, eventdata, handles)
% % hObject    handle to seg_subregion_pushbutton (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% global INSPR4pi
% 
% show_ch1 = [];
% show_ch2 = [];
% 
% for i = 1 : size(INSPR4pi.subregion_ch1,3)
%     tmp1 = INSPR4pi.subregion_ch1(:,:,i);
%     show_ch1(:,:,i) = tmp1./max(tmp1(:));
%     tmp2 = INSPR4pi.subregion_ch2(:,:,i);
%     show_ch2(:,:,i) = tmp2./max(tmp2(:));
% end
% 
% % show PSFs
% col_img = ones(size(INSPR4pi.subregion_ch1,1),2,size(INSPR4pi.subregion_ch1,3));
% 
% PSFs_all = cat(2, show_ch1,col_img,show_ch2);
% h_subregion = dipshow(PSFs_all);
% diptruesize(h_subregion,800)




% --- Executes on button press in seg_process_pushbutton.
function seg_process_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to seg_process_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi
addpath('.\Segmentation\');

if isfield(INSPR4pi,'qd2_registered') 
    display('Segmentation process...');
    boxsz = str2num( get(handles.seg_boxsz_edit, 'String') );
    thresh_dist = str2num( get(handles.seg_thresh_dist_edit, 'String') );
    thresh_low = str2num( get(handles.seg_thresh_low_edit, 'String') );
    thresh_high = str2num( get(handles.seg_thresh_high_edit, 'String') );
    
    %find subregion
    thresh = [thresh_low thresh_high];
    [subregion_chs,seg_display] = crop_subregion_4pi(INSPR4pi.qd1,INSPR4pi.qd2,INSPR4pi.qd3,INSPR4pi.qd4,INSPR4pi.tform_all,boxsz,thresh,thresh_dist,INSPR4pi.setup);

    INSPR4pi.subregion_chs = subregion_chs;
    INSPR4pi.seg_display = seg_display;

    num_subregion = size(subregion_chs,3);
    msgbox({'Finish segmentation!' ['There are ' num2str(num_subregion) ' subregion images']});
    if num_subregion < 1000
         msgbox('Warning! Less number of sub regions!');
    end
    
    display(['Finish segmentation! There are ' num2str(num_subregion) ' subregion images']);
    
%     set(handles.seg_subregion_pushbutton,'Enable','on');
    set(handles.seg_export_pushbutton,'Enable','on');
    set(handles.seg_show_img_pushbutton,'Enable','on');

else
    msgbox('Please process channel registration!');
end


% --- Executes on button press in pupil_data_checkbox.
function pupil_data_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_data_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pupil_data_checkbox


if get(handles.pupil_data_checkbox,'Value')==1
   set(handles.pupil_import_pushbutton,'Enable','on');
else
   set(handles.pupil_import_pushbutton,'Enable','off'); 
end


% --- Executes on button press in pupil_import_pushbutton.
function pupil_import_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_import_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi
[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select the 4Pi subregion images..');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.pupil_import_edit,'String',filename);

   % load file
   tmp = load([pathname filename]);
   
   if isfield(tmp,'subregion_chs')
       if (size(tmp.subregion_chs,4) == 4)
           INSPR4pi.subregion_chs = tmp.subregion_chs;
           
           msgbox('Finish importing subregion data!');
       else
           msgbox('Channel number is not correct!');
       end
   else
       msgbox('Please import the correct data!');
   end
end

function pupil_import_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_import_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_import_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_import_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function pupil_Zrange_low_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_Zrange_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_Zrange_low_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_Zrange_low_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_Zrange_low_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_Zrange_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pupil_Zrange_mid_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_Zrange_mid_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_Zrange_mid_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_Zrange_mid_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_Zrange_mid_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_Zrange_mid_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pupil_Zrange_high_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_Zrange_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_Zrange_high_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_Zrange_high_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_Zrange_high_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_Zrange_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pupil_bin_low_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_bin_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_bin_low_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_bin_low_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_bin_low_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_bin_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pupil_zshift_mode_popupmenu.
function pupil_zshift_mode_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_zshift_mode_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pupil_zshift_mode_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pupil_zshift_mode_popupmenu
global INSPR4pi

contents = get(handles.pupil_zshift_mode_popupmenu,'String'); 

switch contents{get(handles.pupil_zshift_mode_popupmenu,'Value')}
    case 'No shift'
        INSPR4pi.pupil.Zshift_mode = 0;   
    case 'Shift'
        INSPR4pi.pupil.Zshift_mode = 1;
    otherwise
        INSPR4pi.pupil.Zshift_mode = 0;
end

display(['Z shift mode is: ' num2str(INSPR4pi.pupil.Zshift_mode)]);


% --- Executes during object creation, after setting all properties.
function pupil_zshift_mode_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_zshift_mode_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in recon_data_checkbox.
function recon_data_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_data_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_data_checkbox
global INSPR4pi

if get(handles.recon_data_checkbox,'Value')==1
    INSPR4pi.recon.isNewdata = 1;
    set(handles.recon_import_data_pushbutton,'Enable','on');
else
    INSPR4pi.recon.isNewdata = 0; 
    set(handles.recon_import_data_pushbutton,'Enable','off');
end


% --- Executes on button press in recon_import_data_pushbutton.
function recon_import_data_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_import_data_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi
[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select the 4Pi data (you can choose multiple data)...','MultiSelect','on');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.recon_import_data_edit,'String', pathname);
   
   % load one file or muti-files
   if iscell(filename)
       dirN = numel(filename);
       INSPR4pi.recon.datafile_name = filename;
   else
       dirN = 1;
       INSPR4pi.recon.datafile_name{1} = filename; 
   end
   
   %save file path
   INSPR4pi.recon.dirN = dirN;
   INSPR4pi.recon.datapath = pathname;
   
   msgbox({ 'Finish importing the 4Pi data path!' ['You choose ' num2str(dirN) ' data'] });


   
end

function recon_import_data_edit_Callback(hObject, eventdata, handles)
% hObject    handle to recon_import_data_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_import_data_edit as text
%        str2double(get(hObject,'String')) returns contents of recon_import_data_edit as a double


% --- Executes during object creation, after setting all properties.
function recon_import_data_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_import_data_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in recon_tform_checkbox.
function recon_tform_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_tform_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_tform_checkbox
if get(handles.recon_tform_checkbox,'Value')==1
   set(handles.recon_import_tform_pushbutton,'Enable','on');
else
   set(handles.recon_import_tform_pushbutton,'Enable','off'); 
end



% --- Executes on button press in recon_import_tform_pushbutton.
function recon_import_tform_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_import_tform_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi
[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select the transformation model');

if isequal(filename,0)
    disp('User selected Cancel')
else
    disp(['User selected ', fullfile(pathname, filename)])
    
    set(handles.recon_import_tform_edit,'String',filename);
    
    % load file
    tmp = load([pathname filename]);
    
    if isfield(tmp,'tform_all')
        INSPR4pi.tform_all = tmp.tform_all;
        
        msgbox('Finish importing the transformation model!');
    else
        msgbox('Please import the correct data!');
    end
    
end



function recon_import_tform_edit_Callback(hObject, eventdata, handles)
% hObject    handle to recon_import_tform_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_import_tform_edit as text
%        str2double(get(hObject,'String')) returns contents of recon_import_tform_edit as a double


% --- Executes during object creation, after setting all properties.
function recon_import_tform_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_import_tform_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in recon_pupil_checkbox.
function recon_pupil_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_pupil_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_pupil_checkbox
global INSPR4pi

if get(handles.recon_pupil_checkbox,'Value')==1
    INSPR4pi.recon.isNewpupil = 1;
    set(handles.recon_import_pupil_pushbutton,'Enable','on');
else
    INSPR4pi.recon.isNewpupil = 0;

    set(handles.recon_import_pupil_pushbutton,'Enable','off');
end



% --- Executes on button press in recon_import_pupil_pushbutton.
function recon_import_pupil_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_import_pupil_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi
[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select the pupil model (you can choose multiple pupil)...','MultiSelect','on');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.recon_import_pupil_edit,'String', pathname);

   % load one file or muti-files
   probj_all = [];
   loopn = [];
   if iscell(filename)
       loopn = numel(filename);
       for ii = 1 : loopn
           input_file = [pathname filename{ii}];
           tmp = load(input_file);
           if isfield(tmp,'probj')
               probj_all{ii} = tmp.probj;
           else
               msgbox('Please import the correct data!');
               return;
           end
       end
       msgbox({'Finish importing the pupil model!' ['You choose ' num2str(loopn) ' pupil model']});
   else
       loopn = 1;
       input_file = [pathname filename];
       tmp = load(input_file);
       if isfield(tmp,'probj')                 
           probj_all{1} = tmp.probj;
       else
           msgbox('Please import the correct data!');
           return;
       end
       msgbox({'Finish importing the pupil model!' ['You choose ' num2str(loopn) ' pupil model']});
   end
   
   %save all pupil models to INSPR4pi
   INSPR4pi.recon.probj_all = probj_all;
   INSPR4pi.recon.loopn = loopn;
   
end






function recon_import_pupil_edit_Callback(hObject, ~, handles)
% hObject    handle to recon_import_pupil_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_import_pupil_edit as text
%        str2double(get(hObject,'String')) returns contents of recon_import_pupil_edit as a double


% --- Executes during object creation, after setting all properties.
function recon_import_pupil_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_import_pupil_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in recon_seg_checkbox.
function recon_seg_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_seg_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_seg_checkbox
global INSPR4pi

if get(handles.recon_seg_checkbox,'Value')==1
    set(handles.recon_seg_thresh_low_edit,'Enable','on');
    set(handles.recon_seg_thresh_high_edit,'Enable','on');
    
    INSPR4pi.recon.isSeg = 1;
else
    set(handles.recon_seg_thresh_low_edit,'String',num2str(INSPR4pi.recon.seg_thresh_low));
    set(handles.recon_seg_thresh_high_edit,'String',num2str(INSPR4pi.recon.seg_thresh_high));
    set(handles.recon_seg_thresh_low_edit,'Enable','off');
    set(handles.recon_seg_thresh_high_edit,'Enable','off');
    
    INSPR4pi.recon.isSeg = 0;
end




function recon_seg_thresh_low_edit_Callback(hObject, eventdata, handles)
% hObject    handle to recon_seg_thresh_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_seg_thresh_low_edit as text
%        str2double(get(hObject,'String')) returns contents of recon_seg_thresh_low_edit as a double


% --- Executes during object creation, after setting all properties.
function recon_seg_thresh_low_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_seg_thresh_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function recon_seg_thresh_high_edit_Callback(hObject, eventdata, handles)
% hObject    handle to recon_seg_thresh_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_seg_thresh_high_edit as text
%        str2double(get(hObject,'String')) returns contents of recon_seg_thresh_high_edit as a double


% --- Executes during object creation, after setting all properties.
function recon_seg_thresh_high_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_seg_thresh_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in recon_rej_checkbox.
function recon_rej_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_rej_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_rej_checkbox
global INSPR4pi

if get(handles.recon_rej_checkbox,'Value')==1
    INSPR4pi.recon.isRej = 1;
    
    set(handles.recon_rej_change_pushbutton,'Enable','on');
else
    INSPR4pi.recon.isRej = 0;
    
    INSPR4pi.recon.rej.min_photon = 1500;   %rejection parameters
    INSPR4pi.recon.rej.llthreshold = 2500;
    INSPR4pi.recon.rej.loc_uncerxy_max = 15;  % nm
    INSPR4pi.recon.rej.loc_uncerz_max = 10;
    INSPR4pi.recon.rej.zmask_low = -0.6;    % um
    INSPR4pi.recon.rej.zmask_high = 0.6;
    
    set(handles.recon_rej_change_pushbutton,'Enable','off');
end


% --- Executes on button press in recon_rej_change_pushbutton.
function recon_rej_change_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_rej_change_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi 

% Draw rejection figure
f_rej = figure('Position',[600 200 400 450],'Name','Rejection parameters');

%static text
uicontrol(f_rej, 'Style','text','String','Min photon','FontSize',10,'Position',[50 350 100 50],...
    'HorizontalAlignment','left','TooltipString','photon threshold to reject single molecules with photon counts lower than this value');
uicontrol(f_rej, 'Style','text','String','LLR','FontSize',10,'Position',[50 300 100 50],...
    'HorizontalAlignment','left','TooltipString','log-likelihood ratio (LLR) threshold to reject single molecules with LLR higher than this value');
uicontrol(f_rej, 'Style','text','String','Max xy uncertainty (nm)','FontSize',10,'Position',[50 250 100 50],...
    'HorizontalAlignment','left','TooltipString','localization uncertainty in the lateral dimension to reject single molecules with uncertainty higher than this value');
uicontrol(f_rej, 'Style','text','String','Max Z uncertainty (nm)','FontSize',10,'Position',[50 200 100 50],...
    'HorizontalAlignment','left','TooltipString','localization uncertainty in the axial dimension to reject single molecules with uncertainty higher than this value');
uicontrol(f_rej, 'Style','text','String','Z mask (um)','FontSize',10,'Position',[50 150 100 50],...
    'HorizontalAlignment','left','TooltipString','axial position mask to reject single molecules with axial positions beyond this range');

%edit 
h_photon  = uicontrol(f_rej, 'Style','edit','String',num2str(INSPR4pi.recon.rej.min_photon),...
    'FontSize',10,'Position',[180 380 150 25]);
h_ll = uicontrol(f_rej, 'Style','edit','String',num2str(INSPR4pi.recon.rej.llthreshold),...
    'FontSize',10,'Position',[180 330 150 25]);
h_uncerxy = uicontrol(f_rej, 'Style','edit','String',num2str(INSPR4pi.recon.rej.loc_uncerxy_max,'%0.1f'),...
    'FontSize',10,'Position',[180 280 150 25]);
h_uncerz = uicontrol(f_rej, 'Style','edit','String',num2str(INSPR4pi.recon.rej.loc_uncerz_max,'%0.1f'),...
    'FontSize',10,'Position',[180 230 150 25]);
% h_uncer = uicontrol(f_rej, 'Style','edit','String',num2str(INSPR4pi.recon.rej.loc_uncer_max,'%0.3f'),...
%     'FontSize',10,'Position',[180 230 150 25]);
h_zmask_low = uicontrol(f_rej, 'Style','edit','String',num2str(INSPR4pi.recon.rej.zmask_low,'%0.3f'),...
    'FontSize',10,'Position',[180 180 70 25]);
h_zmask_high = uicontrol(f_rej, 'Style','edit','String',num2str(INSPR4pi.recon.rej.zmask_high,'%0.3f'),...
    'FontSize',10,'Position',[260 180 70 25]);

h_rej.h_photon = h_photon;
h_rej.h_ll = h_ll;
% h_rej.h_uncer = h_uncer;
h_rej.h_uncerxy = h_uncerxy;
h_rej.h_uncerz = h_uncerz;
h_rej.h_zmask_low = h_zmask_low;
h_rej.h_zmask_high = h_zmask_high;

%save button
uicontrol(f_rej,'Style','pushbutton','String','Save','Position',[180 80 100 30],'FontSize',10,'Callback',...
    {@recon_rej_save_pushbutton_Callback,h_rej});

%Reset button
uicontrol(f_rej,'Style','pushbutton','String','Reset','Position',[50 80 100 30],'FontSize',10,'Callback',...
    {@recon_rej_reset_pushbutton_Callback,h_rej});


%Reset parameters
function recon_rej_reset_pushbutton_Callback(src,event,t)
global INSPR4pi

default_cfg = load('default_cfg.mat');

%update INSPR4pi
INSPR4pi.recon.rej.min_photon = default_cfg.recon.rej.min_photon;
INSPR4pi.recon.rej.llthreshold = default_cfg.recon.rej.llthreshold;
% INSPR4pi.recon.rej.loc_uncer_max = default_cfg.recon.rej.loc_uncer_max;
INSPR4pi.recon.rej.loc_uncerxy_max = default_cfg.recon.rej.loc_uncerxy_max;
INSPR4pi.recon.rej.loc_uncerz_max = default_cfg.recon.rej.loc_uncerz_max;
INSPR4pi.recon.rej.zmask_low = default_cfg.recon.rej.zmask_low;
INSPR4pi.recon.rej.zmask_high = default_cfg.recon.rej.zmask_high;


%update UI
set(t.h_photon,'String', num2str(INSPR4pi.recon.rej.min_photon));
set(t.h_ll,'String', num2str(INSPR4pi.recon.rej.llthreshold));
% set(t.h_uncer,'String', num2str(INSPR4pi.recon.rej.loc_uncer_max));
set(t.h_uncerxy,'String', num2str(INSPR4pi.recon.rej.loc_uncerxy_max));
set(t.h_uncerz,'String', num2str(INSPR4pi.recon.rej.loc_uncerz_max));
set(t.h_zmask_low,'String', num2str(INSPR4pi.recon.rej.zmask_low));
set(t.h_zmask_high,'String', num2str(INSPR4pi.recon.rej.zmask_high));


%save rejection parameters
function recon_rej_save_pushbutton_Callback(src,event,t)
global INSPR4pi 

min_photon = str2num( get(t.h_photon,'String') ); 
llthreshold = str2num( get(t.h_ll,'String') ); 
% loc_uncer_max = str2num( get(t.h_uncer,'String') ); 
loc_uncerxy_max = str2num( get(t.h_uncerxy,'String') );
loc_uncerz_max = str2num( get(t.h_uncerz,'String') ); 
zmask_low = str2num( get(t.h_zmask_low,'String') ); 
zmask_high = str2num( get(t.h_zmask_high,'String') ); 

INSPR4pi.recon.rej.min_photon = min_photon;
INSPR4pi.recon.rej.llthreshold = llthreshold;
% INSPR4pi.recon.rej.loc_uncer_max = loc_uncer_max;
INSPR4pi.recon.rej.loc_uncerxy_max = loc_uncerxy_max;
INSPR4pi.recon.rej.loc_uncerz_max = loc_uncerz_max;
INSPR4pi.recon.rej.zmask_low = zmask_low;
INSPR4pi.recon.rej.zmask_high = zmask_high;





% --- Executes on button press in recon_dc_checkbox.
function recon_dc_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_dc_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_dc_checkbox

global INSPR4pi

if get(handles.recon_dc_checkbox,'Value')==1
    INSPR4pi.recon.isDC = 1;
    
    set(handles.recon_dc_change_pushbutton,'Enable','on');
else
    INSPR4pi.recon.isDC = 0;
    
    INSPR4pi.recon.dc.frmpfile = 1999;
    INSPR4pi.recon.dc.z_offset = 2000;
    INSPR4pi.recon.dc.step_ini = 400;
    set(handles.recon_dc_change_pushbutton,'Enable','off');
end




% --- Executes on button press in recon_dc_change_pushbutton.
function recon_dc_change_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_dc_change_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi 

% Draw DC figure
f_dc = figure('Position',[600 200 350 300],'Name','Drift correction parameters');

% static text
uicontrol(f_dc, 'Style','text','String','Frame bin','FontSize',10,'Position',[50 200 100 50],...
    'HorizontalAlignment','left','TooltipString','number of frames used to construct individual 3D volumes for 3D drift correction');
uicontrol(f_dc, 'Style','text','String','Initial Z offset (nm)','FontSize',10,'Position',[50 150 150 50],...
    'HorizontalAlignment','left','TooltipString','axial position offset to avoid negative z positions during 3D volume reconstruction');
uicontrol(f_dc, 'Style','text','String','Step interval (nm)','FontSize',10,'Position',[50 100 150 50],...
    'HorizontalAlignment','left','TooltipString','axial step size for multi-section imaging');

% edit 
h_frmp  = uicontrol(f_dc, 'Style','edit','String',num2str(INSPR4pi.recon.dc.frmpfile),...
    'FontSize',10,'Position',[180 230 110 25]);
h_z_offset = uicontrol(f_dc, 'Style','edit','String',num2str(INSPR4pi.recon.dc.z_offset),...
    'FontSize',10,'Position',[180 180 110 25]);
h_step = uicontrol(f_dc, 'Style','edit','String',num2str(INSPR4pi.recon.dc.step_ini),...
    'FontSize',10,'Position',[180 130 110 25]);


h_dc.h_frmp = h_frmp;
h_dc.h_z_offset = h_z_offset;
h_dc.h_step = h_step;

% save button
uicontrol(f_dc,'Style','pushbutton','String','Save','Position',[180 50 100 30],'FontSize',10,'Callback',...
    {@recon_dc_save_pushbutton_Callback,h_dc});
% Reset button
uicontrol(f_dc,'Style','pushbutton','String','Reset','Position',[50 50 100 30],'FontSize',10,'Callback',...
    {@recon_dc_reset_pushbutton_Callback,h_dc});

%Reset parameters
function recon_dc_reset_pushbutton_Callback(src,event,t)
global INSPR4pi

default_cfg = load('default_cfg.mat');

%update INSPR4pi
INSPR4pi.recon.dc.frmpfile = default_cfg.recon.dc.frmpfile;
INSPR4pi.recon.dc.z_offset = default_cfg.recon.dc.z_offset;
INSPR4pi.recon.dc.step_ini = default_cfg.recon.dc.step_ini;


%update UI
set(t.h_frmp,'String', num2str(INSPR4pi.recon.dc.frmpfile));
set(t.h_z_offset,'String', num2str(INSPR4pi.recon.dc.z_offset));
set(t.h_step,'String', num2str(INSPR4pi.recon.dc.step_ini));


%save rejection parameters
function recon_dc_save_pushbutton_Callback(src,event,t)
global INSPR4pi 

frmpfile = str2num( get(t.h_frmp,'String') ); 
z_offset = str2num( get(t.h_z_offset,'String') ); 
step_ini = str2num( get(t.h_step,'String') ); 

INSPR4pi.recon.dc.frmpfile = frmpfile;
INSPR4pi.recon.dc.z_offset = z_offset;
INSPR4pi.recon.dc.step_ini = step_ini;


% --- Executes on button press in recon_process_pushbutton.
function recon_process_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_process_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath('.\4Pi_localization\');
addpath('..\Support\PSF Toolbox_4pi\');

global INSPR4pi
global recon_stop

recon_stop = 0; %initil stop label is none

% check the data
if INSPR4pi.recon.isNewdata == 0
    if ~isfield(INSPR4pi,'qd1') || ~isfield(INSPR4pi,'qd2') || ~isfield(INSPR4pi,'qd3') || ~isfield(INSPR4pi,'qd4')
        msgbox('Please import the data!');
        return;
    end
    INSPR4pi.recon.qd1 = INSPR4pi.qd1;
    INSPR4pi.recon.qd2 = INSPR4pi.qd2;
    INSPR4pi.recon.qd3 = INSPR4pi.qd3;
    INSPR4pi.recon.qd4 = INSPR4pi.qd4;
    INSPR4pi.recon.dirN = 1;
    INSPR4pi.recon.datapath = pwd;
else
    if ~isfield(INSPR4pi.recon,'datapath')
        msgbox('Please import the data!');
        return;
    end
end

if ~isfield(INSPR4pi,'tform_all')
    msgbox('Please import transformation model!');
    return;
end

if INSPR4pi.recon.isNewpupil == 0
    if ~isfield(INSPR4pi,'probj')
        msgbox('Please import pupil model!');
        return;
    end
    INSPR4pi.recon.probj_all{1} = INSPR4pi.probj;
    INSPR4pi.recon.loopn = 1;
else
    if ~isfield(INSPR4pi.recon,'probj_all')
        msgbox('Please import the pupil model!');
        return;
    end
end

if ~isfield(INSPR4pi,'sobj_phi0')
    msgbox('Please import cavity phase!');
    return;
end

if length(INSPR4pi.sobj_phi0) ~= INSPR4pi.recon.dirN
    msgbox('Please make sure the number of estimated cavity phase match with dataset!');
    return;
end

if INSPR4pi.recon.isSeg == 1
    INSPR4pi.recon.seg_thresh_low = str2num( get(handles.recon_seg_thresh_low_edit, 'String') );
    INSPR4pi.recon.seg_thresh_high = str2num( get(handles.recon_seg_thresh_high_edit, 'String') );    
end

if INSPR4pi.setup.is_sCMOS
    disp('Consider independent readout noise in sCMOS camera');
end

if INSPR4pi.recon.isGPU
    %CUDA check
    try, listGPUs;
    catch
        msgbox('Please install CUDA environment! https://developer.nvidia.com/cuda-75-downloads-archive');
        return;
    end  
end

% if INSPR4pi.recon.is_bg == 1    % whether subtract background
%     INSPR4pi.setup.is_bg = 1; 
% else
%     INSPR4pi.setup.is_bg = 0;
% end
    
% 4Pi-SMSN pupil fitting
srobj = analysis3D_from4pi_pupil_v2(INSPR4pi.recon, INSPR4pi.tform_all, INSPR4pi.sobj_phi0, INSPR4pi.setup);

if recon_stop == 1
    msgbox('Stop by user control!');
    return;
end

INSPR4pi.srobj = srobj; 

% set export
msgbox('Finish 4Pi-SMSN reconstruction!');
set(handles.recon_export_pushbutton,'Enable','on');
set(handles.recon_export_csv_pushbutton,'Enable','on');

% make dispaly buttons enable
set(handles.display_show_pushbutton,'Enable','on');
set(handles.display_export_pushbutton,'Enable','off');


% --- Executes on button press in display_set_pushbutton.
function display_set_pushbutton_Callback(~, eventdata, handles)
% hObject    handle to display_set_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

% Draw DC figure
f_display = figure('Position',[600 200 350 260],'Name','Set display parameters');

% static text
uicontrol(f_display, 'Style','text','String','SR image high bound','FontSize',10,'Position',[50 150 150 50],...
    'HorizontalAlignment','left');
% uicontrol(f_display, 'Style','text','String','Image width (pixels)','FontSize',10,'Position',[50 150 150 50],...
%     'HorizontalAlignment','left');
uicontrol(f_display, 'Style','text','String','Zoom in','FontSize',10,'Position',[50 100 100 50],...
    'HorizontalAlignment','left');

% edit 
h_contrast  = uicontrol(f_display, 'Style','edit','String',num2str(INSPR4pi.display.Recon_color_highb),...
    'FontSize',10,'Position',[200 180 100 25]);

h_zm = uicontrol(f_display, 'Style','edit','String',num2str(INSPR4pi.display.zm),...
    'FontSize',10,'Position',[200 130 100 25]);

h_display.h_contrast = h_contrast;
h_display.h_zm = h_zm;

% save button
uicontrol(f_display,'Style','pushbutton','String','save','Position',[180 50 100 30],'FontSize',10,'Callback',...
    {@display_set_save_pushbutton_Callback,h_display});

%save display parameters
function display_set_save_pushbutton_Callback(src,event,t)
global INSPR4pi 

Recon_color_highb = str2num( get(t.h_contrast,'String') ); 
zm = str2num( get(t.h_zm,'String') ); 

INSPR4pi.display.Recon_color_highb = Recon_color_highb;
INSPR4pi.display.zm = zm;


% --- Executes on button press in display_show_pushbutton.
function display_show_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to display_show_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi 

if ~isfield(INSPR4pi,'srobj')
    msgbox('No reconstruction results!');
    return;
end

display('Show 3D reconstrcted Z-color image...');

obj = INSPR4pi.srobj;
sz = INSPR4pi.display.imagesz;
obj.zm = INSPR4pi.display.zm;
obj.Recon_color_highb = INSPR4pi.display.Recon_color_highb;
segnum = 64;
flagstr = [];
if isempty(obj.loc_x_f)||isempty(obj.loc_y_f)||isempty(obj.loc_z_f)
    warning('loc_x(y or z)_f property is empty, reconstructing from raw localization data');
    reconx=obj.loc_x;
    recony=obj.loc_y;
    reconz=obj.loc_z;
    flagstr='raw';
else
    reconx=obj.loc_x_f(:)./obj.Cam_pixelsz;
    recony=obj.loc_y_f(:)./obj.Cam_pixelsz;
    reconz=obj.loc_z_f(:);
    flagstr='dc';
end

if isempty(reconx)||isempty(recony)||isempty(reconz)
    error('Empty matrix detected. Reconstruction will not proceed.');
end

[rch,gch,bch]=srhist_color(sz,obj.zm,reconx,recony,reconz,segnum);
% save colored reconstruction
rchsm = imgaussfilt(rch,1);
gchsm = imgaussfilt(gch,1);
bchsm = imgaussfilt(bch,1);
rchsmst=imstretch_linear(rchsm,0,obj.Recon_color_highb,0,255);
gchsmst=imstretch_linear(gchsm,0,obj.Recon_color_highb,0,255);
bchsmst=imstretch_linear(bchsm,0,obj.Recon_color_highb,0,255);
% colorim=joinchannels('RGB',rchsmst,gchsmst,bchsmst);
colorim = cat(3,rchsmst,gchsmst,bchsmst);
colorim = uint8(colorim);
% dipshow(colorim);

figure; imshow(colorim,[]);
axis tight

INSPR4pi.display.colorim = colorim;
% save backup
imwrite(colorim,fullfile(INSPR4pi.setup.workspace,['SR_Zcol_' flagstr '_backup.tif']));

%enable export button
set(handles.display_export_pushbutton,'Enable','on');



% --- Executes on button press in recon_export_csv_pushbutton.
function recon_export_csv_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_export_csv_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

if ~isfield(INSPR4pi,'srobj')
    msgbox('No reconstruction results!');
    return;
end

obj = INSPR4pi.srobj;
[filename, pathname] = uiputfile(fullfile(INSPR4pi.setup.workspace,'particles.csv'), 'Export to CSV format...');
if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   disp(['User selected ',fullfile(pathname,filename)])
   
   if isempty(obj.loc_x_f)||isempty(obj.loc_y_f)||isempty(obj.loc_z_f)
       warning('loc_x(y or z)_f property is empty, exporting from raw localization data');
       reconx=obj.loc_x;
       recony=obj.loc_y;
       reconz=obj.loc_z;
   else
       reconx=obj.loc_x_f./obj.Cam_pixelsz;
       recony=obj.loc_y_f./obj.Cam_pixelsz;
       reconz=obj.loc_z_f;
   end
   
   [flag]=export2csv(pathname,filename,1,{reconx.*obj.Cam_pixelsz},{recony.*obj.Cam_pixelsz},{reconz},{obj.loc_t});
      
   if flag==1
       msgbox('Export Successfully!');
   else
       msgbox('Export with error!');
   end

end




function pupil_iter_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_iter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_iter_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_iter_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_iter_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_iter_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pupil_min_similarity_edit_Callback(hObject, eventdata, ~)
% hObject    handle to pupil_min_similarity_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_min_similarity_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_min_similarity_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_min_similarity_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_min_similarity_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pupil_zernike_order_popupmenu.
function pupil_zernike_order_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_zernike_order_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pupil_zernike_order_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pupil_zernike_order_popupmenu
global INSPR4pi

contents = get(handles.pupil_zernike_order_popupmenu,'String'); 

switch contents{get(handles.pupil_zernike_order_popupmenu,'Value')}
    case '25'
        INSPR4pi.pupil.ZernikeorderN = 4;   
    case '36'
        INSPR4pi.pupil.ZernikeorderN = 5;
    case '49'
        INSPR4pi.pupil.ZernikeorderN = 6;
    case '64'
        INSPR4pi.pupil.ZernikeorderN = 7;
    otherwise
        INSPR4pi.pupil.ZernikeorderN = 7;
end

display(['shift mode is: ' num2str(INSPR4pi.pupil.ZernikeorderN)]);


% --- Executes during object creation, after setting all properties.
function pupil_zernike_order_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_zernike_order_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pupil_change_pushbutton.
function pupil_change_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_change_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

f_init = figure('Position',[600 400 1000 150],'Name','Initial Zernike value');

cnames = {'Ast','DAst','Coma x','Coma y','1st Sph','Trefoil x','Trefoil y','2nd Ast','2nd DAst',...
    '2nd Coma x','2nd Coma y','2nd Sph','Tetrafoil x','Tetrafoil y','2nd Trefoil x','2nd Trefoil y',...
    '3rd Ast','3rd DAst','3rd Coma x','3rd Coma y','3rd Sph'};
rnames = {'Bottom obj.', 'Top obj.'};
label_bot = INSPR4pi.pupil.init_z_bot;
label_top = INSPR4pi.pupil.init_z_top;
dat_zernike = cat(1,label_bot,label_top); 

t_zernike = uitable(f_init,'Data',dat_zernike,'ColumnName',cnames,'RowName',rnames,'Position',[20 20 800 100]);
t_zernike.ColumnEditable = true;

% INSPR4pi.pupil.init_z

uicontrol(f_init,'Style','pushbutton','String','save','Position',[860 50 100 50],'Callback',{@pupil_zernike_save_pushbutton_Callback,t_zernike});

%save initial zernike value
function pupil_zernike_save_pushbutton_Callback(src,event,t)
global INSPR4pi 

dat_zernike = get(t,'Data');
INSPR4pi.pupil.init_z_bot = dat_zernike(1,:);
INSPR4pi.pupil.init_z_top = dat_zernike(2,:);

display(['initial zernike of bottom objective is: ' num2str(INSPR4pi.pupil.init_z_bot)]);
display(['initial zernike of top objective is: ' num2str(INSPR4pi.pupil.init_z_top)]);


function pupil_blur_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_blur_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pupil_blur_edit as text
%        str2double(get(hObject,'String')) returns contents of pupil_blur_edit as a double


% --- Executes during object creation, after setting all properties.
function pupil_blur_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pupil_blur_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pupil_process_pushbutton.
function pupil_process_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_process_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath('.\4Pi_insitu_PSF_retrieval\');
addpath('.\4Pi_insitu_PSF_retrieval\PR4Pi_4channels\');

global INSPR4pi
global pupil_stop

%check data

if ~isfield(INSPR4pi,'subregion_chs')
    msgbox('Please improt subregion images!');
    return
end

pupil_stop = 0; %initial stop == 0


%import parameters from GUI
display('Import parameters...');

Zrange_low = str2num(get(handles.pupil_Zrange_low_edit,'String'));
Zrange_mid = str2num(get(handles.pupil_Zrange_mid_edit,'String'));
Zrange_high = str2num(get(handles.pupil_Zrange_high_edit,'String'));
bin_lowerBound = str2num(get(handles.pupil_bin_low_edit,'String'));
iter = str2num(get(handles.pupil_iter_edit,'String'));
min_similarity = str2num(get(handles.pupil_min_similarity_edit,'String'));
% blur_sigma = str2num(get(handles.pupil_blur_edit,'String'));

Z_pos = [Zrange_low:Zrange_mid:Zrange_high];

%import to INSPR4pi
INSPR4pi.pupil.Zrange_low = Zrange_low;
INSPR4pi.pupil.Zrange_mid = Zrange_mid;
INSPR4pi.pupil.Zrange_high = Zrange_high;
INSPR4pi.pupil.Z_pos = Z_pos;
INSPR4pi.pupil.bin_lowerBound = bin_lowerBound; %the lower bound of each group images
INSPR4pi.pupil.iter = iter;
INSPR4pi.pupil.min_similarity = min_similarity; %min similarity in NCC calculation


%EMpupil process... unfinished
display('4Pi in situ data process...');

probj = INSPR4Pi_model_generation(INSPR4pi.subregion_chs,INSPR4pi.setup,INSPR4pi.pupil); %this function is initial, need change


if pupil_stop == 1
    msgbox('Stop by user control!'); 
    return;
end

%save probj (pupil)
INSPR4pi.probj = probj;


%enable button
msgbox('Finish pupil estimation!');

set(handles.pupil_export_pushbutton,'Enable','on');
set(handles.pupil_showPSFs_pushbutton,'Enable','on');
set(handles.pupil_show_pupil_pushbutton,'Enable','on');
set(handles.pupil_show_zernike_pushbutton,'Enable','on');
set(handles.pupil_showPSFs_fixed_pushbutton,'Enable','on');


% --- Executes on button press in pupil_export_pushbutton.
function pupil_export_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_export_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi

[filename, pathname] = uiputfile(fullfile(INSPR4pi.setup.workspace,'*.mat'), 'save pupil...');
if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   disp(['User selected ',fullfile(pathname,filename)])

   probj = INSPR4pi.probj;
   save(fullfile(pathname,filename),'probj');
end


% --- Executes on button press in pupil_showPSFs_pushbutton.
function pupil_showPSFs_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_showPSFs_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

disp('Show the Reassembled and retrieved 4Pi-PSFs ...')
% INSPR4pi.probj.genPRfigs('PSF');
% genPupilfigs(INSPR4pi.probj, 'PSF',INSPR4pi.setup.workspace);
genPupilfigs_4ch(INSPR4pi.probj, 'PSF',INSPR4pi.setup.workspace);

% --- Executes on button press in pupil_show_pupil_pushbutton.
function pupil_show_pupil_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_show_pupil_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

disp('Show PR pupil and Zernike pupil ...')
% INSPR4pi.probj.genPRfigs('pupil');
genPupilfigs_4ch(INSPR4pi.probj, 'pupil',INSPR4pi.setup.workspace);

% --- Executes on button press in pupil_show_zernike_pushbutton.
function pupil_show_zernike_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_show_zernike_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi

disp('Show Zernike value in magnitude and phase ...')
% INSPR4pi.probj.genPRfigs('zernike');
% genPupilfigs(INSPR4pi.probj, 'zernike',INSPR4pi.setup.workspace);
genPupilfigs_4ch(INSPR4pi.probj, 'zernike',INSPR4pi.setup.workspace);


% --- Executes on button press in seg_export_pushbutton.
function seg_export_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to seg_export_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi

[filename, pathname] = uiputfile(fullfile(INSPR4pi.setup.workspace,'*.mat'), 'save 4Pi subregion images');
if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   disp(['User selected ',fullfile(pathname,filename)])

   subregion_chs = INSPR4pi.subregion_chs;

   save(fullfile(pathname,filename),'subregion_chs');
end



% --- Executes on button press in pupil_zernike_checkbox.
function pupil_zernike_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_zernike_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pupil_zernike_checkbox
global INSPR4pi

if get(handles.pupil_zernike_checkbox,'Value')==1
    set(handles.pupil_change_pushbutton,'Enable','on');
else
    INSPR4pi.pupil.init_z_bot = zeros(1,21);
    INSPR4pi.pupil.init_z_bot(1) = 1.2;
    INSPR4pi.pupil.init_z_top = zeros(1,21);
    INSPR4pi.pupil.init_z_top(1) = -1.2;

    set(handles.pupil_change_pushbutton,'Enable','off');
end


% --- Executes on button press in registration_process_pushbutton.
function registration_process_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to registration_process_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

if isfield(INSPR4pi,'qd1') && isfield(INSPR4pi,'qd2') && isfield(INSPR4pi,'qd3') && isfield(INSPR4pi,'qd4') && isfield(INSPR4pi,'tform_all')
    display('Registrate 4Pi data...');
        
    data1 = INSPR4pi.qd1;
    for ii = 2 : 4
        data_tmp = eval(['INSPR4pi.qd', num2str(ii)]);    
        rate = mean(mean(mean(data1))) / mean(mean(mean(data_tmp)));
        data_tmp = data_tmp * rate;
        
        data_registered = imwarp(data_tmp,INSPR4pi.tform_all{ii},'cubic','OutputView',imref2d(size(data1)));
        eval(['data', num2str(ii), '= data_registered;']);
    end
    
    INSPR4pi.qd2_registered = data2;
    INSPR4pi.qd3_registered = data3;
    INSPR4pi.qd4_registered = data4;
    
    set(handles.registration_merge_pushbutton,'Enable','on');
    
    %Enable
    set(handles.seg_process_pushbutton,'Enable','on');
    msgbox('Finish 4Pi registration!');
else
    msgbox('Please check the data or transformation model');
end


set(handles.registration_merge_pushbutton,'Enable','on');

% --- Executes on button press in registration_merge_pushbutton.
function registration_merge_pushbutton_Callback(hObject, ~, handles)
% hObject    handle to registration_merge_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

im1 = max(INSPR4pi.qd1,[],3);
imref = im1;
imouts = im1;
for ii = 2 : 4
    im2 = eval(['INSPR4pi.qd', num2str(ii), '_registered']);
    tmp_img = max(im2,[],3);
    
    imref = cat(2,imref,im1);
    imouts = cat(2,imouts,tmp_img);
end

figure; imshowpair(imref, imouts,'ColorChannels',[1 2 0], 'Scaling','joint');
axis tight
title('Registrated images in 4 channels compared to P1 channel');



% --- Executes on button press in recon_export_pushbutton.
function recon_export_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_export_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi

[filename, pathname] = uiputfile(fullfile(INSPR4pi.setup.workspace,'*.mat'), 'save SMLM reconstruction result...');
if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   disp(['User selected ',fullfile(pathname,filename)])

   srobj = INSPR4pi.srobj;
   save(fullfile(pathname,filename),'srobj');
   msgbox('Finish export!');
end


% --- Executes on button press in pupil_stop_pushbutton.
function pupil_stop_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_stop_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global pupil_stop

pupil_stop = 1;


% --- Executes on button press in recon_stop_pushbutton.
function recon_stop_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_stop_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global recon_stop

recon_stop = 1;


% --- Executes on button press in setup_workspace_pushbutton.
function setup_workspace_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to setup_workspace_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi

folder_name = uigetdir(INSPR4pi.setup.workspace,'Select workspace folder');

if isequal(folder_name,0)
   disp('User selected Cancel')
else
   disp(['User selected ', folder_name]);
   set(handles.setup_workspace_edit,'String', folder_name);

   INSPR4pi.setup.workspace = folder_name;
end




function setup_workspace_edit_Callback(hObject, eventdata, handles)
% hObject    handle to setup_workspace_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setup_workspace_edit as text
%        str2double(get(hObject,'String')) returns contents of setup_workspace_edit as a double


% --- Executes during object creation, after setting all properties.
function setup_workspace_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setup_workspace_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in display_export_pushbutton.
function display_export_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to display_export_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

if isfield(INSPR4pi.display,'colorim')
    [filename, pathname] = uiputfile(fullfile(INSPR4pi.setup.workspace,'*.tif'), ...
        'save super resolution image...');
    if isequal(filename,0) || isequal(pathname,0)
        disp('User selected Cancel')
    else
        disp(['User selected ',fullfile(pathname,filename)])
                
        imwrite(INSPR4pi.display.colorim,fullfile(pathname,filename),'TIFF');

        msgbox('Finish export!');
    end

else
    msgbox('No Z projection super-resolution images!');
end



% % --- Executes on selection change in pupil_xyshift_mode_popupmenu.
% function pupil_xyshift_mode_popupmenu_Callback(hObject, eventdata, handles)
% % hObject    handle to pupil_xyshift_mode_popupmenu (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = cellstr(get(hObject,'String')) returns pupil_xyshift_mode_popupmenu contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from pupil_xyshift_mode_popupmenu
% 
% global INSPR4pi
% 
% contents = get(handles.pupil_xyshift_mode_popupmenu,'String'); 
% 
% switch contents{get(handles.pupil_xyshift_mode_popupmenu,'Value')}
%     case 'Separate shift'
%         INSPR4pi.pupil.XYshift_mode = 0;   
%     case 'Together shift'
%         INSPR4pi.pupil.XYshift_mode = 1;
%     otherwise
%         INSPR4pi.pupil.XYshift_mode = 1;
% end
% 
% display(['XY shift mode is: ' num2str(INSPR4pi.pupil.XYshift_mode)]);



% % --- Executes during object creation, after setting all properties.
% function pupil_xyshift_mode_popupmenu_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to pupil_xyshift_mode_popupmenu (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: popupmenu controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end





% --- Executes on button press in pushbutton48.
function pushbutton48_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton48 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton49.
function pushbutton49_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton49 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton50.
function pushbutton50_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function seg_frame_id_edit_Callback(hObject, eventdata, handles)
% hObject    handle to seg_frame_id_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_frame_id_edit as text
%        str2double(get(hObject,'String')) returns contents of seg_frame_id_edit as a double


% --- Executes during object creation, after setting all properties.
function seg_frame_id_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_frame_id_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in seg_show_img_pushbutton.
function seg_show_img_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to seg_show_img_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% show sub-region images 
global INSPR4pi

boxsz = str2num( get(handles.seg_boxsz_edit, 'String') );
num_display = str2num( get(handles.seg_frame_id_edit, 'String') );

if num_display < 1 || num_display > size(INSPR4pi.seg_display.ims_ch1,3)
    msgbox('Please input the correct number!');
end

% P1
display('Show segmentation results in P1 Channel');
raw = INSPR4pi.seg_display.ims_ch1(:,:,num_display)/max(max(INSPR4pi.seg_display.ims_ch1(:,:,num_display)));
index_rec = find(INSPR4pi.seg_display.allcds_mask(:,3) == num_display-1);
num_detect = length(index_rec);

if num_detect ~= 0
    rec_vector = cat(2,INSPR4pi.seg_display.t1(index_rec),INSPR4pi.seg_display.l1(index_rec),repmat(boxsz,length(index_rec),2));
    img_select = insertShape(raw,'Rectangle',rec_vector,'LineWidth',1, 'Color', 'green');
    
    figure; imshow(img_select);
    axis tight
    title('P1 channel');
else
    figure; imshow(raw);
    axis tight
    title('P1 channel');
end

% S2
display('Show segmentation results in S2 Channel');
raw = INSPR4pi.seg_display.ims_ch2(:,:,num_display)/max(max(INSPR4pi.seg_display.ims_ch2(:,:,num_display)));
index_rec = find(INSPR4pi.seg_display.allcds_mask(:,3) == num_display-1);
num_detect = length(index_rec);

if num_detect ~= 0 
    rec_vector = cat(2,INSPR4pi.seg_display.t1(index_rec),INSPR4pi.seg_display.l1(index_rec),repmat(boxsz,length(index_rec),2));
    img_select = insertShape(raw,'Rectangle',rec_vector,'LineWidth',1, 'Color', 'green');
    
    figure; imshow(img_select);
    axis tight
    title('S2 Channel');
else
    figure; imshow(raw);
    axis tight
    title('S2 Channel');
    msgbox('Did not detect sub-region in this frame!');
end

% P2
display('Show segmentation results in P2 Channel');
raw = INSPR4pi.seg_display.ims_ch3(:,:,num_display)/max(max(INSPR4pi.seg_display.ims_ch3(:,:,num_display)));
index_rec = find(INSPR4pi.seg_display.allcds_mask(:,3) == num_display-1);
num_detect = length(index_rec);

if num_detect ~= 0 
    rec_vector = cat(2,INSPR4pi.seg_display.t1(index_rec),INSPR4pi.seg_display.l1(index_rec),repmat(boxsz,length(index_rec),2));
    img_select = insertShape(raw,'Rectangle',rec_vector,'LineWidth',1, 'Color', 'green');
    
    figure; imshow(img_select);
    axis tight
    title('P2 Channel');
else
    figure; imshow(raw);
    axis tight
    title('P2 Channel');
    msgbox('Did not detect sub-region in this frame!');
end

% S1
display('Show segmentation results in S1 Channel');
raw = INSPR4pi.seg_display.ims_ch4(:,:,num_display)/max(max(INSPR4pi.seg_display.ims_ch4(:,:,num_display)));
index_rec = find(INSPR4pi.seg_display.allcds_mask(:,3) == num_display-1);
num_detect = length(index_rec);

if num_detect ~= 0 
    rec_vector = cat(2,INSPR4pi.seg_display.t1(index_rec),INSPR4pi.seg_display.l1(index_rec),repmat(boxsz,length(index_rec),2));
    img_select = insertShape(raw,'Rectangle',rec_vector,'LineWidth',1, 'Color', 'green');
    
    figure; imshow(img_select);
    axis tight
    title('S1 Channel');
else
    figure; imshow(raw);
    axis tight
    title('S1 Channel');
    msgbox('Did not detect sub-region in this frame!');
end

    


% --- Executes on button press in recon_gpu_checkbox.
function recon_gpu_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_gpu_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_gpu_checkbox

global INSPR4pi

if get(handles.recon_gpu_checkbox,'Value')==1
    INSPR4pi.recon.isGPU = 1;
else
    INSPR4pi.recon.isGPU = 0;
end


% --- Executes on button press in display_import_pushbutton.
function display_import_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to display_import_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi
[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select reconstruction result');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   % load file
   tmp = load([pathname filename]);
   
   if isfield(tmp,'srobj')        
       INSPR4pi.srobj = tmp.srobj;
       set(handles.display_show_pushbutton,'Enable','on');
       set(handles.display_export_pushbutton,'Enable','off');
       set(handles.display_import_edit,'String', filename);
       
       msgbox('Finish importing reconstruction result!');
       
   else
       msgbox('Please import the correct data!');
   end
end




function display_import_edit_Callback(hObject, eventdata, handles)
% hObject    handle to display_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of display_import_edit as text
%        str2double(get(hObject,'String')) returns contents of display_import_edit as a double



% --- Executes during object creation, after setting all properties.
function display_import_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display_import_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in recon_bg_checkbox.
function recon_bg_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_bg_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_bg_checkbox

global INSPR4pi

if get(handles.recon_bg_checkbox,'Value')==1
    INSPR4pi.recon.is_bg = 1;
else
    INSPR4pi.recon.is_bg = 0;
end


% --- Executes on button press in data_bg_checkbox.
function data_bg_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to data_bg_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of data_bg_checkbox

global INSPR4pi

if get(handles.data_bg_checkbox,'Value')==1
    INSPR4pi.setup.is_bg = 1;
    set(handles.data_show_bg_pushbutton,'Enable','on');
else
    INSPR4pi.setup.is_bg = 0;
    set(handles.data_show_bg_pushbutton,'Enable','off');
end


% --- Executes on button press in data_show_bg_pushbutton.
function data_show_bg_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to data_show_bg_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi

if isfield(INSPR4pi,'bg_img_1') && isfield(INSPR4pi,'bg_img_2') && isfield(INSPR4pi,'bg_img_3') && isfield(INSPR4pi,'bg_img_4')
    
    disp('Show background in 4 channels...')

    imouts = cat(2,INSPR4pi.bg_img_1,INSPR4pi.bg_img_2,INSPR4pi.bg_img_3,INSPR4pi.bg_img_4);
    figure; imshow(imouts,[]);
    axis tight
    title('Background in 4 channels (P1, S2, P2, S1)');

else
    msgbox('Please import the data again!');
end




% --- Executes on button press in pupil_showPSFs_fixed_pushbutton.
function pupil_showPSFs_fixed_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_showPSFs_fixed_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

disp('Show retrieved 4Pi-PSFs at fixed axial positions')
% genPupilfigs(INSPR4pi.probj, 'PSF',INSPR4pi.setup.workspace);

% INSPR PSFs
z = [-0.8:0.1:0.8];
sz_z = size(z,2);
zind=[1:2:sz_z];
R = 128;


% Set parameters of pupil funcion
oprobj = INSPR4pi.probj;
PRstruct = oprobj.PRstruct;
bxsz = 32;

pxsz = oprobj.Pixelsize;
psfobj = PSF_4pi(PRstruct);  
psfobj.Boxsize = bxsz;
psfobj.Pixelsize = pxsz; 
psfobj.PSFsize = R;
psfobj.nMed = oprobj.PRstruct.RefractiveIndex;
psfobj.Phasediff = oprobj.Phasediff;                         
psfobj.Iratio = oprobj.Iratio;                             
psfobj.ModulationDepth = oprobj.ModulationDepth;                      
psfobj.Phi0 = oprobj.Phi0;                              
psfobj.Zoffset = oprobj.Zoffset;

% add X,Y,Z position
psfobj.Xpos = zeros(1,sz_z);                   
psfobj.Ypos = zeros(1,sz_z);
psfobj.Zpos = z;

% set pupil function
PRstruct1 = oprobj.PRstruct1;
PRstruct2 = oprobj.PRstruct2;

psfobj.precomputeParam();
psfobj.set2Pupil(PRstruct1,PRstruct2); 
% psfobj.setUnifMag();
psfobj.genPupil_4pi_2();     
psfobj.genPSF_4pi_md()

label = {'p1','s2','p2','s1'};                     
psf_4pi = [];                                     
for nn = 1:4
    psfobj.PSFs = psfobj.PSF4pi.(label{nn});
    psfobj.scalePSF();
    psfI = psfobj.ScaledPSFs;
    psfI = psfI *4000 ;
    psf_4pi = cat(4,psf_4pi,psfI);
end

% show PSFs in x-z view
N = length(zind);
qN = 4;
% ind = round(linspace(1,Num,N));
Ri = min([24,psfobj.Boxsize]);
h1 = figure('position',[100,100,100*N,102*4],'Name','4Pi-brainspot PSFs at fixed z positions');
for ii = 1:qN
    psfr = squeeze(psf_4pi(:,:,:,ii));
    Ro = size(psfr,1);
    for jj = 1:N
        ha = axes('position',[(jj-1)/N,(qN-ii)/qN,1/N,1/qN],'parent',h1);
        imagesc(psfr(Ro/2-Ri/2+1:Ro/2+Ri/2,Ro/2-Ri/2+1:Ro/2+Ri/2,zind(jj)));
        axis equal;axis off;
        if ii == 1
        text(2,3, ['z=',num2str(z(zind(jj)),3),'\mum'],'color',[1,1,1],'fontsize',12);
        end
        if jj == 1
        text(2,Ri-3, label{ii},'color',[1,1,1],'fontsize',12);
        end
    end
end
colormap(gray)






% --- Executes on button press in cavity_misobj_checkbox.
function cavity_misobj_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_misobj_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cavity_misobj_checkbox
global INSPR4pi

if get(handles.cavity_misobj_checkbox,'Value')==1
    INSPR4pi.cavity.isObj = 1;
else
    INSPR4pi.cavity.isObj = 0;
end



% --- Executes on button press in pushbutton73.
function pushbutton73_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton73 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cavity_stop_pushbutton.
function cavity_stop_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_stop_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cavity_stop

cavity_stop = 1;

% --- Executes on button press in cavity_export_pushbutton.
function cavity_export_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_export_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi

[filename, pathname] = uiputfile(fullfile(INSPR4pi.setup.workspace,'*.mat'), 'save cavity phase estimation...');
if isequal(filename,0) || isequal(pathname,0)
   disp('User selected Cancel')
else
   disp(['User selected ',fullfile(pathname,filename)])

   sobj_phi0 = INSPR4pi.sobj_phi0;
   save(fullfile(pathname,filename),'sobj_phi0');
   msgbox('Finish export!');
end



% --- Executes on button press in cavity_process_pushbutton.
function cavity_process_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_process_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% addpath('.\CavityPhase_and_objMisalign\');
addpath('.\Dynamic_model_update\');
addpath('..\Support\PSF Toolbox_4pi\');

global INSPR4pi
global cavity_stop

cavity_stop = 0; %initil stop label is none

% check the data
if INSPR4pi.cavity.isNewdata == 0
    if ~isfield(INSPR4pi,'qd1') || ~isfield(INSPR4pi,'qd2') || ~isfield(INSPR4pi,'qd3') || ~isfield(INSPR4pi,'qd4')
        msgbox('Please import the data!');
        return;
    end
    INSPR4pi.cavity.qd1 = INSPR4pi.qd1;
    INSPR4pi.cavity.qd2 = INSPR4pi.qd2;
    INSPR4pi.cavity.qd3 = INSPR4pi.qd3;
    INSPR4pi.cavity.qd4 = INSPR4pi.qd4;
    INSPR4pi.cavity.dirN = 1;
    INSPR4pi.cavity.datapath = pwd;
else
    if ~isfield(INSPR4pi.cavity,'datapath')
        msgbox('Please import the data!');
        return;
    end
end

if ~isfield(INSPR4pi,'tform_all')
    msgbox('Please import transformation model!');
    return;
end

if INSPR4pi.cavity.isNewpupil == 0
    if ~isfield(INSPR4pi,'probj')
        msgbox('Please import pupil model!');
        return;
    end
    INSPR4pi.cavity.probj_all{1} = INSPR4pi.probj;
    INSPR4pi.cavity.loopn = 1;
else
    if ~isfield(INSPR4pi.cavity,'probj_all')
        msgbox('Please import the pupil model!');
        return;
    end
end

if INSPR4pi.cavity.isSeg == 1
    INSPR4pi.cavity.seg_thresh_low = str2num( get(handles.cavity_seg_thresh_low_edit, 'String') );
    INSPR4pi.cavity.seg_thresh_high = str2num( get(handles.cavity_seg_thresh_high_edit, 'String') );    
end


    
% cavity phase and objective misalignment estimation
if INSPR4pi.cavity.isObj == 0   %only consider cavity phase
    
    sobj_phi0 = analysis_cavityPhase(INSPR4pi.cavity, INSPR4pi.tform_all, INSPR4pi.setup);
else    % consider cavity phase + objective misalignment   
    
    sobj_phi0 = analysis_cavityPhase_and_misObj_v2(INSPR4pi.cavity, INSPR4pi.tform_all, INSPR4pi.setup); 
end


if cavity_stop == 1
    msgbox('Stop by user control!');
    return;
end

INSPR4pi.sobj_phi0 = sobj_phi0; 

% set export
msgbox('Finish cavity phase estimation!');
set(handles.cavity_export_pushbutton,'Enable','on');
set(handles.cavity_show_phase_pushbutton,'Enable','on');





function cavity_seg_thresh_high_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_seg_thresh_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cavity_seg_thresh_high_edit as text
%        str2double(get(hObject,'String')) returns contents of cavity_seg_thresh_high_edit as a double


% --- Executes during object creation, after setting all properties.
function cavity_seg_thresh_high_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cavity_seg_thresh_high_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cavity_seg_thresh_low_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_seg_thresh_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cavity_seg_thresh_low_edit as text
%        str2double(get(hObject,'String')) returns contents of cavity_seg_thresh_low_edit as a double


% --- Executes during object creation, after setting all properties.
function cavity_seg_thresh_low_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cavity_seg_thresh_low_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cavity_seg_checkbox.
function cavity_seg_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_seg_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cavity_seg_checkbox

global INSPR4pi

if get(handles.cavity_seg_checkbox,'Value')==1
    set(handles.cavity_seg_thresh_low_edit,'Enable','on');
    set(handles.cavity_seg_thresh_high_edit,'Enable','on');
    
    INSPR4pi.cavity.isSeg = 1;
else
    set(handles.cavity_seg_thresh_low_edit,'String',num2str(INSPR4pi.cavity.seg_thresh_low));
    set(handles.cavity_seg_thresh_high_edit,'String',num2str(INSPR4pi.cavity.seg_thresh_high));
    set(handles.cavity_seg_thresh_low_edit,'Enable','off');
    set(handles.cavity_seg_thresh_high_edit,'Enable','off');
    
    INSPR4pi.cavity.isSeg = 0;
end


function cavity_import_pupil_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_import_pupil_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cavity_import_pupil_edit as text
%        str2double(get(hObject,'String')) returns contents of cavity_import_pupil_edit as a double


% --- Executes during object creation, after setting all properties.
function cavity_import_pupil_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cavity_import_pupil_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cavity_import_pupil_pushbutton.
function cavity_import_pupil_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_import_pupil_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi
[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select the pupil model (you can choose multiple pupil)...','MultiSelect','on');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.cavity_import_pupil_edit,'String', pathname);

   % load one file or muti-files
   probj_all = [];
   loopn = [];
   if iscell(filename)
       loopn = numel(filename);
       for ii = 1 : loopn
           input_file = [pathname filename{ii}];
           tmp = load(input_file);
           if isfield(tmp,'probj')
               probj_all{ii} = tmp.probj;
           else
               msgbox('Please import the correct data!');
               return;
           end
       end
       msgbox({'Finish importing the pupil model!' ['You choose ' num2str(loopn) ' pupil model']});
   else
       loopn = 1;
       input_file = [pathname filename];
       tmp = load(input_file);
       if isfield(tmp,'probj')                 
           probj_all{1} = tmp.probj;
       else
           msgbox('Please import the correct data!');
           return;
       end
       msgbox({'Finish importing the pupil model!' ['You choose ' num2str(loopn) ' pupil model']});
   end
   
   %save all pupil models to INSPR4pi
   INSPR4pi.cavity.probj_all = probj_all;
   INSPR4pi.cavity.loopn = loopn;
   
end



% --- Executes on button press in cavity_pupil_checkbox.
function cavity_pupil_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_pupil_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cavity_pupil_checkbox
global INSPR4pi

if get(handles.cavity_pupil_checkbox,'Value')==1
    INSPR4pi.cavity.isNewpupil = 1;
    set(handles.cavity_import_pupil_pushbutton,'Enable','on');
else
    INSPR4pi.cavity.isNewpupil = 0;

    set(handles.cavity_import_pupil_pushbutton,'Enable','off');
end


function cavity_import_tform_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_import_tform_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cavity_import_tform_edit as text
%        str2double(get(hObject,'String')) returns contents of cavity_import_tform_edit as a double


% --- Executes during object creation, after setting all properties.
function cavity_import_tform_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cavity_import_tform_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cavity_import_tform_pushbutton.
function cavity_import_tform_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_import_tform_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi
[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select the transformation model');

if isequal(filename,0)
    disp('User selected Cancel')
else
    disp(['User selected ', fullfile(pathname, filename)])
    
    set(handles.cavity_import_tform_edit,'String',filename);
    
    % load file
    tmp = load([pathname filename]);
    
    if isfield(tmp,'tform_all')
        INSPR4pi.tform_all = tmp.tform_all;
        
        msgbox('Finish importing the transformation model!');
    else
        msgbox('Please import the correct data!');
    end
    
end




% --- Executes on button press in cavity_tform_checkbox.
function cavity_tform_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_tform_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cavity_tform_checkbox
if get(handles.cavity_tform_checkbox,'Value')==1
   set(handles.cavity_import_tform_pushbutton,'Enable','on');
else
   set(handles.cavity_import_tform_pushbutton,'Enable','off'); 
end


function cavity_import_data_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_import_data_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cavity_import_data_edit as text
%        str2double(get(hObject,'String')) returns contents of cavity_import_data_edit as a double


% --- Executes during object creation, after setting all properties.
function cavity_import_data_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cavity_import_data_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cavity_import_data_pushbutton.
function cavity_import_data_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_import_data_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi
[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select the 4Pi data (you can choose multiple data)...','MultiSelect','on');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.cavity_import_data_edit,'String', pathname);
   
   % load one file or muti-files
   if iscell(filename)
       dirN = numel(filename);
       INSPR4pi.cavity.datafile_name = filename;
   else
       dirN = 1;
       INSPR4pi.cavity.datafile_name{1} = filename; 
   end
   
   %save file path
   INSPR4pi.cavity.dirN = dirN;
   INSPR4pi.cavity.datapath = pathname;
   
   msgbox({ 'Finish importing the 4Pi data path!' ['You choose ' num2str(dirN) ' data'] });
   
end


% --- Executes on button press in cavity_data_checkbox.
function cavity_data_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_data_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cavity_data_checkbox
global INSPR4pi

if get(handles.cavity_data_checkbox,'Value')==1
    INSPR4pi.cavity.isNewdata = 1;
    set(handles.cavity_import_data_pushbutton,'Enable','on');
else
    INSPR4pi.cavity.isNewdata = 0; 
    set(handles.cavity_import_data_pushbutton,'Enable','off');
end


% --- Executes on button press in pupil_interference_checkbox.
function pupil_interference_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_interference_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pupil_interference_checkbox
global INSPR4pi

if get(handles.pupil_interference_checkbox,'Value')==1
    set(handles.pupil_interference_change_pushbutton,'Enable','on');
else
    INSPR4pi.pupil.Phasediff = -1.58;  
    INSPR4pi.pupil.Iratio = 0.92;       
    INSPR4pi.pupil.Phi0 = 0;           
    INSPR4pi.pupil.ModulationDepth = 0.8;
    INSPR4pi.pupil.min_photon = 1000;
    INSPR4pi.pupil.blur_sigma = 2;
    
    set(handles.pupil_interference_change_pushbutton,'Enable','off');
end



% --- Executes on button press in pupil_interference_change_pushbutton.
function pupil_interference_change_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to pupil_interference_change_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global INSPR4pi 

    
% Draw interference figure
f_interference = figure('Position',[600 200 400 470],'Name','Interference parameters');

% static text
uicontrol(f_interference, 'Style','text','String','Phase difference (rad)','FontSize',10,'Position',[50 350 200 50],...
    'HorizontalAlignment','left','TooltipString','phase difference between s- and p-polarizations');
uicontrol(f_interference, 'Style','text','String','Transmission ratio','FontSize',10,'Position',[50 300 200 50],...
    'HorizontalAlignment','left','TooltipString','transmission ratio between bottom and top emission path');
uicontrol(f_interference, 'Style','text','String','Cavity phase (rad)','FontSize',10,'Position',[50 250 200 50],...
    'HorizontalAlignment','left','TooltipString','cavity phase between the two paths when the single emitter is in the common focus of the two objectives');

uicontrol(f_interference, 'Style','text','String','Coherence strength','FontSize',10,'Position',[50 200 200 50],...
    'HorizontalAlignment','left','TooltipString','coherence strength between interferometric PSF and the conventional PSF');
uicontrol(f_interference, 'Style','text','String','Min photon','FontSize',10,'Position',[50 150 200 50],...
    'HorizontalAlignment','left','TooltipString','Min photon of PSF in libray. Reject the molecule with low photon count');
uicontrol(f_interference, 'Style','text','String','Blur parameter','FontSize',10,'Position',[50 100 200 50],...
    'HorizontalAlignment','left','TooltipString','Blur parameter in Fourier space. Smaller value, larger blur');



% edit 
h_Phasediff = uicontrol(f_interference, 'Style','edit','String',num2str(INSPR4pi.pupil.Phasediff),...
    'FontSize',10,'Position',[220 380 110 25]);
h_Iratio = uicontrol(f_interference, 'Style','edit','String',num2str(INSPR4pi.pupil.Iratio),...
    'FontSize',10,'Position',[220 330 110 25]);
h_Phi0 = uicontrol(f_interference, 'Style','edit','String',num2str(INSPR4pi.pupil.Phi0),...
    'FontSize',10,'Position',[220 280 110 25]);

h_ModulationDepth  = uicontrol(f_interference, 'Style','edit','String',num2str(INSPR4pi.pupil.ModulationDepth),...
    'FontSize',10,'Position',[220 230 110 25]);
h_min_photon = uicontrol(f_interference, 'Style','edit','String',num2str(INSPR4pi.pupil.min_photon),...
    'FontSize',10,'Position',[220 180 110 25]);
h_blur_sigma = uicontrol(f_interference, 'Style','edit','String',num2str(INSPR4pi.pupil.blur_sigma),...
    'FontSize',10,'Position',[220 130 110 25]);


h_interference.h_Phasediff = h_Phasediff;
h_interference.h_Iratio = h_Iratio;
h_interference.h_Phi0 = h_Phi0;
h_interference.h_ModulationDepth = h_ModulationDepth;
h_interference.h_min_photon = h_min_photon;
h_interference.h_blur_sigma = h_blur_sigma;

% save button
uicontrol(f_interference,'Style','pushbutton','String','Save','Position',[220 50 100 30],'FontSize',10,'Callback',...
    {@pupil_interference_save_pushbutton_Callback,h_interference});
% Reset button
uicontrol(f_interference,'Style','pushbutton','String','Reset','Position',[50 50 100 30],'FontSize',10,'Callback',...
    {@pupil_interference_reset_pushbutton_Callback,h_interference});

%Reset parameters
function pupil_interference_reset_pushbutton_Callback(src,event,t)
global INSPR4pi

default_cfg = load('default_cfg.mat');

%update INSPR4pi
INSPR4pi.pupil.Phasediff = default_cfg.pupil.Phasediff;
INSPR4pi.pupil.Iratio = default_cfg.pupil.Iratio;
INSPR4pi.pupil.Phi0 = default_cfg.pupil.Phi0;
INSPR4pi.pupil.ModulationDepth = default_cfg.pupil.ModulationDepth;
INSPR4pi.pupil.min_photon = default_cfg.pupil.min_photon;
INSPR4pi.pupil.blur_sigma = default_cfg.pupil.blur_sigma;

%update UI
set(t.h_Phasediff,'String', num2str(INSPR4pi.pupil.Phasediff));
set(t.h_Iratio,'String', num2str(INSPR4pi.pupil.Iratio));
set(t.h_Phi0,'String', num2str(INSPR4pi.pupil.Phi0));
set(t.h_ModulationDepth,'String', num2str(INSPR4pi.pupil.ModulationDepth));
set(t.h_min_photon,'String', num2str(INSPR4pi.pupil.min_photon));
set(t.h_blur_sigma,'String', num2str(INSPR4pi.pupil.blur_sigma));

%save interference parameters
function pupil_interference_save_pushbutton_Callback(src,event,t)
global INSPR4pi 

Phasediff = str2num( get(t.h_Phasediff,'String') ); 
Iratio = str2num( get(t.h_Iratio,'String') ); 
Phi0 = str2num( get(t.h_Phi0,'String') ); 
ModulationDepth = str2num( get(t.h_ModulationDepth,'String') ); 
min_photon = str2num( get(t.h_min_photon,'String') ); 
blur_sigma = str2num( get(t.h_blur_sigma,'String') ); 

INSPR4pi.pupil.Phasediff = Phasediff;
INSPR4pi.pupil.Iratio = Iratio;
INSPR4pi.pupil.Phi0 = Phi0;
INSPR4pi.pupil.ModulationDepth = ModulationDepth;
INSPR4pi.pupil.min_photon = min_photon;
INSPR4pi.pupil.blur_sigma = blur_sigma;



% --- Executes on button press in recon_phi0_checkbox.
function recon_phi0_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_phi0_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_phi0_checkbox
global INSPR4pi

if get(handles.recon_phi0_checkbox,'Value')==1
    INSPR4pi.recon.isNewphi0 = 1;
    set(handles.recon_import_phi0_pushbutton,'Enable','on');
else
    INSPR4pi.recon.isNewphi0 = 0;

    set(handles.recon_import_phi0_pushbutton,'Enable','off');
end




% --- Executes on button press in recon_import_phi0_pushbutton.
function recon_import_phi0_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to recon_import_phi0_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global INSPR4pi
[filename,pathname] = uigetfile(fullfile(INSPR4pi.setup.workspace,'*.mat'),'Select the cavity phase (estimated from Cavity phase and objective misalignment section)');

if isequal(filename,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(pathname, filename)])
   
   set(handles.recon_import_phi0_pushbutton,'String',filename);
   
   % load file
   tmp = load([pathname filename]);
   
   if isfield(tmp,'sobj_phi0')
       INSPR4pi.sobj_phi0 = tmp.sobj_phi0;
       
       msgbox('Finish importing the cavity phase!');
   else
       msgbox('Please import the correct data!');
   end
 
   
end



function recon_import_phi0_edit_Callback(hObject, eventdata, handles)
% hObject    handle to recon_import_phi0_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of recon_import_phi0_edit as text
%        str2double(get(hObject,'String')) returns contents of recon_import_phi0_edit as a double


% --- Executes during object creation, after setting all properties.
function recon_import_phi0_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to recon_import_phi0_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in recon_misobj_checkbox.
function recon_misobj_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to recon_misobj_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of recon_misobj_checkbox

global INSPR4pi

if get(handles.recon_misobj_checkbox,'Value')==1
    INSPR4pi.recon.isObj = 1;
else
    INSPR4pi.recon.isObj = 0;
end



% --- Executes on button press in cavity_show_phase_pushbutton.
function cavity_show_phase_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cavity_show_phase_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global INSPR4pi

if INSPR4pi.cavity.isObj == 0
    phi0_fit_cycle = INSPR4pi.sobj_phi0.phi0_fit_cycle;
    figure; plot(phi0_fit_cycle,'+:','LineWidth',1)
    
    ylim([-3.3,3.3])
    xlabel('number of cycles')
    ylabel('Phi0 estimation (rad)')
    set(gca,'FontSize',12)
else
    if ~isfield(INSPR4pi.sobj_phi0,'misobj_fit_cycle')
        msgbox('Please check whether click Objective misalignment estimation checkbox!');
        return;
    end
    phi0_fit_cycle = INSPR4pi.sobj_phi0.phi0_fit_cycle;
    misobj_fit_cycle = INSPR4pi.sobj_phi0.misobj_fit_cycle;
    
    figure; plot(phi0_fit_cycle,'+:','LineWidth',1)    
    ylim([-3.3,3.3])
    xlabel('number of cycles')
    ylabel('Phi0 estimation (rad)')
    set(gca,'FontSize',12)
    
    figure; plot(misobj_fit_cycle(:,1),'+:','LineWidth',1)
    hold on 
    plot(misobj_fit_cycle(:,2),'+:','LineWidth',1)
    plot(misobj_fit_cycle(:,3),'+:','LineWidth',1)
    hold off
    
    ylim([-2,2])
    legend('x tilt', 'y tilt', 'defocus');
    xlabel('number of cycles')
    ylabel('objective misalignment estimation (Lambda)')
    set(gca,'FontSize',12)
end
