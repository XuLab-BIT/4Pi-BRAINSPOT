%%
% Script for estimating in situ 3D 4Pi-PSF  
% (C) Copyright 2022                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu
%           
function probj = INSPR4Pi_model_generation(subregion_chs,setup,pupil_para)

global pupil_stop

empupil = [];   
%% setup parameters   
empupil.NA = setup.NA;
empupil.Lambda = setup.Lambda;
empupil.nMed = setup.nMed;
empupil.Pixelsize = setup.Pixelsize;   

%% interference parameters
empupil.Phasediff = pupil_para.Phasediff; 
empupil.Iratio = pupil_para.Iratio;
empupil.Phi0 = pupil_para.Phi0;
empupil.Zoffset = 0;
empupil.ModulationDepth = pupil_para.ModulationDepth;

min_photon = pupil_para.min_photon;
empupil.blur_sigma = pupil_para.blur_sigma;


%% pupil parameters
empupil.imsz = size(subregion_chs,1);
empupil.Z_pos = pupil_para.Z_pos;
empupil.bin_lowerBound = pupil_para.bin_lowerBound;
empupil.min_similarity = pupil_para.min_similarity;
empupil.iter = pupil_para.iter;
label_bot = pupil_para.init_z_bot; % initial zernike value in pupil phase at bottom objective
label_top = pupil_para.init_z_top; % initial zernike value in pupil phase at top objective

empupil.ZernikeorderN = pupil_para.ZernikeorderN; 
empupil.Zernike_sz = (empupil.ZernikeorderN+1)^2;
empupil.Zshift_mode = pupil_para.Zshift_mode; % 0 for no XYZ shift; 1 for XYZ shift and update Z shift in each iteration; 
                                              
empupil.iter_mode = 0;  




%% reject low intensity sub-regions
x1 = genIniguess( mean(subregion_chs,4),'median'); 

intmask = x1(:,3) < min_photon;
data_input = subregion_chs(:,:,~intmask,:);


%% Generate initial pupil
empupil.zshift = 0;

probj = gen_init_4PiPupil_v2(empupil,label_top,label_bot);

subregion_chs_norm = subregion_normalization_4pi(data_input);

for iter = 1 : empupil.iter
    display(['Iteration: ' num2str(iter) ' ...']);
    
    
    %% Generate reference Z-postions' PSFs from 4Pi interference patterns
    [ref_planes,empupil] = gen_aber4PiPSF_directlyfromPR(probj, empupil, label_top, label_bot);
 
    
    %% classify single molecules from different channels to giving reference Z-position images
    
    drawnow
    if pupil_stop == 1
        return;
    end
    
    [ims_Zcal_ave_planes, index_record_Zplanes] = classify_4PiPSF_4channels_seperately(subregion_chs_norm,ref_planes, empupil);
    
    
    %% Refine the pupil function and estimate the aberration using average
    drawnow
    if pupil_stop == 1
        return;
    end
    
    % 4 channels
    [probj,empupil] = PR4Pi_fromAveZ_4channels_Zshift(ims_Zcal_ave_planes,index_record_Zplanes,empupil);
    label_top = probj.PRstruct1.Zernike_phase(5:25);
    label_bot = probj.PRstruct2.Zernike_phase(5:25);

%     genPRfigs_4ch(probj, 'PSF'); 
    
end

%%
index_number = length(find(index_record_Zplanes==1));

if index_number < 5
    msgbox('Warning! The range of localization is not enough for reliable model generation!');
end





