%  Script for cavity phase estimation
% (C) Copyright 2022                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu, 
%       
%%
function sobj_phi0 = analysis_cavityPhase(cavity, tform_all, setup_para)

global cavity_stop  %control program stop
% input:
%   cavity: cavity phase parameters from GUI,including data and pupil model
%   tfrom_all: Affine transformation model between 4 channels in 4Pi system
%   setup: setup parameters

% output:
%   sobj_phi0: estimation for cavity phase and objective misalignment


%% setup parameters
setup.workspace = setup_para.workspace;
sobj_phi0 = [];

%% import data, channel registration and cavity phase estimation

phi0_keep_cycle = {};
phi0_fit_cycle = [];
for nn = 1:cavity.dirN
    
    drawnow
    if cavity_stop == 1
        return;
    end
    
    disp(['...Processing Data   ' int2str(nn)]);    
    if cavity.isNewdata == 1
        load([cavity.datapath, cavity.datafile_name{nn}]);
    else
        qd1 = cavity.qd1;
        qd2 = cavity.qd2;
        qd3 = cavity.qd3;
        qd4 = cavity.qd4;
    end
    
    %find subregion
    drawnow
    if cavity_stop == 1
        return;
    end
    
    disp('Image segmentation');
    boxsz = 16;
    thresh = [cavity.seg_thresh_low,cavity.seg_thresh_high];        
    [subregion_chs,frame_num,l,t] = find_4pi_subregion(qd1,qd2,qd3,qd4,tform_all,boxsz,thresh,setup_para,cavity);
   
    
    %Phi0 estimation
    drawnow
    if cavity_stop == 1
        return;
    end
    
    disp('Phi0 estimation');
    step_num = mod(nn-1,cavity.loopn) + 1;
    
    [phi0_guess_keep,phi0_fit] = caliphi0_NCC_v3(subregion_chs,cavity.probj_all{step_num});

    phi0_keep_cycle{nn} = phi0_guess_keep;
    phi0_fit_cycle = [phi0_fit_cycle; phi0_fit];
    
end

sobj_phi0 = struct('phi0_keep_cycle',phi0_keep_cycle,'phi0_fit_cycle',phi0_fit_cycle);
save(fullfile(setup.workspace,'phi0_est'), 'sobj_phi0');


%%

