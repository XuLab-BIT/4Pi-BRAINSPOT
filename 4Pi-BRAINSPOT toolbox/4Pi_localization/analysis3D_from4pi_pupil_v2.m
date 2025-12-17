%  Script for 4Pi-SMSN reconstruction
% (C) Copyright                     The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu
%       
%%
function srobj = analysis3D_from4pi_pupil_v2(recon, tform_all, sobj_phi0, setup_para)

global recon_stop  %control program stop
% input:
%   recon: reconstruction parameters from GUI,including data and pupil
%   model
%   tfrom_all: Affine transformation model in 4Pi system
%   sobj_phi0: cavity phase estimation
%   setup: setup parameters

% output:
%   srobj: reconstruction results


%% setup parameters

setup.pixelsize = setup_para.Pixelsize * 1000;  %nm
setup.n_imm = setup_para.RefractiveIndex;
setup.n_sample = setup_para.nMed;
setup.workspace = setup_para.workspace;

srobj = SRscmos(setup.workspace); 


%% import data, channel registration and localization

loc_x = [];
loc_y = [];
loc_z = [];
loc_t = [];
loc_photons = [];
loc_bg = [];
loc_ll = [];
loc_crlb = [];
loc_step = [];
for nn = 1:recon.dirN
    
    drawnow
    if recon_stop == 1
        return;
    end
    
    disp(['...Processing Data   ' int2str(nn)]);    
    if recon.isNewdata == 1
        load([recon.datapath, recon.datafile_name{nn}]);
    else
        qd1 = recon.qd1;
        qd2 = recon.qd2;
        qd3 = recon.qd3;
        qd4 = recon.qd4;
    end
    
    %find subregion
    drawnow
    if recon_stop == 1
        return;
    end
    
    disp('Image segmentation');
    boxsz = 16;
    thresh = [recon.seg_thresh_low,recon.seg_thresh_high];      
    [subregion_chs,subvar_chs,frame_num,l,t,offset_seg] = crop_4pi_subregion_transModel_sCMOS(qd1,qd2,qd3,qd4,tform_all,boxsz,thresh,setup_para,recon);

    
    %import cavity phase
    if recon.isObj == 1 && isfield(sobj_phi0,'misobj_keep_cycle')
        % consider cavity phase + objective misalignment
        phi0_input = sobj_phi0(1).phi0_fit_cycle(nn);
        misobj_input = sobj_phi0(1).misobj_fit_cycle(nn,:);
    else
        % only consider cavity phase
        phi0_input = sobj_phi0(1).phi0_fit_cycle(nn);  
        misobj_input = [0 0 0];
    end
    
    %localzation
    drawnow
    if recon_stop == 1
        return;
    end
    
    disp('4Pi localization');    
    step_num = mod(nn-1,recon.loopn) + 1;
    
    if recon.isGPU
        [PM,crlbM,errM,maskc] = loc_channel_specific_4Pimodel(subregion_chs,...
            recon.probj_all{step_num},phi0_input,misobj_input,tform_all,offset_seg,subvar_chs);
    
    else 
        % message to GUI
        numcores = feature('numcores');
        msgbox('CPU version is in developing');  
    end
    
    loc_x = [loc_x; PM(:,2)+l(maskc)];
    loc_y = [loc_y; PM(:,1)+t(maskc)];
    
    offset = -(recon.probj_all{step_num}.Zpos(1)+ recon.probj_all{step_num}.Zpos(end))/2;
    loc_z = [loc_z; PM(:,3)+offset];
    loc_t = [loc_t; frame_num(maskc) + (nn-1)*size(qd1,3)];
    loc_photons = [loc_photons; 2*PM(:,4)];
    loc_bg = [loc_bg; PM(:,5)];
    loc_ll = [loc_ll; errM(:,2)];
    loc_crlb = [loc_crlb; crlbM];
    
    loc_step = [loc_step; step_num * ones(size(frame_num(maskc),1),1)];

end

sobj = struct('loc_x',loc_x,'loc_y',loc_y,'loc_z',loc_z,'loc_t',loc_t,'loc_photons',loc_photons,'loc_bg',loc_bg,'loc_ll',loc_ll,'loc_crlb',loc_crlb, 'loc_step', loc_step);
save(fullfile(setup.workspace,'loc_backup'), 'sobj');


%% pre-rejection (LL, z range)
if recon.isRej == 1
    disp('Localization pre-rejection'); 
    
    llmask = sobj.loc_ll > recon.rej.llthreshold;

    uncerxymask = sqrt(sobj.loc_crlb(:,1))*1000 > recon.rej.loc_uncerxy_max | sqrt(sobj.loc_crlb(:,2))*1000 > recon.rej.loc_uncerxy_max;
    uncerzmask = sqrt(sobj.loc_crlb(:,3))*1000 > recon.rej.loc_uncerz_max;
    intmask = sobj.loc_photons < recon.rej.min_photon;

    zmask= sobj.loc_z > recon.rej.zmask_high | sobj.loc_z < recon.rej.zmask_low;
    

    totmask = llmask | uncerzmask | intmask | uncerxymask | zmask;
        
    if recon.is_bg == 1 || (recon.isNewdata == 0 && setup_para.is_bg == 1)
        totmask = uncerzmask | uncerxymask | zmask;
    end
    
    loc_x_keep = sobj.loc_x(~totmask);
    loc_y_keep = sobj.loc_y(~totmask);
    loc_z_keep = sobj.loc_z(~totmask);
    loc_t_keep = sobj.loc_t(~totmask);
    loc_step_keep = sobj.loc_step(~totmask);
    
    loc_ll_keep = sobj.loc_ll(~totmask);
    loc_bg_keep = sobj.loc_bg(~totmask);
    loc_photons_keep = sobj.loc_photons(~totmask);
    loc_crlb_keep = sobj.loc_crlb(~totmask,:);
else 
    loc_x_keep = sobj.loc_x;
    loc_y_keep = sobj.loc_y;
    loc_z_keep = sobj.loc_z;
    loc_t_keep = sobj.loc_t;
    loc_step_keep = sobj.loc_step;
    
    loc_ll_keep = sobj.loc_ll;
    loc_bg_keep = sobj.loc_bg;
    loc_photons_keep = sobj.loc_photons;
    loc_crlb_keep = sobj.loc_crlb;
end


%% 2D or 3D alignment
pixelsize = setup.pixelsize;    %nm

srobj.loc_x = loc_x_keep;% pixel
srobj.loc_y = loc_y_keep;% pixel
srobj.loc_z = loc_z_keep.*1e3 + recon.dc.z_offset;% nm, must be positive value
srobj.loc_t = loc_t_keep;
srobj.loc_step = loc_step_keep;

srobj.loc_ll = loc_ll_keep;
srobj.loc_bg = loc_bg_keep;
srobj.loc_photons = loc_photons_keep;
srobj.loc_crlb = loc_crlb_keep;

if recon.isDC == 1
    disp('Drift correction');

    srobj.frmpfile = recon.dc.frmpfile;
    srobj.loc_cycle = floor(loc_t_keep/srobj.frmpfile);
    srobj.Cam_pixelsz = pixelsize;
    
    srobj.Perform_DriftCorrection3D()
    
    drawnow
    if recon_stop == 1
        return;
    end

    n_imm = setup.n_imm;
    n_sample = setup.n_sample;
    srobj.step_ini = recon.dc.step_ini.*n_sample/n_imm;
    
    srobj.Perform_stackalignment()
end


%%
save(fullfile(setup.workspace,'recon4Pi'),'srobj');


