
%% calculate objective misalignment and phi0 using NCC
             

function [phi0_guess_keep,phi0_fit,misobj_guess_keep,misobj_fit] = caliphi0_and_misObj_NCC_v3(subregion_chs,probj,init_phi0,init_misobj)

if nargin<3
    init_phi0 = 0;
    init_misobj = [0 0 0];
end

%% Set parameters of pupil funcion
oprobj = probj;
PRstruct = oprobj.PRstruct;
bxsz = size(subregion_chs,1);
pxsz = oprobj.Pixelsize;
R = 128;
psfobj = PSF_4pi(PRstruct);  
psfobj.Boxsize = bxsz;
psfobj.Pixelsize = pxsz; % micron
psfobj.PSFsize = R;
psfobj.nMed = oprobj.PRstruct.RefractiveIndex;
psfobj.Phasediff = oprobj.Phasediff;                         % phase difference between s- and p-polarizations
psfobj.Iratio = oprobj.Iratio;                             % transmission ratio between top and bottom emission path, from 0 to 1
psfobj.ModulationDepth = oprobj.ModulationDepth;                      % modulation strength of interferometric PSFs,from 0 to 1 
psfobj.Phi0 = init_phi0;
psfobj.Zoffset = oprobj.Zoffset;

% set pupil function
PRstruct1 = oprobj.PRstruct1;
PRstruct2 = oprobj.PRstruct2;

psfobj.PRstruct.Zernike_phase = PRstruct1.Zernike_phase;
psfobj.precomputeParam();
psfobj.genZernike();
psfobj.set2Pupil(PRstruct1,PRstruct2); 

%% input sub-region images
data = subregion_chs;
data(data<=0) = 1e-6;

x1 = genIniguess(mean(data,4),'median'); 

min_photon = 1000;  
intmask = x1(:,3) < min_photon;
data_input = data(:,:,~intmask,:);

num_input = size(data_input,3);

if num_input < 800
    display('Warning! Less number of molecuels for phi0 estimation!');
elseif num_input > 2000
    num_sample = 2000;
    dist_sample = floor(num_input / num_sample); 
    data_input = data_input(:,:,1:dist_sample:end,:);
end

%% iteratively calculate objective misalignment

psfobj_input = psfobj;
misobj_guess_keep = [];

for ii = 1 : 4
    
    %% set model parameters
    Nz = 41;
    Num_phi0 = 31; 
    Num_misobj = 21;
    z_input = linspace(-0.8,0.8,Nz);
    phi0_input = linspace(-3,3,Num_phi0);
    misobj_input = linspace(-2,2,Num_misobj);
     
    %% esimate x tilt
    disp('estimate x tilt')
    flag_misobj = 1; % 1 for x tilt; 2 for y tilt; 3 for defocus      
    [xtilt_guess,z_guess,mcc_val,~] = genini4pi_misalign_parfor_v4(data_input,psfobj_input,[],init_misobj,z_input,misobj_input,flag_misobj);  %FX changed 
    [xtilt_guess_keep,xtilt_fit,totmask] = est_misobj_from_hist_v2(xtilt_guess,z_guess,mcc_val);
    data_input = data_input(:,:,~totmask,:);
    
    init_misobj(1) = init_misobj(1) + xtilt_fit;
    
    %% esimate y tilt
    disp('estimate y tilt')
    flag_misobj = 2; % 1 for x tilt; 2 for y tilt; 3 for defocus      
    [ytilt_guess,z_guess,mcc_val,~] = genini4pi_misalign_parfor_v4(data_input,psfobj_input,[],init_misobj,z_input,misobj_input,flag_misobj);  %FX changed
    [ytilt_guess_keep,ytilt_fit,totmask] = est_misobj_from_hist_v2(ytilt_guess,z_guess,mcc_val);
    data_input = data_input(:,:,~totmask,:);

    init_misobj(2) = init_misobj(2) + ytilt_fit;
    
    %% esimate defocus
    disp('estimate defocus')
    flag_misobj = 3; % 1 for x tilt; 2 for y tilt; 3 for defocus      
    [defocus_guess,z_guess,mcc_val,~] = genini4pi_misalign_parfor_v4(data_input,psfobj_input,[],init_misobj,z_input,misobj_input,flag_misobj);  %FX changed
    [defocus_guess_keep,defocus_fit,totmask] = est_misobj_from_hist_v2(defocus_guess,z_guess,mcc_val);
    data_input = data_input(:,:,~totmask,:);

    init_misobj(3) = init_misobj(3) + defocus_fit;
    
    %% esimate phi0
    disp('estimate phi0')
    psfobj_input.addMisalign(init_misobj)
    [phi0_guess,z_guess,mcc_val,~] = genini4pi_phi0_parfor_v3(data_input,psfobj_input,[],z_input,phi0_input);  
    [phi0_guess_keep,phi0_fit,totmask] = est_phi0_from_hist_v2(phi0_guess,z_guess,mcc_val);
    data_input = data_input(:,:,~totmask,:);
       
    psfobj_input.Phi0 = phi0_fit;
%     init_misobj
%     phi0_fit


end

misobj_fit = init_misobj;
misobj_guess_keep.xtilt = xtilt_guess_keep;
misobj_guess_keep.ytilt = ytilt_guess_keep;
misobj_guess_keep.defocus = defocus_guess_keep;

display(['Estimated objective misalignment: ' num2str(misobj_fit)]);
display(['Estimated phi0: ' num2str(phi0_fit)]);

disp('finish')




