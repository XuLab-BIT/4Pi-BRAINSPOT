


function [phi0_guess_keep,phi0_fit] = caliphi0_NCC_v3(subregion_chs,probj)


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

psfobj.Zoffset = oprobj.Zoffset;

% set pupil function
PRstruct1 = oprobj.PRstruct1;
PRstruct2 = oprobj.PRstruct2;


psfobj.precomputeParam();

psfobj.set2Pupil(PRstruct1,PRstruct2); 

%% calculate phi0 based on cross correlation between reference model and data

data = subregion_chs;
data(data<=0) = 1e-6;

x1 = genIniguess(mean(data,4),'median'); 

min_photon = 1000; 
intmask = x1(:,3) < min_photon;
data_input = data(:,:,~intmask,:);

num_input = size(data_input,3);


if num_input < 1000
    display('Warning! Less number of molecuels for phi0 estimation!');
elseif num_input > 2000
    num_sample = 2000;
    dist_sample = floor(num_input / num_sample); 
    data_input = data_input(:,:,1:dist_sample:end,:);
end

Nz = 41; 
Num_phi0 = 31;
z_input = linspace(-0.8,0.8,Nz);

phi0_input = linspace(-3,3,Num_phi0);

psfobj_input = psfobj;

[phi0_guess,z_guess,mcc_val,~] = genini4pi_phi0_parfor_v3(data_input,psfobj_input,[],z_input,phi0_input);  %FX changed



%% rejection 

ccmask = mcc_val < 0.6;
zmask = z_guess > 0.4 | z_guess < -0.4;
totmask = ccmask | zmask;
phi0_guess_keep = phi0_guess(~totmask);


%% fit Gaussian curve v3
mean_phi0 = angle(mean(exp(1i.*phi0_guess_keep)));
phi0_input = [phi0_guess_keep-2*pi; phi0_guess_keep; phi0_guess_keep+2*pi];

tmp_phi0_low =  mean_phi0 - pi;
tmp_phi0_high = mean_phi0 + pi;

index_sel = phi0_input > tmp_phi0_low & phi0_input < tmp_phi0_high;
phi0_input_sel = phi0_input(index_sel);

pd = fitdist(phi0_input_sel,'normal');
tmp_mu = pd.mu;

for ii = 1 : 5
    tmp_phi0_low =  tmp_mu - pi;
    tmp_phi0_high = tmp_mu + pi;
    
    index_sel_2 = phi0_input > tmp_phi0_low & phi0_input < tmp_phi0_high;
    phi0_input_sel = phi0_input(index_sel_2);
    
    pd_2 = fitdist(phi0_input_sel,'normal');
    
    if abs(tmp_mu - pd_2.mu) < 0.1
        break
    else
        tmp_mu = pd_2.mu;
    end
end

phi0_fit = pd_2.mu;


%% output 

plotflag = 0;
if plotflag == 1
    mean_phi0 = angle(mean(exp(1i.*phi0_guess)))
    figure; plot(phi0_guess,'o')
    hold on
    plot(phi0_guess*0+mean_phi0,'LineWidth',1)
    
    ylim([-3.8,3.8])
    xlabel('number of sub-regions')
    ylabel('Phi0 estimation (rad)')
    set(gca,'FontSize',12)
    
end
