

function [PcudaM,crlbM,errM,maskc] = loc_channel_specific_4Pimodel(subregion_chs,probj,phi0_input,misobj_input,tform_all,offset_seg,subvar_chs)


%% set parameters from the retrieved pupil function
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


psfobj.PRstruct.Zernike_phase = PRstruct1.Zernike_phase;
psfobj.precomputeParam();
psfobj.genZernike();
psfobj.set2Pupil(PRstruct1,PRstruct2); 
psfobj.addMisalign(misobj_input)


%% generate sample PSF
tic

ref_point = [49, 49];   %calculate channel-specific model
[tform_fitting, ~] = cal_model_affine_4pi(tform_all,ref_point); 

phic = phi0_input;
pixelsize = pxsz;
boxsize = 25;

bin = 4;
psfsize = 128 * bin;
Nzs = 301; 
[samplepsf,startx,starty,startz,dz,dx] = gensamplepsf_4pi_affine(psfobj,pixelsize,psfsize,boxsize,bin,Nzs,phic,tform_fitting);

samplepsf_cuda = [];
samplepsf_4d = [];
for ii = 1:4
    samplepsf_cuda = cat(3,samplepsf_cuda,permute(reshape(samplepsf{ii},boxsize*bin,boxsize*bin,Nzs),[2,1,3,4]));
    samplepsf_4d = cat(4,samplepsf_4d,permute(reshape(samplepsf{ii},boxsize*bin,boxsize*bin,Nzs),[2,1,3,4]));  %For calculating spline parameters     
end
samplepsf_cuda = single(samplepsf_cuda);

toc
%% Calculate spline parameters

[st] = genpsfstruct(samplepsf_4d(:,:,:,1),dx,dz,'matrix');
for ss = 2:4
    [sti] = genpsfstruct(samplepsf_4d(:,:,:,ss),dx,dz,'matrix');
    st = catstruct(st,sti,3);
end

sobj = struct('samplepsf_cuda',samplepsf_cuda,'st',st,'dx',dx,'dz',dz,'startx',startx,'starty',starty,'startz',startz);

%% Estimate initial parameters (pragraming)

tic

data = subregion_chs;
data(data<=0) = 1e-6;


% Get tform_fitting model and affine matrix without translation 
ref_point = [8, 8];   
[tform_fitting, tMatrix_noTranslation] = cal_model_affine_4pi(tform_all,ref_point);


x1 = genIniguess(mean(data,4),'median'); 

Nz = 61;  
z_input = linspace(-1.1,1.1,Nz)+0.0;
[z_guess,mcc_val,psf_model] = genini4pi_z_parfor_affine(data,sobj,z_input,tMatrix_noTranslation);  %FX changed

maskc = mcc_val>0.;

x_next = x1(maskc,1)-bxsz/2;
y_next = x1(maskc,2)-bxsz/2;
z_next = z_guess(maskc);

I_next = x1(maskc,3)*2;
bg_next = x1(maskc,4);
toc



%% Spline localization

data_tmp = data(:,:,maskc,:);
var_tmp = subvar_chs(:,:,maskc,:);

boxsize_data = 16;
All_fit = size(data_tmp, 3);
interval = 20000;
N_loop = ceil(All_fit/interval);

PcudaM = [];
crlbM = [];
errM = [];
for ii = 1 : N_loop
    index_start = (ii-1) * interval + 1;
    if ii == N_loop
        index_end = All_fit;
    else
        index_end = ii * interval;
    end
    
    % input data
    data_selection = data_tmp(:,:,index_start:index_end,:);
    Nfit = size(data_selection,3);
    data_cuda = single(data_selection(:)); %psfsize x Nfit x quadrantN
    
    % input sCMOS noise
    var_selection = var_tmp(:,:,index_start:index_end,:);
    var_input = single(var_selection(:)); %psfsize x Nfit x quadrantN
    
    coords_cuda = cat(2,ones(Nfit,2),zeros(Nfit,1))';
    coords_cuda = single(coords_cuda(:));
    
    offset_seg_tmp = offset_seg(index_start:index_end,:);
    xtmp = x_next(index_start:index_end) + boxsize_data/2;
    ytmp = y_next(index_start:index_end) + boxsize_data/2;
    x0 = cat(2,xtmp,ytmp,z_next(index_start:index_end),I_next(index_start:index_end,:),bg_next(index_start:index_end,:));
    x0i = single(reshape(x0',Nfit*5,1));
    lambda = 0;
    
    iterateN = 50;
                              
    [P_cuda,CG,crlb,err,PSF_cuda] = cuda4pi_sCMOS_loc_transModel(data_cuda,coords_cuda,samplepsf_cuda,...
                                dx,dz,startx,starty,startz,iterateN,Nfit,lambda,x0i,...
                                single(st.Fx),single(st.Fy),single(st.Fz),...
                                single(st.Fxy),single(st.Fxz),single(st.Fyz),single(st.Fxyz),...
                                single(tMatrix_noTranslation),single(offset_seg_tmp),single(var_input)); 


    
    PcudaM = [PcudaM; reshape(P_cuda,5,Nfit)'];
    crlbM = [crlbM; reshape(crlb,5,Nfit)']; 
    errM = [errM; reshape(err,2,Nfit)'];

end

toc 
clear cuda4pi_sCMOS_loc_transModel                         
                     

%%


end