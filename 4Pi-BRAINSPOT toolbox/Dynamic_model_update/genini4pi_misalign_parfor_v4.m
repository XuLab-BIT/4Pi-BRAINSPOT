

%%
function [misobj_guess,z_guess,mcc_val,psf_model] = genini4pi_misalign_parfor_v4(data,psfobj_input,psf_model,init_misobj,z_input,misobj_input,flag_misobj)


R = size(data,1);
Nfit = size(data,3);
Nplane = size(data,4);
Nz = length(z_input);
Num_misobj = length(misobj_input);

tic

if isempty(psf_model)
    psf_model = zeros(R,R,Num_misobj,Nz,Nplane);

    label_ch = {'p1','s2','p2','s1'};
    parfor kk = 1:Num_misobj
        psfobj = psfobj_input;
        psfobj.Boxsize = R;
        psfobj.Xpos = zeros(Nz,1);% pixel
        psfobj.Ypos = zeros(Nz,1);% pixel
        psfobj.Zpos = z_input;% micron
        
        
        % input different objective misalignment parameters  
        aber_misalign = [];
        switch flag_misobj
            case 1
%                 disp('estimate x tilt')
                aber_misalign = init_misobj + [misobj_input(kk) 0 0];
            case 2
%                 disp('estimate y tilt')
                aber_misalign = init_misobj + [0 misobj_input(kk) 0];
            case 3
%                 disp('estimate defocus')
                aber_misalign = init_misobj + [0 0 misobj_input(kk)];
            otherwise
                disp('check the input')
        end
        

        psfobj.addMisalign(aber_misalign)
        
        psfobj.genPupil_4pi_2();
        psfobj.genPSF_4pi_md();
        
        for ss = 1:Nplane
            psfobj.PSFs = psfobj.PSF4pi.(label_ch{ss});
            psfobj.scalePSF();
            psf_model(:,:,kk,:,ss) = psfobj.ScaledPSFs;
        end
    end

end
   
toc


%%
z_guess = zeros(Nfit,1);
misobj_guess = zeros(Nfit,1);
mcc_val = zeros(Nfit,1);
img_fft2 = zeros(R,R,Nfit,Nplane);
ref_fft2 = zeros(R,R,Num_misobj,Nz,Nplane);

% tic
parfor kk = 1:Num_misobj
    for ii = 1: Nz
        for ss = 1:Nplane
            ref = psf_model(:,:,kk,ii,ss);
            ref = ref-mean(ref(:));
            ref = ref/std(ref(:));
            ref_fft2(:,:,kk,ii,ss) = fft2(ref);
        end
    end
end

% toc

%%
% tic
parfor nn = 1:Nfit
    for ss = 1:Nplane
        img = data(:,:,nn,ss);
        img = img-mean(img(:));
        img = img/std(img(:));
        img_fft2(:,:,nn,ss) = fft2(img);
    end
end
% toc


%% calculate cross correlation
tic

parfor nn = 1:Nfit
    mcc_max = 0;
    ind_z = 1;
    ind_misobj = 1;
    mcc_backup = zeros(Num_misobj,Nz);
    for kk = 1:Num_misobj
        for ii = 1:Nz
            mcc = 0;
            cc_value = zeros(R,R);
            for ss = 1:Nplane
                cc_value = abs(ifft2(ref_fft2(:,:,kk,ii,ss) .*conj(img_fft2(:,:,nn,ss))));
                maxa = 1/(R*R)*max(cc_value(:));
                
                mcc = mcc + maxa/Nplane;
                
            end

            mcc_backup(kk,ii) = mcc;
            
        end
    end
        
    [mcc_max, max_idx] = max(mcc_backup(:));

    [ind_misobj,ind_z]=ind2sub(size(mcc_backup),max_idx);
    z_guess(nn) = z_input(ind_z);
    misobj_guess(nn) = misobj_input(ind_misobj);
    mcc_val(nn) = mcc_max;
end
toc



