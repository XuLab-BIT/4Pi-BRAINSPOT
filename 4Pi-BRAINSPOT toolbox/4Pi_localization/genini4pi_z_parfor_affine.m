
function [z_guess,mcc_val,psf_model] = genini4pi_z_parfor_affine(data,sobj,z_input,tMatrix_noTranslation)


R = size(data,1);
Nfit = size(data,3);
Nplane = size(data,4);
Nz = length(z_input);

%% Generate reference PSFs

x = zeros(Nz,1) + R/2;
y = zeros(Nz,1) + R/2;
z = z_input(:);
I = 1;
bg = 0;
coords_cuda = cat(2,ones(Nz,2),zeros(Nz,1))';
coords_cuda = single(coords_cuda(:));

x0 = cat(2,x,y,z,ones(Nz,1).*I,ones(Nz,1).*bg);
x0i = single(reshape(x0',Nz*5,1));
data_cuda = zeros(R,R,Nz*Nplane);
data_cuda = single(data_cuda(:));
offset_int = zeros(Nz,6);


iterateN = 0;
                            
[~,~,~,~,PSF_cuda] = cuda4pi_loc_transModel(data_cuda,coords_cuda,sobj.samplepsf_cuda,...
                                sobj.dx,sobj.dz,sobj.startx,sobj.starty,sobj.startz,iterateN,Nz,0,x0i,...
                                single(sobj.st.Fx),single(sobj.st.Fy),single(sobj.st.Fz),...
                                single(sobj.st.Fxy),single(sobj.st.Fxz),single(sobj.st.Fyz),single(sobj.st.Fxyz),...
                                single(tMatrix_noTranslation),single(offset_int)); 
                            

                            
clear cuda4pi_loc_transModel

psf_model = reshape(PSF_cuda,[R,R,Nz,Nplane]);

%%
z_guess = zeros(Nfit,1);
mcc_val = zeros(Nfit,1);
img_fft2 = zeros(R,R,Nfit,Nplane);
ref_fft2 = zeros(R,R,Nz,Nplane);


%cal ref fft
parfor ii = 1: Nz
    for ss = 1:Nplane
        ref = psf_model(:,:,ii,ss);
        ref = ref-mean(ref(:));
        ref = ref/std(ref(:));
        ref_fft2(:,:,ii,ss) = fft2(ref);
    end
end



%%

parfor nn = 1:Nfit
    %cal img fft
    for ss = 1:Nplane
        img = data(:,:,nn,ss);
        img = img-mean(img(:));
        img = img/std(img(:));
        img_fft2(:,:,nn,ss) = fft2(img);
    end
end


%% calculate cross correlation
% tic
parfor nn = 1:Nfit
    mcc_max = 0;
    ind_z = 1;
    for ii = 1:Nz
        mcc = 0;
        for ss = 1:Nplane
            cc_value = abs(ifft2(ref_fft2(:,:,ii,ss) .*conj(img_fft2(:,:,nn,ss))));
            maxa = 1/(R*R)*max(cc_value(:));
                mcc = mcc + maxa/Nplane;

        end
        if mcc>mcc_max
            ind_z = ii;
            mcc_max = mcc;
        end
        
    end
    z_guess(nn) = z_input(ind_z);
    mcc_val(nn) = mcc_max;
end
% toc



