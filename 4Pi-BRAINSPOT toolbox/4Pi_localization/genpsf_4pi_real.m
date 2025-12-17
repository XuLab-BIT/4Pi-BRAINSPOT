function [psf4pi] = genpsf_4pi_real(x,w,obj)

R = obj.Boxsize;
obj.Xpos = x(:,1).*w(1);
obj.Ypos = x(:,2).*w(2);
obj.Zpos = x(:,3).*w(3);
N = numel(x(:,1));


% generate the 4Pi PSF, with enlarged pupil
obj.precomputeParam();
obj.pad2Pupil();  
obj.genPupil_4pi_2();
obj.genPSF_4pi_md();


label = {'p1','s2','p2','s1'};
psf4pi = zeros(R,R,N,4);
for nn = 1:4
    obj.PSFs = obj.PSF4pi.(label{nn});
    obj.scalePSF();
    psfI = obj.ScaledPSFs;
    psf = psfI;

    I = x(:,nn+3);
    tmp = zeros(1,1,N);
    tmp(1,1,:) = I;
    IL = repmat(tmp,[R,R,1]).*w(4);

    bg = x(:,nn+7);
    tmp = zeros(1,1,N);
    tmp(1,1,:) = bg;
    bgL = repmat(tmp,[R,R,1]).*w(5);
    
    psf4pi(:,:,:,nn) = psf.*IL+bgL;
end

end
