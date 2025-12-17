function [samplepsf,startx,starty,startz,dz,dx] = gensamplepsf_4pi_affine(psfobj,pixelsize,psfsize,boxsize,bin,Nzs,phic,tform)

zpos = linspace(-1.3,1.3,Nzs)';
xpos = zeros(Nzs,1);
ypos = zeros(Nzs,1);
I = ones(Nzs,4);
bg = zeros(Nzs,4);
x1 = cat(2,xpos,ypos,zpos,I,bg);
w = [1,1,1,1,1];

psfobj.Pixelsize = pixelsize/bin;
psfobj.PSFsize = psfsize;
psfobj.Boxsize = boxsize*bin;
psfobj.Phi0 = phic;

[psf2d_fit] = genpsf_4pi_real(x1,w,psfobj);

psf_fit = [];
for ii = 1:4
    if ii ~= 1
        %affine model, fill boundary
        edge_1 = mean(squeeze(psf2d_fit(1,:,:,ii)),1)';
        edge_2 = mean(squeeze(psf2d_fit(:,1,:,ii)),1)';
        edge_3 = mean(squeeze(psf2d_fit(end,:,:,ii)),1)';
        edge_4 = mean(squeeze(psf2d_fit(:,end,:,ii)),1)';
        edge_fill = min([edge_1 edge_2 edge_3 edge_4],[],2);
        
        psf2d_fit(:,:,:,ii) =  imwarp(psf2d_fit(:,:,:,ii),tform{ii-1},'cubic','OutputView',imref2d(size(psf2d_fit(:,:,:,ii))),...
            'FillValues',edge_fill,'SmoothEdges',true);
        
        
    end
    
    psf_fit = cat(2,psf_fit,squeeze(psf2d_fit(:,:,:,ii)));
end
samplepsf = cell(4,1);
for ii = 1:4
    f0 = psf2d_fit(:,:,:,ii);
    N = size(f0,1);
    Nz = size(f0,3);
    F = reshape(permute(f0,[2,1,3]),N*N*Nz,1);

    samplepsf{ii} = F;
end
startx = -0.5*psfobj.Pixelsize*psfobj.Boxsize;
starty = -0.5*psfobj.Pixelsize*psfobj.Boxsize;
startz = zpos(1);
dz = zpos(2)-zpos(1);
dx = psfobj.Pixelsize;
end


