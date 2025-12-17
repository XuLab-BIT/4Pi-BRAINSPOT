% Script for normalizing sub-regions
% (C) Copyright 2020                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu
% 
function subregion_chs_norm = subregion_normalization_4pi(subregion_chs)

imsz = size(subregion_chs,1);
Nfit = size(subregion_chs,3);
Nplane = size(subregion_chs,4);
subregion_chs_norm = zeros(imsz,imsz,Nfit,Nplane);

for nn = 1 : 4
    for ii = 1 : size(subregion_chs,3)
        tmp_ch = subregion_chs(:,:,ii,nn);
        tmp_ch = reshape(zscore(tmp_ch(:)),size(tmp_ch,1),size(tmp_ch,2));
        subregion_chs_norm(:,:,ii,nn) = tmp_ch; 
    end
end

