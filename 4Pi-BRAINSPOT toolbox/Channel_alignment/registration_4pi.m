%
% (C) Copyright 2022                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu
%
%% Script for 4Pi registration
%  input: qd1, qd2, qd3 and qd4
%  output: tform_all (affine parameters in 4 channels)

%%
function tform_all = registration_4pi(qd1,qd2,qd3,qd4)

tform_all = [];
imref = [];
imouts = [];
im1 = qd1;
for ii = 1 : 4
    im2 = eval(['qd', num2str(ii)]);


    [tform,fixed,movingRegistered]= biplane_registration(im1,im2);
    
    tform_all{ii} = tform;
    imref = cat(2,imref,fixed);
    imouts = cat(2,imouts,movingRegistered);
end




