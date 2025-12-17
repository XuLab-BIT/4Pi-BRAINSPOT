% Script for segmentation for 4pi images
% (C) Copyright 2022                The Huang Lab
%
%     All rights reserved           Weldon School of Biomedical Engineering
%                                   Purdue University
%                                   West Lafayette, Indiana
%                                   USA
%
%     Author: Fan Xu
%%  

function [subregion_chs,frame,l,t] = find_4pi_subregion(qd1,qd2,qd3,qd4,tform_all,boxsz,thresh,setup,cavity)

if cavity.isNewdata == 1
    if setup.is_sCMOS   %sCMOS case
        % sCMOS parameters
        offsetim_chs = repmat(setup.sCMOS_input.qd_offset,[1 1 size(qd1,3) 1]);
        varim_chs = repmat(setup.sCMOS_input.qd_var,[1 1 size(qd1,3) 1]);
        gainim_chs = repmat(setup.sCMOS_input.qd_gain,[1 1 size(qd1,3) 1]);
        
        qd1_in = (qd1 - offsetim_chs(:,:,:,1)) ./ gainim_chs(:,:,:,1);
        qd2_in = (qd2 - offsetim_chs(:,:,:,2)) ./ gainim_chs(:,:,:,2);
        qd3_in = (qd3 - offsetim_chs(:,:,:,3)) ./ gainim_chs(:,:,:,3);
        qd4_in = (qd4 - offsetim_chs(:,:,:,4)) ./ gainim_chs(:,:,:,4);
    else    %EMCCD case
        qd1_in = (qd1 - setup.offset) /setup.gain;
        qd2_in = (qd2 - setup.offset) /setup.gain;
        qd3_in = (qd3 - setup.offset) /setup.gain;
        qd4_in = (qd4 - setup.offset) /setup.gain;
    end
    
    qd1_in(qd1_in<=0) = 1e-6;
    qd2_in(qd2_in<=0) = 1e-6;
    qd3_in(qd3_in<=0) = 1e-6;
    qd4_in(qd4_in<=0) = 1e-6;
else
    qd1_in = qd1;
    qd2_in = qd2;
    qd3_in = qd3;
    qd4_in = qd4;
end
 
%% 

data1 = qd1_in;
for ii = 2 : 4
    data_tmp = eval(['qd', num2str(ii), '_in']);
    data_registered = imwarp(data_tmp,tform_all{ii},'cubic','OutputView',imref2d(size(data1)));
    eval(['data', num2str(ii), '= data_registered;']);
end

imsz_original = size(data1);
rangemin=[15, 15];    %5
rangemax=[imsz_original(1)-15, imsz_original(1)-15];

ims_ch1 = single(data1(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));
ims_ch2 = single(data2(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));
ims_ch3 = single(data3(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));
ims_ch4 = single(data4(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));

ims_detect = ims_ch1 + ims_ch2 + ims_ch3 + ims_ch4;

%%  Set parameters
allcds = [];    
if setup.is_imgsz == 1
    thresh_dist = floor(boxsz / sqrt(2))-1;   %Two detected spots must be less than this distance
else 
    thresh_dist = 6;
end

%% Detect single molecules in 4 channels 
display('Use uniform filter and maximum filter to obtain sub_region centers');

imsz = size(ims_detect, 1);
x=[];
y=[];
t=[];

sz = 3; 
[filteredim1] = unif_img(squeeze(ims_detect),sz);
sz = 9;
[filteredim2] = unif_img(squeeze(ims_detect),sz);
im_unif=filteredim1-filteredim2;

sz = 4;
loc_max = (im_unif>=.999999999 * imdilate(im_unif, true(sz)));
im_max_L = loc_max & ((im_unif>thresh(1))&(im_unif<=thresh(2)));
im_max_H = loc_max & (im_unif>thresh(2));


centers_L = findcoord_seg(im_max_L);
centers_L(:,4) = 0;     
centers_H = findcoord_seg(im_max_H);
centers_H(:,4) = 1; 

allcds = cat(1,centers_L,centers_H);




%% Remove too close spots
allcds_keep = [];
fnum = size(ims_detect, 3);

for f = 0 : fnum - 1
    tmp_loc = allcds(allcds(:,3) == f, :);
    tmp_dist = pdist(tmp_loc);
    
    num_tmp = 1;
    tmp_remove = [];
    for ii = 1 : size(tmp_loc, 1) - 1
        for jj = ii + 1 : size(tmp_loc, 1)
            if (tmp_dist(num_tmp) < thresh_dist)
                tmp_remove = [tmp_remove ii jj];
            end
            num_tmp = num_tmp + 1;
        end
    end
    
    tmp_loc(tmp_remove,:) = [];
    allcds_keep = [allcds_keep; tmp_loc];
end

%%
boundmask=(allcds_keep(:,1)<=boxsz/2)|(allcds_keep(:,1)>=imsz-boxsz/2)|(allcds_keep(:,2)<=boxsz/2)|(allcds_keep(:,2)>=imsz-boxsz/2) | allcds_keep(:,4) == 0;
allcds_mask = allcds_keep(~boundmask,:);



%% Segmentation from two channels
[subregion_ch1 l1 t1]=cMakeSubregions(allcds_mask(:,2),allcds_mask(:,1),allcds_mask(:,3),boxsz, ims_ch1); 
[subregion_ch2 l2 t2]=cMakeSubregions(allcds_mask(:,2),allcds_mask(:,1),allcds_mask(:,3),boxsz, ims_ch2); 
[subregion_ch3 l3 t3]=cMakeSubregions(allcds_mask(:,2),allcds_mask(:,1),allcds_mask(:,3),boxsz, ims_ch3); 
[subregion_ch4 l4 t4]=cMakeSubregions(allcds_mask(:,2),allcds_mask(:,1),allcds_mask(:,3),boxsz, ims_ch4); 

%% output
subregion_chs = cat(4,subregion_ch1,subregion_ch2,subregion_ch3,subregion_ch4);
frame = allcds_mask(:,3);
l = l1;
t = t1;

end

