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
%% Script for segmentation from 4Pi single-molecule dataset
%  input: qd1, qd2, qd3, qd4
%  output: subregion_chs

%%
function [subregion_chs,seg_display] = crop_subregion_4pi(qd1,qd2,qd3,qd4,tform_all,boxsz,thresh,thresh_dist,setup)


%%

data1 = qd1;
for ii = 2 : 4
    data_tmp = eval(['qd', num2str(ii)]);
    data_registered = imwarp(data_tmp,tform_all{ii},'cubic','OutputView',imref2d(size(data1)));
    eval(['data', num2str(ii), '= data_registered;']);
end

imsz_original = size(data1);

if setup.is_imgsz == 1
    rangemin=[5, 5];
    rangemax=[imsz_original(1)-5, imsz_original(1)-5];
else 
    rangemin=[1, 1];
    rangemax=[imsz_original(1), imsz_original(1)];
end

ims_ch1 = single(data1(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));
ims_ch2 = single(data2(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));
ims_ch3 = single(data3(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));
ims_ch4 = single(data4(rangemin(1):rangemax(1),rangemin(2):rangemax(2),:));


ims_detect = (ims_ch1 + ims_ch2 + ims_ch3 + ims_ch4);

%% Detect single molecules in summed 4 channels
display('Use uniform filter and maximum filter to obtain sub_region centers');
imsz = size(ims_detect, 1);
allcds = [];   

x=[];
y=[];
t=[];

sz = 3; 
[filteredim1] = unif_img(squeeze(ims_detect),sz);
sz = 9;
[filteredim2] = unif_img(squeeze(ims_detect),sz);
im_unif=filteredim1-filteredim2;


sz = 3;
loc_max=(im_unif>=.999999999 * imdilate(im_unif, true(sz)));
im_max_L = loc_max & ((im_unif>thresh(1))&(im_unif<=thresh(2)));
im_max_H = loc_max & (im_unif>thresh(2));


centers_L = findcoord_seg(im_max_L);
centers_L(:,4) = 0;     %low threshold 
centers_H = findcoord_seg(im_max_H);
centers_H(:,4) = 1;     %high threshold

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

%% selected spots
boundmask=(allcds_keep(:,1)<=boxsz/2)|(allcds_keep(:,1)>=imsz-boxsz/2)|(allcds_keep(:,2)<=boxsz/2)|(allcds_keep(:,2)>=imsz-boxsz/2) | allcds_keep(:,4) == 0;
allcds_mask = allcds_keep(~boundmask,:);


%% Segmentation from 4 channels
[subregion_ch1 l1 t1]=cMakeSubregions(allcds_mask(:,2),allcds_mask(:,1),allcds_mask(:,3),boxsz, ims_ch1); %ims_ch1
[subregion_ch2 l2 t2]=cMakeSubregions(allcds_mask(:,2),allcds_mask(:,1),allcds_mask(:,3),boxsz, ims_ch2); %ims_ch2
[subregion_ch3 l3 t3]=cMakeSubregions(allcds_mask(:,2),allcds_mask(:,1),allcds_mask(:,3),boxsz, ims_ch3); %ims_ch2
[subregion_ch4 l4 t4]=cMakeSubregions(allcds_mask(:,2),allcds_mask(:,1),allcds_mask(:,3),boxsz, ims_ch4); %ims_ch2

subregion_chs = cat(4,subregion_ch1,subregion_ch2,subregion_ch3,subregion_ch4);

%% save display
seg_display.ims_ch1 = ims_ch1;
seg_display.ims_ch2 = ims_ch2;
seg_display.ims_ch3 = ims_ch3;
seg_display.ims_ch4 = ims_ch4;
seg_display.allcds_mask = allcds_mask;
seg_display.t1 = t1;
seg_display.l1 = l1;
