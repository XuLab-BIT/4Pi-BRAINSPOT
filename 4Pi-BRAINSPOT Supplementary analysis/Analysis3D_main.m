close all
clear
clc

addpath('Function\');

path = 'D:\User\datasets\';
file = 'Cell11.mat';

%%  Parameter
para.sz = 144;
para.psz = 129;
para.zm = 20;
para.sigma = 1;
para.cm = jet(64)*0.9+0.1;

para.Recon_color_highb = 0.4;     %%

para.tk = 50;     %%    %nm 
para.wk = 550;     %% nm
para.k = 1;     %%

para.filepath = fullfile(path,file);
para.savepath = path;

load(para.filepath);
obj = struct('x',srobj.loc_x_f,'y',srobj.loc_y_f,'z',srobj.loc_z_f);

%%
close all

subname='NP';
while ~isempty(subname)
    prompt = 'Input subfolder name (Return ''NP'' to stop): ';
    subname = input(prompt,'s');

    if(subname=="NP" || isempty(subname))
        break;
    else
        Analysis_spine_section_3D(obj,para,subname,0);
        Analysis_spine_area_3D(para,subname,0.35);
    end

end