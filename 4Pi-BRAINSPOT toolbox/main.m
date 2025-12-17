%%
% Main script for running INSPR
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
% 
%% set environment paths

support_path = '..\Support'; %check support path

addpath([support_path '\PSF Toolbox_4pi']);
addpath([support_path '\SRsCMOS']);
addpath([support_path '\Helpers']);


%% call GUI

brainspot_4pi_GUI();