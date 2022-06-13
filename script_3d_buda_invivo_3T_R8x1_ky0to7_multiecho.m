%--------------------------------------------------------------------------
%% 3d-epi R8 AP/PA @ 1 mm iso
%--------------------------------------------------------------------------
clear;
close all;

addpath('E:\GitHub_Files\utls\');

load('kspace_ap.mat');load('kspace_pa.mat');
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% BUDA Recon
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
addpath('.\utls');

esp = prot_epi1.iEffectiveEpiEchoSpacing * 1e-6;   % echo spacing;
% Rz = 1;       % Rz factor
% Ry = prot_epi1.lAccelFactPE;     % PAT factor

% 3D-BUDA parameter settings
opt.step_size = 1;
opt.num_iter = 100;   % 100 iterations
opt.tol = 0.1;
opt.winSize = [1,1]*5; % kernel size
opt.lambda_msl = 1.5;  % lambda of hankel low rank matrix
opt.esp = esp;     % echo-spacing

opt.fista = 1;
opt.t_k = 1;opt.fpa = para.fpa;opt.pf = para.pf;

k_ap = kspace_ap;
k_pa = kspace_pa;
ky_idx_ap = find(sq(k_ap(1,:,1,1))~=0);
ky_idx_pa = find(sq(k_pa(1,:,1,2))~=0);

% load B0 maps
load('img_fieldmap.mat')
load('sense_map.mat')

% dt = load_nifti(fpField);
% img_fieldmap = dt.vol;
% img_fieldmap =crop (img_fieldmap, [N(1),N(2),num_slc ]);


tic
    [img_buda] = recon_3d_BUDA_ME(k_ap, k_pa,ky_idx_ap,ky_idx_pa,sens_csm, Ry, Rz,img_fieldmap,opt);
toc
apodization_para = 0.2;

if apodization_para > 0
    for ii = 1:s(img_buda,4)
    img_buda(:,:,:,ii,ndwi) = apodize(img_buda(:,:,:,ii,ndwi),apodization_para);
    end
end

%% show the reconstructed image
imagesc3d2(mean(abs(img_buda),3), N/2, 15, [-1,-1,1]*90, [-0,3e-3], [], 'buda')

imagesc3d2(angle(img_buda(:,:,1,:)), N/2, 12, [-1,-1,1]*90, [-pi,pi], [], 'ap buda')
imagesc3d2(angle(img_buda(:,:,2,:)), N/2, 13, [-1,-1,1]*90, [-pi,pi], [], 'pa buda')

imagesc3d2(img_pocs(:,:,1,:), N/2, 3, [-1,-1,1]*90, [-0,3e-3], [], 'ap pocs')
imagesc3d2(img_pocs(:,:,2,:), N/2, 4, [-1,-1,1]*90, [-0,3e-3], [], 'pa pocs')