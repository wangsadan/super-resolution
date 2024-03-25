clear
clc

addpath(genpath('LTTR_file'))
addpath(genpath('data'))
addpath(genpath('Utilize'))
%% load data
F=create_F();
sf = 8;
sz=[512 512];
s0=1;
psf        =    fspecial('gaussian',8,3);  %%%7,2  %8,3
par.fft_B      =    psf2otf(psf,sz);
par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);
par.P=create_F();

% F=F(:,3:31);
for band = 1:size(F,1)
    div = sum(F(band,:));
    for i = 1:size(F,2)
        F(band,i) = F(band,i)/div;
    end
end

%%
im_structure =load('.\data\face_ms.mat');
S = im_structure.b;
%% im_structure =load('.\data\chart_and_stuffed_toy_msUint8b.mat');
%S0 = im_structure.b;
%  filename=dir('data\face_ms\*.png');
% for k=1:length(filename)
%     I=imread(filename(k).name);
%     S0(:,:,k) =double(I);
% end
% S=rescale(S0,0,1);

% S=S(:,:,3:31);
[M,N,L] = size(S);
S_bar = hyperConvert2D(S);
hyper= par.H(S_bar);
Y_h = hyperConvert3D(hyper, M/sf, N/sf);
Y = hyperConvert3D((F*S_bar), M, N);
%% parameters
% para.lambda_rr=[2e-4 3e-4 4e-4 5e-4 6e-4 7e-4];%4e-4; %2e-4 47.3359///4e-4 1e-3 1e-2 big
% para.mur=[8e-5 1e-4 2e-4 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4];
% para.lambda_tr =[1e-4 1e-3 0.01 0.05 0.1 0.5 1 5 10 100 1e3];
% para.betar = [1e-9 4e-9 8e-9 1e-8 4e-8 8e-8 1e-7 4e-7 8e-7 1e-6 4e-6];
%
% timer=[];
% Zr=[];
% psnr44=[];
% rmse44=[];
% ergas44=[];
% sam44=[];
% uiqi44=[];
% ssim44=[];
% DD44=[];
% CC44=[];

% for rrr=1:length(para.lambda_tr)
    %% low rank
    para.K=160;  %patch number
    
    para.lambda_r= 3e-4; %4e-4
    para.mu=1e-4;%
    
    %% group sparse
    para.lambda_t =1e-4;%1  1.2;
    para.beta =1e-6; %  %1e-7  1e-5 / 1e-3  big
    
    %% unmixing
    par.KK=120;  %120
    
    par.betaG=1e-6; %1e-2
    par.betaS=1e-6; %1e-2
    par.betaS2=1e-6; %1e-2
    
    par.lambda3=2e-3; %2e-4  0.2
    par.lambda5=1e-3; %0.1,1e-3
    
    %tau=0.01;   % 0.001;
    par.lambda4=1e-4;%par.lambda3*tau  2e-6
    par.lambda6=1e-6; %par.lambda5*tau
    
    par.gamma1=5;% 5
    par.lambda_MCP=1;%1
    
    par.vepsilon=0.01; %same with tip2018 cyber2012
    
    %%
    para.outIter=13; %12  10  25  100
    para.error = 1e-4; 
    %%
    t0=clock;
    Z =  LTTRGSWmix(Y_h,Y,F,sf,S,para,par);
    time=etime(clock,t0)
    [psnr4,rmse4, ergas4, sam4, uiqi4,ssim4,DD4,CC4] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z)), 0, 1.0/sf);
    
%     timer(rrr)=time;
%     Zr(:,:,:,rrr)=Z;
%     psnr44(rrr)=psnr4;
%     rmse44(rrr)=rmse4;
%     ergas44(rrr)=ergas4;
%     sam44(rrr)=sam4;
%     uiqi44(rrr)=uiqi4;
%     ssim44(rrr)=ssim4;
%     DD44(rrr)=DD4;
%     CC44(rrr)=CC4;
% end







