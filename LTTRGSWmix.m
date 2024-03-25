function  [Z] =  LTTRGSWmix( HSI, MSI,R,sf,S0,para,par)
%%
K=para.K;
FBm=par.fft_B;

lambda_r=para.lambda_r;
mu=para.mu;

lambda_t=para.lambda_t;
beta = para.beta;

%%
error = para.error;
outIter=para.outIter;

%%
rel_error = [];
rmse2=[];

%%
patchsize=8;
overlap=4;
bparams.block_sz = [patchsize, patchsize];
bparams.overlap_sz=[overlap overlap];
[nr, nc,~]=size(MSI);
L=size(HSI,3);
num1=(nr-patchsize)/(patchsize-overlap)+1;
num2=(nc-patchsize)/(patchsize-overlap)+1;
bparams.block_num=[num1 num2];
fkmeans_omaxpt.careful = 1;
predenoised_blocks = ExtractBlocks(MSI, bparams);
Y2=Unfold(predenoised_blocks,size(predenoised_blocks),4);
[aa ]=fkmeans(Y2,K,fkmeans_omaxpt);

HSI_int=zeros(nr,nc,L);
HSI_int(1:sf:end,1:sf:end,:)=HSI;

FBmC  = conj(FBm);
FBs  = repmat(FBm,[1 1 L]);

FBCs1=repmat(FBmC,[1 1 L]);
HHH=ifft2((fft2(HSI_int).*FBCs1));
HHH1=hyperConvert2D(HHH);
psfY.w=nr/sf;
psfY.h=nc/sf;
psfY.W=nr;
psfY.H=nc;

Zt=HSI_int;
HR_load1=imresize(HSI, sf,'bicubic');
ZE=hyperConvert2D(HR_load1);
MSI3=Unfold(MSI,size(MSI),3);

n_dr=nr/sf;
n_dc=nc/sf;

%%
CCC=R'*MSI3+HHH1;
C1=R'*R+(3*mu+beta)*eye(size(R,2));
[Q1,Lambda]=eig(C1);
Lambda=reshape(diag(Lambda),[1 1 L]);
InvLbd=1./repmat(Lambda,[ sf*n_dr  sf*n_dc 1]);
B2Sum=PPlus(abs(FBs).^2./( sf^2),n_dr,n_dc);
InvDI=1./(B2Sum(1:n_dr,1:n_dc,:)+repmat(Lambda,[n_dr n_dc 1]));
%%
Eny_x_fft   = (abs(psf2otf([+1; -1], [nr,nc,L]))).^2  ;
Eny_y_fft   = (abs(psf2otf([+1, -1], [nr,nc,L]))).^2  ;
Eny_fft  =  Eny_x_fft + Eny_y_fft;

%% Initialization   auxilizary variables
sizeD=[nr,nc,L];
%% for low rank
V1=ZE;                  % V1 : auxiliary variable for Z_<1>
V2=ZE;                  % V2 : auxiliary variable for Z_<1>
V3=ZE;                  % V3 : auxiliary variable for Z_<1>
G1=zeros(size(V1));     % G1 : multiplier for V1-ZE
G2=zeros(size(V2));     % G2 : multiplier for V1-ZE
G3=zeros(size(V3));     % G3 : multiplier for V1-ZE
%% Group sparse parameters
F1 = zeros(sizeD);      % F1 : auxiliary variable for DM_x
F2 = zeros(sizeD);      % F2 : auxiliary variable for DM_y
M3=ZE;                  % M3 : auxiliary variable for Z
W3 = F1;                % W3 : multiplier for DM_x-F1
W4 = F2;                % W4 : multiplier for DM_y-F2
W2 = zeros(sizeD);      % W2 : multiplier for Z - M
W2_3=zeros(size(M3));

%% unmixing
%% Initialization
G=Zt;
betaG=par.betaG;
betaS=par.betaS;
betaS2=par.betaS2;
lambda3=par.lambda3;
lambda4=par.lambda4;
lambda5=par.lambda5;
lambda6=par.lambda6;

gamma1=par.gamma1;
lambda_MCP=par.lambda_MCP;

vepsilon=par.vepsilon;

%% Initialize the dictionary E-D
G33=Unfold(G,sizeD,3);
Q=randperm( size(G33,2) );
D=G33(:,Q(1:par.KK));
D=D./repmat((sqrt(sum(D.^2))+eps), size(D,1), 1);
%% Initialize V1-R
R1=rand(par.KK,size(G33,2));
S=R1;
LambdaS= zeros(size(R1));

Yh3=Unfold(HSI,size(HSI),3);
R2=rand(par.KK,size(Yh3,2));
S2=R2;
LambdaS2= zeros(size(R2));

LambdaG= zeros(size(G));
G_3=ZE;
LambdaG_3=zeros(size(G_3));
%%  main
for iter=1:outIter
    preZt=Zt;
    %% step1 update low rank Z
    
    %% updata Zt
    HR_HSI3=mu*(V1+G1/mu + V2+G2/mu + V3+G3/mu)+beta*(M3-W2_3/beta) + betaG*(G_3- LambdaG_3/betaG);
    C3=CCC+HR_HSI3;
    C30=fft2(reshape((Q1\C3)',[nr nc L   ])).*InvLbd;
    temp  = PPlus_s(C30/( sf^2).*FBs,n_dr,n_dc);
    invQUF = C30-repmat(temp.*InvDI,[ sf  sf 1]).*FBCs1; % The operation: C5bar- temp*(\lambda_j d Im+\Sum_i=1^d Di^2)^{-1}Dv^H)
    VXF    = Q1*reshape(invQUF,[nc*nc L])';
    ZE = reshape(real(ifft2(reshape(VXF',[nr nc L   ]))),[nc*nc L])';
    Zt=hyperConvert3D(ZE,nr, nc );
    
    %% update auxiliary tensors U=Z, V=Z, W =Z
    B1=ZE-G1/(2*mu);
    B1=hyperConvert3D(B1,nr, nc );
    predenoised_blocks1 = ExtractBlocks(B1, bparams);
    Z_block1=zeros( bparams.block_sz(1), bparams.block_sz(2),L, bparams.block_num(1)* bparams.block_num(2));
    
    B2=ZE-G2/(2*mu);
    B2=hyperConvert3D(B2,nr, nc );
    predenoised_blocks2 = ExtractBlocks(B2, bparams);
    Z_block2=zeros( bparams.block_sz(1), bparams.block_sz(2),L, bparams.block_num(1)* bparams.block_num(2));
    
    B3=ZE-G3/(2*mu);
    B3=hyperConvert3D(B3,nr, nc );
    predenoised_blocks3 = ExtractBlocks(B3, bparams);
    Z_block3=zeros( bparams.block_sz(1), bparams.block_sz(2),L, bparams.block_num(1)* bparams.block_num(2));
    
    for mn=1:max(aa)
        gg=find(aa==mn);
        XES=predenoised_blocks1(:,:,:,gg);
        [a, b, c, d]=size(XES);
        a1=min(a*b,c*d);
        a2=min(a*b*c,d);
        a3=a;
        c1=sqrt(a1)/(sqrt(a1)+sqrt(a2)+sqrt(a3));
        c2=sqrt(a2)/(sqrt(a1)+sqrt(a2)+sqrt(a3));
        c3=sqrt(a3)/(sqrt(a1)+sqrt(a2)+sqrt(a3));
        %%
        D1=lambda_r*c1;
        D2=lambda_r*c2;
        D3=lambda_r*c3;
        nsig2=1;
        XES = reshape(XES,[a*b c*d]);
        V11=WNNM( XES, D1/mu, nsig2 );
%         V11new=[];
%        for iii=1:size(XES,2)
%             V11new(:,iii)=softthre( XES(:,iii), D1/2/mu );
%        end
%         WW= max(lambda_MCP-V11new/gamma1,0);
%         for iii=1:size(XES,2)
%             V11new(:,iii)=softthre( XES(:,iii), D1/2/mu*WW(:,iii));
%         end
%         V11=V11new;
        
        V11=reshape(V11,[a b c d]);
        Z_block1(:,:,:,gg)=V11;
        
        XES=predenoised_blocks2(:,:,:,gg);
        [a, b, c, d ]=size(XES);
        XES = reshape(XES,[a*b*c d]);
        V22=WNNM( XES, D2/mu, nsig2 );
%           V22new=[];
%         for iii=1:size(XES,2)
%             V22new(:,iii)=softthre( XES(:,iii), D2/2/mu );
%         end
%         WW= max(lambda_MCP-V22new/gamma1,0);
%         for iii=1:size(XES,2)
%             V22new(:,iii)=softthre( XES(:,iii), D2/2/mu*WW(:,iii));
%         end
%         V22=V22new;
        
        V22=reshape(V22,[a b c d]);
        Z_block2(:,:,:,gg)=V22;
        
        XES=predenoised_blocks3(:,:,:,gg);
        [a, b, c, d ]=size(XES);
        XES = reshape(XES,[a b*c*d]);
        V33=WNNM( XES, D3/mu, nsig2 );
%           V33new=[];
%         for iii=1:size(XES,2)
%             V33new(:,iii)=softthre( XES(:,iii), D3/2/mu );
%         end
%         WW= max(lambda_MCP-V33new/gamma1,0);
%         for iii=1:size(XES,2)
%             V33new(:,iii)=softthre( XES(:,iii), D3/2/mu*WW(:,iii));
%         end
%         V33=V33new;
             
        V33=reshape(V33,[a b c d]);
        Z_block3(:,:,:,gg)=V33;
    end
    
    V1= JointBlocks(Z_block1, bparams);
    V1=hyperConvert2D(V1);
    V2= JointBlocks(Z_block2, bparams);
    V2=hyperConvert2D(V2);
    V3= JointBlocks(Z_block3, bparams);
    V3=hyperConvert2D(V3);
    %% update lagrangian O P Q
    G1=G1+mu*(V1-ZE);
    G2=G2+mu*(V2-ZE);
    G3=G3+mu*(V3-ZE);
    
    %% Step2 %%%%
    %% update group saprse F=DM M=Z
    
    %% - M subproblem update  M=Z
    diffT_p = diffT3(beta * F1 - W3, beta * F2 - W4, sizeD);
    temp1 = reshape(diffT_p + beta*Zt + W2, sizeD);
    z = real(ifftn(fftn(temp1) ./ (beta*Eny_fft + beta)));
    M = reshape(z,sizeD);
    M3=hyperConvert2D(M);
    %% - F1 and F2 subproblem update
    [diff_Mx, diff_My] = diff3(M, sizeD);
    F1 = Thres_21(diff_Mx+ W3/beta, lambda_t/beta);
    F2 = Thres_21(diff_My+ W4/beta, lambda_t/beta);
    
    %% - multiplier  W2 update
    W2 = W2 + beta*(Zt-M);
    W2_3=hyperConvert2D(W2);
    %% - multiplier W3 and W4 update
    W3 = W3+beta*(diff_Mx-F1);
    W4 = W4+beta*(diff_My-F2);
    
    %%  Step3  %%%
    %% update U1-S
    G33=Unfold(G,sizeD,3);
    temp2=D'*G33;
    temp=betaS*R1+LambdaS+lambda3*temp2;
    aa1=gcp;
    parfor jj=1:size(S,2) % solve the linear system column by column
        s=S(:,jj);
        snew= myPCG2(s,betaS,lambda3,D,temp(:,jj));
        Snew(:,jj)=snew;
    end
    S=Snew;
    %% update V1- R1
    for iii=1:size(R1,2)
        Rnew(:,iii)=softthre( S(:,iii)- LambdaS(:,iii)/betaS, lambda4/betaS );
    end
    %%
    WW= max(lambda_MCP-Rnew/gamma1,0);
   
    WW2=1./(abs(Rnew)+vepsilon);
    
    for iii=1:size(R1,2)
        Rnew(:,iii)=softthre( S(:,iii)- LambdaS(:,iii)/betaS, lambda4/betaS*WW(:,iii).*WW2(:,iii));
    end
    
    R1=Rnew;
    
    
    %% update U2-S2
    temp22=D'*Yh3;
    temp=betaS2*R2+LambdaS2+lambda5*temp22;
    aa1=gcp;
    parfor jj=1:size(S2,2) % solve the linear system column by column
        s2=S2(:,jj);
        s2new= myPCG2(s2,betaS2,lambda5,D,temp(:,jj));
        S2new(:,jj)=s2new;
    end
    S2=S2new;
    %% update V2-R2
    for iii=1:size(R2,2)
        R2new(:,iii)=softthre( S2(:,iii)- LambdaS2(:,iii)/betaS2, lambda6/betaS2 );
    end
    WW= max(lambda_MCP-R2new/gamma1,0);
    
    WW2=1./(abs(Rnew)+vepsilon);   %%%Weight
    for iii=1:size(R2,2)
        R2new(:,iii)=softthre( S2(:,iii)- LambdaS2(:,iii)/betaS2, lambda6/betaS2*WW(:,iii).*WW2(:,iii));
    end
    
    R2=R2new;
    
    %% update E-D
    b1     =    sum(S.^2, 2);
    H1     =    G33 - D*S;
    b2     =    sum(S2.^2, 2);
    H2     =    Yh3 - D*S2;
    
    for k = 1 : par.KK
        d_k_pre   =   D(:,k);
        d_k       =   max( D(:,k) + (H1*S(k,:)'+lambda5/lambda3*H2*S2(k,:)')/(b1(k)+lambda5/lambda3*b2(k)), 0 );
        d_k       =   d_k./max(norm(d_k), 1);
        D(:,k)    =   d_k;
        
        H10        =   H1;
        H1         =   H10 - (d_k - d_k_pre)*S(k,:);
        H20        =   H2;
        H2         =   H20 - (d_k - d_k_pre)*S2(k,:);
    end
    
    %% update G
    G=1/(lambda3+betaG)*(lambda3*Fold(D*S,sizeD,3)+betaG*(Zt+LambdaG/betaG));
    G_3=hyperConvert2D(G);
    %% update multipliers
    gamma=1.2;
    LambdaG=LambdaG-min(gamma*betaG,1e5)*(G-Zt);
    LambdaG_3=hyperConvert2D(LambdaG);
    LambdaS=LambdaS-min(gamma*betaS,1e5)*(S-R1);
    LambdaS2=LambdaS2-min(gamma*betaS2,1e5)*(S2-R2);
    
    %% compute the error of iters
    rmse2(iter)=getrmse(double(im2uint8(S0)),double(im2uint8(Zt)));
    rel_error(iter) = norm(Zt(:)-preZt(:),'fro') / norm(preZt(:),'fro');
    fprintf('iter= %d  rmse2=%f rel_error=%f\n', iter,rmse2(iter), rel_error(iter));
    if rel_error(iter) < error
        break;
    end
    
    
end
delete(gcp)  % close the parallel computation
%% output low-rank Z
Z=Zt;
end








