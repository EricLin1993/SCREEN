clear all, clc, close all
load('PureLapalceNMRData.mat')
%%
T1list = logspace(-1,1,40);
T2list = logspace(-1,1,40);
DifCoe = linspace(0,10,40);
D2D = squeeze( D3D(:,1,:) );
s =  D3D(:);
%%
K3 = -(1-2*exp(-taut1.'*(1./T1list)));
K2 = exp(-taut2.'*(1./T2list));
K1 = exp(-b.'*DifCoe);
K = kron(K3,kron(K2,K1)) ;
% ------------ SVD Compress ------------------
[ K,s  ] = IPSVDCompress( K,s, 50 );
%================ SCREEN ================================
sn = s; 
sn_max = max(sn);
sn = sn./sn_max;
ParaIn.KComplexYes = 0;
ParaIn.KMatrixYes = 1;
ParaIn.rel_tol = 1e-8;
ParaIn.K = K;
ParaIn.KSize = size(K);
ParaIn.s = sn;
ParaIn.lambda = 0.01;
ParaIn.nonneg = 1;
tic
[ ParaOut ] = SCREEN( ParaIn );
x = ParaOut.x;
xdiff = ParaOut.xdiff;
toc
    x =x./max(abs(x)); 
    X_tnipm = reshape( x,size(K1,2),size(K2,2),[]);
figure,plot(xdiff,'linewidth',2),title('xdiff')
%% Show Results
Dlable = -log10((fliplr(DifCoe)+eps)*1e-10);
T2lable =  log10(T2list);
T1lable =  log10(T1list);
X3D_tnipm = X_tnipm;
X3D_tnipm(X3D_tnipm < 0.15) = 0; 
figure,hold on,contourslice(T2lable,Dlable,T1lable,X3D_tnipm,T2lable,Dlable,T1lable,50); 
xlabel('log(\bf{\it{T}_1})','fontname','times new roman','fontsize',14);
ylabel('-log_{10}(\bf{\it{D}})','fontname','times new roman','fontsize',14);
zlabel('log(\bf{\it{T}_2})','fontname','times new roman','fontsize',14);
set (gca,'XGrid','on','YGrid','on','ZGrid','on',  'XTick',[-1,0,1],'ZTick',[-1,0,1],'fontsize',14,'fontname','times new roman')
view(3); 
xlim([min(T1lable),max(T1lable)]),ylim([min(Dlable),9.5]),zlim([min(T2lable),max(T2lable)])