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

% ================= 3D-ITAMED =========================================
% Here, for convenience, FISTA, the core algorthm of 3D-ITAMeD, is directly extracted for calculation.    
sn = s; 
sn_max = max(sn);
sn = sn./sn_max;
lambda = 0.01;
tic
[ x,fista_xdiff,fista_obj]=fista(ones(1,size(K,2)), sn, lambda, K,10e4); 
xdiff=fista_xdiff;
figure,plot(fista_obj,'linewidth',2),title('fista obj]')
toc
    x =x./max(abs(x)); 
    X_tnipm = reshape( x,size(K1,2),size(K2,2),[]);
figure, 
figure,plot(xdiff,'linewidth',2),title('xdiff'),
%% Show Results
Dlable = -log10((fliplr(DifCoe)+eps)*1e-10);
T2lable =  log10(T2list);
T1lable =  log10(T1list);
X3D_tnipm = X_tnipm;
X3D_tnipm(X3D_tnipm < 0.15) = 0; 
figure,hold on,contourslice(T2lable,Dlable,T1lable,X3D_tnipm,T2lable,Dlable,T1lable,50); %DifCoe,T2list,T1list,
xlabel('log(\bf{\it{T}_1})','fontname','times new roman','fontsize',14);
ylabel('-log_{10}(\bf{\it{D}})','fontname','times new roman','fontsize',14);
zlabel('log(\bf{\it{T}_2})','fontname','times new roman','fontsize',14);
set (gca,'XGrid','on','YGrid','on','ZGrid','on',  'XTick',[-1,0,1],'ZTick',[-1,0,1],'fontsize',14,'fontname','times new roman')
view(3); %axis tight
xlim([min(T1lable),max(T1lable)]),ylim([min(Dlable),9.5]),zlim([min(T2lable),max(T2lable)])