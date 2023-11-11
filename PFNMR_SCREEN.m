
clear,close all,clc
load('Ny_Nz.mat')
addpath('2DFT_Toolbox')
%%  ---------------  Fully Sampled Data ---------------------
Fid = FID;
[N1,N2]=size(Fid);
Spec2D_Ideal = real(fft2(Fid)); 
Spec2D_Ideal = max(Spec2D_Ideal,0); Spec2D_Ideal = Spec2D_Ideal./max(max(Spec2D_Ideal));
ppm1 = linspace(107.8,131.6,N1) ;ppm2 = linspace(107.8,131.6,N2) ;
Spec2D_Ideal = flip(flip(Spec2D_Ideal,1),2);
figure,contour((ppm1),(ppm2),Spec2D_Ideal,(linspace(0.01,1,50))),title('Full Sampling')
set(gca,'xdir','reverse','ydir','reverse')
%% --------------- 25% NUS randomly chosen with uniform distribution ---------------
mask = zeros(N1,N2);
mask(rand(N1,N2)<0.25)=1;
FidNus = Fid;
FidNus(mask==0) = [];
FidNus = FidNus(:);
%% Fourier NUS Sampling Operator AM
AM = Fu_downsample(mask,N1,N2);
AMh = AM';
%% -------------- SCREEN -----------------------
sn = FidNus; 
sn_max = max(sn);
sn = sn./sn_max;
ParaIn.KComplexYes = 1;
ParaIn.KMatrixYes = 0;
ParaIn.rel_tol = 1e-5;
ParaIn.K = AM;
ParaIn.KSize = [76610,768800];
ParaIn.s = sn;
ParaIn.lambda = 0.01;
ParaIn.nonneg = 0;
tic
[ ParaOut ] = SCREEN( ParaIn );
x = ParaOut.x*sn_max;
xdiff = ParaOut.xdiff;
final_spectrum=reshape(x,N1,N2);
rec_time = toc 
fprintf('Completing, time consumes %0.2f sec \n',rec_time);
%% ================ find diagnol peak Index =========================
Sp = real(final_spectrum);
Sp = max(Sp,0); Sp = Sp./max(max(Sp));
Sp = fftshift(real(Sp));
Sp = flip(flip(Sp,1),2);
figure,contour((ppm1),(ppm2),Sp,(linspace(0.001,1,50))),title('Reconstruted')
set(gca,'xdir','reverse','ydir','reverse')
RecS = real(Sp)./max(max(real(Sp)));
IdeS = real(Spec2D_Ideal)./max(max(real(Spec2D_Ideal)));
rlne = norm(RecS-IdeS)/norm(IdeS)
NUSRate = sum(sum(mask))/(N1*N2);
IdeS_diag = diag(IdeS);RecS_diag = diag(RecS);
[pks_ideal,locs] = findpeaks(IdeS_diag,'minpeakheight',0.02);
pks_rec = RecS_diag(locs); 
CorDiagPeaks= corrcoef(IdeS_diag(locs),RecS_diag(locs))
figure,plot(IdeS_diag),hold on,plot(RecS_diag),title('diagnol extraction'),legend('Full','Rec') 
figure,plot(IdeS_diag(locs),RecS_diag(locs),'bo'),title('DiagnolPeaks Intensity')
%% ================ find cross peak Index =========================
CroPeakFlag = zeros(size(RecS));
for it = 1:length(locs)
   [~,Cro] = findpeaks(RecS(:,locs(it)));
  CroPeakFlag(Cro,locs(it))=1; 
end
CroPeakFlag = CroPeakFlag - diag(ones(1,length(CroPeakFlag)));
RecS_Cross = RecS(CroPeakFlag==1);
IdeS_Cross = IdeS(CroPeakFlag==1);
flag = RecS_Cross>0.004;
RecS_Crossf = RecS_Cross(flag);
% RecS_Crossf = RecS_Crossf./max(RecS_Crossf);
IdeS_Crossf= IdeS_Cross(flag);
% IdeS_Crossf=IdeS_Crossf./max(IdeS_Crossf);
CorCrossPeaks=corrcoef(RecS_Crossf ,IdeS_Crossf)
figure,plot(IdeS_Crossf,RecS_Crossf ,'ro'),
xlim([0,0.07]),ylim([0,0.07]),title('CrossPeaks Intensity')