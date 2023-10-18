close all
clear all
clc

addpath('utils_espirit')

% Split the file due to github file limit
load('K1.mat') 
load('K2.mat')
load('K3.mat')

K = cat(3,K1,K2,K3);

% SMILE
K_SMILE = I2K(reshape(K2I(K),size(K,1),prod(size(K,[2,3])),size(K,4)));

R = 6;
Mask = logical(repmat(GoFIX(size(K_SMILE,2),1,R),[1,size(K_SMILE,1),size(K_SMILE,3)])); % ky kx coil
Mask = permute(Mask,[2,1,3]);

K_SMILE_ob = K_SMILE.*Mask;                                                                     % observed k-space with zero filling

%% SENSE
% CSM
Mask_Calib = false(1,size(K_SMILE,2));
Mask_Calib(end/2-R/2*31:end/2+R/2*31-1) = true;
[~, CSM_SMILE, ~] = walsh(K2I(K_SMILE.*Mask_Calib));

% Forward Operator
fun = @ (x,opt)SENSE(x,opt,CSM_SMILE, Mask,size(K_SMILE));

I_SMILE_SENSE = reshape(lsqr(fun,K_SMILE_ob(Mask),5e-3,1e3,[],[],reshape(sum(K2I(K_SMILE_ob).*conj(CSM_SMILE),3),[],1)),size(K_SMILE,[1,2]));

K_SMILE_SENSE = I2K(I_SMILE_SENSE.*CSM_SMILE);

K_SMILE_SENSE(Mask) = K_SMILE_ob(Mask);

figure; imshow(SSOS(K2I(K_SMILE_SENSE)),[])
disp(['SENSE SNR ',num2str(SNR(K_SMILE_SENSE,K_SMILE))])


%% Functions
function I = K2I(K_POMP)                                                                                                 % k-space to image domain
I = sqrt(size(K_POMP,1)*size(K_POMP,2))*fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(K_POMP,1),2),[],1),[],2),1),2);
end

function K_POMP = I2K(I)                                                                                                 % image domain to k-sapce
K_POMP = 1/sqrt(size(I,1)*size(I,2))*fftshift(fftshift(fft(fft(ifftshift(ifftshift(I,2),1),[],2),[],1),2),1);
end

function I_SSOS = SSOS(I)                                                                                               % image to square root of sum of square
I_SSOS = sum(abs(I).^2,ndims(I)).^0.5;
end

function snr = SNR(X,Ref)                                                                                               % calculate the SNR
snr = -20*log10(norm(X(:)-Ref(:))/norm(Ref(:)));
end

function I_CAIPI = CAIPI (I)
I_CAIPI = I;
for s = 1:size(I,3)
    I_CAIPI(:,:,s,:) = circshift(I(:,:,s,:),round(size(I,2)/size(I,3)*(s-1)),2);
end
end

function I = DeCAIPI (I_CAIPI)
I = I_CAIPI;
for s = 1:size(I,3)
    I(:,:,s,:) = circshift(I_CAIPI(:,:,s,:),-round((size(I_CAIPI,2)/size(I_CAIPI,3)*(s-1))),2);
end
end

function y = SENSE(x,opt,CSM, Mask,SI)
% x image: x y
% csm:     x y  coil

if strcmp(opt,'notransp')
    x = reshape(x,SI(1),SI(2));
    y = I2K(x.*CSM);
    y = y(Mask);
else
    Kdata = zeros(size(Mask));
    Kdata(Mask) = x;
    x = Kdata;
    y = sum(K2I(x).*conj(CSM),3);
    y = y(:);
end
end