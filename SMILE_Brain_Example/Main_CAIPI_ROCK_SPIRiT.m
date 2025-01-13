close all
clear all
clc

addpath('ROCK_SPIRiT')
addpath('utils_espirit')
addpath('utils')

% Split the file due to github file limit
load('K1.mat') 
load('K2.mat')
load('K3.mat')

K = cat(3,K1,K2,K3);

% CAIPI
K_CAIPI = I2K(CAIPI(K2I(K)));


%% ROCK SPIRiT

Mask_Calib = false(1,size(K_CAIPI,2));
Mask_Calib(end/2-15:end/2+15) = true;



% ROCK
[RO_acs,slice_R] = readout_conc_prep(permute(sum(K_CAIPI,3),[1 2 4 3]),permute(K_CAIPI(:,end/2-15:end/2+15,:,:),[1 2 4 3]));

% This part calibrates the ROCK-SPIRiT kernels
[data_kspace,kernel_set,kernel_r,kernel_s] = ROCK_SPIRIT_kernel(permute(sum(K_CAIPI,3),[1 2 4 3]),RO_acs,31,5,slice_R);

% ROCK-SPIRIT Reconstruction
[K_CAIPI_ROCK_SPIRiT,~] = ROCKSPIRIT(data_kspace,[],kernel_set,kernel_r,kernel_s,slice_R,50);

K_CAIPI_ROCK_SPIRiT = K_CAIPI_ROCK_SPIRiT/sqrt(size(K,3));

I_CAIPI_ROCK = K2I(K_CAIPI_ROCK_SPIRiT);

I_CAIPI_ROCK = permute(cat(4,I_CAIPI_ROCK(1:end/6,:,:),I_CAIPI_ROCK(end/6+1:end/6*2,:,:),I_CAIPI_ROCK(end/6*2+1:end/6*3,:,:),...
    I_CAIPI_ROCK(end/6*3+1:end/6*4,:,:),I_CAIPI_ROCK(end/6*4+1:end/6*5,:,:),I_CAIPI_ROCK(end/6*5+1:end,:,:)),[1,2,4,3]);

K_CAIPI_ROCK_SPIRiT = I2K(I_CAIPI_ROCK);

figure; imshow(SSOS(reshape(DeCAIPI(K2I(K_CAIPI_ROCK_SPIRiT)),size(K_CAIPI,1),[],size(K_CAIPI,4))),[])

disp(['ROCK SPIRiT SNR ',num2str(SNR(K_CAIPI_ROCK_SPIRiT,K_CAIPI))])
% save([file_name(1:end-3),'_CAIPI_ROCK_SPIRiT.mat'],'K_CAIPI_ROCK_SPIRiT')


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
% x image: x y slice
% csm:     x y slice coil

if strcmp(opt,'notransp')
    x = reshape(x,SI(1),SI(2),SI(3));
    y = sum(I2K(x.*CSM),3);
    y = y(Mask);
else
    Kdata = zeros(size(Mask));
    Kdata(Mask) = x;
    x = Kdata;
    y = sum(K2I(x).*conj(CSM),4);
    y = y(:);
end
end

function y = A_Fun_CAIPI(x,opt,CSM,Mask,SI)
% x image: x y slice t 1
% csm:     x y slice 1 coil
% Mask:    x y 1     t coil

if strcmp(opt,'notransp')
    x = reshape(x,SI(1),SI(2),SI(3),SI(4),1);
    y = sum(I2K(x.*CSM),3);
    y = y(Mask);
else
    Kdata = zeros(size(Mask));
    Kdata(Mask) = x;
    x = Kdata;
    y = sum( K2I(x).*reshape(ones(SI(4),1),1,1,1,SI(4),1).*conj(CSM),5);
    y = y(:);
end
end