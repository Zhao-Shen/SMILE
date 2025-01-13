close all
clear all
clc

addpath(genpath('grappa-tools'))

% Split the file due to github file limit
load('K1.mat') 
load('K2.mat')
load('K3.mat')

K = cat(3,K1,K2,K3);


% CAIPI
K_CAIPI = I2K(CAIPI(K2I(K)));

figure; imshow(reshape(SSOS(K2I(K)),size(K_CAIPI,1),[]),[]); drawnow


%   Split-slice GRAPPA Multi-band Recon
%   Based on Cauley et al., MRM 2014
%
%   MChiew
%   Nov 2016

%   input is [c,kx,ky,1,t]
%   calib is [c,kx,ky,z]
%   kernel is (kx, ky)

tic

K_CAIPI_SG = single(permute(sg(permute(sum(K_CAIPI,3),[4,1,2,3]), permute(K_CAIPI(:,end/2-15:end/2+15,:,:),[4,1,2,3]), [11,11]),[2 3 4 1]));
SNR(K_CAIPI_SG,K_CAIPI)
figure; imshow(reshape(SSOS(DeCAIPI(K2I(K_CAIPI_SG))),size(K_CAIPI_SG,1),[]),[]); drawnow
toc

tic
K_CAIPI_SPSG = single(permute(spsg(permute(sum(K_CAIPI,3),[4,1,2,3]), permute(K_CAIPI(:,end/2-15:end/2+15,:,:),[4,1,2,3]), [15,15]),[2 3 4 1]));
SNR(K_CAIPI_SPSG,K_CAIPI)
figure; imshow(reshape(SSOS(DeCAIPI(K2I(K_CAIPI_SPSG))),size(K_CAIPI_SPSG,1),[]),[]); drawnow
toc



%% Functions
function I = K2I(K)                                                                                                 % k-space to image domain
I = sqrt(size(K,1)*size(K,2))*fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(K,1),2),[],1),[],2),1),2);
end

function K = I2K(I)                                                                                                 % image domain to k-sapce
K = 1/sqrt(size(I,1)*size(I,2))*fftshift(fftshift(fft(fft(ifftshift(ifftshift(I,2),1),[],2),[],1),2),1);
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
    I_CAIPI(:,:,s,:,:) = circshift(I(:,:,s,:,:),round(size(I,2)/size(I,3)*(s-1)),2);
end
end

function I = DeCAIPI (I_CAIPI)
I = I_CAIPI;
for s = 1:size(I,3)
    I(:,:,s,:,:) = circshift(I_CAIPI(:,:,s,:,:),-round((size(I_CAIPI,2)/size(I_CAIPI,3)*(s-1))),2);
end
end