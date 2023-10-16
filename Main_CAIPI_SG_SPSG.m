close all
clear all
clc

addpath(genpath('grappa-tools'))

for i = 1:20
    switch i
        case 1
            file_name = 'file_brain_AXT2_200_2000020.h5';
        case 2
            file_name = 'file_brain_AXT2_200_2000057.h5';
        case 3
            file_name = 'file_brain_AXT2_200_2000080.h5';
        case 4
            file_name = 'file_brain_AXT2_200_2000092.h5';
        case 5
            file_name = 'file_brain_AXT2_200_2000175.h5';
        case 6
            file_name = 'file_brain_AXT2_200_2000092.h5';
        case 7
            file_name = 'file_brain_AXT2_200_2000290.h5';
        case 8
            file_name = 'file_brain_AXT2_200_2000332.h5';
        case 9
            file_name = 'file_brain_AXT2_200_2000357.h5';
        case 10
            file_name = 'file_brain_AXT2_200_2000360.h5';
        case 11
            file_name = 'file_brain_AXT2_200_2000362.h5';
        case 12
            file_name = 'file_brain_AXT2_200_2000407.h5';
        case 13
            file_name = 'file_brain_AXT2_200_2000417.h5';
        case 14
            file_name = 'file_brain_AXT2_200_2000485.h5';
        case 15
            file_name = 'file_brain_AXT2_200_2000486.h5';
        case 16
            file_name = 'file_brain_AXT2_200_2000488.h5';
        case 17
            file_name = 'file_brain_AXT2_200_2000566.h5';
        case 18
            file_name = 'file_brain_AXT2_200_2000592.h5';
        case 19
            file_name = 'file_brain_AXT2_200_2000600.h5';
        case 20
            file_name = 'file_brain_AXT2_200_2000621.h5';
    end

    K = h5read(file_name,'/kspace');
    K = flip(permute(K.r + K.i.*1i, [2 1 4 3]),1);
    size(K)

    % Coil Compression
    [~,~,V] = svd(reshape(K,[],size(K,4)),'econ');
    K = reshape(reshape(K,[],size(K,4))*V(:,1:12),[size(K,[1 2 3]),12]);

    % choose slice
    K = K(:,:,end/2-3:1:end/2+2,:);
    

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
    

    save([file_name(1:end-3),'_SG_SPSG.mat'],'K_CAIPI_SG','K_CAIPI_SPSG')
end



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