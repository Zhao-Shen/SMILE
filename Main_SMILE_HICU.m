close all
clear all
clc

% Split the file due to github file limit
load('K1.mat') 
load('K2.mat')
load('K3.mat')

K = cat(3,K1,K2,K3);


% SMILE
K_SMILE = I2K(reshape(K2I(K),size(K,1),prod(size(K,[2,3])),size(K,4)));

figure; imshow(SSOS(K2I(K_SMILE)),[])

%% HICU
R = 6;
Mask = logical(repmat(GoFIX(size(K_SMILE,2),1,R),[1,size(K_SMILE,1),size(K_SMILE,3)])); % ky kx coil
Mask = permute(Mask,[2,1,3]);


%% HICU Reconstruction
K_SMILE_ob = K_SMILE.*Mask;                                                                     % observed k-space with zero filling

%% HICU Reconstruction
Kernel_size = [5,5*R,size(K_SMILE,3)];                                                          % [tunable] kernel size: [5,5,5,Nc] is empirically large enough for most 2D+T imaging
Rank = 300;                                                                                     % [tunable] rank
Proj_dim = 1*size(K_SMILE,3);                                                                   % [tunable] projected nullspace dimension: Nc~4*Nc empirically balances between SNR and speed for 2D+T
Denoiser = [];                                                                                  % [tunable] denoising subroutine (optional), no denoiser G = []
Iter_1 = 1;                                                                                     % [tunable] number of iterations: 100 works for R6 and R8
Iter_2 = 1;                                                                                     % [tunable] number of iterations for gradient descent (GD) + exact line search (ELS)
GD_option = 1;                                                                                  % [tunable] options of calculating Grammian and GD, 1: without padding -> accurate & slow, 2. with circular padding approximation applied to Grammain and GD calculation using FFT -> less accurate and fast with large kernels. To reproduce the results in Ref [2], GD_option = 1
Max_time = 1e6;                                                                                 % [tunable] maximum running time (optional) if not empty the algorithm will stop either at maximum number of iterations or maximum running time

% Warm start using center of true k-space!!!
disp('Process the center k-space......');
[~, Null_c, SNR_c, Time_c] = HICUsubroutine_2D(K_SMILE(:,end/2-R/2*31:end/2+R/2*31-1,:), true(size(K_SMILE,1),R*31,size(K_SMILE,3)), K_SMILE(:,end/2-R/2*31:end/2+R/2*31-1,:), [], Kernel_size, Rank, Proj_dim, Denoiser, Iter_1, Iter_2, GD_option, Max_time, squeeze(K_SMILE(:,end/2-R/2*31:end/2+R/2*31-1,:)));



% Process on full k-space array
Iter_1 = 1;                                                                                     % [tunable] number of iterations
Iter_2 = 200;                                                                                    % [tunable] number of iterations for gradient descent + exact line search
Proj_dim = 4*size(K_SMILE,3);                                                                   % [tunable] projected nullspace dimension: Nc~4*Nc empirically balances between SNR and speed for 2D+T
Max_time = 7.2e4- Time_c(end);                                                                  % [tunable] maximum running time (optional) if not empty the algorithm will stop either at maximum number of iterations or maximum running time
Ref = K_SMILE;

disp('Process the full k-space......')
[K_SMILE_HICU, Null, SNR_o, Time_o] = HICUsubroutine_2D(K_SMILE_ob, Mask, K_SMILE_ob, Null_c, Kernel_size, Rank, Proj_dim, Denoiser, Iter_1, Iter_2, GD_option, Max_time, K_SMILE);


SNR(K_SMILE_HICU,K_SMILE)
figure; imshow(SSOS(K2I(K_SMILE_HICU)),[])
% save([file_name(1:end-3),'_SMILE.mat'],'K_SMILE')
% save([file_name(1:end-3),'_SMILE_HICU.mat'],'K_SMILE_HICU')




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