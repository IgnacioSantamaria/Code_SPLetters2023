
% Description: This script reproduces Fig.2 in [1].
% We maximize the sum-rate in a 2-user MAC channel assisted by a fully connected 
% beyond-diagonal RIS (BD-RIS). 
% The direct channel is assumed to be blocked. 
% The solution is obtained through Takagi's factorization of 
% matrix A given in Corollary 1 of [1].
% We include for comparison a diagonal RIS with optimized or random phases
% 
%
% I. Santamaria, March 2023

% [1] I. Santamaria, M. Soleymani, E. Jorswieck, J. Gutierrez, "SNR
% Maximization in Beyind Diagonal RIS-assisted single and multiple antenna
% links, IEEE Signal Processing Letters, vol. 30, pp. 923 - 926, 2023. doi:
% 10.1109/LSP.2023.3296902

format compact
clc; clear;
M = 10:10:300;             % Number of RIS elements
K = 2;                     % Number of users
Ntx = 1;                   % Number of transmit antennas
Nrx = 1;                   % Number of receive antennas

numsim = 100;             % Number of simulations

%% Parameters to plot figures
fs = 12;
lw = 1.5;
ms = 8;

%% Parameter channels
ray_fading = 1;     % Set it "1" for uncorrelated Rayleigh fading channels
blocked = 1;        % Set to "1" if all direct links are blocked
RiceFactor = 3;     % Rician factor (for RIS channels when ray_fading = 0)

% ===== Parameters for Large-scale Path Loss =============
No = -174;          % Noise power (dBm/Hz)
BW = 2*10^6;        % Bandwidth (Hz)
sigma2 = No + 10*log10(BW); % Effective noise power
PdBm = [20 20];          % We take equal powers for the 2 users
P = 10.^((PdBm-sigma2)/10);

pl_0      = -30;     % Path loss at a reference distance (d_0)
alpha_RIS = 3.75;    % Path loss exponent for the RIS links
alpha_U   = 3.75;    % Path loss exponent for the direct linksd

% --- Position of the Tx-Rx/RIS (in meters) -----
PosRx_XYZ = [50, 0, 2];
PosRIS_XYZ = [40, 2, 5];
radius = 10;            % circle of radius 10 to randomly locate the users

% Variables to store the sum rates and their averages
SR_BDRIS = zeros(size(M));
SR_RIS = zeros(size(M));
SR_RISrnd = zeros(size(M));

for nn =1:numsim
    if rem(nn,10) == 0
        %% Plot intermediate results
        figure(1);clf;
        plot(M, SR_BDRIS/nn,'b-s','MarkerSize',ms,'LineWidth',lw);
        hold on
        plot(M, SR_RIS/nn,'r-o','MarkerSize',ms,'LineWidth',lw);
        plot(M, SR_RISrnd/nn,'g-d','MarkerSize',ms,'LineWidth',lw);
        xlabel('Number of RIS elements'); ylabel('Sum Rate (b/s/Hz)');
        legend('Fully Connected BD-RIS ','RIS (opt. phases)','RIS (random phases)', 'Location','best');hold off
        set(findall(gcf,'-property','FontSize'),'FontSize',fs)
        hold off

        disp(['simulation ', int2str(nn), ' out of ' int2str(numsim)])
    end

    for mm = 1:length(M)
        
        % Random positions for the 2 users in a circle centered at (0,0)
        dx = (rand(K,1)-0.5)*sqrt(radius); % xaxis
        dy = (rand(K,1)-0.5)*10 ;          % yaxis
        PosTx_XYZ1 = [dx(1), dy(1), 2];    % Location user 1
        PosTx_XYZ2 = [dx(2), dy(2), 2];    % Location user 2

        %% Generate channels 
        % channel user 1 (for the RIS to Rx we take the realization of this channel)
        [hd1,hT1,hR] = ChannelsMIMO(M(mm),Nrx,Ntx,PosTx_XYZ1, PosRx_XYZ,PosRIS_XYZ,ray_fading, RiceFactor, ...
            pl_0, alpha_RIS,alpha_U,blocked);
        % channel user 2
        [hd2,hT2,~] = ChannelsMIMO(M(mm),Nrx,Ntx,PosTx_XYZ2, PosRx_XYZ,PosRIS_XYZ,ray_fading, RiceFactor, ...
            pl_0, alpha_RIS,alpha_U,blocked);

        HT = [hT1*sqrt(P(1)) hT2*sqrt(P(2))];  % This is the M x 2 channel from users to RIS including the users' powers
        uR = hR/norm(hR);
        [Ftx,Ktx,Gtx] = svd(HT);

        %% Takagi factorization of A (Corollary 1)
        A = (uR*Ftx(:,1)'+ conj(uR*Ftx(:,1)')');    % Symmetric matrix
        [F_SVD,~,G_SVD] = svd(A);
        dummy = F_SVD'*conj(G_SVD);
        % Eigenvectors post-correction
        phases = angle(diag(dummy))/2;
        F_SVD = F_SVD*diag(exp(1i*phases));
        PhiBDRIS = F_SVD*F_SVD.';
        powRIS_Takagi = norm(hR'*PhiBDRIS*HT,'fro').^2;

        heqBDRIS = hR'*PhiBDRIS*HT;     % Equivalent MISO channel

        %% Aproximate solution for a diagonal RIS
        % For a diag RIS a closed-form solution does not exit.
        % The problem can be solved convexified and solved via cvx (Dobres'
        % paper) or via iterative algorithm (one phase at a time). However, for
        % simple (2 user) problem the difference wrt to the approximate
        % solution is negligible.

        hg = Ftx(:,1).*conj(hR);                       
        PhiRIS = diag(exp(-1i*angle(hg)));
        heqRISd = hR'*PhiRIS*HT;             % Equivalent MISO channel

        %% Solution with diagonal RIS random phases
        phasesRIS = 2*pi*rand(M(mm),1);           % Random phases
        PhiRISrnd = diag(exp(1i*phasesRIS)); % diagonal RIS
        heqRISrnd = hR'*PhiRISrnd*HT;        % Equivalent MISO channel

        %% Calculate rates 
        SR_BDRIS(mm) = SR_BDRIS(mm) + log2(real(1 + norm(heqBDRIS)^2));
        SR_RIS(mm) = SR_RIS(mm) +  log2(real(1 + norm(heqRISd)^2));
        SR_RISrnd(mm) = SR_RISrnd(mm) + log2(real(1 + norm(heqRISrnd)^2));
    end

end

%% Plot final results

figure(1);
plot(M, SR_BDRIS/nn,'b-s','MarkerSize',ms,'LineWidth',lw);
hold on
plot(M, SR_RIS/nn,'r-o','MarkerSize',ms,'LineWidth',lw);
plot(M, SR_RISrnd/nn,'g-d','MarkerSize',ms,'LineWidth',lw);
xlabel('Number of RIS elements'); ylabel('Sum Rate (b/s/Hz)');
legend('Fully Connected BD-RIS ','RIS (opt. phases)','RIS (random phases)', 'Location','best');hold off
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
hold off

