
% Description: This script reproduces Fig.1 in [1].
% The Max-SNR BDRIS solution maximizes the SNR in SISO link assisted by a fully connected 
% or group-connected beyond-diagonal RIS (BD-RIS). 
% The direct channel is assumed to be blocked. 
% The solution is obtained through Takagi's factorization of A = (hR*hT'+ conj(hR*hT')')/2.
% 
%
% I. Santamaria, March 2023

% [1] I. Santamaria, M. Soleymani, E. Jorswieck, J. Gutierrez, "SNR
% Maximization in Beyind Diagonal RIS-assisted single and multiple antenna
% links, IEEE Signal Processing Letters, vol. 30, pp. 923 - 926, 2023. doi:
% 10.1109/LSP.2023.3296902

format compact
clc; clear;
M = 2:2:64;           % Number of RIS elements
Ntx = 1;              % Number of transmit antennas
Nrx = 1;              % Number of receive antennas
numsim = 1e4;         % Number of independent simulations (use at least numsim = 1e4 to get smooth curves 1e4

%% Parameters to plot figures
fs = 12;
lw = 1.5;
ms = 8;

%% Parameter channels
ray_fading = 1;     % Set it "1" for uncorrelated Rayleigh fading channels
blocked = 1;        % Set to "1" if all direct links are blocked
RiceFactor = 3;     % Rician factor (for RIS channels when ray_fading = 0)

% ===== Parameters for Large-scale Path Loss =============
No = -174;                    % Noise power (dBm/Hz)
BW = 80*10^5;                 % Bandwidth (Hz)
sigma2 = No + 10*log10(BW);   % Effective noise power

pl_0      = -30;     % Path loss at a reference distance (d_0)
alpha_RIS = 3.75;    % Path loss exponent for the RIS links
alpha_U   = 3.75;    % Path loss exponent for the direct links


% --- Position of the Tx-Rx/RIS (in meters) -----
PosRx_XYZ = [50, 0, 2];
PosRIS_XYZ = [40, 2, 5];
radius = 10;            % circle of radius 10 to randomly locate the Tx


Mg2 = M(M/2==fix(M/2));              % M's multiples of 2   
Mg4 = M(M/4==fix(M/4));              % M's multiples of 4    
Mg8 = M(M/8==fix(M/8));              % M's multiples of 8    

PowerBDRISfc = zeros(size(M));       % Power fully-conected BDRIS
PowerRISfc = zeros(size(M));         % Power RIS
PowerBDRISg2 = zeros(size(Mg2));       % Power BDRIS with group size 2
PowerRISg2 = zeros(size(Mg2));         % Power RIS
PowerBDRISg4 = zeros(size(Mg4));       % Power BDRIS with group size 4
PowerRISg4 = zeros(size(Mg4));         % Power RIS
PowerBDRISg8 = zeros(size(Mg8));       % Power BDRIS with group size 8
PowerRISg8 = zeros(size(Mg8));         % Power RIS

for nn = 1:numsim   % Monte Carlo simulations
    if rem(nn,100) == 0
        disp(['simulation ', int2str(nn), ' out of ' int2str(numsim)])
        %% Plot intermediate results
        figure(1);clf;
        plot(M, 10*log10(PowerBDRISfc./PowerRISfc),'b-s','MarkerSize',ms,'LineWidth',lw);
        hold on
        plot(Mg8, 10*log10(PowerBDRISg8./PowerRISg8),'m-*','MarkerSize',ms,'LineWidth',lw);
        plot(Mg4, 10*log10(PowerBDRISg4./PowerRISg4),'k-d','MarkerSize',ms,'LineWidth',lw);
        plot(Mg2, 10*log10(PowerBDRISg2./PowerRISg2),'r-o','MarkerSize',ms,'LineWidth',lw);
        xlabel('Number of RIS elements'); ylabel('SNR Gain (dB)');
        title(['simulation ', int2str(nn), ' out of ' int2str(numsim)])
        legend('Fully-Connected','Group-Size = 8','Group-Size = 4','Group-Size = 2', 'Location','best');hold off
        set(findall(gcf,'-property','FontSize'),'FontSize',fs)
        hold off
    end
    %% fully-connected BDRIS
    for mm = 1:length(M)
        %mm
        %% Generate channels
        PosTx_XYZ = [(rand-0.5)*sqrt(radius), (rand-0.5)*sqrt(radius), 2];
        %PosTx_XYZ = [0, 0, 2];
        [hd,hT,hR] = ChannelsMIMO(M(mm),Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,ray_fading, RiceFactor, ...
            pl_0, alpha_RIS,alpha_U,blocked);

        %% Fully-connected BDRIS
        [PhiBDRISfc,PowBDRISfc] = BDRIS_SISO(hd,hR,hT,M(mm));
        %% Diagonal RIS
        [PhiRIS,PowRIS] = RIS_SISO(hd,hR,hT);
        % SNR Gain
        PowerBDRISfc(mm) = PowerBDRISfc(mm) + PowBDRISfc;
        PowerRISfc(mm) = PowerRISfc(mm) + PowRIS;
    end
    %% BDRIS with group size 2
    for mm = 1:length(Mg2)  % for BDRIS with group size 2

        %% Generate channels
         PosTx_XYZ = [sign(randn)*rand*sqrt(radius), sign(randn)*rand*sqrt(radius), 2];
        %PosTx_XYZ = [0, 0, 2];
        [hd,hT,hR] = ChannelsMIMO(Mg2(mm),Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,ray_fading, RiceFactor, ...
            pl_0, alpha_RIS,alpha_U,blocked);

        %% Fully-connected BDRIS
        [PhiBDRISg2,PowBDRISg2] = BDRIS_SISO(hd,hR,hT,2);
        %% Diagonal RIS
        [PhiRISg2,PowRISg2] = RIS_SISO(hd,hR,hT);
        % SNR Gain
        PowerBDRISg2(mm) = PowerBDRISg2(mm) + PowBDRISg2;
        PowerRISg2(mm) = PowerRISg2(mm) + PowRISg2;
    end
    %% BDRIS with group size 4
    for mm = 1:length(Mg4)  % for BDRIS with group size 2

        %% Generate channels
        PosTx_XYZ = [sign(randn)*rand*sqrt(radius), sign(randn)*rand*sqrt(radius), 2];
        %PosTx_XYZ = [0, 0, 2];
        [hd,hT,hR] = ChannelsMIMO(Mg4(mm),Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,ray_fading, RiceFactor, ...
            pl_0, alpha_RIS,alpha_U,blocked);

        %% Fully-connected BDRIS
        [PhiBDRISg4,PowBDRISg4] = BDRIS_SISO(hd,hR,hT,4);
        %% Diagonal RIS
        [PhiRISg4,PowRISg4] = RIS_SISO(hd,hR,hT);
        % SNR Gain
        PowerBDRISg4(mm) = PowerBDRISg4(mm) + PowBDRISg4;
        PowerRISg4(mm) = PowerRISg4(mm) + PowRISg4;
    end

    %% BDRIS with group size 8
    for mm = 1:length(Mg8)  % for BDRIS with group size 2

        %% Generate channels
         PosTx_XYZ = [sign(randn)*rand*sqrt(radius), sign(randn)*rand*sqrt(radius), 2];
        %PosTx_XYZ = [0, 0, 2];
        [hd,hT,hR] = ChannelsMIMO(Mg8(mm),Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,ray_fading, RiceFactor, ...
            pl_0, alpha_RIS,alpha_U,blocked);

        %% Fully-connected BDRIS
        [PhiBDRISg8,PowBDRISg8] = BDRIS_SISO(hd,hR,hT,8);
        %% Diagonal RIS
        [PhiRISg8,PowRISg8] = RIS_SISO(hd,hR,hT);
        % SNR Gain
        PowerBDRISg8(mm) = PowerBDRISg8(mm) + PowBDRISg8;
        PowerRISg8(mm) = PowerRISg8(mm) + PowRISg8;
    end
    
end

%% Plot final results
figure(1);clf;
plot(M, 10*log10(PowerBDRISfc./PowerRISfc),'b-s','MarkerSize',ms,'LineWidth',lw);
hold on
plot(Mg8, 10*log10(PowerBDRISg8./PowerRISg8),'m-*','MarkerSize',ms,'LineWidth',lw);
plot(Mg4, 10*log10(PowerBDRISg4./PowerRISg4),'k-d','MarkerSize',ms,'LineWidth',lw);
plot(Mg2, 10*log10(PowerBDRISg2./PowerRISg2),'r-o','MarkerSize',ms,'LineWidth',lw);
xlabel('Number of RIS elements'); ylabel('SNR Gain (dB)');
 legend('Fully-Connected','Group-Size = 8','Group-Size = 4','Group-Size = 2', 'Location','best');hold off
set(findall(gcf,'-property','FontSize'),'FontSize',fs)
hold off

