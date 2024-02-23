function [H,G,F] = ChannelsMIMO(M,Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,ray_fading, RiceFactor, ...
    pl_0, alpha_RIS,alpha_U,blocked)

% Description: Generate channels for a RIS-assisted MIMO link
%
% Input parameters:
% M: Number of RIS elements
% Nrx, Ntx: number of transmit and receive antennas
% PosTx_XYZ, PosRx_XYZ, PosRIS_XYZ: positions of Tx, Rx and RIS (x,y,z)
% ray_fading: 0 if all channels are Rayleigh
% RiceFactor, pl_0, alpha_RIS, alpha_U, blocked: Ricena factor, path loss
% exponents, etc
%
% Output parameters:
% H, G, F channels. H is Nrx x Ntx , G is M x Ntx, F is M x Nrx  

%% Coordinates of Tx, Rx and RIS
x_tx = PosTx_XYZ(1);
y_tx = PosTx_XYZ(2);
z_tx = PosTx_XYZ(3);
x_rx = PosRx_XYZ(1);
y_rx = PosRx_XYZ(2);
z_rx = PosRx_XYZ(3);
x_ris = PosRIS_XYZ(1);
y_ris = PosRIS_XYZ(2);
z_ris = PosRIS_XYZ(3);

d_RIS_rx  = sqrt((x_ris-x_rx).^2+(y_ris-y_rx).^2+(z_ris-z_rx).^2);
d_tx_RIS  = sqrt((x_ris-x_tx).^2+(y_ris-y_tx).^2+(z_ris-z_tx).^2);


%% Links Tx-RIS & RIS-Rx
pl_tx_ris_db  = pl_0-10*alpha_RIS*log10(d_tx_RIS);
pl_tx_ris_eff = 10^((pl_tx_ris_db)/20);

pl_ris_rx_db  = pl_0-10*alpha_RIS*log10(d_RIS_rx);
pl_ris_rx_eff = 10^((pl_ris_rx_db)/20);

if ray_fading == 1
    G = pl_tx_ris_eff* ...
        (1/sqrt(2)*randn(M,Ntx)+1i*1/sqrt(2)*randn(M,Ntx));

    F = pl_ris_rx_eff * ...
        (1/sqrt(2)*randn(M,Nrx)+1i*1/sqrt(2)*randn(M,Nrx));
else
    % ======================================================
    % =============== Modeling Rician Fading ===============
    % ======================================================
    % ---- Modeling the LOS links Tx-RIS-Rx ---
    phi_AoD1 = 2*pi*rand;  % Random angle of arrival and departure
    phi_AoA1 = 2*pi*rand;
   
    a_D_r = exp(1i*pi*(0:M-1)'*sin(phi_AoD1));
    a_D_t = exp(1i*pi*(0:Ntx-1)'*sin(phi_AoA1));
    
    G = pl_tx_ris_eff* ...
        ((sqrt(RiceFactor)/sqrt(RiceFactor+1))*a_D_r*a_D_t' ...
        + (1/sqrt(RiceFactor+1))*(1/sqrt(2)*randn(M,Ntx)+1i*1/sqrt(2)*randn(M,Ntx)));
    % ----------------------------------------
    clear a_D_t a_D_r
    phi_AoD1 = 2*pi*rand;  % New angles for F
    phi_AoA1 = 2*pi*rand;

    a_D_r = exp(1i*pi*(0:M-1)'*sin(phi_AoD1));
    a_D_t = exp(1i*pi*(0:Nrx-1)'*sin(phi_AoA1));
   
    F = pl_ris_rx_eff * ...
        ((sqrt(RiceFactor)/sqrt(RiceFactor+1))*a_D_r*a_D_t' ...
        +(1/sqrt(RiceFactor+1))*(1/sqrt(2)*randn(M,Nrx)+1i*1/sqrt(2)*randn(M,Nrx)));
end

%% Direct MIMO channel H
d_tx_rx  = sqrt((x_tx-x_rx)^2+(y_tx-y_rx)^2);
pl_tx_rx_db  = pl_0-10*alpha_U*log10(d_tx_rx);
pl_tx_rx_eff = 10^((pl_tx_rx_db)/20);
H = pl_tx_rx_eff*(1/sqrt(2)*randn(Nrx,Ntx)+...
    1i*1/sqrt(2)*randn(Nrx,Ntx));
if blocked == 1
    H = 0* H;
end

