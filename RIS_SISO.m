function [Phi,Pow] = RIS_SISO(hd,hR,hT)

% Description: Finds the diagonal Reconfigurable Intelligent Surface
% (RIS) that maximizes the SNR computed as |hd + hR'*Phi*hT|^2.
% Phi = diag(e^{1i*theta1},\ldots, e^{1i*thetaM}) is an MxM diagonal matrix

% Input parameters: 
%
% h: direct link
% hR: Mx1 column vector from receiver to BDRIS
% ht: Mx1 column vector from transmitter to BDRIS

%
% Output parameters:
% Phi: MxM RIS matrix
% Pow: power achieved by the solution computed as  |\h_R^H\Phi\h_T|^2


Phi = diag(exp(-1i*angle(conj(hR).*hT)));      % RIS solution that maximizes SNR (or capacity)
if hd ~=0                                % rotation to align with the direct channel
    alpha = angle((hR'*Phi*hT)*hd);    
    Phi = exp(1i*alpha)*Phi;
end
Channel = hd + hR'*Phi*hT;
Pow = abs(Channel).^2;