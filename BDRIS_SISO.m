function [Phi,Pow] = BDRIS_SISO(hd,hR,hT,Mg)

% Description: Finds the Beyond-diagonal Reconfigurable Intelligent Surface
% (BDRIS) with a group-connected architecture. For each group, the solution
% is given by the Takagi's factorization of a certain matrix as described
% in [1].

% Input parameters: 
%
% h: direct link
% hR: Mx1 column vector from receiver to BDRIS
% ht: Mx1 column vector from transmitter to BDRIS
% Mg: we divide the M elements in G groups of Mg = M/G elements each
%
% Output parameters:
% Phi: MxM BDRIS matrix, with Mg=M/G it is a block-diagonal matrix where each block is unitary and symmetric 
% Pow: power achieved by the solution computed as  |\sum_g \h_{R,g}^H\Phi\h_{T,g}|^2


% [1] I. Santamaria, M. Soleymani, E. Jorswieck, J. Gutierrez, "SNR
% Maximization in Beyind Diagonal RIS-assietd single and multiple anetnna
% links, IEEE Signal Porcessing Letters, vol. 30, pp. 923 - 926, 2023. doi:
% 10.1109/LSP.2023.3296902

M = length(hT);   % Number of BDRIS elements 
G = fix(M/Mg);    % number of groups of Mg elements each
Channel = 0;      % equivalent SISO channel
Phi = [];         % to store the BDRIS matrix
for gg = 1:G
    % find gth block matrix
    hTg = hT((gg-1)*Mg+1: gg*Mg);  % coefficients Tx channel
    hRg = hR((gg-1)*Mg+1: gg*Mg);  % coefficients Rx channel

    %% Here we compute the Takagi's decomposition of the MgxMg matrix Ag
    % to this end we use the SVD  as explained in Remarks 2 of [1]
    Ag = (hRg*hTg'+ conj(hRg*hTg')')/2;   % complex and symmetric
    [F_SVD,~,G_SVD] = svd(Ag);
    dummy = F_SVD'*conj(G_SVD);
    % Eigenvectors post-correction
    phases = angle(diag(dummy))/2;
    F_SVD = F_SVD*diag(exp(1i*phases));
    Phig = F_SVD*F_SVD.';                % MgxMg unitary and symmetric block
    Phi = blkdiag(Phi,Phig);             % we add the block to the BDRIS matrix
    Channel = Channel +  hRg'*Phig*hTg;  % equivalent SISO channel
end
if hd ~=0                                % rotation to align with the direct channel
    alpha = angle((hR'*Phi*hT)*hd);    
    Phi = exp(1i*alpha)*Phi;
end
Channel = hd + hR'*Phi*hT;
Pow = abs(Channel).^2;
