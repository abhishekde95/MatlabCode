function [L_norm1,M_norm1,S_norm1,L_norm2,M_norm2,S_norm2] = normalize_resp(L_exc1,M_exc1,S_exc1,L_exc2,M_exc2,S_exc2,mode)
% This function normalizes the responses from L, M and S cones
% 0- don't normalize, 1- do normalize
if mode == 0
    L_norm1 = L_exc1;
    M_norm1 = M_exc1;
    S_norm1 = S_exc1;
    L_norm2 = L_exc2;
    M_norm2 = M_exc2;
    S_norm2 = S_exc2;
else
    L_norm1 = L_exc1/norm([L_exc1;M_exc1;S_exc1]);
    M_norm1 = M_exc1/norm([L_exc1;M_exc1;S_exc1]);
    S_norm1 = S_exc1/norm([L_exc1;M_exc1;S_exc1]);
    L_norm2 = L_exc2/norm([L_exc2;M_exc2;S_exc2]);
    M_norm2 = M_exc2/norm([L_exc2;M_exc2;S_exc2]);
    S_norm2 = S_exc2/norm([L_exc2;M_exc2;S_exc2]);
end
end

