function B = func_opt_stem_tx(V_bar)

Y0 = 1 / 50;
[NT,NS] = size(V_bar);

R = real(V_bar)'; R1 = R(:,1:NS-1); R2 = R(:,NS:end);
J = imag(V_bar)'; J1 = J(:,1:NS-1); J2 = J(:,NS:end);

% Optimize B22
B22_22 = zeros(NT-NS+1,NT-NS+1);
for i = 1:NT-NS+1
    B22_22_tmp = [J1, J(:,NS-1+i)] \ R(:,NS-1+i);
    B22_22(i,i) = Y0 * B22_22_tmp(NS);
end
B22_12 = nan(NS-1,NT-NS+1);
for j = 1:NS-1
    for i = 1:NT-NS+1
        B22_12_tmp = [J1, J(:,NS-1+i)] \ R(:,NS-1+i);
        B22_12(j,i) = Y0 * B22_12_tmp(j);
    end
end
B22_21 = B22_12';
[UJ,SJ,VJ] = svd(J1);
B22_11 = VJ * inv(SJ(1:NS-1,:)) * UJ(:,1:NS-1)' * (Y0*R1 - J2*B22_21);
B22 = [[B22_11, B22_12]; [B22_21, B22_22]];

% Optimize other blocks of B
B12 = -Y0*J - R*B22;
B21 = B12';
B11 = -R * B21;

% Compute B, T and check capacity-achieving condition
B = [[B11, B12]; [B21, B22]];

end