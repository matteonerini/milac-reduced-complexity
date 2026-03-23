function B = func_opt_stem_rx(U_bar)

Y0 = 1 / 50;
[NR,NS] = size(U_bar);

R = real(U_bar)'; R1 = R(:,1:NS-1); R2 = R(:,NS:end);
J = imag(U_bar)'; J1 = J(:,1:NS-1); J2 = J(:,NS:end);

% Optimize B22
B11_22 = zeros(NR-NS+1,NR-NS+1);
for i = 1:NR-NS+1
    B11_22_tmp = [J1, J(:,NS-1+i)] \ R(:,NS-1+i);
    B11_22(i,i) = -Y0 * B11_22_tmp(NS);
end
B11_12 = nan(NS-1,NR-NS+1);
for j = 1:NS-1
    for i = 1:NR-NS+1
        B11_12_tmp = [J1, J(:,NS-1+i)] \ R(:,NS-1+i);
        B11_12(j,i) = -Y0 * B11_12_tmp(j);
    end
end
B11_21 = B11_12';
[UJ,SJ,VJ] = svd(J1);
B11_11 = VJ * inv(SJ(1:NS-1,:)) * UJ(:,1:NS-1)' * (-Y0*R1 - J2*B11_21);
B11 = [[B11_11, B11_12]; [B11_21, B11_22]];

% Optimize other blocks of B
B21 = Y0*J - R*B11;
B12 = B21';
B22 = -R * B12;

% Compute B, T and check capacity-achieving condition
B = [[B11, B12]; [B21, B22]];

end