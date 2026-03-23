function B = func_opt_fully_rx(U,NS)

Y0 = 1 / 50;

reU = real(U);
imU = imag(U);
imU_inv = inv(imU);
imU_reU = imU\reU;
reU_imU = reU/imU;
B = Y0 * [[-reU_imU, imU_inv(1:NS,:)'];
          [imU_inv(1:NS,:), -imU_reU(1:NS,1:NS)]];

end