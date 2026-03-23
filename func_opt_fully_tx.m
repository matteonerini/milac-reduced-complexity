function B = func_opt_fully_tx(V,NS)

Y0 = 1 / 50;

reV = real(V);
imV = imag(V);
imV_inv = inv(imV);
imV_reV = imV\reV;
reV_imV = reV/imV;
B = Y0 * [[imV_reV(1:NS,1:NS), -imV_inv(1:NS,:)];
          [-imV_inv(1:NS,:)', reV_imV]];

end