close all
clear ;
clc;
pkg load image

O = imread('original.tif');
img = imread('4.tif');
M = medfilt2(img);
I = imsharpen(M,'Radius',0.5,'Amount',3);

ft=fftshift(fft2(I));
[m,n]=size(ft);

% define some functions
norm_img = @(img) (img - min(img(:))) / (max(img(:)) - min(img(:)));
show_spec = @(img) imshow(norm_img(log(abs(img)-min(abs(img(:)))+1.0001)));
gNotch = @(v,mu,cov) 1-exp(-0.5*sum((bsxfun(@minus,v,mu).*(cov\bsxfun(@minus,v,mu)))));

% by inspection
cx = 251;
cy = 188;

% distance of noise from center
wx1 = 240-cx;
wx2 = 262-cx;
wy1 = 158-cy;
wy2 = 218-cy;

% create notch filter
filt = ones(m,n);

% use gaussian notch with standard deviation of 5
sigma = 5;
[y,x] = meshgrid(1:n, 1:m);
X = [y(:) x(:)].';
filt = filt .* reshape(gNotch(X,[cx+wx1;cy+wy1],eye(2)*sigma^2),[m,n]);
filt = filt .* reshape(gNotch(X,[cx+wx2;cy+wy2],eye(2)*sigma^2),[m,n]);
filt = filt .* reshape(gNotch(X,[cx-wx1;cy-wy1],eye(2)*sigma^2),[m,n]);
filt = filt .* reshape(gNotch(X,[cx-wx2;cy-wy2],eye(2)*sigma^2),[m,n]);

% apply filter
ft = ft .* filt;

% compute inverse
ifft_ = ifft2(ifftshift(ft));
input_noise = norm_img(ifft_);


#geometric mean filter
[m,n] = size(input_noise);                  

output = zeros(m,n);                        %output image set with placeholder values of all zeros

val = 1;                                    %variable to hold new pixel value


for i = 2:m-2                               %loop through each pixel in original image

    for j = 2:n-2                           %compute geometric mean of 3x3 window around pixel

        p = input_noise(i-1,j-1);

        q = input_noise(i-1,j);

        r = input_noise(i-1,j+1);

        s = input_noise(i,j-1);

        t = input_noise(i,j);

        u = input_noise(i,j+1);

        v = input_noise(i+1,j-1);

        w = input_noise(i+1,j);

        x = input_noise(i+1,j+1);

        val = (p*q*r*s*t*u*v*w*x)^(1/9);

        output(i,j) = val;                  %set output pixel to computed geometric mean

        val = 1;                            %reset val for next pixel

    end

end

J = input_noise;

figure(1);imshow(O);title('Original Picture');

figure(2);imshow(img);title('Picture with Noise');

figure(3);imshow(log(1.25+abs(J./0.52)),'Radius',1.21,'Amount',2.5);
title('Picture without Noise');

