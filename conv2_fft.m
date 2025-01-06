function c_fft = conv2_fft(a,b,shape)
% function c_fft = conv2_fft(a,b) use c_fft=ifft2(fft2(a_zp).*fft2(b_zp)) to
% speed-up conv2 where a_zp and b_zp are a and b patched with zeros.
%
% If shape == 'full', returns the full two-dimensional convolution (default).
% If shape == 'same', return the center part of the convolution of the same
% size as a (same with conv2). In this case, dimensions of a and b need to
% be odd numbers. Since fft is fastest when the size is 2^N, dimensions are
% best set at 2^N-1
%
% Example: a=rand(1023);b=rand(1023);tic;c=conv2_fft(a,b,'same');toc
%
% First created on 2009/01/07, by Kuen-Yu Tsai, National Taiwan University,
% all rights reserved. 
% 01 2009/01/07 creation KT
% 02 2009/01/07 debug& Cast variable to "single" data type CHLiu 
%
%
a=[1 1 1 1 1;
   1 0 0 0 1;
   1 0 0 0 1;
   1 0 0 0 1;
   1 1 1 1 1];
b=[0 0 0 0 0;
   0 0 0 0 0;
   0 0 1 0 0;
   0 0 0 0 0;
   0 0 0 0 0];
% b=[0 0 0;
%    0 1 0;
%    0 0 0;];
 c=conv2(a,b,'full');
 cs=conv2(a,b,'same');

[ma,na]=size(a);[mb,nb]=size(b);
a_zp=cast(zeros(ma+mb-1, na+nb-1),'single');a_zp(1:ma, 1:na)=a;
b_zp=cast(zeros(ma+mb-1, na+nb-1),'single');b_zp(1:mb, 1:nb)=b;
c_fft=ifft2(fft2(a_zp).*fft2(b_zp));
%為何要fft

if(shape == 'full');
    c_fft = c_fft;
end 

if(shape == 'same');
    ma_half=(ma-1)/2;na_half=(na-1)/2;
    c_same_m_1st = (ma+mb)/2 - ma_half;
    c_same_n_1st = (na+nb)/2 - na_half;
    
    c_fft_same=c_fft( c_same_m_1st:c_same_m_1st+ma-1, c_same_n_1st:c_same_n_1st+na-1);
    c_fft=c_fft_same;
end

