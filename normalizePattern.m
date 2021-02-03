function In = normalizePattern(I)

%We normalize the interferograms for the 2 step methods which
%require normalized patters

[NR, NC]=size(I);
[u,v]=meshgrid(1:NC, 1:NR);

%Temporal variables
u0=floor(NC/2)+1; v0=floor(NR/2)+1;
u=u-u0; v=v-v0;

%High pass Fourier filtering by Gaussian filter with sigma=freq
freq = 0.5;
H=1-exp(-(u.^2+v.^2)/(2*(freq)^2));
C=fft2(I);
CH=C.*ifftshift(H);
I=real(ifft2(CH));

C=fft2(I);
CH=C.*ifftshift(H);
In=real(ifft2(CH));