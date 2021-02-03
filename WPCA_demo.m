% This Matlab script shows fringe pattern demodulation examples using 
% weigthed PCA [1], PCA [2] and two-step demodulation methods. The results obtained 
% by these methods are compared.

%This script reproduces all the results shown in the paper [1]. To visualize 
%the different results change the variable medida (values between 1 to 7).

%   Javier Vargas
%   29/11/20 
%   Copyright 2020
%   Universidad Complutense de Madrid
%   $ Revision: 1.0.0.0 $
%   $ Date: 22/06/19 $
%
%NOTE: If you use or redistribute these functions, please cite the papers
%[1], [2] and/or [3]
%
%REFERENCES 
%
%[1] "Robust weighted principal components analysis demodulation algorithm 
% for phase-shifting interferometry", send to publication in Optics Express
%
%[2] J Vargas, JA Quiroga, T Belenguer, "Phase shifting interferometry
%based on principal component analysis", Optics Letters 36(8) (2011)
%
%[3] T. M. Kreis, P. P Jueptner, "Fourier transform evaluation of 
%   interference patterns: demodulation and sign ambiguity", Proc of SPIE,
%   vol. 1553 (1992)

%%
warning off;
clear all;
close all;
clc;

%Cases shown in the paper [1]

%measure  =1; (simulation)
%measure = 2; (simulation)
%measure = 4; (simulation)
%measure = 8; (simulation)
%measure = 5; (experimental)
%measure = 7; (experimental)

medida = 7;

%% Generate and load data for the different cases
switch medida    
         
    case 1
                
        %Patterns to two step demodulation:
        num1 = 1;
        num2 = 2;

        %Number of rows and columns of the simulated interferograms
        N = 450;     
        
        %Number of simulated fringe patterns
        num  = 300;
        
        %Noise level
        noise = 0.8;
               
        %alpha
        alpha = 1;       
        
        %Actual phase map
        P = 4*peaks(N);
        
        %Phase-steps
        ps =  rand(1,num)*2*pi;
                
        %Amplitude of the tilt aberration
        f1=15*(randn(1,num));
        
        % The first two frames are not affected by tilt aberration
        f1(1)=0;
        f1(2)=0;
        ps(1)=0;
        ps(2)=pi/2;
               
        %We define the background (a) and modulation (b) maps
        [X Y]=meshgrid(1:N);        
        X = X-mean2(X);
        X = X/max(X(:));
        Y = Y-mean2(Y);
        Y = Y/max(Y(:));

        %amplitude term
        a = (X./50);
        %modulation term
        b = exp(-0.5.*(X.^2+Y.^2)./1e4);

        %Processing Mask
        Mask=(ones(size(X))>0.5);
        
        S = 0;
        %We compose the different simulated fringe patterns       
        for i=1:num
            Iref = b.*(cos(P+X*f1(i)+ps(i)));
            In = noise*randn(N).*b;
            p = (X./50)+Iref+In;
            S = S + mean( abs(Iref(:))./ std(In(:)));
            I(:,:,i) = (p.*Mask);
        end
        
        S = S/num;
        
        pwRef = mod(P,2*pi);
                
    case 2
                
        %Patterns to two step demodulation:
        num1 = 1;
        num2 = 2;

        %Number of rows and columns of the simulated interferograms
        N = 450;     
        
        %Number of simulated fringe patterns
        num  = 300;
        
        %Noise level
        noise = 0.8;
        
        %alpha
        alpha = 1;       
        
        %Actual phase map
        P = 4*peaks(N);
        
        %Phase-steps
        ps =  rand(1,num)*2*pi;
                
        %Amplitude of the tilt aberration
        f1=2*(randn(1,num));
        f1(1)=0;
        f1(2)=0;
        ps(1) = 0;
        ps(2) = pi/2;
               
        %We define the background (a) and modulation (b) maps
        [X Y]=meshgrid(1:N);        
        X = X-mean2(X);
        X = X/max(X(:));
        Y = Y-mean2(Y);
        Y = Y/max(Y(:));

        %amplitude term
        a = (X./50);
        %modulation term
        b = exp(-0.5.*(X.^2+Y.^2)./1e4);

        %Processing Mask
        Mask=(ones(size(X))>0.5);
        
        %We compose the different simulated fringe patterns       
        for i=1:num
            p = (X./50)+b.*(cos(P+X*f1(i)+ps(i)))+noise*randn(N).*b;
            I(:,:,i) = (p.*Mask);
        end
        
        pwRef = mod(P,2*pi);

    case 3
                
        %Patterns to two step demodulation:
        num1 = 1;
        num2 = 2;

        %Number of rows and columns of the simulated interferograms
        N = 450;        
        
        %Number of simulated fringe patterns
        num  = 300;

        %Noise level
        noise = 0.8;
        
        %alpha
        alpha = 1.5;
                
        %Actual phase map
        P = 4*peaks(N);
        
        %Phase-steps
        ps =  rand(1,num)*2*pi;
        ps(1) = 0;
        ps(2) = pi/2;
        
        %Varying phase factor
        f1=(linspace(1,1.5,num));
        
        %We define the background (a) and modulation (b) maps
        [X Y]=meshgrid(1:N);
        X = X-mean2(X);
        Y = Y-mean2(Y);
        
        a = (X./50);
        b = exp(-0.5.*(X.^2+Y.^2)./1e4);

        X = X/max(X(:));
        Y = Y/max(Y(:));
        
        %Processing Mask
        Mask=(ones(size(X))>0.5);
        
        %We compose the different simulated fringe patterns
        for i=1:num
            p = (X./50)+b.*(cos(P*f1(i)+ps(i)))+noise*randn(N).*b;
            I(:,:,i) = p.*Mask;
        end
        
        pwRef = mod(P,2*pi);
                
    case 4
                
        %Patterns to two step demodulation:
        num1 = 1;
        num2 = 2;

        %Number of rows and columns of the simulated interferograms
        N = 450;        
        
        %Number of simulated fringe patterns
        num  = 300;

        %Noise level
        noise = 0.8;
        
        %alpha
        alpha = 1.5;
                
        %Actual phase map
        P = 4*peaks(N);
        
        %Phase-steps
        ps =  rand(1,num)*2*pi;
        ps(1) = 0;
        ps(2) = pi/2;
        
        %We define the background (a) and modulation (b) maps
        [X Y]=meshgrid(1:N);
        X = X-mean2(X);
        Y = Y-mean2(Y);
        X = X/max(X(:));
        Y = Y/max(Y(:));
        
        a = (X./50);
        b = exp(-0.5.*(X.^2+Y.^2)./1e4);
       
        %Temporal aberration by Zernike Polynomials
        n = 3;
        [Theta,R]=cart2pol(X,Y);
        R = R(:);
        Theta = Theta(:);
        R=R./(max(R)+eps);
        Zpol = reshape(Zernike(3, R, Theta),450,450);
        clear Theta R s;

        %Varying amplitude
        f1=10*randn(1,N);
        f1(1)=0;
        f1(2)=0;
        
        %Processing Mask
        Mask=(ones(size(X))>0.5);
        
        %We compose the different simulated fringe patterns
        for i=1:num
            p = (X./50)+b.*(cos(P+ps(i)+f1(i)*Zpol))+noise*randn(N).*b;
            I(:,:,i) = p.*Mask;
        end
        
        pwRef = mod(P,2*pi);

            
    case 8
                
        %Patterns to two step demodulation:
        num1 = 1;
        num2 = 2;

        %Number of rows and columns of the simulated interferograms
        N = 450;        
        
        %Number of simulated fringe patterns
        num  = 300;

        %Noise level
        noise = 0.8;
        
        %alpha
        alpha = 1.5;
               
        %Phase-steps
        ps =  rand(1,num)*2*pi;
        ps(1) = 0;
        ps(2) = pi/2;
        
        %We define the background (a) and modulation (b) maps
        [X Y]=meshgrid(1:N);
        X = X-mean2(X);
        Y = Y-mean2(Y);
        X = X/max(X(:));
        Y = Y/max(Y(:));
        
        %Calculating actual phase map
        [Theta R]=cart2pol(X,Y);
        Z1 = Zernike(5, R(:), Theta(:));
        Z2 = Zernike(8, R(:), Theta(:));
        Z3 = Zernike(13, R(:), Theta(:));
        P = reshape(6*Z1+5*Z2+2*Z3,N,N);
    
        a = (X./50);
        b = exp(-0.5.*(X.^2+Y.^2)./1e4);
       
        %Temporal aberration by Zernike Polynomials
        n = 3;
        [Theta,R]=cart2pol(X,Y);
        R = R(:);
        Theta = Theta(:);
        R=R./(max(R)+eps);
        Zpol = reshape(Zernike(3, R, Theta),450,450);
        clear Theta R s;

        %Varying amplitude
        f1=10*randn(1,N);
        f1(1)=0;
        f1(2)=0;
        
        %Processing Mask
        Mask=(ones(size(X))>0.5);
        
        %We compose the different simulated fringe patterns
        for i=1:num
            p = (X./50)+b.*(cos(P+ps(i)+f1(i)*Zpol))+noise*randn(N).*b;
            I(:,:,i) = p.*Mask;
        end
        
        pwRef = mod(P,2*pi);

        
    case 5
        
        %Patterns to two step demodulation:
        num1 = 3;
        num2 = 4;
                
        %In file RealInt.mat there are 19 real phase-shifted       
        load ShouyuExp2
               
        %alpha
        alpha = 5;
        
        %Mask
        Mask = ones(1024,1024)>0.5;
 
        %In the cases medida 5, 6 and 7, we do not have a reference phase so we asume
        %that the reference phase is the PCA phase to correct the possible
        %piston only.
        K = 2;
        pwRef = pcaDemod(I,K,Mask);
        pwRef = mod(pwRef,2*pi);
        P = Unwrap_TIE_DCT_Iter(pwRef);
       
    case 6
          
        %Patterns to two step demodulation:
        num1 = 3;
        num2 = 4;
        
        %Experimental patterns    
        load ShouyuExp
        I = Iexp;
        clear Iexp;        
        
        %alpha
        alpha = 10;
        
        %Mask
        Mask = ones(1024,1024)>0.5;
 
        K = 2;
        
        %In the cases medida 5, 6 and 7, we do not have a reference phase so we asume
        %that the reference phase is the PCA phase to correct the possible
        %piston only.
        pwRef=pcaDemod(I,K,Mask);
        pwRef = mod(pwRef,2*pi);
        P = Unwrap_TIE_DCT_Iter(pwRef); 
        
            
    case 7
          
        %Patterns to two step demodulation:
        num1 = 2;
        num2 = 3;
        
        %Experimental patterns    
        load ShouyuExp3
        
        %alpha
        alpha = 3;
        
        %Mask
        Mask = ones(582,782)>0.5;
 
        K = 2;
        
        %In the cases medida 5, 6 and 7, we do not have a reference phase so we asume
        %that the reference phase is the PCA phase to correct the possible piston
        pwRef=pcaDemod(I,K,Mask);
        pwRef = mod(pwRef,2*pi);
        P = Unwrap_TIE_DCT_Iter(pwRef);   
end

%% Analysis

%Two step phase demodulation
tic
[pw,Mod,PC] = kreisDemod(I(:,:,num1),I(:,:,num2),Mask);
t1=toc;
fprintf('elapsed time Kreis is: %.2f seconds. \n',t1')

seed = PC;

%Number of Principal Components
K = 2;

%Weighted PCA results
tic
[pw2,Mod,PC2,V2,W2]=pcaWDemod(I,K,Mask,seed,alpha);
t2=toc;
fprintf('elapsed time PCA is: %.2f seconds. \n',t2')

%PCA results
tic
[pw3,Mod,PC3,V3,Xm3]=pcaDemod(I,K,Mask);
t3=toc;
fprintf('elapsed time WPCA is: %.2f seconds. \n',t3')

%GS results
tic
% We normalize the patterns for the two-step methods that require
% normalized patters
In1 = normalizePattern(I(:,:,num1));
In2 = normalizePattern(I(:,:,num2));
[pw4,Mod]=gsDemod(In1,In2,Mask);
t4=toc;
fprintf('elapsed time GS is: %.2f seconds. \n',t4')

%OF results
tic
[pw5,Mod]=rofDemod(I(:,:,num1),I(:,:,num2),5,0.3,35);
t5=toc;
fprintf('elapsed time OF is: %.2f seconds. \n',t5')

pw = mod(pw,2*pi);
pw = pw-min(pw(:));
pw = 2*pi*pw./max(pw(:));

pw2 = mod(pw2,2*pi);
pw3 = mod(pw3,2*pi);
pw4 = mod(pw4,2*pi);
pw5 = mod(pw5,2*pi);

%Correct the piston
pw    = adjustPiston(pwRef,pw);
pw2   = adjustPiston(pwRef,pw2);
pw3   = adjustPiston(pwRef,pw3);
pw4   = adjustPiston(pwRef,pw4);
pw5   = adjustPiston(pwRef,pw5);

%Unwrap the phases
up=Unwrap_TIE_DCT_Iter(pw);
up2=Unwrap_TIE_DCT_Iter(pw2);
up3=Unwrap_TIE_DCT_Iter(pw3);
up4=Unwrap_TIE_DCT_Iter(pw4);
up5=Unwrap_TIE_DCT_Iter(pw5);

%Remove a possible remaining piston term
up  = up-mean2(up);
up2 = up2-mean2(up2);
up3 = up3-mean2(up3);
up4 = up4-mean2(up4);
up5 = up5-mean2(up5);

%% Visualization
%In the cases medida 5, 6 and 7, we do not have a reference phase so we cannot
%obtain rms errors

if (medida ~= 5) && (medida ~= 6) && (medida ~= 7)
    P   = P -mean2(P);
    rms  = sqrt(mean2((P-up).^2));
    rms2 = sqrt(mean2((P-up2).^2));
    rms3 = sqrt(mean2((P-up3).^2));
    rms4 = sqrt(mean2((P-up4).^2));
    rms5 = sqrt(mean2((P-up5).^2));

    figure,
    multi = cat(3,pwRef,pw,pw2,pw3);
    montage(multi,'Size', [2 2],'DisplayRange', [0 6]);
    colorbar
        
    figure,
    multi2 = cat(3,P,up,up2,up3);
    %montage(multi2,'Size', [2 2],'DisplayRange',[-20 35]);
    montage(multi2,'Size', [2 2],'DisplayRange',[-15 60]);
    colorbar
    colormap jet
   
    %UP - Kreis
    %UP2 - WPCA
    %UP4 - GS
    %UP5 - OF
    figure,
    %multi3 = cat(3,up,up2,up4,up5);
    %montage(multi3,'Size', [1 4],'DisplayRange',[-20 35]);
    multi3 = cat(3,pw,pw2,pw4,pw5);
    montage(multi3,'Size', [1 4],'DisplayRange',[0 2*pi]);
    colorbar
    colormap jet
    
    fprintf("\n");
    fprintf("RMS Kreis: %4.2f - RMS WPCA: %4.2f - RMW PCA: %4.2f - RMW GS: %4.2f - RMW OF: %4.2f\n",rms,rms2,rms3,rms4,rms5);
    fprintf("------------- \n");
    
    figure,
    montage(I(:,:,1:20),'Size', [4 5],'DisplayRange', [-1 1]);
    
    figure,
    plot(W2,'ko')
    grid on
    
else
    
    figure,
    multi = cat(3,pw,pw2,pw3);
    montage(multi,'Size', [1 3],'DisplayRange', [0 5]);
    colorbar
      
        
    if (medida == 7)
        dr = [-60 160];
    else
        dr  = [-20 35];
    end
    
    figure,
    multi2 = cat(3,up,up2,up3);
    montage(multi2,'Size', [1 3],'DisplayRange', dr);
    colorbar
    colormap jet
    
    %UP - Kreis
    %UP2 - WPCA
    %UP4 - GS
    %UP5 - OF
    figure,
    multi3 = cat(3,up,up2,up4,up5);
    %montage(multi3,'Size', [1 4],'DisplayRange',[-20 35]);
    montage(multi3,'Size', [1 4],'DisplayRange',[-80 120]);
    %multi3 = cat(3,pw,pw2,pw4,pw5);
    %montage(multi3,'Size', [1 4],'DisplayRange',[0 2*pi]);
    colorbar
    colormap jet
    
    figure,
    montage(I,'Size', [4 10],'DisplayRange', [0 220]);
    
    figure,
    plot(W2,'ko')
    grid on
    
    if ((medida ~= 7))
        [X Y]=meshgrid(1:1024);
        [S, Zzer, Zpol] = AjusteZernike(X, Y, up3, Mask, 1:15);
        [S2, Zzer2, Zpol2] = AjusteZernike(X, Y, up2, Mask, 1:15);
        figure, bar((abs(S)-abs(S2)))
    end
    
end