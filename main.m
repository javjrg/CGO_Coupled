%% CGO solutions (vectorial case)
clear all;
%% parameters
% choose the type of example
% example1: A=I*sigma.
% example2: A=[2 0;0 3]*sigma.
% example3: A=[2 1;1 3]*sigma
% example4: A=[sigma1 sigma2;sigma3 sigma4]
example = 'example0';

s   =   2.1;  %square [-s,s)
m   =   8;
M   =   2^m;
h   =   (2*s)/M;
k   =   2;
% Create Grid
[x1,x2] =   crea_grid(m,h);
z       =   complex(x1,x2);
%% functions eta, g and b (Cauchy and Beurling Transforms)
eta=crea_eta(x1,x2,M,2);
g=eta./(pi*z);
g(abs(z)==0)=0;
b=eta./(pi*z.^2);
b(abs(z)==0)=0;
% Fourier transform
Fg=fft2(fftshift(g));
Fb=fft2(fftshift(b));
%% Sigma y mu
if example == 'example0'
    disp('Example 0.')
    sigma_1   =   0.5*sigma5(z);
    sigma_2   =   0.2*sigma6(z);
    sigma_3   =   0.3*sigma7(z);
    sigma_4   =   0.5*sigma8(z);
end 
if example == 'example1'
    disp('Example 1.')
    sigma_1   =   sigma1(z);
    sigma_2   =   0*sigma1(z);
    sigma_3   =   0*sigma1(z);
    sigma_4   =   sigma1(z);
end 
if example == 'example2'
    disp('Example 2.')
    sigma_1   =   2*sigma1(z);
    sigma_2   =   0*sigma1(z);
    sigma_3   =   0*sigma1(z);
    sigma_4   =   3*sigma1(z);
end 
if example == 'example3'
    disp('Example 3.')
    sigma_1   =   2*sigma1(z);
    sigma_2   =   1*sigma1(z);
    sigma_3   =   1*sigma1(z);
    sigma_4   =   3*sigma1(z);
end 
if example == 'example4'
    disp('Example 4.')
    sigma_1   =   sigma1(z);
    sigma_2   =   sigma2(z);
    sigma_3   =   sigma3(z);
    sigma_4   =   sigma4(z);
end 
sigma_2(abs(z)>1)    =   0;
sigma_3(abs(z)>1)    =   0;
% SIGMA MATRIX
SIGMA   =   [sigma_1 sigma_2; sigma_3 sigma_4];
% Create MU matrix  (MU=(I-SIGMA)(I+SIGMA)^{-1})
% Determinant of Sigma
detS   =   ((1+sigma_1).*(1+sigma_4))-(sigma_2.*sigma_3);
mu1     =   ((1-sigma_1).*(1+sigma_4)-(sigma_2.*sigma_3))./detS;
mu2     =   (2*(sigma_1.*sigma_2))./detS;
mu3     =   (2*(sigma_3.*sigma_4))./detS;
mu4     =   ((1+sigma_1).*(1-sigma_4)-(sigma_2.*sigma_3))./detS;
MU      =   [mu1 mu2;mu3 mu4];
% Plot Sigma and MU
sigma_1tmp    =   sigma_1;
sigma_1tmp(abs(z)>1)  =   nan;
f1=figure(1);
surf(x1,x2,sigma_1tmp);
saveas(f1,'example5_sigma_1','png')
sigma_2tmp    =   sigma_2;
sigma_2tmp(abs(z)>1)  =   nan;
f2=figure(2);
surf(x1,x2,sigma_2tmp);
saveas(f2,'example5_sigma_2','png')
sigma_3tmp    =   sigma_3;
sigma_3tmp(abs(z)>1)  =   nan;
f3=figure(3);
surf(x1,x2,sigma_3tmp);
saveas(f3,'example5_sigma_3','png')
sigma_4tmp    =   sigma_4;
sigma_4tmp(abs(z)>1)  =   nan;
f4=figure(4);
surf(x1,x2,sigma_4tmp);
saveas(f4,'example5_sigma_4','png')
%% alpha and nu
E=exp(-1i*((k*z) + (conj(k)*conj(z))));
a1  =   -1i*conj(k)*E.*mu1;
a2  =   -1i*conj(k)*E.*mu2;
a3  =   -1i*conj(k)*E.*mu3;
a4  =   -1i*conj(k)*E.*mu4;
nu1 =   E.*mu1;
nu2 =   E.*mu2;
nu3 =   E.*mu3;
nu4 =   E.*mu4;
alpha   =   [a1 a2;a3 a4];
nu      =   [nu1 nu2;nu3 nu4];
A=[a1(:)+a2(:);a3(:)+a4(:)];
%% Solve for V from   (I-A*conj(A))V=-conj(alpha) using GMRES
v = gmres('operator', -conj(A), 50, 1e-6, 500, [], [], -conj(A), Fg, Fb,h,nu1,nu2,nu3,nu4,a1,a2,a3,a4,M);
%% Calculate U=(I-A*rho)V
vbar    =   conj(v);
v1      =   v(1:M^2);
v2      =   v((M^2)+1:2*M^2);
v1tmp   =   reshape(v1,M,M);
v2tmp   =   reshape(v2,M,M);
vtmp    =   [v1;v2];
Sv1   =   (h^2*(ifft2(Fb.*fft2(conj(v1tmp)))));
Sv2   =   (h^2*(ifft2(Fb.*fft2(conj(v2tmp)))));
Pv1   =   (h^2*(ifft2(Fg.*fft2(conj(v1tmp)))));
Pv2   =   (h^2*(ifft2(Fg.*fft2(conj(v2tmp)))));
Av1    =  (-conj(nu1).*Sv1-conj(nu2).*Sv2)-(conj(a1).*Pv1+conj(a2).*Pv2);
Av2    =  (-conj(nu3).*Sv1-conj(nu4).*Sv2)-(conj(a3).*Pv1+conj(a4).*Pv2);
Av      =   [Av1(:);Av2(:)];
u       =   vtmp-Av;
%%  Compute N = - P*conj(U)
ubar    =   conj(u);
u1      =   u(1:M^2);
u2      =   u((M^2)+1:2*M^2);
u1tmp   =   reshape(u1,M,M);
u2tmp   =   reshape(u2,M,M);
Pu1   =   (h^2*(ifft2(Fg.*fft2(conj(u1tmp)))));
Pu2   =   (h^2*(ifft2(Fg.*fft2(conj(u2tmp)))));
N1  =   -Pu1;
N2  =   -Pu2;
%% Plot real(N)
N1tmp    =   real(N1);
N1tmp(abs(z)>1)  =   nan;
f5=figure(5);
surf(x1,x2,N1tmp);
saveas(f5,'example5_N_1_m8','png')
N2tmp    =   real(N2);
N2tmp(abs(z)>1)  =   nan;
f6=figure(6);
surf(x1,x2,N2tmp);
saveas(f6,'example5_N_2_m8','png')
%% functions
function [x1,x2] = crea_grid(m,h)
j1=-2^(m-1):2^(m-1)-1;
j2=-2^(m-1):2^(m-1)-1;
hj1=h*j1;
hj2=h*j2;
[x1,x2]=meshgrid(hj1,hj2);
end
function eta = crea_eta(x1,x2,M,s)
etatemp=zeros(M,M);
for i=1:M
for j=1:M      
        if sqrt(x1(1,i)^2+x2(j,1)^2)<2
        etatemp(i,j)=1;
        end
        if 2<=sqrt(x1(1,i)^2+x2(j,1)^2) && sqrt(x1(1,i)^2+x2(j,1)^2)<2+((s-2)/2)
        etatemp(i,j)=(-2/(s-2))*(sqrt(x1(1,i)^2+x2(j,1)^2)-2)+1;
        end
        if 2+((s-2)/2)<=sqrt(x1(1,i)^2+x2(j,1)^2)
        etatemp(i,j)=0;
        end      
end
end
 eta = etatemp;
end
%% Sigma
function result = sigma1(z)
% Conductivities of heart and lung (background is 1)
heart = 2;	
lung  = 0.7;
% Initialize
[zfila,zcol] = size(z);
z           = z(:);
result         = ones(size(z));
x1          = real(z);
x2          = imag(z);
% Build coarse representation of heart. Planar point (hc1,hc2) is the center of the ellipse
% describing the heart; numbers he1 and he2 give the eccentrities with respect to radius hR.
hc1 = -.1;
hc2 = .4;
he1 = .8;
he2 = 1;
hR  = .2;
% Compute elliptical "distance" of the evaluation points from heart
hd  = sqrt(he1*(x1-hc1).^2 + he2*(x2-hc2).^2);
% Set value of conductivity inside the heart 
result(hd <= hR) = heart;
% Build coarse representation of two lungs
l1c1  = .5;
l1c2  = 0;
l1e1  = 3;
l1e2  = 1;
l1R   = .5;
fii   = -pi/7;
rot11 = cos(fii);
rot12 = sin(fii);
rot21 = -sin(fii);
rot22 = cos(fii);
l1d   = sqrt(l1e1*((rot11*x1+rot12*x2)-l1c1).^2 + l1e2*((rot21*x1+rot22*x2)-l1c2).^2);
result(l1d <= l1R) = lung;
l2c1 = -.6;
l2c2 = 0;
l2e1 = 3;
l2e2 = 1;
l2R  = .4;
fii   = pi/7;
rot11 = cos(fii);
rot12 = sin(fii);
rot21 = -sin(fii);
rot22 = cos(fii);
l2d  = sqrt(l2e1*((rot11*x1+rot12*x2)-l2c1).^2 + l2e2*((rot21*x1+rot22*x2)-l2c2).^2);
result(l2d <= l2R) = lung;
 result = reshape(result,[zfila,zcol]);
end
