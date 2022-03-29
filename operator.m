% Implements the operator
% v |->  v-Aconj(A)v,
% where A=(-conj(nu)S-conj(alpha)P).
function result = operator(v,Fg,Fb,h,nu1,nu2,nu3,nu4,a1,a2,a3,a4,M)
v       =   v(:);
v1      =   v(1:M^2);
v2      =   v((M^2)+1:2*M^2);
vtmp    =   [v1;v2];
v1tmp   =   reshape(v1,M,M);
v2tmp   =   reshape(v2,M,M);
Sconjv1   =   (h^2*(ifft2(Fb.*fft2(conj(v1tmp)))));
Sconjv2   =   (h^2*(ifft2(Fb.*fft2(conj(v2tmp)))));
Pconjv1   =   (h^2*(ifft2(Fg.*fft2(conj(v1tmp)))));
Pconjv2   =   (h^2*(ifft2(Fg.*fft2(conj(v2tmp)))));
A1    =  (-conj(nu1).*Sconjv1-conj(nu2).*Sconjv2)-(conj(a1).*Pconjv1+conj(a2).*Pconjv2);
A2    =  (-conj(nu3).*Sconjv1-conj(nu4).*Sconjv2)-(conj(a3).*Pconjv1+conj(a4).*Pconjv2);
A1bar   =   conj(A1);
A2bar   =   conj(A2);
SA1   =   (h^2*(ifft2(Fb.*fft2(A1bar))));
SA2   =   (h^2*(ifft2(Fb.*fft2(A2bar))));
PA1   =   (h^2*(ifft2(Fg.*fft2(A1bar))));
PA2   =   (h^2*(ifft2(Fg.*fft2(A2bar))));
AA1    =  (-conj(nu1).*SA1-conj(nu2).*SA2)-(conj(a1).*PA1+conj(a2).*PA2);
AA2    =  (-conj(nu3).*SA1-conj(nu4).*SA2)-(conj(a3).*PA1+conj(a4).*PA2);
AA1    =    AA1(:);
AA2    =    AA2(:);
A      =    [AA1;AA2];
result  =	vtmp-A;
