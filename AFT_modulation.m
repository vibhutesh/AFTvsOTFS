function x_with_CP = AFT_modulation(N_AFT,Num_OFDM_sym, N_CP, c1, c2, X)
%% AFT Modulation:
Gamma_c1 = zeros(N_AFT);
Gamma_c2 = zeros(N_AFT);
x_with_CP = zeros(Num_OFDM_sym, N_AFT + N_CP);
for n = 0:N_AFT-1
    Gamma_c1(n+1, n+1) = exp(-1i*2*pi*(c1*n^2));
    Gamma_c2(n+1, n+1) = exp(-1i*2*pi*(c2*n^2));
end
for i=1:Num_OFDM_sym
    x                   = transpose(sqrt(N_AFT)*transpose(diag(Gamma_c1)').*ifft(transpose(diag(Gamma_c2)'.*(X((i-1)*N_AFT+1:i*(N_AFT)))), N_AFT));
    v                   = N_CP:-1:1;
    x_with_CP(i , :)    = [x(N_AFT-N_CP+1:N_AFT).*exp(-1i*2*pi*c1.*(N_AFT^2 - 2*N_AFT*(v))) x];%adding CP
end
end