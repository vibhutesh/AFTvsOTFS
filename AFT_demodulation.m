function y = AFT_demodulation(N_AFT,Num_OFDM_sym, c0, c1, c2,r_AFT)
%% AFT demodulation:
Gamma_c0 = zeros(N_AFT);
Gamma_c1 = zeros(N_AFT);
Gamma_c2 = zeros(N_AFT);
y = zeros(size(r_AFT));
for n = 0:N_AFT-1
    Gamma_c0(n+1, n+1) = exp(1i*2*pi*(c0*n));
    Gamma_c1(n+1, n+1) = exp(-1i*2*pi*(c1*n^2));
    Gamma_c2(n+1, n+1) = exp(-1i*2*pi*(c2*n^2));
end
for i=1:Num_OFDM_sym
    y(i, :) = (1/sqrt(N_AFT))*transpose((diag(Gamma_c2)).*fft(transpose(transpose(diag(Gamma_c1*Gamma_c0)).*(r_AFT(i, :))), N_AFT));
end
end