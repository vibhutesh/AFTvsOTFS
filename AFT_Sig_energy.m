function sig_energy = AFT_Sig_energy(N_AFT, Num_OFDM_sym, taps, delay_taps, Doppler_taps, chan_coef,x_with_CP)
%% wireless channel and noise
N_CP = max(delay_taps);
s_chan_mat = zeros(size(x_with_CP));
for i = 1:Num_OFDM_sym
    for itao = 1:taps
        s_chan_mat(i, :) = s_chan_mat(i, :)+exp(1i*2*pi*Doppler_taps(itao)*delay_taps(itao))*chan_coef(itao)*exp(-1j*2*pi ...
            *(-N_CP:N_AFT-1)*Doppler_taps(itao)).*circshift(x_with_CP(i, :) ,delay_taps(itao), 2);
    end
end
s_chan = reshape(s_chan_mat, 1, []);
sig_energy = sum(abs(s_chan).^2)/length(s_chan);
end