function [c0, c1, c2] = ComputeC0_C1_for2path(Doppler_taps, delay_taps)
f0 = Doppler_taps(1);
f1 = Doppler_taps(2);
l0 = delay_taps(1);
l1 = delay_taps(2);
c0 = (f0*l1 - f1*l0)/(l1 - l0);
c1 = (f0-f1)/(2*(l1 - l0));
c2 = 0;
end