function [signal, received] = signals(fc, slope, R, v, t)

c = 3e8;
delta_t = t - (R + t * v) / c;

signal = cos(2 * pi * (fc * t + slope * t^2 / 2));
received = cos(2 * pi * (fc * delta_t + slope * delta_t^2 / 2));
end