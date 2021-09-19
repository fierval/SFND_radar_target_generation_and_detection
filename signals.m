function [signal, received] = signals(fc, slope, R, v, t)

c = 3e8;
% delta_t = t(curr) - tau
% tau = send_t + receive_t
% tau = 2 * sent_t
% tau = (range_0 + v*t)/c
delta_t = t - 2 * (R + t * v) / c;

signal = cos(2 * pi * (fc * t + slope * t^2 / 2));
received = cos(2 * pi * (fc * delta_t + slope * delta_t^2 / 2));
end