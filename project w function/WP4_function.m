function [Bit_rate,channel_gain_avg,L] = WP4_function(all_channel_gains,closest_sat_dists, CarrierFrequency)


f = CarrierFrequency;
c = 3e8;
L = zeros (1, length(closest_sat_dists));
G_gs = 1;
G_sat = 10^5;
N0 = 3.98e-21;
BW = 5e6;
N = length(closest_sat_dists);
P_Tx = 10e-3; % 10 mW

%free-space path loss
for i = 1: length(closest_sat_dists)

    L(i) = ((c)/(4*pi*closest_sat_dists(i,1)*f))^2; 

end

channel_gain_avg = zeros (1, length(all_channel_gains));

for j = 1 : length(all_channel_gains)

    channel_gain_avg(j) = mean(abs(all_channel_gains{j}))^2;
    
end

PRX = zeros (1, length(closest_sat_dists));

for z = 1 : length(closest_sat_dists) 

    PRX(z) = P_Tx * G_sat * G_gs * channel_gain_avg (z) * L(z);

end

                     
Bit_rate = 0;


for k = 1 : length(all_channel_gains)

    Bit_rate = Bit_rate + ( BW * (log2(1 + PRX(k)/(N0 * BW))));

end

Bit_rate = Bit_rate / N;

disp(Bit_rate);



end

