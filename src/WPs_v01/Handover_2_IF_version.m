function [closest_sat_elevations_discrete_2] = Handover_2_IF_version(closest_sat_elevations_discrete,closest_sat_indices, Handover_Type)

closest_sat_elevations_discrete_2 = zeros (length(closest_sat_elevations_discrete),1);

sat_idx = closest_sat_indices(1,1);
sat_idx_optimal = closest_sat_indices(1,1);
closest_sat_elevations_discrete_2 (1) = closest_sat_elevations_discrete(1,1);
Handover_counter = 0;
HO_opt_counter = 0;

if Handover_Type == 0

    for i = 2 : length(closest_sat_elevations_discrete)
    if sat_idx_optimal ~= closest_sat_indices(i,1)
        
        HO_opt_counter = HO_opt_counter + 1;
        sat_idx_optimal = closest_sat_indices(i,1);    
    
    end
    end
    disp(HO_opt_counter);
    disp("type 0");

else

    for i = 2 : length(closest_sat_elevations_discrete)
            
        col = find (closest_sat_indices(i,:) == sat_idx);
    
        if ~isempty(col) && closest_sat_elevations_discrete(i,col) ~= 0
    
            closest_sat_elevations_discrete_2 (i) = closest_sat_elevations_discrete(i,col);
    
        else
            closest_sat_elevations_discrete_2 (i) = closest_sat_elevations_discrete(i,1);
            sat_idx = closest_sat_indices(i,1);
            Handover_counter = Handover_counter + 1;
    
        end
    end

    disp(Handover_counter);
    disp("type 1")
end

end