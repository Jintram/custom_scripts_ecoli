function dts = general_check_delta_t(times)

dts = [];

for idx = 1:length(times)-1
    dts(end+1) = times(idx+1)-times(idx);    
end

for idx = 1:length(times)-1
    if times(idx+1)~=times(idx)
        disp('Warning: timesteps are not equal!');
        break;
    end
end

end