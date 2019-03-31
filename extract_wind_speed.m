function v_wind = extract_wind_speed(t_n,data)

if t_n == data(1,1)
    v_wind = data(2,1);
elseif t_n == data(1,end)
    v_wind = data(2,end);
else
    index1 = max(find(data(1,:) < t_n));
    index2 = min(find(data(1,:) > t_n));
    x = [data(1,index1) data(1,index2)];
    v = [data(2,index1) data(2,index2)];
    
    v_wind = interp1(x,v,t_n);
end

end