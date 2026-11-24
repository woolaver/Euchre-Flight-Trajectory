%test nose cone heating functions


M_inf = 8;
alt = 1;
alpha = 0;
T_nose = 297; %K

vehicleProperties = half_cone_properties();

Q_cone_heating = cone_heating(vehicleProperties.nose_STL, alpha, alt, M_inf, T_nose, vehicleProperties);
Q_fayRiddell = fayRiddellMethod(vehicleProperties, M_inf, alt, T_nose);

disp("Cone heating: " + Q_cone_heating)
disp("Fay Riddell: " + Q_fayRiddell)