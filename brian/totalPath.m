totalTime = 1200;
totPath = 0;

for i = 1:numel(ratData_pos_timeStamps)-1
    currPath = sqrt((ratData_pos_x(i) - ratData_pos_x(i+1))^2 + ...
        (ratData_pos_y(i) - ratData_pos_y(i+1))^2);
    totPath = totPath + currPath;
end

totPath
vel = totPath/totalTime