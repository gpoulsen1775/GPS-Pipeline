% Reads in the vehicle data and uses the satellite and receiver classes to
% translate the information and output hopefully the same vehicle
% coordinates (or at least close). To use, put all vehicle data in a
% vehicle.dat file in the same folder as this pipeline, as well as the
% satellite and receiver classes. Then simply run this script and the
% satellite.log and receiver.log files will be created automatically. 

% Note: if a satellite.log file already exists this will not overwrite it
% but instead add onto it so either change the file name or delete it
% before attempting to run this again

fidd = fopen('vehicle.dat', 'r');
C = fscanf(fidd, '%f');
fclose(fidd);
C = C';
j = 1;
vehicleDataStreams = length(C) / 10;
s = Satellite();

for k = 1:vehicleDataStreams
    vehValues = C(j:j+9);
    %Send the vehicle values to the satellite and interpret them
    satTempLog = s.sendData(vehValues);

    %Make a receiver from the satellite data and send the receiver data
    r = Receiver();
    r.sendData();
    %Process the satellite log and receiver log so we
    %can write them into their own files neatly
    command = 'cat satellite_temp.log >> satellite.log';
    status = system(command);
    command = 'cat receiver_temp.log >> receiver.log';
    status = system(command);
    
    %Move onto the next line of vehicle data
    j = j+10;
end

disp("Finished.");

            