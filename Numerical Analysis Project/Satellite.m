classdef Satellite
    properties
        %Constant values to be read in by data.dat
        pi_data;
        c_light;
        R_earth;
        sidereal;
        
        %The values of the satellites to be read in by data.dat
        sat_u;
        sat_v;
        sat_per;
        sat_alt;
        sat_pha;
    end
    methods
        %Creates the Satellite class object with the constants and
        %satellite values
        function obj = Satellite()
            fid = fopen('data.dat','r');
            formatSpec = '%f %*[^\n]';
            A = textscan(fid,formatSpec);
            fclose(fid);
            vdata = A{1};

            obj.sat_u = zeros(24,3);
            obj.sat_v = zeros(24,3);
            obj.sat_per = zeros(24,1);
            obj.sat_alt = zeros(24,1);
            obj.sat_pha = zeros(24,1);

            obj.pi_data = vdata(1);
            obj.c_light = vdata(2);
            obj.R_earth = vdata(3);
            obj.sidereal = vdata(4);

            j = 4;
            for k=1:24
                obj.sat_u(k,:) = vdata(j+1:j+3);
                obj.sat_v(k,:) = vdata(j+4:j+6);
                obj.sat_per(k) = vdata(j+7);
                obj.sat_alt(k) = vdata(j+8);
                obj.sat_pha(k) = vdata(j+9);
                j = j+9;
            end
        end
        
        %A function to convert DMS to Radians
        function radians = convertToRadians(Satellite,degrees, minutes, seconds, sign)
            deg = degrees + (minutes/60) + (seconds/3600);
            radians = deg*(Satellite.pi_data)*sign/180;
        end
        
        %A function to convert from radians to DMS
        function DMS = convertFromRadians(Satellite,radians)
            DMS(4) = 1;
            if (radians < 0)
                DMS(4) = -1;
                radians = -radians;
            end
            decimalDeg = radians*180/Satellite.pi_data;
            DMS(1) = floor(decimalDeg);
            DMS(2) = floor((decimalDeg - DMS(1)) * 60);
            DMS(3) = round((decimalDeg - DMS(1) - (DMS(2) / 60)) * 3600, 2);
        end
        
        %A function to convert DMS coordinates in latitude and longitude to
        %cartesian coordinates
        function [x,y,z] = convertGivenToCartesian(Satellite, values)
            latRad = convertToRadians(Satellite, values(2), values(3), values(4), values(5));
            longRad = convertToRadians(Satellite, values(6), values(7), values (8), values(9));
            
            x = cos(latRad) * (Satellite.R_earth + values(10)) * cos(2*Satellite.pi_data/Satellite.sidereal*values(1) + longRad);
            y = cos(latRad) * (Satellite.R_earth + values(10)) * sin(2*Satellite.pi_data/Satellite.sidereal*values(1) + longRad);
            z = (Satellite.R_earth + values(10)) * sin(latRad);
        end
        
        %A function to convert cartesian coordinates at some time t to DMS
        %coordinates in longitude and latitude
        function values = convertCartesianToDMS(Satellite, xyzt)
            longRad = wrapTo2Pi(atan2(xyzt(2), xyzt(1)) - 2 * Satellite.pi_data / Satellite.sidereal * xyzt(4)); 
            xyztRadius = sqrt(xyzt(1)^2 + xyzt(2)^2 + xyzt(3)^2);
            altitude = round(xyztRadius - Satellite.R_earth, 2);
            latRad = asin(xyzt(3)/xyztRadius);
            
            latDMS = convertFromRadians(Satellite,latRad);
            longDMS = convertFromRadians(Satellite,longRad);
            
            values(1) = xyzt(4);
            values(2) = latDMS(1);
            values(3) = latDMS(2);
            values(4) = latDMS(3);
            values(5) = latDMS(4);
            values(6) = longDMS(1);
            values(7) = longDMS(2);
            values(8) = longDMS(3);
            values(9) = longDMS(4);
            values(10) = altitude;
        end
        
        %A function to find all the satellites above the horizon
        function sats = findSatellitesAboveTheHorizon(Satellite, values)
            sats = [];
            [vx, vy, vz] = convertGivenToCartesian(Satellite, values);
            vehPos = [vx, vy, vz];
            tv = values(1);
            for k=1:24
                satAtTv = calculateSatellitePosition(Satellite, k, tv);
                %Thinking about this as a plane
                %plane = vx*(satAtTv(1)-vx) + vy*(satAtTv(2)-vy) + vz*(satAtTv(3)-vz);
                
                %An alternate way of thinking about this
                plane = dot(vehPos, satAtTv) - dot(vehPos, vehPos);
                if (plane >= 1)
                    sats = [sats k];
                end 
            end
        end
        
        %A function to find when the satellite should send the signal to
        %the vehicle
        function ts = findTimeToSendSignal(Satellite, values, num)
            ts = newtonsMethodForTs(Satellite, values, num);
        end
        
        %A function for calculating the position of a satellite at some
        %time
        function xyz = calculateSatellitePosition(Satellite, num, time)
            xyz = (Satellite.R_earth + Satellite.sat_alt(num)) * ((Satellite.sat_u(num,:)*cos((2*Satellite.pi_data*time/(Satellite.sidereal/2)) + Satellite.sat_pha(num)))+ ((Satellite.sat_v(num,:) * sin(2*Satellite.pi_data*time/(Satellite.sidereal/2) + Satellite.sat_pha(num)))));
        end
        
        function xyz = calculateSatelliteDerivative(Satellite, num, time)
            p = Satellite.sidereal / 2;
            xyz = (Satellite.R_earth + Satellite.sat_alt(num)) * ((Satellite.sat_u(num,:)*-1*2*Satellite.pi_data/p*sin((2*Satellite.pi_data*time/(Satellite.sidereal/2)) + Satellite.sat_pha(num)))+ ((Satellite.sat_v(num,:) * 2*Satellite.pi_data/p*cos(2*Satellite.pi_data*time/(Satellite.sidereal/2) + Satellite.sat_pha(num)))));
        end
        
        %This just gives us an approximation, it's not as accurate as we
        %may want so we'll develop the acutal Newton's method
        function ts = alternateMethodForTs(Satellite, vehValues, num)
            tv = vehValues(1);
            [xv, yv, zv] = convertGivenToCartesian(Satellite, vehValues);
            veh = [xv yv zv];
            s_tv = calculateSatellitePosition(Satellite, num, tv);
            s_tv_der = calculateSatelliteDerivative(Satellite, num, tv);
            
            %For quadratic formula ax^2 + bx + c
            %Note that these values were found by expanding out
            %c^2*deltat^2 = norm(xs(tv)-xv(tv)-xs'(tv)*deltat)^2
            %and then proceeding to use the quadratic formula to solve for
            %delta t and then ts
            a = norm(s_tv_der)^2 - Satellite.c_light^2;
            b_x = 2*s_tv_der(1)*(xv-s_tv(1));
            b_y = 2*s_tv_der(2)*(yv-s_tv(2));
            b_z = 2*s_tv_der(3)*(zv-s_tv(3));
            b = b_x + b_y + b_z;
            c = norm(s_tv - veh)^2;
            
            %Finding deltaT using quadratic formula
            deltaT = zeros(2,1);
            d = sqrt(b^2 - 4*a*c);
            deltaT(1) = ( -b + d ) / (2*a);
            deltaT(2) = ( -b - d ) / (2*a);
            
            %We get two possible values for deltaT and thus two possible
            %values for ts, how 
            ts1 = tv - deltaT(1);
            ts2 = tv - deltaT(2);
            if ts1 > ts2
                ts = ts2;
            else
                ts = ts1;
            end
        end
        
        function value = fToFindRoot(Satellite, values, ts, num)
            [xv, yv, zv] = convertGivenToCartesian(Satellite, values);
            veh = [xv yv zv];
            xs = calculateSatellitePosition(Satellite, num, ts);
            value = norm(veh - xs) - Satellite.c_light*(values(1) - ts);
        end
        
        function deriv = approxDerivForF(Satellite, values, ts, num)
            delta  = 1.0 * 10^-4;
            fx = fToFindRoot(Satellite, values, ts, num);
            fdelta = fToFindRoot(Satellite, values, ts+delta, num);
            deriv = (fdelta - fx) / delta;
        end
        
        function approx = newtonsMethodForTs(Satellite, values, num)
            t0 = alternateMethodForTs(Satellite, values, num);
            tolerance = 1.0*10^(-11);
            Max_iterations = 1000;
            for i=1:Max_iterations
                funcT = fToFindRoot(Satellite, values, t0, num);
                funcDeriv = approxDerivForF(Satellite, values, t0, num);
                t = t0 - funcT/funcDeriv;
                if abs(t - t0) < tolerance
                    approx = t;
                    return;
                else
                    t0 = t;
                    continue;
                end
            end
            %Failed to find a root
            approx = 0;
        end
        
        function loggedData = sendData(Satellite, values)
            sats = findSatellitesAboveTheHorizon(Satellite, values);
            fileID = fopen('satellite_temp.log','w');
            fprintf(fileID, '%.11f %d %.2f %.2f %d %d %.2f %.2f %d %.2f\n', values);
            loggedData = num2str(values, '%.11f %d %.2f %.2f %d %d %.2f %.2f %d %.2f\n');
            loggedData = split(loggedData)'; 
            for sat = 1:length(sats)
                ts = findTimeToSendSignal(Satellite, values, sats(sat));
                xs = calculateSatellitePosition(Satellite, sats(sat), ts);
                satData = [sats(sat); ts; xs';]';
                satString = num2str(satData, '%d %.11f %.2f %.2f %.2f\n');
                satString = split(satString)';
                loggedData = [loggedData satString];
                fprintf(fileID, '%d %.11f %.2f %.2f %.2f\n', satData);
            end
            fclose(fileID);
        end
        
        
    end
end

