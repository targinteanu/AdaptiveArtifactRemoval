function chanObj = buildChannelObj(label, pos1, pos2, pos3, coordtype)
% Compatible with EEGlab 

chanObj.labels = label; % name

if nargin < 2
    % use EEGlab 
    % check whether EEGlab is there? Try starting it? Warn user to check path? 
    chanObj = pop_chanedit(chanObj)
else

    if nargin < 5
        % assume cartesian XYZ
        coordtype = 'Cartesian'; 
    end

    if strcmpi(coordtype, 'Spherical') | strcmpi(coordtype, 'Sphere')
        % EEGlab convention: 
        %   sph_theta about +z from +x, in degrees 
        %   phi from xy plane to +z, in degrees 
        %   theta in xy plane clockwise [about -z] from +x, in degrees
        rho = pos1; 
        sph_theta = pos2 * pi/180; % radians
        phi = pos3 * pi/180; % radians
        x = rho.*cos(phi).*cos(sph_theta); 
        y = rho.*cos(phi).*sin(sph_theta);
        z = rho.*sin(phi); 
        theta = -atan(y./x); % radians
        r = sqrt(x.^2 + y.^2);
        % radians to degrees: 
        sph_theta = sph_theta * 180/pi; 
        phi = phi * 180/pi;
        theta = theta * 180/pi;
    end
    if strcmpi(coordtype, 'Cartesian') | strcmpi(coordtype, 'Rectangular') | strcmpi(coordtype, 'XYZ')
        x = pos1; 
        y = pos2; 
        z = pos3;
        rho = sqrt(x.^2 + y.^2 + z.^2);
        r = sqrt(x.^2 + y.^2);
        theta = -atan(y./x); % EEGlab convention, radians
        % mathematical convention: 
        phi = acos(z./rho); % radians
        sph_theta = sign(y)*acos(x./r); % radians
        % radians to degrees: 
        sph_theta = sph_theta * 180/pi;
        phi = phi * 180/pi;
        phi = 90 - phi; % mathematical to EEGlab convention
        theta = theta * 180/pi; 
    end
    if strcmpi(coordtype, 'Polar') | strcmpi(coordtype, 'Cylindrical') | strcmpi(coordtype, 'Cylinder')
        r = pos1; 
        theta = pos2;
        z = pos3; 
        theta = -theta; % math convention 
        theta = pi/180 * theta; % radians 
        x = r.*cos(theta); 
        y = r.*sin(theta); 
        % mathematical convention: 
        phi = acos(z./rho); % radians
        sph_theta = sign(y)*acos(x./r); % radians
        % radians to degrees: 
        sph_theta = sph_theta * 180/pi;
        phi = phi * 180/pi;
        phi = 90 - phi; % mathematical to EEGlab convention
        theta = -theta * 180/pi; % mathematical ti EEGlab convention
    end

% ** polar radius (r) is not how EEGlab works. Unclear why. **

chanObj.X = x; % cartesian 
chanObj.Y = y; % cartesian 
chanObj.Z = z; % cartesian
chanObj.theta = theta; % polar angle 
chanObj.radius = r; % polar radius
chanObj.sph_theta = sph_theta; % spherical horizon angle
chanObj.sph_phi = phi; % spherical azimuth angle
chanObj.sph_radius = rho; % spherical radius

end
end