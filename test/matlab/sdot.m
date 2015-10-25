function [SDOT] = sdot(theta)

    SDOT(1,1) = -sin(theta(3));
    SDOT(2,1) = -cos(theta(3));
    SDOT(3,1) = 0.0;

    SDOT(1,2) = -sin(theta(1))*sin(theta(3)) + cos(theta(1))*cos(theta(3));
    SDOT(2,2) = -sin(theta(1))*cos(theta(3)) - cos(theta(1))*sin(theta(3));
    SDOT(3,2) = -cos(theta(1)) ;

    SDOT(1,3) = 0.0;
    SDOT(2,3) = 0.0;
    SDOT(3,3) = 0.0;
    
end