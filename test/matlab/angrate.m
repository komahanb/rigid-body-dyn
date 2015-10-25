function [SIB] = angrate(theta)

    SIB(1,1) = cos(theta(3));
    SIB(2,1) = -sin(theta(3));
    SIB(3,1) = 0.0;

    SIB(1,2) = cos(theta(1))*sin(theta(3));
    SIB(2,2) = cos(theta(1))*cos(theta(3));
    SIB(3,2) = -sin(theta(1)) ;

    SIB(1,3) = 0.0;
    SIB(2,3) = 0.0;
    SIB(3,3) = 1.0;
    
end