model Tank
  Real x(start = 10);
  parameter Real a = 0.2;
  input Real u;
equation
  	der(x) = u-a*x;
end Tank;
