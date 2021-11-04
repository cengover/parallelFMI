class Tank
  Real x(start = 10);
  parameter Real a = 0.2;
  input Real u;
equation
  	der(x) = u-a*x;
end Tank;

class TanksInSeries
	parameter Integer n;
	Tank Tanks[n];  
equation
	for i in 1:n loop
		if i == 1 then 
		Tanks[i].u = 1.0; 
		else 
		Tanks[i].u = Tanks[i-1].x*Tanks[i-1].a;
		end if;
	end for;
end TanksInSeries;

model SerialTanks
	TanksInSeries TT(n=1);
end SerialTanks;
