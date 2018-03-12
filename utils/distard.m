
function R = distard(X1,X2,ells)

	[N1,D] = size(X1);
	N2 = size(X2,1);

	R = zeros(N1,N2);
	for i=1:D
		R = R + (X1(:,i) - X2(:,i)').^2 / ells(i)^2;
	end
	R = sqrt(R);
end



