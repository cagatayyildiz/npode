function xN = normalizeExp(x)
xN = exp(x - max(x));
xN = xN / sum(xN);
end

