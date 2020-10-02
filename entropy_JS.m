function h = entropy_JS(prob)
prob = prob(:);
probnew = prob(find(prob~=0));

h = -sum(probnew.*log(probnew));