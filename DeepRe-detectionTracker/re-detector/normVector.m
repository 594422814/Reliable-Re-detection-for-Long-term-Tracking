function  VV = normVector(V)
% normalize each column for the input V

n = size(V,2);
VV =zeros(size(V));
for i = 1 : n
    if norm(V(:,i))~=0
        k = norm(V(:,i));
        VV(:,i) = V(:,i)/k;
    else
        VV(:,i) = V(:,i);
    end
end