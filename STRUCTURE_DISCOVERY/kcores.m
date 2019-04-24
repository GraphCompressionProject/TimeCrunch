function [A,fullDelSet] = kcores(A, k)
%A = spconvert(inp_list);
%sizeofA = size(A);
%nnodes = sizeofA(1);
%if sizeofA(1) ~= sizeofA(2)
%	nnodes = max(sizeofA);
%	A(nnodes,nnodes) = 0;
%end

%A = A + A';
%[i,j,k] = find(A);
%A(i(find(k==2)),j(find(k==2))) = 1;
%A = A - diag(diag(A));

nnodes = size(A);
nnodes = nnodes(1);
changed_flag = 1;
fullDelSet = find(sum(A,2) == 0)';
while changed_flag == 1
	changed_flag = 0;
	delSet = [];
    	degs = sum(A,2);
	for i = 1:nnodes
		if degs(i) < k && degs(i) > 0
			delSet = [delSet i];
			changed_flag = 1;
		end
	end
	fullDelSet = [fullDelSet delSet];
	nDel = numel(delSet);
	for j = 1:nDel
		A(delSet(j),:) = 0;
		A(:,delSet(j)) = 0;
	end
end
fullDelSet = sort(unique(fullDelSet));

%[src, dest, val] = find(A);
%decomp = [src dest val];
end
