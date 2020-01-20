function [edgePot,edgeStruct]=CreateGridUGMModel(NumFils, NumCols, K, lambda)
%
%
% NumFils, NumCols: image dimension
% K: number of states
% lambda: smoothing factor
tic

%Vector counter
idx=1;
%Fill in sparse matrix
for i=1:NumFils
    for j=1:NumCols
        p = (j-1)*NumFils+i;
        if i~=1
            %Adjacency with north pixel
            idx_Ai(idx)= p; 
            idx_Aj(idx) = p-1; 
            a_ij(idx) = 1;
            idx=idx+1; 
        end
        if i~=NumFils
            %Adjacency with south pixel
            idx_Ai(idx)= p; 
            idx_Aj(idx) = p+1; 
            a_ij(idx) = 1;
            idx=idx+1; 
        end
        if j~=1
            %Adjacency with west pixel
            idx_Ai(idx)= p; 
            idx_Aj(idx) = p-NumFils; 
            a_ij(idx) = 1;
            idx=idx+1;
        end
        if j~=NumCols
            %Adjacency with east pixel
            idx_Ai(idx)= p; 
            idx_Aj(idx) = p+NumFils; 
            a_ij(idx) = 1;
            idx=idx+1;
        end
        
    end
end

%Create edge structure
nNodes = NumFils * NumCols;
adj = sparse(idx_Ai, idx_Aj, a_ij, nNodes, nNodes); % 
edgeStruct = UGM_makeEdgeStruct(adj,K);

%Create the edge potential structure
edgePot = zeros(K,K,edgeStruct.nEdges);
for e = 1:edgeStruct.nEdges
   n1 = edgeStruct.edgeEnds(e,1);
   n2 = edgeStruct.edgeEnds(e,2);
   
   pot_same = exp(lambda(2));
   pot_dif = exp(lambda(1));
   
   arr_diag = eye(K);
   arr_outer = ones(K)-arr_diag;
   arr_diag = arr_diag*pot_same;
   arr_outer = arr_outer*pot_dif;
   edgePot(:,:,e) = arr_diag+arr_outer;   
end

toc;