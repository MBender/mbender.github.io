%%%%%%
%% SylvesterMEP : Computes the eigenvalues of a multiparameter eigenvalue problem
%%
%% Input:
%%   + A : cell array of size alpha x (alpha+1) of matrices A{i,j}.
%%         Every matrix sqaure, and at every row, every matrix has the same dimension.
%%    
%% Assume: The last coordinate of every eigenvalue is different to zero.
%%         The input MEP is regular and non-singular.
%%    
%% Output:
%%   + sols : Coordinates of the eigenvalues of a multiparameter eigenvalue problem A
%%
%%%%%%
%% The following code generates a valid random example
%
% sizeMatrix = [2;2;2;2];
% dimEign = size(sizeMatrix,1);
% A = cell(dimEign,dimEign+1);
% for i=1:dimEign
%     for j=1:dimEign+1
%         A{i,j} = rand(sizeMatrix(i));
%     end
% end
% sols = sylvesterMEP(A);
%
%%%%%%
%% This file is a toy implementation of the techniques developed in the 2020 paper
%% "Koszul-type determinantal formulas for families of mixed multilinear systems"
%% By MR Bender, J-C Faug`ere, A Mantzaflaris, and E Tsigaridas.
%%%%%%


function sols = sylvesterMEP(A)

alpha = length(A)-1;
beta = zeros(alpha,1);

for i=1:alpha
    beta(i) = size(A{i,1},1)-1;
end

maxDims = cat(1,arrayfun(@(x) x+1, beta), alpha+1);

prodsMaxDims = cat(1,1,fold(@(lst,el) cat(1,lst,lst(end)*el), maxDims));

nsols = prodsMaxDims(alpha+1);

macMatMEP = zeros(prodsMaxDims(alpha+2) - prodsMaxDims(alpha+1),prodsMaxDims(alpha+2));
size(macMatMEP)
indRow = 1;
for i=1:alpha
    maxLwOff = prodsMaxDims(i)-1;
    maxUpOff = prodsMaxDims(alpha+1)/prodsMaxDims(i+1)-1;
    for indInBeta=1:beta(i)+1
        for lwOff=0:maxLwOff
            for upOff=0:maxUpOff
                %for each monomial in f[betaa[i], indBeta]
                offset = lwOff + upOff*prodsMaxDims(i+1);
                for yvar=0:beta(i)
                    for xvar=0:alpha
                        indCol = offset + prodsMaxDims(i)*yvar + prodsMaxDims(alpha+1)*xvar + 1; 
                        
                        macMatMEP(indRow,indCol) = A{i,xvar+1}(indInBeta,yvar+1);
                    end
                end
                indRow = indRow + 1;
            end
        end
    end
end


f0s = rand(alpha+1);
% f0s = eye(alpha+1); % change following for generic computation

macMatF0s = cell(alpha+1,1);

for indF0=1:alpha+1
  macMatF0s{indF0} = zeros(prodsMaxDims(alpha+1),prodsMaxDims(alpha+2));
for offCol=0:prodsMaxDims(alpha+1)-1
    for xvar=0:alpha
        indCol = offCol + xvar*prodsMaxDims(alpha+1) + 1;
        macMatF0s{indF0}(offCol+1,indCol) = f0s(indF0,xvar+1);
    end
end
end


AmB = (macMatMEP(:,1:prodsMaxDims(alpha+2) - prodsMaxDims(alpha+1))^-1)*macMatMEP(:,prodsMaxDims(alpha+2) - prodsMaxDims(alpha+1)+1:end);

multMap = cell(alpha+1,1);

for indF0=1:alpha+1
  multMap{indF0} =  macMatF0s{indF0}(:,prodsMaxDims(alpha+2) - prodsMaxDims(alpha+1)+1:end) -  macMatF0s{indF0}(:,1:prodsMaxDims(alpha+2) - prodsMaxDims(alpha+1))*AmB;
end


eigns = zeros(nsols,alpha+1);
[U,T] = schur(multMap{1},'complex');
eigns(:,1) = diag(T);
for indF0=2:alpha+1
    eigns(:,indF0) = diag(ctranspose(U)*multMap{indF0}*U);
end

sols = eigns*(transpose(f0s)^-1);

end 

