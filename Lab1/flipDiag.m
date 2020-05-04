function flipped = flipDiag(matrix)
    flipped = matrix;
    %extracting diagonal and flipping
    d = flipud(diag(flipped));
    %replacing flipped diagonal
    for i = 1:length(flipped)
        flipped(i,i) = d(i);
    end
end