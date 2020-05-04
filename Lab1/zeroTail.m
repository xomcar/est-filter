function out = zeroTail(array)
    %selecting the number of item to set to zero
    itemsToZero = 3;
    %resetting
    array((end - itemsToZero + 1) : end) = 0;
    out = array;
end