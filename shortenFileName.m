function shortFileName = shortenFileName(longFileName, termchar)

if nargin < 2
    termchar = '_';
end

lastIdx = find(longFileName == termchar); 

if isempty(lastIdx)
    shortFileName = longFileName;
else
    lastIdx = max(1, lastIdx(1));
    shortFileName = longFileName(1:(lastIdx-1));
end

end