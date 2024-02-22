function filtObj = buildFilterObj(preFilt, postFilt, Num_Taps, Step_Size, ...
    preFiltDir, postFiltDir, useDLMS, useNLMS, ...
    Starting_Weights, Current_Weights)

if nargin < 5
    preFiltDir = ones(size(preFilt)); 
    postFiltDir = ones(size(postFilt));
end
if nargin < 7
    useDLMS = false;
end
if nargin < 8
    useNLMS = false;
end
if nargin < 9
    Starting_Weights = [];
end
if nargin < 10
    Current_Weights = Starting_Weights;
end

filtObj.preFiltDir = preFiltDir;
filtObj.postFiltDir = postFiltDir;
filtObj.useDLMS = useDLMS;
filtObj.useNLMS = useNLMS;

filtObj.preFilt = preFilt; 
filtObj.postFilt = postFilt; 

filtObj.Num_Taps = Num_Taps;
filtObj.Step_Size = Step_Size;

filtObj.Starting_Weights = Starting_Weights; 
filtObj.Current_Weights = Current_Weights;

end