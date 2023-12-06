function filtObj = buildFilterObj(HPF, LPF, Num_Taps, Step_Size, useFiltFilt, ...
    Starting_Weights, Current_Weights)

if nargin < 5
    useFiltFilt = false;
end
if nargin < 6
    Starting_Weights = [];
end
if nargin < 7
    Current_Weights = Starting_Weights;
end

filtObj.useFiltFilt = useFiltFilt;

filtObj.HPF = HPF; 
filtObj.LPF = LPF; 

filtObj.Num_Taps = Num_Taps;
filtObj.Step_Size = Step_Size;

filtObj.Starting_Weights = Starting_Weights; 
filtObj.Current_Weights = Current_Weights;

end