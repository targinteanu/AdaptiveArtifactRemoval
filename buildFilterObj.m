function filtObj = buildFilterObj(HPF, LPF, Num_Taps, Step_Size)

filtObj.HPF = HPF; 
filtObj.LPF = LPF; 

filtObj.Num_Taps = Num_Taps;
filtObj.Step_Size = Step_Size;

filtObj.Starting_Weights = Starting_Weights; 
filtObj.Current_Weights = Current_Weights;

end