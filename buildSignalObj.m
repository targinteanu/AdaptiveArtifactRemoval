function sigObj = buildSignalObj(Data_Unfiltered, Times, Noise_Reference, ...
                                 SampleRate, Channels, ...
                                 Train_Time_Bounds, Preview_Time_Bounds, Analyze_Time_Bounds, ...
                                 Data_HPF, Data_HPF_LMS, Data_HPF_LMS_LPF, Data_BPF)

sigObj.Data_Unfiltered = Data_Unfiltered;
sigObj.Data_HPF = Data_HPF; 
sigObj.Data_BPF = Data_BPF; 
sigObj.Data_HPF_LMS = Data_HPF_LMS; 
sigObj.Data_HPF_LMS_LPF = Data_HPF_LMS_LPF; 

sigObj.Noise_Reference = Noise_Reference; 
sigObj.Times = Times;
sigObj.SampleRate = SampleRate; 
sigObj.Chanels = Channels;

sigObj.Train_Time_Bounds = Train_Time_Bounds; 
sigObj.Preview_Time_Bounds = Preview_Time_Bounds; 
sigObj.Analyze_Time_Bounds = Analyze_Time_Bounds;

end