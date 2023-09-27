function sigObj = buildSignalObj(Data_Unfiltered, Times, Noise_Reference, ...
                                 SampleRate, Channels, ...
                                 Train_Time_Bounds, Preview_Time_Bounds, Analyze_Time_Bounds, ...
                                 Data_HPF, Data_LMS, Data_LMS_LPF, Data_BPF, ...
                                 Times_Train, Times_Test, Noise_Reference_Train, Noise_Reference_Test, ...
                                 Data_Unfiltered_Train, Data_Unfiltered_Test, ...
                                 Data_HPF_Train, Data_HPF_Test, Data_LMS_Train, Data_LMS_Test, ...
                                 Data_LMS_LPF_Train, Data_LMS_LPF_Test, Data_BPF_Train, Data_BPF_Test)

sigObj.Data_Ground_Truth = Data_Ground_Truth;

sigObj.Data_Unfiltered = Data_Unfiltered;
sigObj.Data_HPF = Data_HPF; 
sigObj.Data_BPF = Data_BPF; 
sigObj.Data_LMS = Data_LMS; 
sigObj.Data_LMS_LPF = Data_LMS_LPF; 

sigObj.Noise_Reference = Noise_Reference; 
sigObj.Times = Times;
sigObj.SampleRate = SampleRate; 
sigObj.Chanels = Channels;

sigObj.Train_Time_Bounds = Train_Time_Bounds; 
sigObj.Preview_Time_Bounds = Preview_Time_Bounds; 
sigObj.Analyze_Time_Bounds = Analyze_Time_Bounds;

sigObj.Times_Train = Times_Train; 
sigObj.Times_Test = Times_Test; 
sigObj.Noise_Reference_Train = Noise_Reference_Train; 
sigObj.Noise_Reference_Test = Noise_Reference_Test; 
sigObj.Data_Unfiltered_Train = Data_Unfiltered_Train; 
sigObj.Data_Unfiltered_Test = Data_Unfiltered_Test; 
sigObj.Data_HPF_Train = Data_HPF_Train; 
sigObj.Data_HPF_Test = Data_HPF_Test; 
sigObj.Data_BPF_Train = Data_BPF_Train; 
sigObj.Data_BPF_Test = Data_BPF_Test; 
sigObj.Data_LMS_Train = Data_LMS_Train; 
sigObj.Data_LMS_Test = Data_LMS_Test; 
sigObj.Data_LMS_LPF_Train = Data_LMS_LPF_Train; 
sigObj.Data_LMS_LPF_Test = Data_LMS_LPF_Test; 

end