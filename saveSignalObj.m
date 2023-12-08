function saveSignalObj(filename, sigObj)
% saves signal object to sigObj.mat file, which can later be read by the
% script loadSignalObj

Data_Ground_Truth = sigObj.Data_Ground_Truth;

Data_Unfiltered = sigObj.Data_Unfiltered;
Data_HPF = sigObj.Data_HPF; 
Data_BPF = sigObj.Data_BPF; 
Data_LMS = sigObj.Data_LMS; 
Data_LMS_LPF = sigObj.Data_LMS_LPF; 

Noise_Reference = sigObj.Noise_Reference; 
Times = sigObj.Times;
SampleRate = sigObj.SampleRate; 
Channels = sigObj.Channels;

Train_Time_Bounds = sigObj.Train_Time_Bounds; 
Preview_Time_Bounds = sigObj.Preview_Time_Bounds; 
%Analyze_Time_Bounds = sigObj.Analyze_Time_Bounds;
Preview_Channel_Index = sigObj.Preview_Channel_Index;

Times_Train = sigObj.Times_Train; 
Times_Test = sigObj.Times_Test; 
Noise_Reference_Train = sigObj.Noise_Reference_Train; 
Noise_Reference_Test = sigObj.Noise_Reference_Test; 
Data_Unfiltered_Train = sigObj.Data_Unfiltered_Train; 
Data_Unfiltered_Test = sigObj.Data_Unfiltered_Test; 
Data_HPF_Train = sigObj.Data_HPF_Train; 
Data_HPF_Test = sigObj.Data_HPF_Test; 
Data_BPF_Train = sigObj.Data_BPF_Train; 
Data_BPF_Test = sigObj.Data_BPF_Test; 
Data_LMS_Train = sigObj.Data_LMS_Train; 
Data_LMS_Test = sigObj.Data_LMS_Test; 
Data_LMS_LPF_Train = sigObj.Data_LMS_LPF_Train; 
Data_LMS_LPF_Test = sigObj.Data_LMS_LPF_Test; 

clear sigObj;
save([filename,'_sigObj.mat']);

end