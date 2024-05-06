#To run the code for energy calibration from the ASPET .bin

#Go to the directory

#Create folders: EBOX ; EBOX/Spectrum ; EnergyCal ; EnergyCal/fitting

**********First part- Analyzing peak from the energy spetrum of each channel
#complie ReadBinary file

>g++ BinaryRead9-LAB221024_HighestPeak_PickUp_v1.cpp -o BinaryRead9-LAB221024_HighestPeak_PickUp_v1 `root-config --cflags --glibs` -lSpectrum

#This program will read the spectrum of each channel, pick up the highest and clearest energy peak (not 1.2MeV of Na22 source) and give information of 
#gmslID / sticID / Channel / Peak_area / ToT_ADCvalues / Resolution_ADC / and real Energy value --> in a [input]_HighestPeak.txt file, in the same directory with the [input].bin file 
# The program also dispaly the energy spectra of all channel and the corresponding peak pickup and fitting
# the 2D map of the peak area for each channel is also displayed, total count of each channel, peak area from fitting, peak area of good channels only, and resolution of peak fitting for each channel.

#To execute the program run; there are examples below with 6 arguments 

>./BinaryRead9-LAB221024_HighestPeak_PickUp_v1 Data_pointsource/101740.889_C1BG 100 800 0.6 7 306.82
>./BinaryRead9-LAB221024_HighestPeak_PickUp_v1 Data_pointsource/173228.492_SET3_HV1500-1600_T38-36-37-37C_BG_60S 100 800 0.6 7 306.82
>./BinaryRead9-LAB221024_HighestPeak_PickUp_v1 Data_pointsource/100104.298_A3-2MINS 100 800 0.4 7 511.0
>./BinaryRead9-LAB221024_HighestPeak_PickUp_v1 Data_pointsource/091448.906_A1-2MINS 100 800 0.4 7 511.0
>./BinaryRead9-LAB221024_HighestPeak_PickUp_v1 Data_pointsource/163920.920_SET3_NA22_D21CM_CENTER_5CMPHANTOM_60S 100 800 0.6 7 511.0
>./BinaryRead9-LAB221024_HighestPeak_PickUp_v1 Data_pointsource/174133.387_SET3_HV1500-1600_T38-39-37-38C_NA22-D6CM-CENTER_L-52.5MM_60S 100 800 0.6 7 511.0
>./BinaryRead9-LAB221024_HighestPeak_PickUp_v1 Data_pointsource/182452.671_SET3_HV1500-1600_T38-39-37-38C_CS137-D6CM-CENTER_L-52.5MM_60S 100 800 0.6 7 661.657
>./BinaryRead9-LAB221024_HighestPeak_PickUp_v1 Data_pointsource/181821.792_SET3_HV1500-1600_T38-39-37-38C_CO60-D6CM-CENTER_L-52.5MM_300S 100 800 0.1 3 1332.492

# argv[1]: file name of the input binary file (.bin) from ASPET measurement
# argv[2], argv[3]: minimum and maximum ToT_ADC values of the spectrum for all channels
# argv[4]: peak finding parameter; recommendation:  0.6 for 306keV, 661keV and 511keV from Na22; 0.4 for 511 keV from PAG measurement; 0.1 for 1.3MeV of Co60
# argv[5]: Background finding parameter; recommendation: 7 for 306keV, 661 keV and 511keV; 3 for 1.3MeV 
# argv[6]: Real energy value (highest and cleanest) of the calibration source (keV) 
# argv[7], argv[8] (OPTIONAL): minimum and maximum Peak_area values of the spectrum for all channels to filter the peak fitting (default is 100 & 10e+5 for each channel) 
#{If the .bin is not exceptionally large (>10e+5 count in the peak of each channel); only 6 arguments are nessessary}

#An alternative program to read 2 peaks (2 highest peaks) from the spectrum 

>g++ BinaryRead9-LAB221024_HighestPeak_PickUp_2BG_v2.cpp -o BinaryRead9-LAB221024_HighestPeak_PickUp_2BG_v2 `root-config --cflags --glibs` -lSpectrum

#This program perform similar task as the first program but give information of 2 peaks (the 2 hightest peaks)
#The second program take one more argument as the real energy value of the second peak
./BinaryRead9-LAB221024_HighestPeak_PickUp_2BG_v2 Data_pointsource/101740.889_C1BG 100 800 0.4 7 306.82 201.83
./BinaryRead9-LAB221024_HighestPeak_PickUp_2BG_v2 Data_pointsource/173228.492_SET3_HV1500-1600_T38-36-37-37C_BG_60S 100 800 0.4 7 306.82 201.83

**********Second part- Analyzing peak from the energy spetrum of each channel

# The program below reading the .txt output from the ReadBinary program above and execute energy calibration for the ASPET
# [P0]*(exp(x/[P1]))+[P2]  : is the fitting cunction applied
# The output file is "EnergyCal/Energy_Calibration_v1.txt"
# gmslID / sticID / Channel / P0 / error0 / P1 / error1 / P2 / error2 / Chi2
# only the fitting with Chi2<1 is in the output otherwise all parameters are set as 0

# To compile the energy calibration program

>g++ Read_BinOut_Calibrate_v1.cpp -o Read_BinOut_Calibrate_v1 `root-config --cflags --glibs` -lSpectrum

# To execute the program, thers is a example below with the arguments are the _HighestPeak (filename) output from the BinaryRead

>./Read_BinOut_Calibrate_v1 Data_pointsource/173228.492_SET3_HV1500-1600_T38-36-37-37C_BG_60S_HighestPeak Data_pointsource/174133.387_SET3_HV1500-1600_T38-39-37-38C_NA22-D6CM-CENTER_L-52.5MM_60S_HighestPeak Data_pointsource/182452.671_SET3_HV1500-1600_T38-39-37-38C_CS137-D6CM-CENTER_L-52.5MM_60S_HighestPeak Data_pointsource/181821.792_SET3_HV1500-1600_T38-39-37-38C_CO60-D6CM-CENTER_L-52.5MM_300S_HighestPeak
>./Read_BinOut_Calibrate_v1 Data_pointsource/173228.492_SET3_HV1500-1600_T38-36-37-37C_BG_60S_HighestPeak Data_pointsource/174133.387_SET3_HV1500-1600_T38-39-37-38C_NA22-D6CM-CENTER_L-52.5MM_60S_HighestPeak Data_pointsource/182452.671_SET3_HV1500-1600_T38-39-37-38C_CS137-D6CM-CENTER_L-52.5MM_60S_HighestPeak Data_pointsource/181821.792_SET3_HV1500-1600_T38-39-37-38C_CO60-D6CM-CENTER_L-52.5MM_300S_HighestPeak Data_pointsource/101740.889_C1BG_SecondHighestPeak_v2
>./Read_BinOut_Calibrate_v1 Data_pointsource/173228.492_SET3_HV1500-1600_T38-36-37-37C_BG_60S_HighestPeak Data_pointsource/174133.387_SET3_HV1500-1600_T38-39-37-38C_NA22-D6CM-CENTER_L-52.5MM_60S_HighestPeak Data_pointsource/182452.671_SET3_HV1500-1600_T38-39-37-38C_CS137-D6CM-CENTER_L-52.5MM_60S_HighestPeak Data_pointsource/181821.792_SET3_HV1500-1600_T38-39-37-38C_CO60-D6CM-CENTER_L-52.5MM_300S_HighestPeak Data_pointsource/173228.492_SET3_HV1500-1600_T38-36-37-37C_BG_60S_SecondHighestPeak_v2
# To run this program, at least 3 input files are required (corresponding to 3 peaks for 3 parameter from the energy fitting function) (number of argument are not limited)
#  This program will give the 2Dmap of Chi2 for all channel to evaluate the fitting performance
# Fitting functions are also displayed for all channels in EnergyCal/fitting/ directory

##################################### Update correction for energy calibration ###############################################################
#To run the code for energy calibration update from the ASPET .bin (this program should read the measurement of BG 306.82 keV or PAG 511 keV as dominator)

#Go to the directory

#complie ReadBinary file

g++ BinaryRead9-LAB240506_EnergyCalCorrection.cpp -o BinaryRead9-LAB240506_EnergyCalCorrection `root-config --cflags --glibs` -lSpectrum


#This program will read the spectrum of each channel, pick up the highest and clearest energy peak and estimate the deviation from the calibration file in keV

# The program also displays the energy spectra of all channels and the corresponding peak pickup and fitting

#To execute the program; there are examples below with 6 arguments 

./BinaryRead9-LAB240506_EnergyCalCorrection Data/101740.889_C1BG Energy_Calibration_221226_v1.txt 100 800 0.6 7 306.82

argv[1]: The measurement of the investigated date (BG or PAG 511 keV spectra): .bin file

argv[2]: the original energy calibration file

argv[3-4]: boundary in ToT of the energy spectra

argv[5-6]: peak finding parameters for the spectrum

argv[7]: The known dominating gamma line in the measurement spectrum

#The program will create folders: CalUpdate ; CalUpdate/Spectrum 


