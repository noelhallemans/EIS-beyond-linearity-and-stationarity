# EIS-beyond-linearity-and-stationarity

Here we share data and code discussed in the paper "Electrochemical impedance spectroscopy beyond linearity and stationarity---A critical review" (https://www.sciencedirect.com/science/article/pii/S0013468623011143). These are free to use! However, if you use this data/code for publications, I would ask you to cite the paper. 

The MATLAB script `writeORPmultisine.m` generates a sampled odd random phase multisine. Such a multisine can be written in a text file and applied through a potentiostat that can apply user-defined excitations (for instance the Gamry Interface 5000E).

EDIT: The phase angle of the impedance data above 10 Hz is affected by the potentiostat. The applied current has to be oversampled more to obtain correct impedance data with the Gamry 5000E potentiostat. This is not possible yet.

Please only use the impedance data for frequencies lower than 10 Hz.

In the folder `Data Samsung 48X`, odd random phase multisines are applied to a Li-ion battery using the Gamry Interface 5000E. The battery we use is the commercially available Samsung INR21700-48X cell, a 4.8 Ah 21 700 cell format with cathodes based on lithiated metal oxide (Co, Ni, Al) and anodes based on intercalation graphite and blended Si. Measurements are performed in a thermal chamber at 5°C and 25°C.

For classical EIS experiments (`classicalEIS_25degrees.mat` and `classicalEIS_5degrees.mat`), the multisine is applied in steady state at different SOC levels (10,20,...,80% at 5°C and 10,20,...,90% at 25°C). 10 periods of the multisine are measured. The MATLAB script `classicalEIS.m` plots the measured data at a particular temperature and SOC, and estimates the impedance. Using this data and code one can recreate Fig. 3, 11 and 14 of the paper.

For operando EIS experiments (`operandoEIS_25degrees.mat` and `operandoEIS_5degrees.mat`), the multisine is applied superimposed on a C/2 charging current (2.4 A). At 5°C and 25°C, 29 and 31 periods of the multisine are measured, respectively. The MATLAB script `operandoEIS.m` plots the measured data. We also include estimated time-varying impedance data using operando EIS (https://www.sciencedirect.com/science/article/pii/S037877532200982X). Using this data and code one can recreate Fig. 9, 15 and 16 of the paper.

Noël
