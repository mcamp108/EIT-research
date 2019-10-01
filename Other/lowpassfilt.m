function dataOut= lowpassfilt(N, dataIn, FR)
    lpfilt= designfilt('lowpassiir', 'FilterOrder', N, 'PassbandFrequency', 6, 'SampleRate', FR);
    dataOut= filter(lpfilt, dataIn);
end