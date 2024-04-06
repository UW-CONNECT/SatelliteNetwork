 def checkXCORR3(self, possible_samples):               
        detected = False 
        peakdx = 0
        
        xcorr_arr = [];
        xm_arr= [] 
        window_size = int(len(self.REF_PREAMBLE)/2)
        for i in range(0,len(possible_samples) - 2*(window_size)) :     
            window_1 = possible_samples[ i : ( i + window_size ) ]
            window_2 = possible_samples[( i +window_size) : ( i +window_size+window_size)]
            
   
            # this scales the correlation between 0,1
            fft_window1 = (window_1 - np.mean(window_1)) / (np.std(window_1) * len(window_1))
            fft_window2 = (window_2 - np.mean(window_2)) / (np.std(window_2))
            
            xcorr_val = np.squeeze(abs(np.correlate(fft_window1, fft_window2)))
            xcorr_arr.append(xcorr_val) 
                    
        xcorr_arr = np.array(xcorr_arr, dtype="object")
        argm = np.argmax(xcorr_arr)
        if ( (argm != len(xcorr_arr)-1) ):
            # detected = True
            # print("Packet detected.")
            detected=True
                
        #imax_peak = np.argmax(xcorr_arr)
        max_xcorr = max(xcorr_arr)
        peakdx = argm           
                   
        #[peaks, properties] = signal.find_peaks(xcorr_arr)
        
        #if peaks is not None: 
     
        # if detected: 
            # plt.figure(2)
            # plt.plot(xcorr_arr)
            # plt.axvline(x = peakdx, color = 'r')  
            # print(detected)
            # #plt.figure(3)
            # #plt.plot(xm_arr)        
            # plt.show()       
        
        return detected, peakdx    
        
    def checkXCORR2(self, possible_samples):               
        detected = False 
        peakdx = 0
        
        xcorr_arr = [];
        xm_arr= [] 
        # t = np.linspace(0, len(self.REF_PREAMBLE)/self.FS, len(self.REF_PREAMBLE))
        window_2 = self.REF_PREAMBLE
        # for i in range(0,len(possible_samples) - 3*len(self.REF_PREAMBLE)) :     
        for i in range(0,len(possible_samples) - len(self.REF_PREAMBLE)) :
            window_1 = possible_samples[ i : ( i + len(self.REF_PREAMBLE)) ]         
            # window_2 = 
            
            
            # this scales the correlation between 0,1
            fft_window1 = (window_1 - np.mean(window_1)) / (np.std(window_1) * len(window_1))
            fft_window2 = (window_2 - np.mean(window_2)) / (np.std(window_2))
            
            xcorr_val = np.squeeze(abs(np.correlate(fft_window1, fft_window2)))     
            
            
            xcorr_arr.append(xcorr_val)
            #xcorr_arr.append(xcorr_val) 
                    
        xcorr_arr = np.array(xcorr_arr, dtype="object")
        
        if len(xcorr_arr) == 0:
            return detected, peakdx
        argm = np.argmax(xcorr_arr)
        # # reject samples without any correlation 
        #print(max(xcorr_arr))
        # detected = True
        print(max(xcorr_arr))
        # if (max(xcorr_arr) > .1 and (argm != len(xcorr_arr)-1) and (argm != 0) ):
            # detected = True
            # print("Packet detected.")
        detected = True
                
        #imax_peak = np.argmax(xcorr_arr)
        max_xcorr = max(xcorr_arr)
        peakdx = argm           
                   
        #[peaks, properties] = signal.find_peaks(xcorr_arr)
        
        #if peaks is not None: 
     
        if detected: 
            plt.figure(2)
            plt.plot(xcorr_arr)
            plt.axvline(x = peakdx, color = 'r')  
            #print(detected)
            #plt.figure(3)
            #plt.plot(xm_arr)        
            plt.show()       
        
        return detected, peakdx       
    
    
    def checkXCORR(self, possible_samples):               
        detected = False 
        peakdx = 0
        
        xcorr_arr = [];
        xm_arr= [] 
        for i in range(0,len(possible_samples) - 3*self.PREAMBLE_SIZE) :     
            window_1 = possible_samples[ i : ( i + self.PREAMBLE_SIZE ) ]
            window_3 = possible_samples[( i +self.PREAMBLE_SIZE) : ( i+2*self.PREAMBLE_SIZE)]
            window_2 = possible_samples[( i +2*self.PREAMBLE_SIZE) : ( i+3*self.PREAMBLE_SIZE)]
            # this scales the correlation between 0,1
            fft_window1 = (window_1 - np.mean(window_1)) / (np.std(window_1) * len(window_1))
            fft_window2 = (window_2 - np.mean(window_2)) / (np.std(window_2))
            
            xcorr_val = np.squeeze(abs(np.correlate(fft_window1, fft_window2)))
            
            # fft_window3 = (window_3 - np.mean(window_3)) / (np.std(window_3))
            # xcorr_val2 = np.squeeze(abs(np.correlate(fft_window1, fft_window3)))            
            
            
            # fft_window2 = (window_2 - np.mean(window_2)) / (np.std(window_2) * len(window_2))
            # fft_window3 = (window_3 - np.mean(window_3)) / (np.std(window_3))
            
            # xcorr_val3 = np.squeeze(abs(np.correlate(fft_window2, fft_window3)))
            # print(xcorr_val3)
            # xcorr_arr.append(xcorr_val * xcorr_val2 * xcorr_val3 * 1000)
            # xcorr_arr.append((xcorr_val + xcorr_val2 + xcorr_val3)**2 )
            xcorr_arr.append(xcorr_val) 
                    
        xcorr_arr = np.array(xcorr_arr, dtype="object")
        
        # xcorr_arr = xcorr_arr / max(xcorr_arr)
        
        if len(xcorr_arr) == 0:
            return detected, peakdx
        argm = np.argmax(xcorr_arr)
        # # reject samples without any correlation 
        #print(max(xcorr_arr))
        # detected = True
        # if (max(xcorr_arr) > .04 and (argm != len(xcorr_arr)-1) and (argm != 0) ):
        if (max(xcorr_arr) > .15 and (argm != len(xcorr_arr)-1) and (argm != 0) ):
            detected = True
            print("Packet detected.")
        # detected=True
                
        #imax_peak = np.argmax(xcorr_arr)
        max_xcorr = max(xcorr_arr)
        peakdx = argm           
                   
        #[peaks, properties] = signal.find_peaks(xcorr_arr)
        
        #if peaks is not None: 
     
        # if detected: 
            # plt.figure(2)
            # plt.plot(xcorr_arr)
            # plt.axvline(x = peakdx, color = 'r')  
            # print(detected)
            # #plt.figure(3)
            # #plt.plot(xm_arr)        
            # plt.show()       
        
        return detected, peakdx    