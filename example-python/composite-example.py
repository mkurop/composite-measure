from composite.composite_measure import composite
import numpy as np
from scipy.io import wavfile
from typing import Tuple
import os

def read_wav_and_normalize(file_ : str) -> Tuple[np.ndarray,int]:

    if not os.path.exists(file_):
        raise ValueError(f"The file {file_} does not exist.")

    sr, s = wavfile.read(file_)

    s = s/2.**15
    

    np.random.seed(1000) # causes the dither is same on each run

    s += np.random.randn(*s.shape)*1.e-6 # add dither to improve numerical behaviour

    return np.float32(s), int(sr)

if __name__ == "__main__":

    clean_speech, clean_speech_sampling_frequency = read_wav_and_normalize("../data/clear8.wav")

    enhanced_speech, enhanced_speech_sampling_frequency = read_wav_and_normalize("../data/enhanced8.wav")

    sig, bak, ovl = composite(clean_speech, float(clean_speech_sampling_frequency), enhanced_speech, float(enhanced_speech_sampling_frequency))

    print(f"{sig}, {bak}, {ovl}")
