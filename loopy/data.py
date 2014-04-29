
import numpy

def smooth(x,winlen=5,window='flat'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        winlen: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming',
            'bartlett', 'blackman' flat window will produce a moving average
            smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead
    of a string

    NOTE: length(output) != length(input), to correct this: return
    y[(winlen/2-1):-(winlen/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < winlen:
        raise ValueError, "Input vector needs to be bigger than window size."

    if winlen % 2 == 0:
        winlen += 1

    if winlen<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    winh = (winlen-1)/2
    s=numpy.r_[x[winh:0:-1],x,x[-2:-(winh+2):-1]]

    if window == 'flat': #moving average
        w=numpy.ones(winlen,'d')
    else:
        w=eval('numpy.'+window+'(winlen)')

    w /= w.sum()
    y=numpy.convolve(w,s,mode='valid')

    return y
