from __future__ import division, print_function
import numpy as np
import ezca
import signal, sys

z = ezca.Ezca(prefix='C1:')

def main():
    '''
    Write a sine wave to an EPICS channel 
    until a SIGINT is received.
    Only a single channel is supported for now.
    Example usage:
        python epicsSineWave.py C1:SUS-ETMY_LRBiasSet 10 20
    will write a sine wave at 20 Hz of amplitude 10 to that channel. 
    '''
    def interrupt_handler(*args):
        print('Terminating excitation...')
        z.write(chan, initValue)
        sys.exit(0)
    #Get the arguments from the command line
    if len(sys.argv) < 4:
        print(main.__doc__)
        return 
    chan = sys.argv[1].strip('C1:')
    amp = float(sys.argv[2])
    f = float(sys.argv[3])
    initValue = z.read(chan)
    print('Writing a sine wave of amplitude {} to {}. Use CTRL+C to stop the excitation.'.format(amp, chan))
    t = np.linspace(0, 1/f, 1000)
    y = amp * np.sin(2*np.pi*f*t)
    signal.signal(signal.SIGINT, interrupt_handler)
    while True:
        for yy in y:
            z.write(chan, yy)
    return



if __name__ == '__main__':
    main()
