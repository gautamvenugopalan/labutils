def getOLRIN(opts, tAvg):
    chans = []
    for jj in opts:
        chans.append('C1:SUS-'+jj+'_OPLEV_SUM')
    avgs=cds.avg(tAvg, chans)
    for ii, jj in zip(opts,avgs):
        print('Factor for {} is {}'.format(ii,1./jj))
    return(avgs)

