import numpy                as np
import matplotlib.pyplot    as plt
import matplotlib.patches   as patches
import math
import numpy.matlib         as matlib
import seaborn              as sns
import pandas               as pd
import warnings
import os
import pdb
from IPython.core.debugger import Tracer

class nrn_network(object):

    def __init__(self, neurons):
        print(np.shape(neurons))
        self.neurons = [neurons[i] for i in range(0, len(neurons))]

    def meanFR(self,norm=1,normrange=(36*3600, 60*3600),binsz=3600,trange=(36, 180*3600),plotmean=1,plotindiv=0,saveflag=None,savedir=None):

        edges   = np.arange(trange[0],(trange[1]+(binsz)),binsz)
        aaa     = np.empty( (np.size(self.neurons),np.size(edges)-1) )
        aaa[:]  = np.nan

        # figure out the dep condition passed to this script:
        deplist = [e.deprived for e in self.neurons]
        if np.unique(deplist) == 0:
            netcond = 'control'
        elif np.unique(deplist) == 1:
            netcond = 'deprived'
        else:
            netcond = 'mixed'

        print('Detected condition is {}'.format(netcond))

        ncount = -1
        for i in self.neurons:
            # Check to see if the cell is online at baseline and on for more than 70% of all recording
            basecheck           = np.logical_and( np.min(i.onTime[0])<(36*3600), np.max(i.offTime[0])>(60*3600) )
            try:
                seventycheck    = sum([b for b in i.offTime-i.onTime ])>(0.7*(trange[1] - trange[0]))
            except:
                Tracer()()

            if np.any([i.onTime>trange[1], i.offTime>trange[1]]):
                i.onTime    = np.delete(i.onTime, np.where(i.onTime>trange[1]))
                i.offTime   = np.delete(i.offTime, np.where(i.offTime>trange[1]))
                i.offTime   = np.append(i.offTime, trange[1])


            if basecheck and seventycheck:
                ncount += 1
                tmp         = np.histogram(i.time,edges)
                spkcounts   = tmp[0].astype(float) # this is the raw histogram count of spikes per bin. you need to operate on this and clear out/edit bins that intersect with on and off time.

                if np.size(i.onTime) == 1 and i.onTime[0]<trange[1]:
                    g = next(t for t in edges if t>i.onTime[0]) #Find the first bin edge that is larger than the onTime
                    g = np.where(edges==g)[0][0]
                    spkcounts[0:g] = np.nan

                    try:
                        h = next(t for t in edges if t>i.offTime[0]) #Find the first bin edge that is larger than the offTime
                        h = np.where(edges==h)[0][0]-1 # Set this to the last bin that is smaller than the offtime
                        spkcounts[h:] = np.nan
                    except:
                        print('Neuron online past time range! Nice.')

                    aaa[ncount,:] = spkcounts/binsz # write this in Hz

                elif np.size(i.onTime)>1 and i.onTime[0]<trange[1]: # case in which there are multiple on and off times

                    # first delete spikes leading the OG on time
                    g = next(t for t in edges if t>i.onTime[0]) #Find the first bin edge that is larger than the onTime
                    g = np.where(edges==g)[0][0]
                    spkcounts[0:g] = np.nan

                    # then delete spikes trailing the final offtime
                    try:
                        h = next(t for t in edges if t>i.offTime[-1]) # Find the first bin edge that is larger than the offTime
                        h = np.where(edges==h)[0][0]-1 # Set this to the last bin that is smaller than the offtime
                        spkcounts[h:] = np.nan
                    except:
                        print('Neuron online past time range! Nice.')

                    aaa[ncount,:] = spkcounts/binsz # write this in Hz
                    # now go deal with off periods in bins in the middle of the recording
                    otcount = 0
                    for k in i.offTime[0:-1]: # cycle through all on times until the final one
                        otcount += 1
                        try:
                            gt = next(t for t in edges if t>k) # Find the last bin edge that is smaller than the mid recording off time.
                        except:
                            print('Caught line 74')
                            pdb.set_trace()

                        gt = np.where(edges==gt)[0][0]-1
                        et = next(t for t in edges if t>i.onTime[otcount])
                        et = np.where(edges == et)[0][0] # find the first bin that is larger than the resumed recording time

                        aaa[ncount,gt:et] = np.nan # write affected bins as nan and then see if you can recover any of them

                        timepre     = k - edges[gt] # time btwn the start of bin and the off time.
                        timepost    = edges[et] - i.onTime[otcount] # time btwn on time and end of bin.

                        # deal with instances in which the on and off times occur in different bins:
                        if gt != et-1:

                            if timepre>0.5*binsz:
                                presubbin           = np.sum(np.logical_and(i.time>edges[gt], i.time<k))
                                aaa[ncount,gt]      = presubbin/binsz # write the estimated FR value in the bin containing the offT
                            if timepost>0.5*binsz:
                                postsubbin          = np.sum(np.logical_and(i.time>i.onTime[otcount], i.time<edges[et]))
                                aaa[ncount,et-1]    = postsubbin/binsz # write the estimated FR value in the bin containing the onT

                        elif gt == et-1:
                            # here deal with the circumstance in which you have aaa chunk in the middle of aaa single block missing and then calculate the spikes and on times remaining for that one
                            tottime = np.sum((timepre,timepost)) # total amount of time that cell is ON in the affected bin

                            if tottime>0.5*binsz:
                                presub1bin      = np.sum(np.logical_and(i.time<k, i.time>edges[gt]))
                                postsub1bin     = np.sum(np.logical_and(i.time>i.onTime[otcount], i.time<edges[et]))
                                estimatedrate   = (presub1bin + postsub1bin)/tottime
                                aaa[ncount,gt]  = estimatedrate


        aaa         = np.delete(aaa,np.where(np.nansum(aaa,1) == 0),0) # delete rows that are empty (nan)
        b0          = np.int((normrange[0]/binsz)-1) # in bin number, where does baseline start?
        b1          = np.int((normrange[1]/binsz)-1)

        # Plot each of the traces individually here.
        #T = range(aaa.shape[1])

        #for i in range(aaa.shape[0]):
        #    plt.plot(T, aaa[i,:])

        #plt.show()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            basemean    = np.nanmean(aaa[:,b0:b1],1)
            basemean    = np.rot90(np.matlib.repmat(basemean,np.shape(aaa)[1],1),3)
            a2          = np.divide(aaa,basemean)   # baseline normalized
            a3          = np.nanmean(a2,0)          # mean normalized trace
            SEMs        = np.divide(np.nanstd(a2,0),np.sqrt(np.sum(a2[:,:]>0,0)))

            # set up for aaa pandas dataframe to play extra friendly with the seaborn plotting package
            a2panda     = np.reshape(a2,np.size(a2))
            bincenters  = (edges[0:-1]+1800)/3600
            tmpanda     = np.ndarray.flatten(np.matlib.repmat(bincenters,np.shape(a2)[0],1))
            cellnum     = np.repeat(range(0,np.shape(a2)[0]), np.shape(a2)[1] )
            condition   = np.repeat(netcond,np.shape(a2panda)[0])


            d = {'Time (hours)' : pd.Series(tmpanda), 'Firing Rate (norm)' : pd.Series(a2panda), 'Condition' : pd.Series(condition), 'neuron' : pd.Series(cellnum)}
            df = pd.DataFrame(d)

            if plotmean:
                plt.ion()
                #plt.rcParams['font.family'] = 'serif'
                fig     = plt.figure()
                ax      = fig.add_subplot(111, frame_on=False)
                #plt.plot(edges[0:-1]/3600,a3)
                sns.tsplot(data=df, time="Time (hours)", value="Firing Rate (norm)", condition="Condition", unit="neuron",estimator=np.nanmean)

                lt_off      = np.arange(12*3600,trange[1]+12*3600,24*3600)
                ylims       = ax.get_ylim()[1]
                xlims       = ax.get_xlim()
                # cycle through the lights-off times and plot aaa transparent grey bar to indicate the dark hours
                for p in lt_off/3600:
                    print(p)
                    ax.add_patch(patches.Rectangle((p, 0), 12, ylims, facecolor="grey",alpha=0.2, edgecolor="none"))

                ax.plot(xlims,(1, 1),'--r', alpha = 0.3 )

                ncellsinfig = np.shape(aaa)[0]
                ylm         = ax.get_ylim()[1]
                xlm         = ax.get_xlim()
                xtxt        = xlm[0] + 0.2*(xlm[1]-xlm[0] )
                ytxt        = ylm*0.7
                txt         = 'ncells: ' + str( ncellsinfig )
                ax.text(xtxt,ytxt,txt)


                # plot the arrow to indicate lid suture:
                ax.arrow(61,1.6, 0, -0.4, head_width=5, head_length=0.14, linewidth=2.5, fc='r', ec='k')
                ax.set_ylim([0,2])

                if savedir:
                    ftitle = input('Entirefigure title/FN.')
                    plt.title(ftitle)
                    fig.show()
                    plt.savefig(savedir + os.sep + ftitle + '.pdf')
                else:
                    ftitle = input('Entirefigure title/FN.')
                    plt.title(ftitle)
                    fig.show()


            if plotindiv:
                plt.ion()
                # run this to plot all cells individually
                for ee in np.arange(0,np.shape(aaa)[0]):
                    fig = plt.figure(ee)
                    fig.add_subplot(111,frame_on=False)
                    fig.canvas.set_window_title('Neuron{}'.format(ee))
                    plt.plot(aaa[ee,:])
                    fig.show() # use plt.close('all') to quickly clear figures

    def epocheffect(self, rebound = (125, 180), savedir = None):
        import re
        import os
        import csv
        import pandas as pd

        temp = [e.deprived for e in self.neurons]
        if np.unique(temp)[0] == 1:
            dep = 'Depr'
        elif np.unique(temp)[0] == 0:
            dep = 'Control'

        celltype = input('What type of cell are you working on?  ')

        # Find the state codes
        temp = self.neurons[0].HDF5_tag[0]
        slashes = [m.start() for m in re.finditer('/', temp)][-1]
        temp = temp[0:slashes+1] + 'STATETIMES'

        if not os.path.exists(temp):
            y = input('Where are the statetimes for these data? Enter directory: ')
            os.chdir(y)

        d = {} # this will be your statetimes list in local memory
        for root,dirs,files in os.walk(temp):
            for file in files:
                if file.endswith(".csv"):
                    nm      = file[0:file.find('_')]
                    f       = open(root + '/' + file, 'r')
                    reader  = csv.reader(f)
                    row     = list(reader)

                    stm     = np.zeros(np.shape(row))
                    count   = 0
                    for e in row:
                        stm[count, 0] = np.int(e[0])
                        stm[count, 1] = np.float(e[1])
                        count+=1

                    d[nm]  = stm

        # this will hold the average output from each cell
        nrnlevelones    = np.zeros([np.size(self.neurons), 100])
        nrnleveltwos    = np.zeros([np.size(self.neurons), 100])
        nrnlevelfours   = np.zeros([np.size(self.neurons), 100])
        nrnlevelfives   = np.zeros([np.size(self.neurons), 100])

        nrncount = 0;
        for ee in self.neurons:

            if ee.animal in d:
                print('Processing {} for sleep/wake state FR changes.'.format(ee.animal))
                thesetimes = d[ee.animal]

                kills1 = np.argmin( thesetimes[:,1] < (rebound[0] * 3600) )
                kills2 = np.argmax( thesetimes[:,1] > (rebound[1] * 3600) )

                thesetimes  = thesetimes[kills1+1 : kills2,:] # this is the list of times relevant to your analysis window

                delta100s   = np.zeros([np.shape(thesetimes)[0],100])

                for ss in np.arange(np.shape(thesetimes)[0]-1):
                    t0  = thesetimes[ss,1]
                    t1  = thesetimes[ss+1,1]

                    if t1-t0>60: # time minimum of 60 sec
                        spk = ee.time[(ee.time>t0) * (ee.time<=t1)]

                        ltemp           = np.linspace( t0, t1,101 )
                        bar100          = np.histogram(spk, ltemp)
                        first20         = np.nanmean(bar100[0][0:19])
                        if np.sum(bar100[0][0:10])>5: # if NOT 0
                            bar100norm      = bar100[0]/first20
                            delta100s[ss,:] = bar100norm
                        else:
                            delta100s[ss,:] = np.repeat(999,100)
                    else:
                        delta100s[ss,:] = np.repeat(999,100)

                # get rid of zero rows:
                skills      = np.where( np.sum(delta100s,axis=1) == 0)[0]
                delta100s   = np.delete(delta100s,skills,axis = 0)
                thesetimes  = np.delete(thesetimes,skills,axis = 0)

                # and 999 flags
                skills9      = np.where( np.sum(delta100s,axis=1) == 999*100)[0]
                delta100s   = np.delete(delta100s,skills9,axis = 0)
                thesetimes  = np.delete(thesetimes,skills9,axis = 0)

                nankills    = np.unique(np.where(np.isnan(delta100s))[0])
                infkills    = np.unique(np.where(np.isinf(delta100s))[0])
                kills       = np.append(nankills,infkills)
                kills       = np.unique(kills)

                # index the state entries by specific state (1,2,4,5)
                ones    = np.where(thesetimes[:,0] == 1)
                twos    = np.where(thesetimes[:,0] == 2)
                fours   = np.where(thesetimes[:,0] == 4)
                fives   = np.where(thesetimes[:,0] == 5)

                # find overlap of state subets and the "kills" variable
                onekovrlp = np.where(np.in1d(ones,kills))
                twokovrlp = np.where(np.in1d(twos,kills))
                fourkovrlp = np.where(np.in1d(fours,kills))
                fivekovrlp = np.where(np.in1d(fives,kills))

                # remove the junk entries found above
                ones = np.delete(ones,onekovrlp)
                twos = np.delete(twos,twokovrlp)
                fours = np.delete(fours,fourkovrlp)
                fives = np.delete(fives,fivekovrlp)

                # make the means
                onemeans = np.mean(delta100s[ones,:],0)
                twomeans = np.mean(delta100s[twos,:],0)
                fourmeans = np.mean(delta100s[fours,:],0)
                fivemeans = np.mean(delta100s[fives,:],0)

                # assign the neuron means to the network level variable
                nrnlevelones[nrncount,:]    = onemeans
                nrnleveltwos[nrncount,:]    = twomeans
                nrnlevelfours[nrncount,:]   = fourmeans
                nrnlevelfives[nrncount,:]   = fivemeans

                nrncount += 1


        # get rid of zero rows:
        zkills1 = np.where( np.sum(nrnlevelones,axis=1) == 0)[0]
        nrnlevelones = np.delete(nrnlevelones,zkills1,axis = 0)

        zkills2 = np.where( np.sum(nrnleveltwos,axis=1) == 0)[0]
        nrnleveltwos = np.delete(nrnleveltwos,zkills2,axis = 0)

        zkills4 = np.where( np.sum(nrnlevelfours,axis=1) == 0)[0]
        nrnlevelfours = np.delete(nrnlevelfours,zkills4,axis = 0)

        zkills5 = np.where( np.sum(nrnlevelfives,axis=1) == 0)[0]
        nrnlevelfives = np.delete(nrnlevelfives,zkills5,axis = 0)



        sns.set(style="darkgrid", palette="Set2")

        plt.ion()
        fp = plt.subplot()
        sns.tsplot(nrnlevelones,color = '#069af3',ax = fp, condition='REM',legend = True)
        sns.tsplot(nrnleveltwos,color = '#d6b4fc',ax = fp, condition='NREM',legend = True)
        sns.tsplot(nrnlevelfours,color = '#fe2c54',ax = fp, condition='Active',legend = True)
        sns.tsplot(nrnlevelfives,color = '#fcb001',ax = fp, condition='Quiet',legend = True)

        fp.set(xlabel='% epoch time', ylabel='normalized firing rate')
        fp.set_title('{} {} Rebound Epoch Effect'.format(dep,celltype))

        plt.show()

        if savedir:
            g = plt.gcf()
            g.savefig( savedir + os.sep  + '{}_{}_KDE_SW_EpochsEffect.pdf'.format(celltype,dep) )



    def swbasics(self,binmin = 60,savedir = None):
        '''binmin is the minimum epoch duration (in seconds) required to be considered in these analyses.'''
        import re
        import os
        import csv
        import pandas as pd

        if savedir:
            ftitle = input('What Cell Type is this?')


        temp = [e.deprived for e in self.neurons]
        if np.unique(temp)[0] == 1:
            dep = 'Depr'
        elif np.unique(temp)[0] == 0:
            dep = 'Control'


        # Find the state codes
        temp = self.neurons[0].HDF5_tag[0]
        slashes = [m.start() for m in re.finditer('/', temp)][-1]
        temp = temp[0:slashes+1] + 'STATETIMES'

        if not os.path.exists(temp):
            y = input('Where are the statetimes for these data? Enter directory: ')
            os.chdir(y)

        d = {} # this will be your statetimes list in local memory
        for root,dirs,files in os.walk(temp):
            for file in files:
                if file.endswith(".csv"):
                    nm      = file[0:file.find('_')]
                    f       = open(root + '/' + file, 'r')
                    reader  = csv.reader(f)
                    row     = list(reader)

                    stm     = np.zeros(np.shape(row))
                    count   = 0
                    for e in row:
                        stm[count, 0] = np.int(e[0])
                        stm[count, 1] = np.float(e[1])
                        count+=1

                    d[nm]  = stm

        fr4st   = np.zeros([np.size(self.neurons),4])
        cv4st   = np.zeros([np.size(self.neurons),4])
        frac    = np.zeros([np.size(self.neurons)])


        nrncount = 0
        for ee in self.neurons:

            if ee.animal in d:
                print('Processing neuron {} from {} for sleep/wake descriptives.'.format(nrncount,ee.animal))
                thesetimes = d[ee.animal]

                # Select only the relevant times for analyses:
                kills1 = []
                kills2 = []
                if ee.deprived == 1:
                    print('Deprived')
                    kills1  = np.argmin( thesetimes[:,1] < (36 * 3600) )
                    kills2  = np.argmax( thesetimes[:,1] > (60 * 3600) )
                else:
                    print('Control')
                    kills1  = np.argmin( thesetimes[:,1] < (36 * 3600) )
                    kills2  = next((i for i, x in enumerate(thesetimes[:,1] > ee.offTime[-1] ) if x), None)
                if kills2 is None:
                    kills2  = np.shape(thesetimes)[0]

                thesetimes  = thesetimes[kills1+1 : kills2,:] # this is the list of times relevant to your analysis window
                frates      = np.zeros([np.shape(thesetimes)[0],2])
                tempCV      = np.zeros([np.shape(thesetimes)[0],2])
                FRauto      = np.zeros([np.shape(thesetimes)[0],4])

                ecount = 0 # cycle through epochs
                for ss in np.arange(np.shape(thesetimes)[0]-2):
                    spk     = []
                    t0      = thesetimes[ss,1]
                    t1      = thesetimes[ss+1,1]
                    t2      = thesetimes[ss+2,1]
                    edur    = t1-t0
                    spk     = ee.time[(ee.time>t0) * (ee.time<=t1)]

                    # THIS IS BASIC STATS:
                    if edur>=binmin and np.size(spk)>5:
                        # Calculate the firing rate in this epoch
                        frates[ecount,0] = np.size(spk)/edur
                        frates[ecount,1] = thesetimes[ss,0]

                        # calculate the CV of the spiking in each of the epochs
                        tmpisi           = np.diff(spk)
                        tempCV[ecount,0] = np.std(tmpisi)/np.mean(tmpisi)
                        tempCV[ecount,1] = thesetimes[ss,0]

                        # Check to see if the next epoch is also greater than binmin and has more than five spikes. if it does, then you can compare the FR of the cell in epoch n and epoch n + 1.
                        if (t2-t1)>=binmin and np.size(ee.time[(ee.time>t0) * (ee.time<=t1)])>5:
                            FRauto[ecount,0] = np.size(spk)/edur # Calculate FR for epoch n
                            FRauto[ecount,1] = np.size(ee.time[(ee.time>t0) * (ee.time<=t1)]) / (t2-t1) # FR for epoch n + 1
                            FRauto[ecount,2] = thesetimes[ss,0] # Behavioral state code for epoch n
                            FRauto[ecount,3] = thesetimes[ss+1,0] # Behavioral state code for epoch n + 1

                    ecount +=1


                # Find the empty rows in FR and CV matrices. Clear those rows.
                killzeros   = np.where(frates[:,1] == 0)
                frates      = np.delete(frates, killzeros, axis = 0)
                tempCV      = np.delete(tempCV, killzeros, axis = 0)

                # Find the empty rows in the autocorrelation matrix. Clear those rows.
                killautoz   = np.where(FRauto[:,1] == 0)
                FRauto      = np.delete(FRauto, killautoz, axis = 0)

                # set up firing rate by epoch matrix
                fr4st[nrncount,0]   = np.mean( frates[ np.where(frates[:,1] == 1), 0 ] )
                fr4st[nrncount,1]   = np.mean( frates[ np.where(frates[:,1] == 2), 0 ] )
                fr4st[nrncount,2]   = np.mean( frates[ np.where(frates[:,1] == 4), 0 ] )
                fr4st[nrncount,3]   = np.mean( frates[ np.where(frates[:,1] == 5), 0 ] )

                #set up CV by epoch matrix
                cv4st[nrncount,0]   = np.nanmean( tempCV[ np.where(tempCV[:,1] == 1), 0 ] )
                cv4st[nrncount,1]   = np.nanmean( tempCV[ np.where(tempCV[:,1] == 2), 0 ] )
                cv4st[nrncount,2]   = np.nanmean( tempCV[ np.where(tempCV[:,1] == 4), 0 ] )
                cv4st[nrncount,3]   = np.nanmean( tempCV[ np.where(tempCV[:,1] == 5), 0 ] )

                #set up the FR/state autocorrelation matrix
                frac[nrncount]      = np.corrcoef(FRauto[:,0],FRauto[:,1])[0,1]

                nrncount +=1

        killnrn = np.where(np.sum(fr4st,axis = 1) == 0)
        fr4st   = np.delete(fr4st, killnrn, axis = 0)
        cv4st   = np.delete(cv4st,killnrn, axis = 0)
        frac    = np.delete(frac, killnrn)

        perm    =  np.argsort([2,3,1,0])
        AQNR    = fr4st[:,perm]
        AQNRcv  = cv4st[:,perm]

        # calculate the basic stats?
        normAQNR = AQNR/ AQNR[:,0][:, np.newaxis]

        normFR   = pd.DataFrame(data = normAQNR, columns = ['Active','Quiet','NREM','REM'], index=np.arange(0,np.shape(AQNR)[0]),)
        CVdf     = pd.DataFrame(data = AQNRcv, columns = ['Active','Quiet','NREM','REM'], index=np.arange(0,np.shape(AQNRcv)[0]),)
        AQNRclrs = ['#fe2c54','#fcb001','#d6b4fc','#069af3']

        # plot the FR and ISI CV together on a single figure
        fig1, (ax1, ax2) = plt.subplots(ncols=2)
        sns.barplot(data = normFR, palette = AQNRclrs, ax = ax1)
        ax1.set_title('Normalized FR by State')
        ax1.set_ylabel('Firing Rate (norm)')
        sns.barplot(data = CVdf, palette = AQNRclrs, ax = ax2)
        ax2.set_title('ISI CV by state')
        ax2.set_ylabel('CV')

        if savedir:
            fig1.savefig(savedir + os.sep + '{}_{}_SW_Basics'.format(dep, ftitle) + '.pdf')


        # make a dot plot thingie.
        mxfr    = np.ceil(np.max(AQNR)*1.1)
        df      = pd.DataFrame(data = AQNR, columns = ['Active','Quiet','NREM','REM'], index=np.arange(0,np.shape(AQNR)[0]) )

        df['cellnum'] = [chr(i) for i in (np.arange(0,np.shape(AQNR)[0]))]

        # Make the PairGrid
        plt.ion()
        g = sns.PairGrid(df.sort_values("Active", ascending=False), x_vars=df.columns[:-1], y_vars=["cellnum"],  size=10, aspect=.25)
        g.map(sns.stripplot, size=10, orient="h", palette="Blues_d", edgecolor="gray")

        g.set(xlim=(0, mxfr), xlabel="Hz", ylabel="")

        # Use semantically meaningful titles for the columns
        titles = ["Active Wake", "Quiet Wake", "NREM", "REM"]

        for ax, title in zip(g.axes.flat, titles):

            # Set a different title for each axes
            ax.set(title=title)

            # Make the grid horizontal instead of vertical
            ax.xaxis.grid(False)
            ax.yaxis.grid(True)

        sns.despine(left=True, bottom=True)

        if savedir:
            plt.savefig(savedir + os.sep + '{}_{}_SW_SortedDots'.format(dep, ftitle) + '.pdf')

    def _pullautocorrstates(self,d,binmin):
        '''This is an internal function used to examine a single nrn network and return the correlation coefficient of the relationship between firing rate in behavior epoch "n" and epoch "n+1". Cells with very stable firing rates will have a high correlation coefficient. This is a way to ask whether or not there are systematic shifts in FR between states that may not retain directionality (i.e. the FR may decrease or increase, but it changes nonetheless). '''
        frac    = np.zeros([np.size(self.neurons)])

        nrncount = 0
        for ee in self.neurons:

            if ee.animal in d:
                print('Processing neuron {} from {} for sleep/wake descriptives.'.format(nrncount,ee.animal))
                thesetimes = d[ee.animal]

                # Select only the relevant times for analyses:
                kills1 = []
                kills2 = []
                if ee.deprived == 1:
                    print('Deprived')
                    kills1  = np.argmin( thesetimes[:,1] < (36 * 3600) )
                    kills2  = np.argmax( thesetimes[:,1] > (60 * 3600) )
                else:
                    print('Control')
                    kills1  = np.argmin( thesetimes[:,1] < (36 * 3600) )
                    kills2  = next((i for i, x in enumerate(thesetimes[:,1] > ee.offTime[-1] ) if x), None)
                if kills2 is None:
                    kills2  = np.shape(thesetimes)[0]

                thesetimes  = thesetimes[kills1+1 : kills2,:] # this is the list of times relevant to your analysis window
                FRauto      = np.zeros([np.shape(thesetimes)[0],4])

                ecount = 0 # cycle through epochs
                for ss in np.arange(np.shape(thesetimes)[0]-2):
                    spk     = []
                    t0      = thesetimes[ss,1]
                    t1      = thesetimes[ss+1,1]
                    t2      = thesetimes[ss+2,1]
                    edur    = t1-t0
                    edur2   = t2-t1
                    spk     = ee.time[(ee.time>t0) * (ee.time<=t1)]
                    spk2    = ee.time[(ee.time>t0) * (ee.time<=t1)]

                    # THIS IS BASIC STATS:
                    if edur>=binmin and np.size(spk)>5 and edur2>=binmin and np.size(spk2)>5:

                        FRauto[ecount,0] = np.size(spk)/edur # Calculate FR for epoch n
                        FRauto[ecount,1] = np.size(spk2)/edur2 # FR for epoch n + 1
                        FRauto[ecount,2] = thesetimes[ss,0] # Behavioral state code for epoch n
                        FRauto[ecount,3] = thesetimes[ss+1,0] # Behavioral state code for epoch n + 1

                    ecount +=1

                # Find the empty rows in the autocorrelation matrix. Clear those rows.
                killautoz   = np.where(FRauto[:,1] == 0)
                FRauto      = np.delete(FRauto, killautoz, axis = 0)

                #set up the FR/state autocorrelation matrix
                frac[nrncount]      = np.corrcoef(FRauto[:,0],FRauto[:,1])[0,1]

                nrncount +=1

        return frac


    def two_net_stateautocor(net1,net2,binmin=60,savedir=None):
        '''SCript looks at how well the FR in arousal state epoch "n" is correlated with epoch "n+1". The idea here is to check whether this changes arcross LD, SW, cell type, and MD. You'll need to keep building this out, but the basic framework is in place.

            INPUTS:
                binmin: the minimum epoch duration (in seconds) required to be considered in these analyses.
                net1:       The first network for analysis
                net2:       The second network for analysis
                savedir:    The directory for saving figures and any other output  '''


        import re
        import os
        import csv
        import pandas as pd
        # Find the state codes
        temp    = net1.neurons[0].HDF5_tag[0]
        slashes = [m.start() for m in re.finditer('/', temp)][-1]
        temp    = temp[0:slashes+1] + 'STATETIMES'

        if not os.path.exists(temp):
            y = input('Where are the statetimes for these data? Enter directory: ')
            os.chdir(y)

        d = {} # this will be your statetimes list in local memory
        for root,dirs,files in os.walk(temp):
            for file in files:
                if file.endswith(".csv"):
                    nm      = file[0:file.find('_')]
                    f       = open(root + '/' + file, 'r')
                    reader  = csv.reader(f)
                    row     = list(reader)

                    stm     = np.zeros(np.shape(row))
                    count   = 0
                    for e in row:
                        stm[count, 0] = np.int(e[0])
                        stm[count, 1] = np.float(e[1])
                        count+=1

                    d[nm]  = stm

        list_of_nets = [net1,net2]

        TEST    = [[i for i in list(mbt.nrn_network._pullautocorrstates(net,d,binmin))] for net in list_of_nets]

        for i in np.arange(0,np.size(TEST)):
            j = np.asarray(TEST[i])
            j = np.delete(j,np.where(j == 0) )
            TEST[i] = j

        automean = np.zeros(np.size(TEST))
        autosem = np.zeros(np.size(TEST))
        for i in np.arange(0,np.size(TEST)):
            automean[i] = np.mean(TEST[i])
            autosem[i]  = np.std(TEST[i])/np.sqrt(np.size(TEST[i]))

        clrs = sns.xkcd_palette(['sky blue', 'pinky red'])
        figx, ax1 = plt.subplots(ncols =1)
        g = plt.bar([1,2],automean, color=clrs,axes = ax1)
        ax1.set_xticks([1,2])
        ax1.set_xticklabels(['RSU','FS'])
        ax1.errorbar([1,2],automean,autosem,linestyle = 'None')
        ax1.set_ylabel('Mean FR Corr Coeff (R) btwn States')
        figx.show()
