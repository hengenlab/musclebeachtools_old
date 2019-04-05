import os
import numpy as np
import h5py
import time
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches
import statsmodels.stats.api as sms
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
from colour import Color
import pdb


def spikeshapedist(datfile,xsort,ysort,qlim,zsort=None,filt=0,filtval=None,saveflag=None,savedir=None):
    """ Plot the bivariate waveform shape distribution using two or three sorting variable (xsort and ysort, optionally zsort), typicaly neg_pos_time and tailslope, but others (i.e. halfwidth) may be appropriate. Limit the quality <= qlim (e.g. only include single units by setting qlim to 2), and use a filtering variable (a single integer or binary flag such as contstat or deprived) to only include certain cells that meet this condition. This will generate a series of plots if in two dimensions, or a single 3D scatter plot if zsort is included. """
    from scipy.stats import kendalltau

    xvals   = []
    yvals   = []

    if zsort:
        zvals = []
    else:
        zvals = None

    # go through the data file and check to see how many of your cells match the initial sorting variables - collect each of the x and y sorting variable values for each of the cells that pass the quality and filt tests
    for dat in datfile:

        f       = h5py.File(dat,'r')
        l1      = list(f.keys())

        if l1[0] == 'neurons':
            neurons = f['neurons']
                    # if you need to check what variables are available in the members of neurons, run this line:
                    # datfields = list(neurons[list(neurons.keys())[0]].keys())
            ncells  = len(list(neurons))
            xvtmp   = np.zeros((ncells))
            yvtmp   = np.zeros((ncells))
            if zsort:
                zvtmp   = np.zeros((ncells))

            for ee in range(0, ncells):

                thiscell    = neurons['neuron_'+str(ee)]

                #check quality of the neuron if filt is ON
                if filt:
                    if thiscell['quality'][0]<=qlim and thiscell[filt][0]==filtval:
                        print('Including neuron that meets {} filtering'.format(filt))
                        xvtmp[ee]   = thiscell[xsort][0]
                        yvtmp[ee]   = thiscell[ysort][0]
                        if zsort:
                            zvtmp[ee] = thiscell[zsort][0]
                # check quality of the neuron if filt is OFF:
                else:
                    if thiscell['quality'][0]<=qlim:
                        xvtmp[ee]   = thiscell[xsort][0]
                        yvtmp[ee]   = thiscell[ysort][0]
                        if zsort:
                            zvtmp[ee] = thiscell[zsort][0]

            zkills  = np.where(xvtmp == 0)
            xvtmp   = np.delete(xvtmp,zkills)
            yvtmp   = np.delete(yvtmp,zkills)
            if zsort:
                zvtmp = np.delete(zvtmp,zkills)
            xvals.append(xvtmp)
            yvals.append(yvtmp)
            if zsort:
                zvals.append(zvtmp)

    # turn the array of arrays into a single array:
    xvals   = np.asarray([j for i in xvals for j in i])
    yvals   = np.asarray([j for i in yvals for j in i])
    if zsort:
        zvals   = np.asarray([j for i in zvals for j in i])

    sns.set(style="ticks")
    plt.ion()

    xlabl   = xsort.replace("_"," ")
    ylabl   = ysort.replace("_"," ")
    if zsort:
        zlabl  = zsort.replace("_"," ")

    if not zvals:
        ax = sns.jointplot(xvals, yvals, kind="hex", stat_func=kendalltau, color="#4CB391")
        plt.show()

        kdeplot = sns.jointplot(xvals, yvals, kind="kde")
        scat    = sns.jointplot(xvals, yvals)

        g = sns.jointplot(xvals, yvals, kind="kde",alpha=0.7) #color="#4CB391"
        g.plot_joint(plt.scatter, c='m', s=30, linewidth=1, marker="+", alpha=0.2)
        g.ax_joint.collections[0].set_alpha(0)
        g.set_axis_labels("${}$".format(xlabl), "${}$".format(ylabl))
    else:
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(xvals, yvals, zvals, c = 'b', marker='o')
        ax.set_xlabel(xlabl)
        ax.set_ylabel(ylabl)
        ax.set_zlabel(zlabl)
        figtext = (xlabl+ '_'+ ylabl+ '_'+ zlabl)
        fig.savefig(savedir + figtext + '_3D_scatter.pdf')

    plt.show()

    if saveflag:
        figtext = input('ID text please!! ')
        kdeplot.savefig(savedir+ figtext + '_SpikeShape_KDEplot.pdf')
        scat.savefig(savedir+figtext+'_SpikeShape_Scatterplot.pdf')
        ax.savefig(savedir+figtext+'_SpikeShape_hexPlot.pdf')
        g.savefig(savedir+figtext+'_SpikeShape_Scat_KDEplot.pdf')


def FSRSUdescript(datfile,xsort=None,xlims=None,ysort=None,ylims=None,filt1=None,filt1val=None,filt2=None,filt2val=None,qlim=2,name1=None,name2=None,preoff=None,saveflag=None,savedir=None,verbose=False):
    '''This will produce plots of  ISI CV, mean FR, and peak FR for FS and RSU populations. Filtering variable, 'filt1' and 'filt2' can either be 'cont_stat' or 'deprived', respectively. Do NOT switch these assignments. If you choose neither of those, you'll look at all cells, but this will automatically skip deprived times and off periods. '''

# ##------------------------------------------------------------------------------
# ## to run this independenty for hard core debugging, uncomment this and run it (you'll still have to set up some of your own variables though)
#     filt1 = None
#     filt2 = None
#     preoff= None
#     xsort = 'neg_pos_time'
#     ysort = 'tailSlope'
#     qlim  = 2
#     # negpos time
#     RSUxlims 	= (0.5, 1)
#     FSxlims 	= (0, 0.4)
#
#     # tailslope
#     RSUylims 	= (20, 200)
#     FSylims 	= (-300, 0.0)
#
#     xlims       = [RSUxlims,FSxlims]
#     ylims       = [RSUylims,FSylims]
#     saveflag = None
#     import time
#     import pandas as pd
# ##------------------------------------------------------------------------------

    import scipy.stats as stats

    if saveflag and not savedir:
        saveflag = 0


    columns = ['FR mean', 'FR max', 'ISI cv']
    FS      = pd.DataFrame(columns=columns)
    RSU     = pd.DataFrame(columns=columns)

    for dat in datfile:

        f       = h5py.File(dat,'r')
        l1      = list(f.keys())

        if l1[0] == 'neurons':
            neurons = f['neurons']
                    # if you need to check what variables are available in the members of neurons, run this line:
                    # datfields = list(neurons[list(neurons.keys())[0]].keys())

            ncells  = np.size(neurons)
            group1  = np.zeros((ncells,3))
            group2  = np.zeros((ncells,3))

            count1  = 0
            count2  = 0
            for ee in np.arange(0, ncells):
                print('ee is {}'.format(ee))

                thiscell    = neurons['neuron_'+str(ee)]

                # deal with weird packaging of on/off times
                if np.shape(thiscell['onTime']): # only proceed if this isn't empty

                    onts = thiscell['onTime'][:]
                    offts = thiscell['offTime'][:]

                    if np.size(np.shape(offts)) == 2:
                        offts = np.squeeze(offts)
                    else:
                        offts = offts

                    if np.size(np.shape(onts)) == 2:
                        onts = np.squeeze(onts)
                    else:
                        onts = onts
                else: # if it is empty
                    onts = []

                # don't bother processing cells from the deprived hemisphere that weren't on for at least 6h of baseline (but this should only apply to cells if you're interested in the continuous cells). Also skip cells of quality "4":
                # first line checks to see that cell has an ON TIME and that it's quality is at least 3.
                if thiscell['quality'][0] < 3 and np.size(onts) != 0:

                    try:
                        # if the first on time is over 48h AND it's a deprived cell, skip the cell.
                        quickskip = filt1 == 'cont_stat' and onts[0] > (48*3600) and thiscell['deprived'][0] == 1
                        if quickskip:
                            if verbose: print('failed ontime/continuous/deprived test')
                    except:
                        print('@ Line 198')
                        pdb.set_trace()
                # if those filters fail, set this thing to "skip"
                else:
                    if verbose: print('Failed quality and/or on times test.')
                    quickskip = True

                if xlims[0][1]>thiscell[xsort][0]>xlims[0][0] and ylims[0][1]>thiscell[ysort][0]>ylims[0][0]:
                    celltype = 1 # this is an RSU
                elif xlims[1][1]>thiscell[xsort][0]>xlims[1][0] and ylims[1][1]>thiscell[ysort][0]>ylims[1][0]:
                    celltype = 2 # this is an FSs cell
                else:
                    if verbose: print('No cell type assignment.')
                    celltype = 999 # this doesn't fit in either category. Outlier.
                    quickskip = True

                # tell me who the cell is for skipping...
                if quickskip:
                    print('Skipping a quality {} cell, neuron {}'.format(thiscell['quality'][0], ee))

                if not quickskip:

                    filt1dat    = 999
                    #  on/off times to reflect an arbitrary endpoint before the end of the recording:
                    if preoff:
                        try:
                            newontimes  = onts[onts<preoff]
                            newofftimes = offts[offts<preoff]
                        except ValueError:
                            print('Caught at line 197! Value error')
                            pdb.set_trace()

                        if np.size(newontimes) - np.size(newofftimes) == 1:
                            newofftimes = np.append(newofftimes,preoff)

                    elif preoff is None:
                        # this is the case in which there isn't a "pre" ending time:
                        newontimes  = onts
                        newofftimes = offts

                    # Update the cont_stat for "pre" ending:
                    if filt1 is 'cont_stat' and preoff:
                        ender    = np.min([preoff, offts[-1]])
                        expttott = ender - onts[0]

                        # Check to see if the cell is continuous with updated "end" time:
                        if expttott/preoff>=0.8:
                            # This updates fields on cells that are continuous
                            filt1dat    = 1
                        else:
                            # NOT continuous
                            filt1dat    = 0

                    elif filt1 is 'cont_stat' and not preoff:
                        filt1dat = thiscell['cont_stat'][0]

                    print('at line 226')
                    # decide whether to keep neurons based on the various combinations of filt1 and filt 2 (none, one, or both)
                    keeps = None
                    if not (not filt1 or filt2):
                        if thiscell['quality'][0]<=qlim and filt1dat==filt1val:
                            if verbose: print('here1')
                            keeps = 1
                            print('Neuron {} from {} that meets {} filtering'.format( ee, np.array(thiscell['animal'][:], dtype=np.int8).tostring().decode("ascii"), filt1) )
                        else:
                            keeps = None

                    elif filt1 and filt2:
                        if thiscell['quality'][0]<=qlim and filt1dat==filt1val and thiscell[filt2][0]==filt2val:
                            if verbose: print('here2')
                            keeps = 1
                            print('Including neuron that meets {} and {} filtering'.format(filt1, filt2))
                        else:
                            keeps = None

                    elif not filt1 and not filt2:
                        if thiscell['quality'][0]<=qlim:
                            if verbose: print('here3')
                            keeps = 1
                        else:
                            keeps = None

                    # use this to eliminate bogus cells that aren't worth fixing at this time. Make whatever rule you need so that the junk cells trigger "skipbs" assignment. Cells will then have to be "keeps" as well as NOT "skipbs"
                    skipbs = 0
                    if (np.shape(newofftimes)[0] == np.shape(newontimes)[0] == 1) and newofftimes[0] < newontimes[0]:
                        skipbs = 1

                    if keeps and not skipbs:
                        if verbose: print('here4')
                        print('Keeping a {} quality cell, cont stat is {}'.format(thiscell['quality'][0], thiscell['cont_stat'][0]))
                            # Find the peak firing rate in bins of binwidth:
                        frmax       = []
                        frmean      = []
                        binwidth    = 1
                        edges       = []
                        try:
                            edges   = np.arange(newontimes[0],newofftimes[-1],binwidth)
                        except IndexError:
                            print('Line 144 in FScellproject.py')
                            pdb.set_trace()
                        print('at 261')

                        # If the cell is in the deprived hemisphere, only calculate the mean firing rate and ISI CV based on the "on" time before lid suture (this assumes lid suture occuring at approximately 60h of recording), otherwise use the 'rate_CELL' field of the saved data structure and all "on" time ISIs.
                        tcheck0     = time.time()
                        if thiscell['deprived'][0] == 1:
                            if verbose: print('here5')
                            # add 60h as an off time
                            temp    = newontimes[newontimes<(60*3600)]
                            temp2   = newofftimes[newofftimes<(60*3600)]
                            temp2   = np.append(temp2,60*3600)

                            thiscv  = []
                            tcount  = 0
                            #dcount  = 0
                            times   = thiscell['time'][:]
                            #dtemp   = np.zeros(np.size(times))
                            ns      = np.zeros(np.size(temp))
                            tpass   = np.zeros(np.size(temp))

                            # go through the on times
                            for gg in temp: # gg is the on time
                                ot             = temp2[tcount] # ot is the off time in order of gg
                                # Mean FR:
                                ns[tcount]     = np.sum((thiscell['time']>gg)*(thiscell['time']<ot )) # n spikes in this window
                                tpass[tcount]  = ot-gg # elapsed time

                                # CV:
                                #timevar = times[ (times<ot)*(times>gg) ]
                                #diffvar = np.diff(timevar)
                                #dtemp[dcount:dcount+np.size(diffvar)] = diffvar

                                bins = np.arange(gg,ot,600) # divide this block up by ten minute segments

                                if verbose: print('here5.1')

                                for i in np.arange(0,np.size(bins)-1): # go through all of the ten minute bins above to find the cv
                                    t1      = np.argmax( times > bins[i] )
                                    t2      = np.argmin( times < bins[i+1] )-1
                                    isi     = np.diff( times[t1:t2] )
                                    cv      = np.nanstd( isi )/np.nanmean( isi )
                                    thiscv  = np.append( thiscv,cv )
                                    if verbose: print('here5.2')

                                # #use this if you want to look at individual cells that make it into the data you're working with:
                                # #could add in a programmatic elimination of ISIs over a certain value
                                #print('line 264')
                                #pdb.set_trace()
                                # #add to the counters tracking the on/off times and the number of ISIs
                                #dcount += np.size(diffvar)
                                tcount +=1

                            if verbose: print('here5.3')

                            thiscv = np.nanmean(thiscv)
                            frmean = np.sum(ns)/np.sum(tpass)
                            if verbose: print('This cv is {} \nThis frmean is {}'.format(thiscv, frmean))
                            if verbose: print('here5.4')

                        # and now for control hemisphere cells:
                        elif thiscell['deprived'][0] == 0:

                            if verbose: print('here6')
                            #frmean      = thiscell['rate_CELL'][0] # this is a cheap way to get these data...

                            # calculate the CV of the ISI distribution for this cell:
                            thiscv      = []
                            times       = thiscell['time'][:]

                            # Calculate firing rates and CV for control hemisphere:
                            xcount  = 0
                            ns_c      = np.zeros(np.size(newontimes))
                            et   = np.zeros(np.size(newontimes))

                            for xx in np.arange(0,np.size(newontimes)):

                                et[xcount] = newofftimes[xx] - newontimes[xx] # elapsed time
                                ns_c[xcount] = np.sum((thiscell['time']>newontimes[xx])*(thiscell['time']<newofftimes[xx] )) # number of spikes in this window

                                bins = np.arange(newontimes[xx],newofftimes[xx],600)
                                xcount += 1

                                for i in np.arange(0,np.size(bins)-1):

                                    t1      = np.argmax( times > bins[i] )
                                    t2      = np.argmin(times < bins[i+1] )-1
                                    isi     = np.diff(times[t1:t2])
                                    cv      = np.nanstd(isi)/np.nanmean(isi)
                                    thiscv  = np.append(thiscv,cv)


                            frmean = np.sum(ns_c)/np.sum(et)
                            thiscv = np.nanmean(thiscv)
                            if verbose: print('here7')
                            if verbose: print('This cv is {} \nThis frmean is {}'.format(thiscv, frmean))

                        tcheck1 = time.time() - tcheck0

                        if verbose: print('here8')
                        print('Processing step 1 took {} seconds'.format(tcheck1))

                        maxtemp = (np.histogram(thiscell['time'],edges)[0])/binwidth
                        mxmean  = np.mean(maxtemp)
                        mxstd   = np.std(maxtemp)
                        frmax   = mxmean + 3*mxstd
                        #frmax   = np.max(np.histogram(thiscell['time'],edges)[0])/binwidth

                        # C.V. is standard deviation divided by the mean
                        print('This CV is {}. '.format(thiscv))

                        if verbose: print('here9')
                        if celltype == 1:

                            if verbose: print('here10')
                            # This is an RSU
                            group1[count1,0] = frmean
                            group1[count1,1] = frmax
                            group1[count1,2] = thiscv
                            count1 += 1

                        elif celltype == 2:

                            if verbose: print('here11')
                            # This is an FS cell
                            group2[count2,0] = frmean
                            group2[count2,1] = frmax
                            group2[count2,2] = thiscv
                            count2 += 1
                    else:
                        if verbose: print('here12')
                        # NOT 'keeps'
                        pass




            if np.sum(group1[:,1])>0:
                if verbose: print('here13')
                RSUend      = np.max(np.where(group1[:,1] == next(bb for bb in reversed(group1[:,1]) if bb>0 )  )[0])  +1

                rsutempdf   = pd.DataFrame(data = group1[0:RSUend,:],columns=columns)
                RSU         = RSU.append(rsutempdf)

            else:
                if verbose: print('here14')
                RSUend      = None

            if np.sum(group2[:,1])>0:
                if verbose: print('here15')
                FSend       = np.max(np.where(group2[:,1] == next(bb for bb in reversed(group2[:,1]) if bb>0 )  ) [0]) +1

                fstempdf    = pd.DataFrame(data = group2[0:FSend,:],columns=columns)
                FS          = FS.append(fstempdf)

            else:
                if verbose: print('here16')
                FSend       = None

    # Plot everything...
    print('line 367')

    # Make violin plots first. So far, these look cool but they're not very informative in these datasets
    plt.ion()
    sns.set(style="ticks")
    clrs = sns.xkcd_palette(['sky blue', 'pinky red'])

    fig10 = plt.figure(10,figsize=(12, 6), dpi=100)
    ax1 = plt.subplot(1,3,1)
    ax1.set_xlabel('Cell type')
    ax1.set_ylabel('Firing Rate (Hz)')
    ax1.set_title('Mean Firing Rate')
    sns.violinplot(data=[RSU['FR mean'],FS['FR mean']], palette="Set1", bw=.3, cut=1, linewidth=1, inner="stick")
    sns.despine(left=True, bottom=True)
    ylim = ax1.get_ylim()[1]
    ax1.set_ylim([0,ylim])
    ax1.set_xticklabels(['RSU','FS'])

    ax2 = plt.subplot(1,3,2)
    ax2.set_xlabel('Cell type')
    ax2.set_ylabel('Firing Rate (Hz)')
    ax2.set_title('Max Firing Rate')
    sns.violinplot(data=[RSU['FR max'],FS['FR max']], palette="Set1", bw=.3, cut=1, linewidth=1, inner="stick")
    sns.despine(left=True, bottom=True)
    ylim = ax2.get_ylim()[1]
    ax2.set_ylim([0,ylim])
    ax2.set_xticklabels(['RSU','FS'])

    ax3 = plt.subplot(1,3,3)
    ax3.set_xlabel('Cell type')
    ax3.set_ylabel('Coefficient of Variation')
    ax3.set_title('ISI CV')
    sns.violinplot(data=[RSU['ISI cv'],FS['ISI cv']], palette="Set1", bw=.3, cut=1, linewidth=1, inner="stick")
    sns.despine(left=True, bottom=True)
    ylim = ax3.get_ylim()[1]
    ax3.set_ylim([0,ylim])
    ax3.set_xticklabels(['RSU','FS'])

    plt.show()

    # Make a two color comparison (bar plot) of the mean/max FR and CV of FS and RSU
    fig11 = plt.figure(11,figsize=(12, 6), dpi=100)
    sns.set(style="ticks")

    # This is the mean firing rate bar plot
    ax1     = plt.subplot(1,3,1)
    x       = [1,2]
    y       = [np.mean(RSU['FR mean']), np.mean(FS['FR mean'])]
    yerr    = [ stats.sem(RSU['FR mean'],nan_policy='omit'), stats.sem(FS['FR mean'],nan_policy = 'omit') ]# omit the NaN entries so that the SEM stats runs smoothly
    ax1.set_ylim( [ 0, 1.5 * np.max( np.add( y, yerr ) ) ] )
    ax1.bar(x,y,color=clrs)
    ax1.errorbar(x,y,yerr=yerr, ls='none')
    plt.xticks(x, ['RSU','FS'],rotation='vertical')
    sns.despine(offset=10, trim=True)
    ax1.set_title('Mean Firing Rate')
    ax1.set_ylabel('Hz')


    frstats = stats.ttest_ind(RSU['FR mean'], FS['FR mean'], nan_policy = 'omit')
    ypos = 1.1 * np.ceil(np.max( np.add( y, yerr ) ) )
    s = ('p = {}'.format('%.5f' % frstats.pvalue))
    ax1.text(1, ypos, s, fontdict=None, withdash=False)

    # this is the max firing rate bar plot
    ax2     = plt.subplot(1,3,2)
    x       = [1,2]
    y       = [np.mean(RSU['FR max']), np.mean(FS['FR max'])]
    yerr    = [ stats.sem(RSU['FR max'], nan_policy = 'omit'), stats.sem(FS['FR max'],nan_policy = 'omit') ]  # mask the NaN entries so that the SEM stats runs smoothly
    ax2.set_ylim( [ 0, 1.5 * np.ceil(np.max( np.add( y, yerr ) ) )] )
    ax2.bar(x,y,color=clrs)
    ax2.errorbar(x,y,yerr=yerr, ls='none')
    plt.xticks(x, ['RSU','FS'],rotation='vertical')
    sns.despine(offset=10, trim=True)
    ax2.set_title('Max Firing Rate in {} sec bins'.format(binwidth))
    ax2.set_ylabel('Hz')

    frmaxstats = stats.ttest_ind(RSU['FR max'], FS['FR max'], nan_policy = 'omit')
    ypos = 1.1 * np.ceil(np.max( np.add( y, yerr ) ) )
    s = ('p = {}'.format('%.5f' % frmaxstats.pvalue))
    ax2.text(1, ypos, s, fontdict=None, withdash=False)

    # this is the CV bar plot
    ax3     = plt.subplot(1,3,3)
    x       = [1,2]
    y       = [np.mean(RSU['ISI cv']), np.mean(FS['ISI cv'])]
    yerr    = [ stats.sem(RSU['ISI cv'], nan_policy = 'omit'), stats.sem( FS['ISI cv'], nan_policy = 'omit' ) ]# mask the NaN entries so that the SEM stats runs smoothly
    ax3.set_ylim( [ 0, 1.5 * np.max( np.add( y, yerr ) ) ] )
    ax3.bar(x,y,color=clrs)
    ax3.errorbar(x,y,yerr=yerr, ls='none')
    plt.xticks(x, ['RSU','FS'],rotation='vertical')
    sns.despine(offset=10, trim=True)
    ax3.set_title('ISI CV'.format(binwidth))
    ax3.set_ylabel('Coef Var')

    cvstats = stats.ttest_ind(RSU['ISI cv'], FS['ISI cv'], nan_policy = 'omit')
    ypos = 1.1 * np.ceil(np.max( np.add( y, yerr ) ) )
    s = ('p = {}'.format('%.5f' % cvstats.pvalue))
    ax3.text(1, ypos, s, fontdict=None, withdash=False)


    plt.show()

    if saveflag:
        fig10.savefig(savedir+'FSRSU_violin.pdf')
        fig11.savefig(savedir+'FSRSU_meanbar.pdf')

def unitylineplot(nnet,time0=(36*3600, 42*3600),time05=(54*3600, 60*3600), time1=(40*3600, 46*3600),time2=(96*3600, 108*3600),time3=(180*3600,198*3600),logscale=1, savedir=None):

    # figure out if you're dealing with deprived or control cells:
    temp = [e.deprived for e in nnet.neurons]
    if np.unique(temp)[0] == 1:
        dep = 'Depr'
    elif np.unique(temp)[0] == 0:
        dep = 'Control'

    if savedir:
        figtitle    = input('Figure and file unique ID? ')

    celltype = input('What type of cell are you working on?    ')

    foldmtx     = np.zeros((np.size(nnet.neurons),5))
    count       = 0
    nrncount    = 0
    for ee in nnet.neurons:

        if np.size(np.shape(ee.offTime)) == 2:
            offts = np.squeeze(ee.offTime)
        else:
            offts = ee.offTime

        if np.size(np.shape(ee.onTime)) == 2:
            onts = np.squeeze(ee.onTime)
        else:
            onts = ee.onTime

        # Check that the neuron is online?
        if onts[0]<time1[1] and offts[-1]>time3[0]:

            # Time bin 0 - baseline startpoint
            s0 = np.sum((time0[0] < ee.time) & (ee.time < time0[1]))

            #if any((time1[0] < ee.onTime) & (ee.onTime < time1[1])):
            if np.any( np.logical_and( time0[0] < onts, onts < time0[1] ) ):
                thison  = onts[np.where((time0[0] < onts) & (onts < time0[1]))]
                thisoff = offts[np.where((time0[0] < offts) & (offts < time0[1]))]

                if np.size(thison) == np.size(thisoff)+1:
                    thisoff = np.append(time0[0],thisoff)

                elif np.size(thison) == np.size(thisoff)-1:
                    thison = np.append(thison,time0[1])

                tlost = np.zeros(np.size(thison))
                for bads in np.arange(0,np.size(thison)):
                    #print(bads)
                    tlost[bads] = thison[bads] - thisoff[bads]
                tlost = sum(tlost)

                unfixs0 = s0/(time0[1]-time0[0])
                s0 /= (time0[1]-time0[0])-tlost
                #print('Fixing s0 from {} to {}'.format(unfixs0,s0))

            else:
                s0 /= (time0[1]-time0[0])

            # Time bin 05 - baseline drift
            s05 = np.sum((time05[0] < ee.time) & (ee.time < time05[1]))

            if np.any(np.logical_and(time05[0] < onts, onts < time05[1])):
                thison  = onts[np.where((time05[0] < onts) & (onts < time05[1]))]
                thisoff = offts[np.where((time05[0] < offts) & (offts < time05[1]))]

                if np.size(thison) == np.size(thisoff)+1:
                    thisoff = np.append(time1[0],thisoff)

                elif np.size(thison) == np.size(thisoff)-1:
                    thison = np.append(thison,time1[1])

                tlost = np.zeros(np.size(thison))
                for bads in np.arange(0,np.size(thison)):
                    #print(bads)
                    tlost[bads] = thison[bads] - thisoff[bads]
                tlost = sum(tlost)

                s05 /= (time05[1]-time05[0])-tlost

            else:
                s05 /= (time05[1]-time05[0])

            # Time bin 1 - Full Baseline Sample
            s1 = np.sum((time1[0] < ee.time) & (ee.time < time1[1]))

            #if any((time1[0] < ee.onTime) & (ee.onTime < time1[1])):
            if np.any(np.logical_and(time1[0] < onts, onts < time1[1])):
                thison  = onts[np.where((time1[0] < onts) & (onts < time1[1]))]
                thisoff = offts[np.where((time1[0] < offts) & (offts < time1[1]))]

                if np.size(thison) == np.size(thisoff)+1:
                    thisoff = np.append(time1[0],thisoff)

                elif np.size(thison) == np.size(thisoff)-1:
                    thison = np.append(thison,time1[1])

                tlost = np.zeros(np.size(thison))
                for bads in np.arange(0,np.size(thison)):
                    #print(bads)
                    tlost[bads] = thison[bads] - thisoff[bads]
                tlost = sum(tlost)

                unfixs0 = s1/(time1[1]-time1[0])
                s1 /= (time1[1]-time1[0])-tlost
                #print('Fixing s0 from {} to {}'.format(unfixs0,s0))

            else:
                s1 /= (time1[1]-time1[0])

            if np.log10(s1)<-18:
                print('gotcha {}'.format(nrncount))
                pdb.set_trace()
            else:
                pass

            # Time bin 2 - early MD
            s2 = np.sum((time2[0] < ee.time) & (ee.time < time2[1]))

            if np.any( np.logical_and( time2[0] < onts, onts < time2[1] ) ):

                thison  = onts[np.where((time2[0] < onts) & (onts < time2[1]))]
                thisoff = offts[np.where((time2[0] < offts) & (offts < time2[1]))]

                if np.size(thison) == np.size(thisoff)+1:
                    thisoff = np.append(time2[0],thisoff)

                elif np.size(thison) == np.size(thisoff)-1:
                    thison = np.append(thison,time2[1])


                tlost = np.zeros(np.size(thison))
                for bads in np.arange(0,np.size(thison)):
                    #print(bads)
                    tlost[bads] = thison[bads] - thisoff[bads]
                tlost = sum(tlost)

                unfixs2 = s2/(time2[1]-time2[0])
                s2 /= (time2[1]-time2[0])-tlost
                #print('Fixing s2 from {} to {}'.format(unfixs2,s2))

            else:
                s2 /= (time2[1]-time2[0])

            # Time bin 3 - late MD
            s3 = np.sum((time3[0] < ee.time) & (ee.time < time3[1]))

            if (time3[0] < offts).any() & (offts < time3[1]).any():

                thison  = onts[np.where((time3[0] < onts) & (onts < time3[1]))]
                thisoff = offts[np.where((time3[0] < offts) & (offts < time3[1]))]

                if np.size(thison) == np.size(thisoff)+1:
                    thisoff = np.append(time3[0],thisoff)

                elif np.size(thison) == np.size(thisoff)-1:
                    thison = np.append(thison,time3[1])


                tlost = np.zeros(np.size(thison))
                for bads in np.arange(0,np.size(thison)):
                    #print(bads)
                    tlost[bads] = thison[bads] - thisoff[bads]
                tlost = sum(tlost)

                unfixs3 = s3/(time3[1]-time3[0])
                s3 /= (time3[1]-time3[0])-tlost
                #print('Fixing s3 from {} to {}'.format(unfixs3,s3))

            else:

                if offts[-1]>time3[1]:
                    s3 /= (time3[1]-time3[0])

                else:
                    dt  = ee.offTime[-1] - time3[0]
                    s3 /= dt

            foldmtx[count,:] = [s0, s05, s1, s2, s3]
            count+=1
        nrncount+=1

    # get rid of zero entries
    fldkill = np.where(np.sum(foldmtx,axis=1)==0)
    foldmtx = np.delete(foldmtx,fldkill,axis = 0)

    if logscale:
        logfoldmtx = np.log10(foldmtx)

    # load into pandas
    df = pd.DataFrame(logfoldmtx, columns=['EarlyBL','LateBL','FullBaseLine', 'MD1', 'Rebound'])

    # PLOTTING -----------------------------------------------------------------

    plt.ion()
        # tangerine #ff9408 - early MD, cornflower #6a79f7 - late MD, aqua green #12e193 - baseline2
    sns.set(style="darkgrid", color_codes=True)
    try:
        ax0 = sns.jointplot("EarlyBL", "LateBL", data=df, kind="reg", color='#12e193', size=10)
        ax1 = sns.jointplot("FullBaseLine", "MD1", data=df, kind="reg", color="#ff9408", size=10)
        ax2 = sns.jointplot("FullBaseLine", "Rebound", data=df, kind="reg", color="#6a79f7", size=10)
    except:
        print('Stuck at 635 in FScellsproject.py. Solve this problem, genius! ')
        pdb.set_trace()


    # if the user doesn't specify the limits for plotting, autoscale it
    ax0.ax_joint.autoscale(enable=True, axis='both', tight=None)
    ax1.ax_joint.autoscale(enable=True, axis='both', tight=None)
    ax2.ax_joint.autoscale(enable=True, axis='both', tight=None)

    # given the above scaling, go figure out where the unity line should sit
    axmax = np.round(np.max(logfoldmtx)*2)/2
    axmin = np.round(np.floor(np.min(logfoldmtx)*2))/2


    ax0.ax_joint.plot([axmin,axmax], np.poly1d(np.polyfit([axmin,axmax], [axmin,axmax], 1))([axmin,axmax]), color = '#738595',linewidth = 3, ls = '--', alpha = 0.7)
    ax1.ax_joint.plot([axmin,axmax], np.poly1d(np.polyfit([axmin,axmax], [axmin,axmax], 1))([axmin,axmax]), color = '#738595',linewidth = 3, ls = '--', alpha = 0.7)
    ax2.ax_joint.plot([axmin,axmax], np.poly1d(np.polyfit([axmin,axmax], [axmin,axmax], 1))([axmin,axmax]), color = '#738595',linewidth = 3, ls = '--', alpha = 0.7)

    ax0.ax_joint.set_xlabel('Early Base')
    ax0.ax_joint.set_ylabel('Late Base')
    ax1.ax_joint.set_xlabel('Baseline')
    ax1.ax_joint.set_ylabel('Early MD')
    ax2.ax_joint.set_xlabel('Baseline')
    ax2.ax_joint.set_ylabel('Late MD' )

    ax0.ax_joint.set_title('{} {} Baseline Drift'.format(celltype, dep))
    ax1.ax_joint.set_title('{} {} Early MD Drift'.format(celltype, dep))
    ax2.ax_joint.set_title('{} {} Late MD Drift'.format(celltype, dep))

    # Write in the time periods:
    ax0.ax_joint.text(ax0.ax_joint.get_xlim()[1]*0.001,ax0.ax_joint.get_ylim()[1]*0.8, 'BL Early {} to {}'.format(time0[0]/3600, time0[1]/3600) )
    ax0.ax_joint.text(ax0.ax_joint.get_xlim()[1]*0.001,ax0.ax_joint.get_ylim()[1]*0.7, 'BL Late {} to {}'.format(time05[0]/3600, time05[1]/3600) )

    ax1.ax_joint.text(ax0.ax_joint.get_xlim()[1]*0.001,ax0.ax_joint.get_ylim()[1]*0.8, 'BL {} to {}'.format(time1[0]/3600, time1[1]/3600) )
    ax1.ax_joint.text(ax0.ax_joint.get_xlim()[1]*0.001,ax0.ax_joint.get_ylim()[1]*0.7, 'MD Early {} to {}'.format(time2[0]/3600, time2[1]/3600) )

    ax2.ax_joint.text(ax0.ax_joint.get_xlim()[1]*0.001,ax0.ax_joint.get_ylim()[1]*0.8, 'BL {} to {}'.format(time1[0]/3600, time1[1]/3600) )
    ax2.ax_joint.text(ax0.ax_joint.get_xlim()[1]*0.001,ax0.ax_joint.get_ylim()[1]*0.7, 'MD Late {} to {}'.format(time3[0]/3600, time3[1]/3600) )

    if savedir:
        ax0.savefig(savedir + os.sep + figtitle + '{}_{}_Baselineunity.pdf'.format(celltype,dep))
        ax1.savefig(savedir + os.sep + figtitle + '{}_{}_earlyMDunity.pdf'.format(celltype,dep))
        ax2.savefig(savedir + os.sep + figtitle + '{}_{}_lateMDunity.pdf'.format(celltype,dep))
    else:
        pass

    # Calculate the percent drift across the time periods (from baseline)
    drift0      = 100*np.abs((foldmtx[:,1]-foldmtx[:,0])/foldmtx[:,0]) # Baseline drift from ealy to late BL
    drift1      = 100*np.abs((foldmtx[:,3]-foldmtx[:,2])/foldmtx[:,2]) # Early MD drift from full baseline
    drift2      = 100*np.abs((foldmtx[:,4]-foldmtx[:,2])/foldmtx[:,2]) # Late MD drift from full baseline
    semdrift0   = np.std(drift0)/np.sqrt(np.size(drift0))
    semdrift1   = np.std(drift1)/np.sqrt(np.size(drift1))
    semdrift2   = np.std(drift1)/np.sqrt(np.size(drift2))

    #clrs = ['#ff9408', "#6a79f7"]
    # tangerine #ff9408 - early MD, cornflower #6a79f7 - late MD, aqua green #12e193 - baseline2
    clrs = ['#12e193', '#ff9408','#6a79f7']
    fig100, ax100 = plt.subplots(1)
    sns.set(style="dark")
    sns.barplot([1,2,3],[np.mean(drift0), np.mean(drift1), np.mean(drift2)],palette=clrs,yerr=[semdrift0, semdrift1, semdrift2], ax = ax100)
    ax100.set_xticklabels(['Baseline','EarlyMD', 'LateMD'])
    ax100.set_ylabel('Mean % Change from Early Baseline')
    ax100.text(0,ax100.get_ylim()[1]*0.8, 'BL Early {} : {}\nBL Late {} : {}\nFull BL {} : {}\nMD Early {} : {}\nMD Late {} : {}'.format(time0[0]/3600, time0[1]/3600,time05[0]/3600, time05[1]/3600,time1[0]/3600,time1[1]/3600,time2[0]/3600,time2[1]/3600,time3[0]/3600,time3[1]/3600) )
    ax100.set_title('{} {} Mean abs fold change from Early BL'.format(celltype, dep))

    if savedir:
        print('saving')
        fig100.savefig(savedir + os.sep + figtitle + '{}_{}_AbsFold_bar.pdf'.format(celltype,dep))
    else:
        pass

    # Plot a filled kernel density estimate

    pv0 = stats.ks_2samp(drift0,drift1)[1]*3
    pv1 = stats.ks_2samp(drift0,drift2)[1]*3
    pv2 = stats.ks_2samp(drift1,drift2)[1]*3


    fig101, ax101 = plt.subplots(1)
    sns.distplot(drift0, hist=False, color='#12e193', kde_kws={"shade": True}, ax = ax101 ) # Fill
    sns.distplot(drift0, hist=False, color='#738595', kde_kws={"shade": False}, ax = ax101 ) # Edge

    sns.distplot(drift1, hist=False, color="#ff9408", kde_kws={"shade": True}, ax = ax101) # Fill
    sns.distplot(drift1, hist=False, color="#738595", kde_kws={"shade": False}, ax = ax101) # Edge

    sns.distplot(drift2, hist=False, color="#6a79f7", kde_kws={"shade": True}, ax = ax101) # Fill
    sns.distplot(drift2, hist=False, color="#738595", kde_kws={"shade": False}, ax = ax101) # Edge
    ax101.set_xlabel('Percent Change from Baseline')
    ax101.set_ylabel('Kernel Density Estimate')
    ax101.set_title('{} {} % Change KDE'.format(celltype, dep))
    ax101.set_xlim([-75, 275])

    ax101.text(ax101.get_xlim()[1]*0.5,ax101.get_ylim()[1]*0.7, 'BL Early {} : {}\nBL Late {} : {}\nFull BL {} : {}\nMD Early {} : {}\nMD Late {} : {}'.format(time0[0]/3600, time0[1]/3600,time05[0]/3600, time05[1]/3600,time1[0]/3600,time1[1]/3600,time2[0]/3600,time2[1]/3600,time3[0]/3600,time3[1]/3600) )

    # add stats text information:
    if pv0<0.05:
        ax101.text(ax101.get_xlim()[1]*0.5,ax101.get_ylim()[1]*0.5, 'KS base v eMD: p = {}'.format(np.around(pv0,4)))
    else:
        ax101.text(ax101.get_xlim()[1]*0.5,ax101.get_ylim()[1]*0.5, 'KS base v eMD: ns')
    if pv1<0.05:
        ax101.text(ax101.get_xlim()[1]*0.5,ax101.get_ylim()[1]*0.45, 'KS base v lMD: p = {}'.format(np.around(pv1,4)))
    else:
        ax101.text(ax101.get_xlim()[1]*0.5,ax101.get_ylim()[1]*0.45, 'KS base v lMD: ns')
    if pv2<0.05:
        ax101.text(ax101.get_xlim()[1]*0.5,ax101.get_ylim()[1]*0.4, 'KS eMD v lMD: p = {}'.format(np.around(pv2,4)))
    else:
        ax101.text(ax101.get_xlim()[1]*0.5,ax101.get_ylim()[1]*0.4, 'KS eMD v lMD: ns')


    if savedir:
        fig101.savefig(savedir + os.sep + figtitle + '{}_{}_AbsFold_KDE.pdf'.format(celltype,dep))
    else:
        pass

    return (ax0, ax1, ax2, ax100, ax101, drift0, drift1, drift2)

def unityblobs(nnet,time1=(50*3600, 60*3600), time2=(110*3600, 140*3600),time3=(190*3600,200*3600),logscale=1,savedir=None,lilbin = 300, dispersion = 'quartiles'):

    if savedir:
        figtitle    = input('Figure and file title? ')

    t1mtx = np.zeros( [np.size(nnet.neurons), 3 ])
    t2mtx = np.zeros( [np.size(nnet.neurons), 3 ])
    t3mtx = np.zeros( [np.size(nnet.neurons), 3 ])

    count = 0
    for ee in nnet.neurons:

        if np.size(np.shape(ee.offTime)) == 2:
            offts = np.squeeze(ee.offTime)
        else:
            offts = ee.offTime

        if np.size(np.shape(ee.onTime)) == 2:
            onts = np.squeeze(ee.onTime)
        else:
            onts = ee.onTime

        wholedge    = np.arange(0,ee.time[-1],lilbin)
        wholex      = np.histogram(ee.time,wholedge)[0]/lilbin

        # Convert bins before initial ON time and after final OFF time to NAN:
        prenans             = np.where(wholedge<=onts[0])[0][-1] + 1 # Find the last bin edge BEFORE the initial on time. Add 1 to account for the difference in the number of edges and the number of bins
        postnans            = np.where(wholedge<offts[-1])[0][-1]
        wholex[0:prenans]   = np.nan
        wholex[postnans:]   = np.nan

        # if there are intermediate on/off times, deal with converting those bins to nan here:
        midonts = None
        midoffs = None
        if np.size(onts) > 1:
            midonts = onts[1:]
            midoffs = offts[0:-1]

            for mo in np.arange(0,np.size(midonts)):
                offn   = np.where(wholedge<midoffs[mo])[0][-1]  # off occurs first mid experiment.
                onn    = np.where(wholedge<midonts[mo])[0][-1]  # then on resumes recording.
                wholex[offn:onn+1] = np.nan # need to add 1 here because of how python indexes

        # Get BIG BIN stats here for making plots:
        t1on            = np.where(wholedge<time1[0])[0][-1]
        t1off           = np.where(wholedge<time1[1])[0][-1] + 1
        test1           = wholex[t1on:t1off]
        t1mtx[count,0]  = np.nanmedian(test1) # mean in column 1
        t1nonan         = test1
        t1nonan         = np.delete( t1nonan, np.where(np.isnan(t1nonan) ) )
        if dispersion == 'stdev':
            t1mtx[count,1]  = np.nanmedian(test1) - np.nanstd(test1) # get stdev
            t1mtx[count,2]  = np.nanmedian(test1) - np.nanstd(test1) # get stdev
        elif dispersion == 'quartiles':
            t1mtx[count,1]  = np.percentile(t1nonan,25) # get first quartile of FR range
            t1mtx[count,2]  = np.percentile(t1nonan,75) # get third quartile of FR range
        elif dispersion == 'confint':
            t1mtx[count,1:] = sms.DescrStatsW(t1nonan).tconfint_mean() # get the 95% confidence interval

        t2on            = np.where(wholedge<time2[0])[0][-1]
        t2off           = np.where(wholedge<time2[1])[0][-1] + 1
        test2           = wholex[t2on:t2off]
        t2mtx[count,0]  = np.nanmedian(test2)
        t2nonan         = test2
        t2nonan         = np.delete( t2nonan, np.where(np.isnan(t2nonan) ) )
        if dispersion == 'stdev':
            t2mtx[count,1]  = np.nanmedian(test2) - np.nanstd(test2) # get stdev
            t2mtx[count,2]  = np.nanmedian(test2) - np.nanstd(test2) # get stdev
        elif dispersion == 'quartiles':
            t2mtx[count,1]  = np.percentile(t2nonan,25) # get first quartile of FR range
            t2mtx[count,2]  = np.percentile(t2nonan,75) # get third quartile of FR range
        elif dispersion == 'confint':
            t2mtx[count,1:] = sms.DescrStatsW(t2nonan).tconfint_mean() # get the 95% confidence interval

        # If the cell isn't online for at least 6h prior to the end of the time3 big bin, take the offtime of the cell and then use 1/2 of the time3 big bin duration to look at the firing rate mean and distribution at the end of the cell's recorded period. Use that as a stand in for the time3 big bin.
        if offts[-1]>(time3[0]+(3600*6)):
            t3on            = np.where(wholedge<time3[0])[0][-1]
            t3off           = np.where(wholedge<time3[1])[0][-1] + 1
            test3           = wholex[t3on:t3off]
            t3mtx[count,0]  = np.nanmedian(test3)
            t3nonan         = test3
            t3nonan         = np.delete( t3nonan, np.where(np.isnan(t3nonan) ) )
            if dispersion == 'stdev':
                t3mtx[count,1]  = np.nanmedian(test3) - np.nanstd(test3) # get stdev
                t3mtx[count,2]  = np.nanmedian(test3) - np.nanstd(test3) # get stdev
            elif dispersion == 'quartiles':
                t3mtx[count,1]  = np.percentile(t3nonan,25) # get first quartile of FR range
                t3mtx[count,2]  = np.percentile(t3nonan,75) # get third quartile of FR range
            elif dispersion == 'confint':
                t3mtx[count,1:] = sms.DescrStatsW(t3nonan).tconfint_mean() # get the 95% confidence interval
        else:
            t3dur           = time3[1] - time3[0]
            t3tempon        = offts[-1] - 0.5*t3dur
            t3on            = np.where(wholedge<t3tempon)[0][-1]
            t3off           = np.where(wholedge<offts[-1])[0][-1] + 1
            test3           = wholex[t3on:t3off]
            t3mtx[count,0]  = np.nanmedian(test3)
            t3nonan         = test3
            t3nonan         = np.delete( t3nonan, np.where(np.isnan(t3nonan) ) )
            if dispersion == 'stdev':
                t3mtx[count,1]  = np.nanmedian(test3) - np.nanstd(test3) # get stdev
                t3mtx[count,2]  = np.nanmedian(test3) - np.nanstd(test3) # get stdev
            elif dispersion == 'quartiles':
                t3mtx[count,1]  = np.percentile(t3nonan,25) # get first quartile of FR range
                t3mtx[count,2]  = np.percentile(t3nonan,75) # get third quartile of FR range
            elif dispersion == 'confint':
                t3mtx[count,1:] = sms.DescrStatsW(t3nonan).tconfint_mean()

        count+=1

    # SORT BY FIRING RATE....
    srt     = np.argsort(t1mtx[:,0])

    # Plotting
    fig     = plt.figure()
    ax      = fig.add_subplot(111,aspect='equal')

    plt.ion()
    count1 = 0
    for e in srt:

        cent1    = t1mtx[e,0]
        rad1     = t1mtx[e,2] - cent1
        circ1    = plt.Circle( (count1,cent1), color = 'xkcd:seafoam', radius=rad1, fill=True, alpha = 0.7)

        cent2    = t2mtx[e,0]
        rad2     = t2mtx[e,2] - cent2
        circ2    = plt.Circle( (count1,cent2), color = 'xkcd:tangerine', radius=rad2, fill=True, alpha = 0.7)

        cent3    = t3mtx[e,0]
        rad3     = t3mtx[e,2] - cent3
        circ3    = plt.Circle( (count1,cent3), color = 'xkcd:cornflower', radius=rad3, fill=True, alpha = 0.7)

        ax.add_patch(circ1)
        ax.add_patch(circ2)
        ax.add_patch(circ3)

        count1 +=1


    baseline_patch  = mpatches.Patch(color='xkcd:seafoam', label='Baseline')
    earlymd_patch   = mpatches.Patch(color='xkcd:tangerine', label='Early MD')
    latemd_patch    = mpatches.Patch(color='xkcd:cornflower', label='Late MD')
    plt.legend(handles=[baseline_patch,earlymd_patch,latemd_patch])

    if logscale == 1:
        ax.set_yscale('log')
    ax.autoscale(enable=True, axis='both', tight=None)
    ax.set_xlabel('Cell rank by baseline FR')
    ax.set_ylabel('Firing Rate (Hz)')
    ax.set_title(figtitle)
    plt.show()

    if savedir:
        fig.savefig(savedir + os.sep + figtitle + 'BLOBDRIFT.pdf')

def kdeblobdrift(nnet, time1=(50*3600, 60*3600), time2=(110*3600, 140*3600), time3=(190*3600,200*3600), logscale=1, savedir=None, lilbin = 300, dispersion = 'KDE'):

    if savedir:
        figtitle    = input('Figure and file title? ')

    nbins1 = np.size(np.arange(time1[0],time1[1],lilbin))
    nbins2 = np.size(np.arange(time2[0],time2[1],lilbin))
    nbins3 = np.size(np.arange(time3[0],time3[1],lilbin))

    t1mtx = np.zeros( [np.size(nnet.neurons), nbins1 ])
    t2mtx = np.zeros( [np.size(nnet.neurons), nbins2 ])
    t3mtx = np.zeros( [np.size(nnet.neurons), nbins3 ])


    count = 0
    for ee in nnet.neurons:

        if np.size(np.shape(ee.offTime)) == 2:
            offts = np.squeeze(ee.offTime)
        else:
            offts = ee.offTime

        if np.size(np.shape(ee.onTime)) == 2:
            onts = np.squeeze(ee.onTime)
        else:
            onts = ee.onTime

        wholedge    = np.arange(0,ee.time[-1],lilbin) # make lil bin edges through whole experiment
        wholex      = np.histogram(ee.time,wholedge)[0]/lilbin # make bin FR counts through the whole experiment

        # Convert bins before initial ON time and after final OFF time to NAN:
        prenans             = np.where(wholedge<=onts[0])[0][-1] + 1 # Find the last bin edge BEFORE the initial on time. Add 1 to account for the difference in the number of edges and the number of bins
        postnans            = np.where(wholedge<offts[-1])[0][-1]
        wholex[0:prenans]   = np.nan
        wholex[postnans:]   = np.nan

        # if there are intermediate on/off times, deal with converting those bins to nan here:
        midonts = None
        midoffs = None
        if np.size(onts) > 1:
            midonts = onts[1:]
            midoffs = offts[0:-1]

            for mo in np.arange(0,np.size(midonts)):
                offn   = np.where(wholedge<midoffs[mo])[0][-1]  # off occurs first mid experiment.
                onn    = np.where(wholedge<midonts[mo])[0][-1]  # then on resumes recording.
                wholex[offn:onn+1] = np.nan # need to add 1 here because of how python indexes

        # Get BIG BIN stats here for making plots:
        t1on            = np.where(wholedge<time1[0])[0][-1] +1
        t1off           = np.where(wholedge<time1[1])[0][-1] +1
        t1mtx[count,:]  = np.nan
        tmp1            = wholex[t1on:t1off][~np.isnan(wholex[t1on:t1off])]
        t1mtx[count,0:np.size(tmp1)]  = tmp1

        t2on            = np.where(wholedge<time2[0])[0][-1] +1
        t2off           = np.where(wholedge<time2[1])[0][-1] +1
        t2mtx[count,:]  = np.nan
        tmp2            = wholex[t2on:t2off][~np.isnan(wholex[t2on:t2off])]
        t2mtx[count,0:np.size(tmp2)]  = tmp2

        # If the cell isn't online for at least 6h prior to the end of the time3 big bin, take the offtime of the cell and then use 1/2 of the time3 big bin duration to look at the firing rate mean and distribution at the end of the cell's recorded period. Use that as a stand in for the time3 big bin.
        if offts[-1]>time3[0]:
            t3on            = np.where(wholedge<time3[0])[0][-1] +1
            t3off           = np.where(wholedge<time3[1])[0][-1] +1

            if np.size(wholedge[t3on:t3off]/3600) == nbins3:
                t3mtx[count,:]  = wholex[t3on:t3off]
            else: # just trust me, this was necessary for an obnoxious bullshit thing.
                t3mtx[count,:]  = np.nan # by hitting this else statement, you don't have enough data to fill the entire array. start with all nans, then overwrite the parts that have data. (see next line)
                t3mtx[count,0:np.size(wholex[t3on:t3off])]  = wholex[t3on:t3off]

        else:
            t3dur           = time3[1] - time3[0]
            t3tempon        = offts[-1] - 0.5*t3dur
            t3on            = np.where(wholedge<t3tempon)[0][-1] +1
            t3off           = np.where(wholedge<offts[-1])[0][-1] +1
            t3mtx[count,:]  = np.nan # by hitting this else statement, you don't have enough data to fill the entire array. start with all nans, then overwrite the parts that have data. (see next line)
            t3mtx[count,0:np.size(wholex[t3on:t3off])]  = wholex[t3on:t3off]

        count+=1

    # eliminate cells whose drop is "offline":
    t2nans = np.where(np.isnan(np.mean(t2mtx,axis = 1)))

    t1mtx = np.delete(t1mtx, t2nans, axis = 0)
    t2mtx = np.delete(t2mtx, t2nans, axis = 0)
    t3mtx = np.delete(t3mtx, t2nans, axis = 0)


    mean1   = np.mean(t1mtx,axis = 1)
    srt1    = np.argsort(mean1)

    # select a random subset, still in order, of the FR sorted array in srt1:
    grdtot      = 25
    srtsubset   = np.sort(np.random.choice(np.arange(0,np.shape(t1mtx)[0]),grdtot,replace=False))

    # Plotting
    sqside  = np.int(np.ceil(np.sqrt(grdtot)))
    fig, ax = plt.subplots(nrows=sqside,ncols=sqside)

    axcols = np.squeeze(np.matlib.repmat(np.arange(0,sqside),1,sqside))
    axrows = np.squeeze(np.repeat(np.arange(0,sqside),sqside))

    # Make colormaps out of individual colors you want to use:
    sf = sns.light_palette("xkcd:seafoam", as_cmap=True)
    tn = sns.light_palette("xkcd:tangerine", as_cmap=True)
    cf = sns.light_palette("xkcd:cornflower", as_cmap=True)

    plt.ion()
    xval    = 10
    rng1    = 0.5
    sns.set_style(style='white')
    sns.despine()
    axcount = 0
    for e in srtsubset:

        # e is an FR ordered random subset of the cells. Tag pulls the proper cell from the sorted matrix.
        tag = srt1[e]

        plt.sca(ax[axrows[axcount],axcols[axcount]])

        #nan fix:
        t1tmp   = np.delete(t1mtx[tag,:], np.where(np.isnan(t1mtx[tag,:])))
        np.random.shuffle(t1tmp) # make more circular
        xvals0  = np.linspace(xv0,xv1,np.size(t1tmp) )

        t2tmp   = np.delete(t2mtx[tag,:], np.where(np.isnan(t2mtx[tag,:])))
        np.random.shuffle(t2tmp) # make more circular
        xvals1  = np.linspace(xv0,xv1,np.size(t2tmp) )

        t3tmp   = np.delete(t3mtx[tag,:], np.where( np.isnan(t3mtx[tag,:]) ) )
        np.random.shuffle(t3tmp) # make more circular
        xvals2  = np.linspace(xv0,xv1,np.size(t3tmp) )


        # set the xvals for the plotting:
        xv00    = xval -0.2*np.max(t1tmp)
        xv01    = xval +0.2*np.max(t1tmp)
        xvals0  = np.linspace(xv00,xv01,np.size(t1tmp) )

        xv10    = xval -0.2*np.max(t2tmp)
        xv11    = xval +0.2*np.max(t2tmp)
        xvals1  = np.linspace(xv10,xv11,np.size(t2tmp) )

        xv20    = xval -0.2*np.max(t3tmp)
        xv21    = xval +0.2*np.max(t3tmp)
        xvals2  = np.linspace(xv20,xv21,np.size(t3tmp) )

        sns.kdeplot(xvals0, t1tmp, kind="kde",cmap = sf, shade_lowest=False, alpha=0.9)
        sns.kdeplot(xvals1, t2tmp, kind="kde",cmap = tn, shade_lowest=False, alpha=0.7)
        sns.kdeplot(xvals2, t3tmp, kind="kde",cmap = cf, shade_lowest=False, alpha=0.8)

        ax[axrows[axcount],axcols[axcount]].autoscale(enable=True, axis='both', tight=True)

        axcount +=1

def kdedrift(nnet,tstart = 36*3600, tstop = 200*3600, savedir=None, windowsize = 12*3600, stepsize = 3600, BLT = (36*3600, 42*3600)):
    '''Consider drift from baseline firing rates in control hemisphere cells. Use this to plot the kde for each position as you slide a window over the dataset. this will give a sense of whether or not there's a time dependent shift in the distribution of drift from baseline firing rates.   '''

    temp = [e.deprived for e in nnet.neurons]
    if np.unique(temp)[0] == 1:
        dep = 'Depr'
    elif np.unique(temp)[0] == 0:
        dep = 'Control'

    if savedir:
        figtitle    = input('Figure and file unique ID? ')

    celltype = input('What type of cell are you working on?    ')

    edges   = np.arange(tstart,tstop-windowsize+stepsize, stepsize)
    mtrx    = np.zeros([np.size(nnet.neurons),np.size(edges) -1 ]) # take one off of the edges because the first edge is the starting position for your baseline calculation and you're not including it in the change matrix.
    df      = pd.DataFrame()

    count = 0
    for ee in nnet.neurons:
        if np.size(np.shape(ee.offTime)) == 2:
            offts = np.squeeze(ee.offTime)
        else:
            offts = ee.offTime

        if np.size(np.shape(ee.onTime)) == 2:
            onts = np.squeeze(ee.onTime)
        else:
            onts = ee.onTime

        if onts[0]<BLT[1]: # make sure cell is online here

            basetimes = ee.time[np.logical_and(ee.time>=BLT[0], ee.time<BLT[1])]

            # Set up the baseline rate:
            if np.any( np.logical_and( BLT[0] < onts, onts < BLT[1] ) ):

                thison  = onts[ np.where( (BLT[0] < onts) & (onts < BLT[1]) ) ]
                thisoff = offts[ np.where( (BLT[0] < offts) & (offts < BLT[1]) ) ]

                if np.size(thison)  == np.size(thisoff)+1:
                    thisoff = np.append(BLT[0],thisoff)

                elif np.size(thison) == np.size(thisoff)-1:
                    thison = np.append(thison,BLT[1])

                tlost = np.zeros(np.size(thison))
                for bads in np.arange(0,np.size(thison)):
                    # Delete any spikes that are recorded in off times and also sum up the time transpired in off times.
                    basetimes   = np.delete( basetimes, np.where( np.logical_and( basetimes>thisoff[bads], basetimes<thison[bads] ) ) )
                    tlost[bads] = thison[bads] - thisoff[bads]

                tlost       = sum(tlost)
                baserate    = np.size(basetimes)/( (BLT[1] - BLT[0]) -tlost)

            else:
                baserate = np.size(basetimes)/(BLT[1] - BLT[0])

            # now loop through all of the windows and calculate the difference here:
            edgcount = 0
            for t0 in edges[1:]:
                t1 = t0+windowsize

                if t1<= offts[-1]: # don't operate past the final offtime!!!
                    tmptimes = ee.time[ np.logical_and( ee.time>=t0, ee.time<t1 ) ]

                    # Set up the window rate:
                    if np.any( np.logical_and( t0 < onts, onts < t1 ) ) or np.any(np.logical_and (t0 < offts, offts < t1 )):

                        thison  = onts[ np.where( (t0 < onts) & (onts < t1) ) ]
                        thisoff = offts[ np.where( (t0 < offts) & (offts < t1) ) ]

                        if np.size(thison)  == np.size(thisoff)+1:
                            thisoff = np.append(t0,thisoff)

                        elif np.size(thison) == np.size(thisoff)-1:
                            thison = np.append(thison,t1)

                        tlost = np.zeros(np.size(thison))
                        for bads in np.arange(0,np.size(thison)):
                            # Delete any spikes that are recorded in off times and also sum up the time transpired in off times.
                            tmptimes    = np.delete( tmptimes, np.where( np.logical_and( tmptimes>thisoff[bads], tmptimes<thison[bads] ) ) )
                            tlost[bads] = thison[bads] - thisoff[bads]

                        tlost   = sum(tlost)
                        tmprate = np.size(tmptimes)/(windowsize-tlost)

                        if tlost<(0.5*windowsize):
                            mtrx[count, edgcount] = 100*np.abs((tmprate-baserate)/baserate)
                        else:
                            mtrx[count, edgcount] = 999 # go find 999's later and delete them, or replace with mean?

                    else:
                        tmprate = np.size(tmptimes)/windowsize
                        mtrx[count, edgcount] = 100*np.abs((tmprate-baserate)/baserate)

                elif t1>=offts[-1]: # how to handle the bins that slide past the final off time of a given unit
                    mtrx[count, edgcount] = 999

                if mtrx[count, edgcount]>1000:
                    pdb.set_trace()

                edgcount +=1

            count +=1

    # get rid of zero rows:
    zkills = np.where( np.sum(mtrx,axis=1) == 0)[0]
    mtrx = np.delete(mtrx,zkills,axis = 0)

    # Get rid of 999 flagged entries (convert to nan)
    mtrx[np.where(mtrx == 999)] = np.nan

    df = pd.DataFrame(mtrx)

    clr1    = Color('#75bbfd') # sky blue
    colors = None
    colors  = list(clr1.range_to(Color('#be03fd'),df.shape[1])) # bright purple

    sns.set_style('dark')
    fig103, ax103 = plt.subplots(1)
    for a in np.arange(1, df.shape[1]):
        temp    = df[a]
        temp    = temp.drop(np.where(np.isnan(temp))[0] )
        tmpclr  = colors[a]
        ccode   = tmpclr.get_hex()
        sns.distplot(temp, hist=False, color = ccode, kde_kws={"shade": False}, ax = ax103 )

    ax103.set_xlabel('Percent FR change from baseline')
    ax103.set_ylabel('KDE')
    ax103.set_title('{} {} fold change dist. in {}h bins'.format(dep,celltype,windowsize/3600))

    ax103.set_xlim([-75,375])

    if savedir:
        fig103.savefig( savedir + os.sep + figtitle + '{}_{}_KDE_drift.pdf'.format(celltype,dep) )
