#!/usr/bin/python3
# ----------------------------------------------------------------------------
# Copyright (c) 2017 Massachusetts Institute of Technology (MIT)
# All rights reserved.
#
# Distributed under the terms of the BSD 3-clause license.
#
# The full license is in the LICENSE file, distributed with this software.
# ----------------------------------------------------------------------------
"""Create a spectral time intensity summary plot for a data set."""


import datetime
import optparse
import os
import string
import sys
import time
import traceback

import dateutil
import digital_rf as drf
import matplotlib.gridspec
import matplotlib.mlab
import matplotlib.pyplot
import matplotlib.pyplot as plt

import numpy
import numpy.fft
import pytz
import scipy
import scipy.signal

from collections import OrderedDict

fDict = OrderedDict()
fDict['14670000Hz'] = {'cfreq':14670000, 'label':'CHU\n 14670.0 kHz'}
fDict['14095600Hz'] = {'cfreq':14095600, 'label':'Ham\n 14095.6 kHz'}
fDict[ '7850000Hz'] = {'cfreq': 7850000, 'label':'CHU\n 7850.0 kHz'}
fDict[ '7038600Hz'] = {'cfreq': 7038600, 'label':'Ham\n 7038.6 kHz'}
fDict[ '3330000Hz'] = {'cfreq': 3330000, 'label':'CHU\n 3330.0 kHz'}
fDict[ '3568600Hz'] = {'cfreq': 3568600, 'label':'Ham\n 3568.6 kHz'}

class DataPlotter(object):
    def __init__(self, control):
        """Initialize a data plotter for STI plotting."""
        self.control = control

        # open digital RF path
        self.dio    = drf.DigitalRFReader(self.control.path)

        self.sTime  = dateutil.parser.parse(self.control.start)
        self.eTime  = dateutil.parser.parse(self.control.end)

        # Figure setup
        scale = 1.5
        self.f = matplotlib.pyplot.figure(figsize=(scale*7, scale*numpy.min([numpy.max([4, self.control.frames]), 7])), dpi=128)
        self.gridspec = matplotlib.gridspec.GridSpec(self.control.frames, 1)

        self.subplots = []
        """ Setup the subplots for this display """
        for n in numpy.arange(self.control.frames):
            ax = self.f.add_subplot(self.gridspec[n])
            self.subplots.append(ax)

        self.frame = 0

    def samp2date(self,samp):
        return datetime.datetime.utcfromtimestamp(samp/self.sr)

    def plot(self,control):
        """Iterate over the data set and plot the STI into the subplot panels.

        Each panel is divided into a provided number of bins of a given
        integration length. Strides between the panels are made between
        integrations.

        """
        self.control        = control

        ch                  = self.control.channel.split(':')
        self.channel        = ch[0]
        self.sub_channel    = int(ch[1])

        # initialize outside the loop to avoid memory leak
        # initial plotting scales
        vmin = 0
        vmax = 0

        sr      = self.dio.get_properties(self.channel)['samples_per_second']
        self.sr = sr

        # initial time info
        b = self.dio.get_bounds(self.channel)
        if self.control.verbose:
            print('Channel:        ' ,self.control.channel)
            print('Sample Rate:    ', sr)
            print('Channel Bounds: ', b)
            print('Channel Bounds: ', self.samp2date(b[0]),self.samp2date(b[1]))

        if self.control.start:
            dtst0 = dateutil.parser.parse(self.control.start)
            st0 = (dtst0 - datetime.datetime(1970, 1, 1, tzinfo=pytz.utc)).total_seconds()
            st0 = int(st0 * sr)
        else:
            st0 = int(b[0])

        if self.control.end:
            dtet0 = dateutil.parser.parse(self.control.end)
            et0 = (dtet0 - datetime.datetime(1970, 1, 1, tzinfo=pytz.utc)).total_seconds()
            et0 = int(et0 * sr)
        else:
            et0 = int(b[1])

        if self.control.verbose:
            # Samples since Unix Epoch (1970 Jan 1)
            print('Plot Sample Start st0: ', st0,self.samp2date(st0))
            print('Plot Sample End   et0: ', et0,self.samp2date(et0))

        blocks              = self.control.bins # Number of time bins
        samples_per_stripe  = self.control.num_fft * self.control.integration * self.control.decimation
        total_samples       = blocks * samples_per_stripe

        if total_samples > (et0 - st0):
            print('Insufficient samples for %d samples per stripe and %d blocks between %ld and %ld' % (samples_per_stripe, blocks, st0, et0))
            return

        stripe_stride       = (et0 - st0) / blocks
        bin_stride          = stripe_stride / self.control.bins

        start_sample        = st0

        # get metadata
        # this could be done better to ensure we catch frequency or sample rate
        # changes
#        mdt = self.dio.read_metadata(st0, et0, self.channel)
#        print(t1-t0)
#        try:
#            md = mdt[list(mdt.keys())[0]]
#            cfreq = md['center_frequencies'].ravel()[self.sub_channel]
#        except (IndexError, KeyError):
#            cfreq = 0.0

        cfreq   = fDict[self.channel].get('cfreq')

        if self.control.verbose:
            print('Processing Info: Frame: {!s}/{!s} Bins: {!s} samples_per_stripe: {!s} bin_stride: {!s}'.format(
                self.frame,self.control.frames, self.control.bins, samples_per_stripe, bin_stride))

        sti_psd_data    = numpy.zeros([self.control.num_fft, self.control.bins], numpy.float)
        sti_times       = numpy.zeros([self.control.bins], numpy.complex128)

        good_data    = False
        for b in numpy.arange(self.control.bins):
            if self.control.verbose:
                print('Read Vector :', self.channel, self.samp2date(start_sample), start_sample, samples_per_stripe)
            sti_times[b] = start_sample / sr

            try:
                data = self.dio.read_vector(start_sample, samples_per_stripe, self.channel, self.sub_channel)
            except:
                start_sample += stripe_stride
                continue

            good_data = True

            if self.control.decimation > 1:
                data = scipy.signal.decimate(data, self.control.decimation)
                sample_freq = sr / self.control.decimation
            else:
                sample_freq = sr

            if self.control.mean:
                detrend_fn = matplotlib.mlab.detrend_mean
            else:
                detrend_fn = matplotlib.mlab.detrend_none

            try:
                psd_data, freq_axis = matplotlib.mlab.psd(
                    data, NFFT=self.control.num_fft, Fs=float(sample_freq), detrend=detrend_fn, scale_by_freq=False)
            except:
                traceback.print_exc(file=sys.stdout)

            sti_psd_data[:, b] = numpy.real( 10.0 * numpy.log10(numpy.abs(psd_data) + 1E-12))


            start_sample += stripe_stride

        # Now Plot the Data
        ax = self.subplots[self.frame]

        if good_data:
            # determine image x-y extent
            extent = (
                0,
                self.control.bins,
                numpy.min(freq_axis) / 1e3,
                numpy.max(freq_axis) / 1e3,
            )

            # determine image color extent in log scale units
            Pss = sti_psd_data
            if self.control.zaxis:
                vmin = int(string.split(self.control.zaxis, ':')[0])
                vmax = int(string.split(self.control.zaxis, ':')[1])
            else:
                vmin = numpy.real(numpy.median(Pss) - 6.0)
                vmax = numpy.real(numpy.median(Pss) + (numpy.max(Pss) - numpy.median(Pss)) * 0.61803398875 + 50.0)

            samp_t0 = matplotlib.dates.date2num(datetime.datetime.utcfromtimestamp(numpy.real(sti_times[0])))
            samp_t1 = matplotlib.dates.date2num(datetime.datetime.utcfromtimestamp(numpy.real(sti_times[-1])))
            _extent  = [samp_t0, samp_t1, extent[2], extent[3]]
            im = ax.imshow(sti_psd_data, cmap='jet', origin='lower', extent=_extent,interpolation='nearest', vmin=vmin, vmax=vmax, aspect='auto')

            plt.sca(ax)
            plt.colorbar(im,orientation='vertical')
            self.im     = im

        ylabel  = fDict[self.channel]['label']
        ax.set_ylabel(ylabel)
        ax.set_xlim(self.sTime,self.eTime)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H%M"))

        if self.frame != self.control.frames-1:
            xtls = ax.get_xticklabels()
            for xtl in xtls:
                xtl.set_visible(False)

        if self.frame == self.control.frames-1:
            ax.set_xlabel('Time [UT]')

        self.frame += 1

    def save_figure(self):
        ll = []
        ll.append('{!s} - {!s}'.format(self.sTime.strftime('%Y %b %d %H%M UT'),self.eTime.strftime('%Y %b %d %H%M UT')))
        ll.append('{!s} ({!s})'.format(self.control.title, self.control.path))
#        self.f.suptitle('\n'.join(ll),va='bottom')
        self.f.text(0.5,1.0,'\n'.join(ll),ha='center',fontdict={'size':'large'})
    
        self.gridspec.update()

        self.f.tight_layout()

#        self.f.subplots_adjust(top=0.95, right=0.88)
#        cax = self.f.add_axes([0.9, 0.12, 0.02, 0.80])
#        self.f.colorbar(self.im, cax=cax)

        print("Save plot as {}".format(self.control.outname))
        matplotlib.pyplot.savefig(self.control.outname,bbox_inches='tight')
        matplotlib.pyplot.close(self.f)


def parse_command_line(str_input=None):
    parser = optparse.OptionParser()
    parser.add_option("-t", "--title",      dest="title",       default='Arrival Heights / McMurdo',    help="Use title provided for the data.")
    parser.add_option("-s", "--start",      dest="start",       default=None,                           help="Use the provided start time instead of the first time in the data. format is ISO8601: 2015-11-01T15:24:00Z")
    parser.add_option("-e", "--end",        dest="end",         default=None,                           help="Use the provided end time for the plot. format is ISO8601: 2015-11-01T15:24:00Z")
    parser.add_option("-p", "--path",       dest="path",                                                help="Use data from the provided digital RF data <path>.")
    parser.add_option("-c", "--channel",    dest="channel",     default="ch0:0",                        help="Use data from the provided digital RF channel <channel>:<subchannel>.")
    parser.add_option("-l", "--length",     dest="length",      default=0.04,       type="float",           help="The default data length in seconds for unframed data.")
    parser.add_option("-b", "--bins",       dest="bins",        default=128,        type="int",             help="The number of time bins for the STI.")
    parser.add_option("-f", "--frames",     dest="frames",      default=1,          type="int",             help="The number of sub-panel frames in the plot.")
    parser.add_option("-n", "--num_fft",    dest="num_fft",     default=1024,       type="int",             help="The number of FFT bints for the STI.")
    parser.add_option("-i", "--integration",dest="integration", default=1,          type="int",             help="The number of rasters to integrate for each plot.")
    parser.add_option("-d", "--decimation", dest="decimation",  default=1,          type="int",             help="The decimation factor for the data (integer).")
    parser.add_option("-z", "--zaxis",      dest="zaxis",       default=None,       type="string",          help="zaxis colorbar setting e.g. -50:50")
    parser.add_option("-o", "--outname",    dest="outname",     default='sti.png',  type=str,               help="Name of file that figure will be saved under.")
    parser.add_option("-m", "--mean",       dest="mean",        default=False,      action="store_true",    help="Remove the mean from the data at the PSD processing step.")
    parser.add_option("-v", "--verbose",    dest="verbose",     default=True,       action="store_true",    help="Print status messages to stdout.")
    parser.add_option("-a", "--appear",     dest="appear",      default=False,      action="store_true",    help="Makes the plot appear through pyplot show.")

    if str_input is None:
        (options, args) = parser.parse_args()
    else:
        (options, args) = parser.parse_args(str_input)

    return (options, args)


if __name__ == "__main__":
    """
        Needed to add main function to use outside functions outside of module.
    """
#    sDate       = datetime.datetime(2019,1,5,tzinfo=pytz.utc)
#    eDate       = datetime.datetime(2019,1,9,tzinfo=pytz.utc)

    sDate       = datetime.datetime(2019,1,3,tzinfo=pytz.utc)
    eDate       = datetime.datetime(2019,1,22,tzinfo=pytz.utc)

    dates       = [sDate]
    while dates[-1] < eDate:
        dates.append(dates[-1]+datetime.timedelta(days=1))

    cfreqs      = []
    cfreqs.append('14670000')
    cfreqs.append('14095600')
    cfreqs.append('7850000')
    cfreqs.append('7038600')
    cfreqs.append('3330000')
    cfreqs.append('3568600')


    str_input       = '-p /home/icerx-vm/ICERX/hf_data/'.split()

    # Parse the Command Line for configuration
    (options, args) = parse_command_line(str_input)

#    dates   = dates[0:1]
    for sDate in dates:
        options.start   = sDate.isoformat()
        options.end     = (sDate + datetime.timedelta(days=1)).isoformat()
        options.outname = '{!s}.png'.format(sDate.strftime('%Y%m%d'))
        options.frames  = len(fDict)

        # Activate the DataPlotter
        dpc = DataPlotter(options)
        for ch in fDict.keys():
            options.channel = '{!s}:0'.format(ch)
            dpc.plot(options)

        dpc.save_figure()
