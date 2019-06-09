#!/usr/bin/env python3
import pytz
import os
import datetime

path        = '/home/icerx-vm/ICERX/hf_data'
samp_rate   = 48000
center_freq = 14095600
channel     = 'WSPR_14095'
sTime       = datetime.datetime(2019,1,5,1,0,tzinfo=pytz.utc)
secs        = 60
eTime       = sTime + datetime.timedelta(seconds=secs)
samp_1      = samp_rate*secs

#cmd = 'drf_plot.py -i /home/icerx-vm/ICERX/hf_data -a 2018-12-04T20:08:00 -o 14095000 -p specgram -c WSPR_14095 -r 0:8000000 -b 1024 -l -s specgram.png'
#cmd = 'drf_plot.py -i {!s} -a {!s} -o {:.0f} -p specgram -c {!s} -r 0:{:.0f} -b 1024 -l -s specgram.png'.format(path,sTime.isoformat(),center_freq,channel,samp_1)
#cmd = 'drf_plot.py -i {!s} -a {!s} -o {:.0f} -p specgram -r 0:{:.0f} -b 1024 -l -s specgram.png'.format(path,sTime.isoformat(),center_freq,samp_1)

dct = {}
dct['sTime']    = sTime.isoformat()
dct['eTime']    = eTime.isoformat()
dct['path']     = '/home/icerx-vm/ICERX/hf_data/'
dct['channel']  = '{!s}:0'.format(center_freq)
#dct['channel']  = 'mooses:1'

#cmd = 'drf_sti.py  -p {path} -c {channel} -s {sTime} -e {eTime} -n 1024 sti.png'.format(**dct)
cmd = './drf_sti.py  -p {path} -c {channel} -s {sTime} -e {eTime} -o sti.png'.format(**dct)
print(cmd)
os.system(cmd)
