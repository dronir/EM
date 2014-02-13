#!/usr/bin/python

import tempfile as tf
import sys
import os

from string import Template

if(len(sys.argv) < 2):
  print "Usage: previewMedium options mediumFileName"
  sys.exit()

if '--quick' in sys.argv:
  transmission = 'transparent'
else:
  transmission = 'opaque'


templateFile = open('/home/hannu/work/Code/em/python/previewTemplate.rib','r')
ribTemplate  = Template(templateFile.read())
ribString    = ribTemplate.substitute(mediumFile = os.getcwd() + "/" + sys.argv[-1], nSamples=64.0,transmission=transmission)

#print ribString

tFiled = tf.mkstemp(suffix=".rib")
tFile  = open(tFiled[1],'w')
tFile.write(ribString)
tFile.close()

os.system("renderdl -d -progress " + tFile.name)

os.remove(tFiled[1])
