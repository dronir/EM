#!/usr/bin/python

import sys
import os

xrf = "/home/hannu/work/Projects/BepiColombo/XR/bin/xrfiphs"
fIn = "/home/hannu/work/Projects/BepiColombo/XR/tests/xrfiphs/TiTube.in"

os.system("echo blah")
os.spawnl(os.P_WAIT, xrf+" "+ fIn)
