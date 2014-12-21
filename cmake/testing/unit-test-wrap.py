#!/usr/bin/env python

import sys
import subprocess

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--build-dir',
                  type='string',
                  action='store',
                  default=None)
parser.add_option('--name',
                  type='string',
                  action='store',
                  default=None)
(options, args) = parser.parse_args()

test_dir = '%s/unit-tests/%s' % (options.build_dir, options.name)

# run test
command = 'cd %s; ./%s.x' % (test_dir, options.name)
p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(out, err) = p.communicate()

# store stdout and stderr
f = open('%s/%s.log' % (test_dir, options.name), 'w')
f.write(out)
f.close()
f = open('%s/%s.err' % (test_dir, options.name), 'w')
f.write(err)
f.close()

# print output to screen
f = open('%s/%s.log' % (test_dir, options.name), 'r')
for line in f.readlines():
  print line.strip('\n')
print '\n'
f.close()

# check if unit test passed
last_line = file('%s/%s.log' % (test_dir, options.name), 'r').readlines()[-1]

if 'OK' in last_line:
    sys.exit(0)
else:
    sys.exit(1)
