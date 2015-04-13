#!/usr/bin/env python

import subprocess
import os
import random

path = os.getcwd();
anchoring = {}

f_anchoring = open('anch.in')

j = 0
for line in f_anchoring:
	j = j + 1
	anchoring[j] = line

for i in range(j):
	cmd = 'rm -r ' + str(i)
	subprocess.call(cmd, shell=True)

	#Make separate directories
	cmd = 'mkdir ' + str(i)
	subprocess.call(cmd, shell=True)

	sub = open('submit','w')
	sub_copy = open('submit.copy','r')

	for line in sub_copy:
		line = line.replace("JOB_NAME",anchoring[i+1].rstrip())
		sub.write(line)
	sub_copy.close()
	sub.close()

	cmd3 = 'cp submit ' + str(i)
	cmd5 = 'cp channel ' + str(i)
	subprocess.call(cmd3, shell=True)
	subprocess.call(cmd5, shell=True)

	read = open('param.copy','r')
	copy = open('param.in','w')
	rand_seed = random.randint(111111,9999999)

	for line in read:
		anchoring_value = anchoring[i+1].rstrip()  
		line = line.replace("TOP",anchoring_value)	
		line = line.replace("RANDOM",str(rand_seed))	
		copy.write(line)
	copy.close()
	read.close()

	cmd2 = 'cp param.in ' +str(i)
	subprocess.call(cmd2, shell=True)

	submit_path = str(i)
	os.chdir(submit_path)

	cmd7 = 'sbatch submit'
	subprocess.call(cmd7, shell=True)

	os.chdir(path)

	cmd = 'rm param.in'
	subprocess.call(cmd, shell=True)
	cmd = 'rm submit'
	subprocess.call(cmd, shell=True)

print 'DONE'
