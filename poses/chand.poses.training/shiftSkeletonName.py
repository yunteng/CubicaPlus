import os
import sys
import string
if len(sys.argv) >= 3:
	# inputName = sys.argv[1]
	# outputName = sys.argv[2]
	startFrame = int(sys.argv[1])
	endFrame = int(sys.argv[2])
	shiftAmount = int(sys.argv[3])
	print "args: startFrame ", startFrame, "endFrame ", endFrame, "shiftAmount ", shiftAmount

	frames = endFrame - startFrame + 1

	if shiftAmount > 0:
		frame = endFrame
		for i in xrange(frames):
			inputName = 'ode.motion.' + format(frame, '04d') + '.skeleton'
			outputName = 'ode.motion.' + format(frame + shiftAmount, '04d') + '.skeleton'
			command = 'mv ' + inputName + ' ' + outputName
			print command
			os.system(command)
			frame = frame - 1
	else:
		frame = startFrame
		for i in xrange(frames):
			inputName = 'ode.motion.' + format(frame, '04d') + '.skeleton'
			outputName = 'ode.motion.' + format(frame + shiftAmount, '04d') + '.skeleton'
			command = 'mv ' + inputName + ' ' + outputName
			print command
			os.system(command)
			frame = frame + 1
	# inputDataDir = 'chand.medium.3/'
	# outputDataDir = 'chand.medium.training.short/'
	# inputData = os.popen('ls ' + inputDataDir + inputName + '*').read().split('\n')
	# for inputF in inputData:
	# 	if len(inputF) > 4:
	# 		# inputF = inputDataDir + f
	# 		outputF = string.replace(inputF, inputName, outputName)
	# 		outputF = string.replace(outputF, inputDataDir, outputDataDir)
	# 		os.system('cp ' + inputF + ' ' + outputF)
			# print 'cp ' + inputF + ' ' + outputF

	
	# outputMotionDir = 'chand.poses.new/'
	# inputF = os.popen('ls ' + inputMotionDir + '*' + inputName + '*').read().split('\n')[0]
	# # print inputF
	# outputF = string.replace(inputF, inputName, outputName)
	# outputF = string.replace(outputF, inputMotionDir, outputMotionDir)
	# os.system('cp ' + inputF + ' ' + outputF)
	# print outputF
			# print 'cp ' + inputF + ' ' + outputF
	# print inputData

