import imageio
frameNo = 50

filenames = []
for i in range(frameNo):
	filenames.append("ClippedAttempt2_" + str(i+1) + ".png")


images = []
for filename in filenames:
	print filename
	images.append(imageio.imread(filename))
imageio.mimsave('movie.gif', images)