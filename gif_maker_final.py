import imageio

def gifMaker(frameNo, baseName, gifName):

	filenames = []
	for i in range(frameNo):
		if i < 10:
		filenames.append(baseName + str(i) + ".png")

	images = []
	for filename in filenames:
		print filename
		images.append(imageio.imread(filename))

	imageio.mimsave(gifName, images)

