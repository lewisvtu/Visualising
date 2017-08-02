import imageio
frameNo = 360

filenames = []
for i in range(frameNo):
	filenames.append("spline_tangential_" + str(i+1) + ".png")


images = []
for filename in filenames:
	print filename
	images.append(imageio.imread(filename))
imageio.mimsave('movie.gif', images)