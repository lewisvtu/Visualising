import imageio
frameNo = 30

filenames = []
for i in range(frameNo):
	if i < 10:
		filenames.append("gas_00000" + str(i) + ".png")
	else:
		filenames.append("gas_0000" + str(i) + ".png")



images = []
for filename in filenames:
	print filename
	images.append(imageio.imread(filename))
imageio.mimsave('small_time.gif', images)