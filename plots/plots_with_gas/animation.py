import imageio
images=[]

plot_path="./" #creating a path to store the plots only if it does not exist


p=[600,630,660,690]
for img in p:
    images.append(imageio.imread(plot_path+"snap"+str(img) + '.png'))
    #os.remove(plot_path+"snap"+str(k) + ".png")
imageio.mimsave(plot_path+"animation.gif", images, duration = 1/1)
     