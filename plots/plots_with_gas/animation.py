import imageio
images=[]

plot_path="./" #creating a path to store the plots only if it does not exist


snapshot_start=596
snapshot_end=690
for img in range(snapshot_start,snapshot_end+1):
    print(img)
    images.append(imageio.imread(plot_path+"17clusters_snap"+str(img) + '.png'))
    #os.remove(plot_path+"snap"+str(k) + ".png")
imageio.mimsave(plot_path+"17clusters_animation.gif", images, duration = 1/5)
     